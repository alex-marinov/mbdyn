/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_SOCKET

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <strings.h>
#include <signal.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/un.h>
#include <arpa/inet.h>

#include "mbc.h"
#include "sock.h"

/* private flags for internal use */
enum sock_flags_t {
	MBC_SF_VALID = 0x1U
};

/* validate command
 *
 * command needs to be set in mbc->cmd
 */
int
mbc_check_cmd(mbc_t *mbc)
{
	switch (mbc->cmd) {
	case ES_REGULAR_DATA:
	case ES_GOTO_NEXT_STEP:
	case ES_ABORT:
	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP:
	case ES_NEGOTIATION:
	case ES_OK:
		return 0;
	}

	fprintf(stderr, "unknown cmd from peer\n");
	return -1;
}

/* get command from peer
 *
 * command is stored in mbc->cmd
 */
int
mbc_get_cmd(mbc_t *mbc)
{
	ssize_t rc;

	rc = recv(mbc->sock, (void *)&mbc->cmd, sizeof(mbc->cmd),
		mbc->recv_flags);
	if (rc != sizeof(mbc->cmd)) {
		fprintf(stderr, "recv(cmd=%u) failed\n", mbc->cmd);
		return -1;
	}

	if (mbc_check_cmd(mbc)) {
		return -1;
	}

	return 0;
}

/* put command to peer
 *
 * command needs to be set in mbc->cmd
 */
int
mbc_put_cmd(mbc_t *mbc)
{
	ssize_t	rc;

	if (mbc_check_cmd(mbc)) {
		return -1;
	}

	if (mbc->verbose) {
		fprintf(stdout, "cmd to peer: %u\n", mbc->cmd);
	}

	rc = send(mbc->sock, (const void *)&mbc->cmd, sizeof(mbc->cmd),
		mbc->send_flags);
	if (rc != sizeof(mbc->cmd)) {
		fprintf(stderr, "send(cmd=%u) failed (%d)\n", mbc->cmd, rc);
		return -1;
	}

	return 0;
}

/* initialize communication - do not call directly */
static int
mbc_init(mbc_t *mbc, struct sockaddr *addr, socklen_t socklen)
{
	unsigned long timeout = mbc->timeout*1000000;
	unsigned long useconds = 100000;

	if (mbc->sock < 0) {
		fprintf(stderr, "unable to create socket\n");
		return -1;
	}

	for ( ; ; ) {
		if (connect(mbc->sock, addr, socklen) < 0) {
			int save_errno = errno;
			const char *msg;

			if (timeout != 0) {
				switch (save_errno) {
				case ECONNREFUSED:	// inet
				case ENOENT:		// unix
					/* Socket does not exist yet; retry */
					usleep(useconds);
					if (mbc->timeout > 0) {
						timeout -= useconds;
					}
					continue;
				}
			}

			/* Connect failed */
			msg = strerror(save_errno);
			fprintf(stderr, "unable to connect to peer (%d: %s)\n",
				save_errno, msg);
			return -1;
		}

		break;
	}

	/* MSG_NOSIGNAL disables SIGPIPE
	 * MSG_WAITALL blocks recv() until all data is available
	 */
	mbc->recv_flags |= (MSG_NOSIGNAL | MSG_WAITALL);
	mbc->send_flags |= MSG_NOSIGNAL;

	mbc->sock_flags = MBC_SF_VALID;

	return 0;
}

/* initialize communication using inet socket
 *
 * mbc must be a pointer to a valid mbc_t structure
 * host and port must be defined
 */
int
mbc_inet_init(mbc_t *mbc, const char *host, short unsigned port)
{
	struct sockaddr_in addr;

	if (host == NULL) {
		fprintf(stderr, "host must be defined\n");
		return -1;
	}

	if (port == 0) {
		fprintf(stderr, "port must be defined\n");
		return -1;
	}

	mbc->sock = mbdyn_make_inet_socket(&addr, host, port, 0, NULL);

	return mbc_init(mbc, (struct sockaddr *)&addr, sizeof(addr));
}

/* initialize communication using unix socket
 *
 * mbc must be a pointer to a valid mbc_t structure
 * path must be defined
 */
int
mbc_unix_init(mbc_t *mbc, const char *path)
{
	struct sockaddr_un addr;

	if (path == NULL) {
		fprintf(stderr, "path must be defined\n");
		return -1;
	}

	mbc->sock = mbdyn_make_named_socket(&addr, path, 0, NULL);

	return mbc_init(mbc, (struct sockaddr *)&addr, sizeof(addr));
}

/* destroy communication
 *
 * does NOT free the mbc structure
 */
int
mbc_destroy(mbc_t *mbc)
{
	/* TODO: send "abort"? */
	if (mbc->sock >= 0) {
		close(mbc->sock);
	}

	return 0;
}

/*
 * nodal stuff
 */

/* get nodal motion from peer
 *
 * if mbc->rigid, access rigid motion using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->nodes > 0, access nodal motion using macros MBC_N_*
 */
int
mbc_nodal_get_motion(mbc_nodal_t *mbc)
{
	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	if (mbc->mbc.cmd == ES_ABORT) {
		fprintf(stdout, "got ABORT from peer\n");
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		if (mbc->rigid) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_R_KINEMATICS(mbc),
				MBC_R_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_R_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%u) rigid failed (%d)\n",
					MBC_R_KINEMATICS_SIZE(mbc), rc);
				return -1;
			}
		}

		if (mbc->nodes > 0) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_N_KINEMATICS(mbc),
				MBC_N_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_N_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%u) x, theta, xP, omega failed (%d)\n",
					MBC_N_KINEMATICS_SIZE(mbc), rc);
				return -1;
			}
		}
	}

	return 0;
}

/* put forces to peer
 *
 * if mbc->rigid, force and moment must be set in storage pointed to
 *	by macros MBC_F, MBC_M
 * if mbc->nodes > 0, nodal forces must be set in storage pointed to
 *	by macro MBC_N_F, MBC_N_M
 */
int
mbc_nodal_put_forces(mbc_nodal_t *mbc, int last)
{
	if (last) {
		if (mbc->mbc.data_and_next) {
			mbc->mbc.cmd = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;

		} else {
			mbc->mbc.cmd = ES_GOTO_NEXT_STEP;
		}

	} else {
		mbc->mbc.cmd = ES_REGULAR_DATA;
	}

	if (mbc_put_cmd((mbc_t *)mbc)) {
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		/* rigid */
		if (mbc->rigid) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_R_DYNAMICS(mbc),
				MBC_R_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_R_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%u) rigid failed (%d)\n",
					MBC_R_DYNAMICS_SIZE(mbc), rc);
				return -1;
			}
		}

		/* nodal */
		if (mbc->nodes > 0) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_N_DYNAMICS(mbc),
				MBC_N_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_N_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%u) nodes failed (%d)\n",
					MBC_N_DYNAMICS_SIZE(mbc), rc);
				return -1;
			}
		}
	}

	return 0;
}

/* initialize nodal data
 *
 * mbc must be a pointer to a valid mbc_nodal_t structure
 *
 * at least rigid body motion must be defined (mbc->rigid != 0),
 * or nodes must be > 0
 *
 * if nodes > 0, mallocs memory that needs to be freed calling
 * mbc_nodal_destroy()
 */
int
mbc_nodal_init(mbc_nodal_t *mbc, unsigned nodes, unsigned flags)
{
	if (!mbc->rigid && nodes == 0) {
		fprintf(stderr, "need at least 1 node or rigid body data\n");
		return -1;
	}

	mbc->nodes = nodes;
	mbc->flags = flags;

	if (mbc->nodes > 0) {
		unsigned rot = (mbc->flags & MBC_ROT_MASK);
		double *ptr;

		if (rot == 0) {
			fprintf(stderr, "need orientation parametrization in flags\n");
			return -1;
		}

		mbc->k_size = (3 + 3 + 3); /* x, xp, omega */

		switch (rot) {
		case MBC_ROT_MAT:
			mbc->k_size += 9;
			break;

		case MBC_ROT_THETA:
		case MBC_ROT_EULER_123:
			mbc->k_size += 3;
			break;

		default:
			fprintf(stderr, "unknown orientation parametrization %u in flags\n", rot);
			return -1;
		}

		if (mbc->flags & MBC_ACCELS) {
			mbc->k_size += 6;
		}

		mbc->k_size *= mbc->nodes;

		mbc->n_x = NULL;
		mbc->n_theta = NULL;
		mbc->n_r = NULL;
		mbc->n_euler_123 = NULL;
		mbc->n_xp = NULL;
		mbc->n_omega = NULL;
		mbc->n_xpp = NULL;
		mbc->n_omegap = NULL;
		mbc->n_f = NULL;
		mbc->n_m = NULL;

		ptr = (double *)malloc(MBC_N_SIZE(mbc));
		if (mbc->n_x == NULL) {
			fprintf(stderr, "nodal data malloc failed\n");
			return -1;
		}

		mbc->n_x = ptr;
		ptr += 3*nodes;

		switch (rot) {
		case MBC_ROT_MAT:
			mbc->n_r = ptr;
			ptr += 9*nodes;
			break;

		case MBC_ROT_THETA:
			mbc->n_theta = ptr;
			ptr += 3*nodes;
			break;

		case MBC_ROT_EULER_123:
			mbc->n_euler_123 = ptr;
			ptr += 3*nodes;
			break;
		}

		mbc->n_xp = ptr;
		ptr += 3*nodes;

		mbc->n_omega = ptr;
		ptr += 3*nodes;

		if (mbc->flags & MBC_ACCELS) {
			mbc->n_xpp = ptr;
			ptr += 3*nodes;

			mbc->n_omegap = ptr;
			ptr += 3*nodes;
		}

		mbc->n_f = ptr;
		ptr += 3*nodes;

		mbc->n_m = ptr;
		ptr += 3*nodes;
	}

	return 0;
}

/* negotiate nodal data
 *
 * mbc must be a pointer to a valid mbc_nodal_t structure
 *
 * at least rigid body motion must be defined (mbc->rigid != 0),
 * or nodes must be > 0
 *
 * the socket must be initialized and connected
 * sends a negotiation request to the master
 */
int
mbc_nodal_negotiate_request(mbc_nodal_t *mbc)
{
	int rc;
	uint8_t *uint8_ptr;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint8_t) + sizeof(uint32_t)];

	if (!mbc->rigid && mbc->nodes == 0) {
		fprintf(stderr, "need at least 1 node or rigid body data\n");
		return -1;
	}

	if (!(mbc->mbc.sock_flags & MBC_SF_VALID)) {
		fprintf(stderr, "socket is not valid\n");
		return -1;
	}

	mbc->mbc.cmd = ES_NEGOTIATION;
	mbc_put_cmd((mbc_t *)mbc);

	uint8_ptr = (uint8_t *)&buf[0];
	uint8_ptr[0] = (uint8_t)MBC_NODAL;
	if (mbc->rigid) {
		uint8_ptr[0] |= MBC_REF_NODE;
	}

	uint32_ptr = (uint32_t *)&uint8_ptr[1];
	uint32_ptr[0] = mbc->nodes;

	rc = send(mbc->mbc.sock, (const void *)buf, sizeof(buf),
		mbc->mbc.send_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "send negotiate request failed (%d)\n", rc);
		return -1;
	}

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	switch (mbc->mbc.cmd) {
	case ES_ABORT:
		fprintf(stdout, "got ABORT from peer\n");
		return -1;

	case ES_OK:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%u from peer\n", mbc->mbc.cmd);
		return -1;
	}

	return 0;
}

int
mbc_nodal_negotiate_response(mbc_nodal_t *mbc)
{
	int rc;
	uint8_t *uint8_ptr;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint8_t) + sizeof(uint8_t) + sizeof(uint32_t)];
	unsigned uNodal;
	unsigned uRef;

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	switch (mbc->mbc.cmd) {
	case ES_NEGOTIATION:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%u from peer\n", mbc->mbc.cmd);
		return -1;
	}

	rc = recv(mbc->mbc.sock, (void *)buf, sizeof(buf), mbc->mbc.recv_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "recv negotiate request failed\n");
		return -1;
	}

	rc = 0;

	uint8_ptr = (uint8_t *)&buf[0];
	uNodal = uint8_ptr[0] & MBC_MODAL_NODAL_MASK;
	uRef = uint8_ptr[0] & MBC_REF_NODE;
	if (uNodal != MBC_NODAL) {
		rc++;
	}

	if ((uRef && !mbc->rigid) || (!uRef && mbc->rigid)) {
		rc++;
	}

	uint32_ptr = (uint32_t *)&uint8_ptr[1];
	if (uint32_ptr[0] != mbc->nodes) {
		rc++;
	}

	if (rc) {
		mbc->mbc.cmd = ES_ABORT;

	} else {
		mbc->mbc.cmd = ES_OK;
	}

	mbc_put_cmd((mbc_t *)mbc);

	return 0;
}

/* destroy nodal data
 *
 * does NOT free the mbc structure
 */
int
mbc_nodal_destroy(mbc_nodal_t *mbc)
{
	if (mbc->n_x) {
		free(mbc->n_x);
		mbc->n_x = NULL;
		mbc->n_theta = NULL;
		mbc->n_r = NULL;
		mbc->n_euler_123 = NULL;
		mbc->n_xp = NULL;
		mbc->n_omega = NULL;
		mbc->n_xpp = NULL;
		mbc->n_omegap = NULL;
		mbc->n_f = NULL;
		mbc->n_m = NULL;
	}

	return 0;
}



/*
 * modal stuff
 */

/* get modal motion from peer
 *
 * if mbc->rigid, access rigid motion using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->modes > 0, access modal motion using macros MBC_Q, MBC_QP
 */
int
mbc_modal_get_motion(mbc_modal_t *mbc)
{
	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	if (mbc->mbc.cmd == ES_ABORT) {
		fprintf(stdout, "got ABORT from peer\n");
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		if (mbc->rigid) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_R_KINEMATICS(mbc),
				MBC_R_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_R_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%u) rigid failed (%d)\n",
					MBC_R_KINEMATICS_SIZE(mbc), rc);
				return -1;
			}
		}

		if (mbc->modes > 0) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_M_KINEMATICS(mbc),
				MBC_M_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_M_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%u) q, qP failed (%d)\n",
					MBC_M_KINEMATICS_SIZE(mbc), rc);
				return -1;
			}
		}
	}

	return 0;
}

/* put forces to peer
 *
 * if mbc->rigid, force and moment must be set in storage pointed to
 *	by macros MBC_F, MBC_M
 * if mbc->modes > 0, modal forces must be set in storage pointed to
 *	by macro MBC_P
 */
int
mbc_modal_put_forces(mbc_modal_t *mbc, int last)
{
	if (last) {
		if (mbc->mbc.data_and_next) {
			mbc->mbc.cmd = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;

		} else {
			mbc->mbc.cmd = ES_GOTO_NEXT_STEP;
		}

	} else {
		mbc->mbc.cmd = ES_REGULAR_DATA;
	}

	if (mbc_put_cmd((mbc_t *)mbc)) {
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		/* rigid */
		if (mbc->rigid) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_R_DYNAMICS(mbc),
				MBC_R_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_R_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%u) rigid failed (%d)\n",
					MBC_R_DYNAMICS_SIZE(mbc), rc);
				return -1;
			}
		}

		/* modal */
		if (mbc->modes > 0) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_M_DYNAMICS(mbc),
				MBC_M_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_M_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%u) modes failed (%d)\n",
					MBC_M_DYNAMICS_SIZE(mbc), rc);
				return -1;
			}
		}
	}

	return 0;
}

/* initialize modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least rigid body motion must be defined (mbc->rigid != 0),
 * or modes must be > 0
 *
 * if modes > 0, mallocs memory that needs to be freed calling
 * mbc_modal_destroy()
 */
int
mbc_modal_init(mbc_modal_t *mbc, unsigned modes)
{
	if (!mbc->rigid && modes == 0) {
		fprintf(stderr, "need at least 1 mode or rigid body data\n");
		return -1;
	}

	mbc->modes = modes;

	if (mbc->modes > 0) {
		mbc->m = (double *)malloc(MBC_M_SIZE(mbc));
	}

	return 0;
}

/* negotiate modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least rigid body motion must be defined (mbc->rigid != 0),
 * or modes must be > 0
 *
 * the socket must be initialized and connected
 * sends a negotiation request to the master
 */
int
mbc_modal_negotiate_request(mbc_modal_t *mbc)
{
	int rc;
	uint8_t *uint8_ptr;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint8_t) + sizeof(uint32_t)];

	if (!mbc->rigid && mbc->modes == 0) {
		fprintf(stderr, "need at least 1 mode or rigid body data\n");
		return -1;
	}

	if (!(mbc->mbc.sock_flags & MBC_SF_VALID)) {
		fprintf(stderr, "socket is not valid\n");
		return -1;
	}

	mbc->mbc.cmd = ES_NEGOTIATION;
	mbc_put_cmd((mbc_t *)mbc);

	uint8_ptr = (uint8_t *)&buf[0];
	uint8_ptr[0] = (uint8_t)MBC_MODAL;
	if (mbc->rigid) {
		uint8_ptr[0] |= MBC_REF_NODE;
	}

	uint32_ptr = (uint32_t *)&uint8_ptr[1];
	uint32_ptr[0] = mbc->modes;

	rc = send(mbc->mbc.sock, (const void *)buf, sizeof(buf),
		mbc->mbc.send_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "send negotiate request failed (%d)\n", rc);
		return -1;
	}

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	switch (mbc->mbc.cmd) {
	case ES_ABORT:
		fprintf(stdout, "got ABORT from peer\n");
		return -1;

	case ES_OK:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%u from peer\n", mbc->mbc.cmd);
		return -1;
	}

	return 0;
}

int
mbc_modal_negotiate_response(mbc_modal_t *mbc)
{
	int rc;
	uint8_t *uint8_ptr;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint8_t) + sizeof(uint8_t) + sizeof(uint32_t)];
	unsigned uModal;
	unsigned uRef;

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %u\n", mbc->mbc.cmd);
	}

	switch (mbc->mbc.cmd) {
	case ES_NEGOTIATION:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%u from peer\n", mbc->mbc.cmd);
		return -1;
	}

	rc = recv(mbc->mbc.sock, (void *)buf, sizeof(buf), mbc->mbc.recv_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "recv negotiate request failed\n");
		return -1;
	}

	rc = 0;

	uint8_ptr = (uint8_t *)&buf[0];
	uModal = uint8_ptr[0] & MBC_MODAL_NODAL_MASK;
	uRef = uint8_ptr[0] & MBC_REF_NODE;
	if (uModal != MBC_MODAL) {
		rc++;
	}

	if ((uRef && !mbc->rigid) || (!uRef && mbc->rigid)) {
		rc++;
	}

	uint32_ptr = (uint32_t *)&uint8_ptr[1];
	if (uint32_ptr[0] != mbc->modes) {
		rc++;
	}

	if (rc) {
		mbc->mbc.cmd = ES_ABORT;

	} else {
		mbc->mbc.cmd = ES_OK;
	}

	mbc_put_cmd((mbc_t *)mbc);

	return 0;
}

/* destroy modal data
 *
 * does NOT free the mbc structure
 */
int
mbc_modal_destroy(mbc_modal_t *mbc)
{
	if (mbc->m) {
		free(mbc->m);
		mbc->m = NULL;
	}

	return 0;
}

#endif /* USE_SOCKET */
