/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

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

static const char *
mbc_cmd2str(uint8_t cmd)
{
	switch (cmd) {
	case ES_REGULAR_DATA: return "REGULAR_DATA";
	case ES_GOTO_NEXT_STEP: return "GOTO_NEXT_STEP";
	case ES_ABORT: return "ABORT";
	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP: return "REGULAR_DATA_AND_GOTO_NEXT_STEP";
	case ES_NEGOTIATION: return "NEGOTIATION";
	case ES_OK: return "OK";
	default:
		break;
	}

	return "UNKNOWN";
}

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

	fprintf(stderr, "unknown cmd (%lu) from peer\n", (unsigned long)mbc->cmd);
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
		fprintf(stderr, "recv(cmd=%lu) failed\n", (unsigned long)mbc->cmd);
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
		fprintf(stdout, "cmd to peer: %lu (%s)\n",
			(unsigned long)mbc->cmd, mbc_cmd2str(mbc->cmd));
	}

	rc = send(mbc->sock, (const void *)&mbc->cmd, sizeof(mbc->cmd),
		mbc->send_flags);
	if (rc != sizeof(mbc->cmd)) {
		fprintf(stderr, "send(cmd=%lu) failed (%ld)\n", (unsigned long)mbc->cmd, (long)rc);
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
				case ECONNREFUSED:	/* inet */
				case ENOENT:		/* unix */
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
			fprintf(stderr, "unable to connect to peer (%ld: %s)\n",
				(long)save_errno, msg);
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
static int
mbc_destroy(mbc_t *mbc)
{
	/* TODO: send "abort"? */
	if (mbc->sock >= 0) {
		close(mbc->sock);
		mbc->sock = -1;
	}

	return 0;
}

/*
 * reference node stuff
 */

static int
mbc_rigid_init(mbc_rigid_t *mbc, unsigned refnode,
	unsigned labels, unsigned *rotp, unsigned accels)
{
	uint32_t offset = 0;
	int do_retry = 1;

	mbc->k_size = 0;
	mbc->r_k_label = -1;
	mbc->r_k_x = -1;
	mbc->r_k_theta = -1;
	mbc->r_k_r = -1;
	mbc->r_k_euler_123 = -1;
	mbc->r_k_xp = -1;
	mbc->r_k_omega = -1;
	mbc->r_k_xpp = -1;
	mbc->r_k_omegap = -1;
	mbc->d_size = 0;
	mbc->r_d_label = -1;
	mbc->r_d_f = -1;
	mbc->r_d_m = -1;
	
	if (!refnode) {
		return 0;
	}

	if (labels) {
		mbc->r_k_label = offset/sizeof(uint32_t);
		/* skip 2 uint32_t to ensure alignment of double */
		offset += 2*sizeof(uint32_t);
	}

	mbc->r_k_x = offset/sizeof(double);
	offset += 3*sizeof(double);

retry:;
	switch (MBC_U_REF_NODE_ROT_2_ROT(*rotp)) {
	case MBC_ROT_NONE:
		if (do_retry) {
			*rotp = MBC_U_ROT_2_REF_NODE_ROT(*rotp);
			do_retry = 0;
			goto retry;
		}
		fprintf(stderr, "rotation must be defined for reference node\n");
		return -1;

	case MBC_ROT_THETA:
		mbc->r_k_theta = offset/sizeof(double);
		offset += 3*sizeof(double);
		break;

	case MBC_ROT_MAT:
		mbc->r_k_r = offset/sizeof(double);
		offset += 9*sizeof(double);
		break;

	case MBC_ROT_EULER_123:
		mbc->r_k_euler_123 = offset/sizeof(double);
		offset += 3*sizeof(double);
		break;

	default:
		fprintf(stderr, "unknown rotation mode 0x%lx\n", (unsigned long)*rotp);
		return -1;
	}

	mbc->r_k_xp = offset/sizeof(double);
	offset += 3*sizeof(double);

	mbc->r_k_omega = offset/sizeof(double);
	offset += 3*sizeof(double);

	if (accels) {
		mbc->r_k_xpp = offset/sizeof(double);
		offset += 3*sizeof(double);

		mbc->r_k_omegap = offset/sizeof(double);
		offset += 3*sizeof(double);
	}

	mbc->k_size = offset;

	if (labels) {
		mbc->r_d_label = offset/sizeof(uint32_t);
		/* skip 2 uint32_t to ensure alignment of double */
		offset += 2*sizeof(uint32_t);
		mbc->d_size += 2*sizeof(uint32_t);
	}

	mbc->r_d_f = offset/sizeof(double);
	offset += 3*sizeof(double);
	mbc->d_size += 3*sizeof(double);

	mbc->r_d_m = offset/sizeof(double);
	offset += 3*sizeof(double);
	mbc->d_size += 3*sizeof(double);

	return 0;
}

/*
 * nodal stuff
 */

/* get nodal motion from peer
 *
 * if "rigid", access reference node motion
 * using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->nodes > 0, access nodal motion using macros MBC_N_*
 */
int
mbc_nodal_get_motion(mbc_nodal_t *mbc)
{
	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	if (mbc->mbc.cmd == ES_ABORT) {
		fprintf(stdout, "got ABORT from peer\n");
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		if (MBC_F_REF_NODE(mbc)) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_R_KINEMATICS(mbc),
				MBC_R_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_R_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%lu) reference node failed (%ld)\n",
					(unsigned long)MBC_R_KINEMATICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}

		if (mbc->nodes > 0) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_N_KINEMATICS(mbc),
				MBC_N_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc != MBC_N_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%lu) x, theta, xP, omega failed (%ld)\n",
					(unsigned long)MBC_N_KINEMATICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}
	}

	return 0;
}

/* put forces to peer
 *
 * if "rigid", force and moment must be set in storage pointed to
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
		/* reference node */
		if (MBC_F_REF_NODE(mbc)) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_R_DYNAMICS(mbc),
				MBC_R_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_R_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%lu) reference node failed (%ld)\n",
					(unsigned long)MBC_R_DYNAMICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}

		/* nodal */
		if (mbc->nodes > 0) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_N_DYNAMICS(mbc),
				MBC_N_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc != MBC_N_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%lu) nodes failed (%ld)\n",
					(unsigned long)MBC_N_DYNAMICS_SIZE(mbc), (long)rc);
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
 * at least reference node motion must be defined (refnode != 0),
 * or nodes must be > 0
 *
 * if nodes > 0, mallocs memory that needs to be freed calling
 * mbc_nodal_destroy()
 *
 * rot must be one of MBC_ROT_*
 *
 * if accelerations != 0 accelerations are also output
 */
int
mbc_nodal_init(mbc_nodal_t *mbc, unsigned refnode, unsigned nodes,
	unsigned labels, unsigned rot, unsigned accels)
{
	mbc->nodes = nodes;
	MBC_F_SET(mbc, MBC_NODAL);

	if (refnode) {
		MBC_F_SET_REF_NODE(mbc);
	}

	if (!MBC_F_REF_NODE(mbc) && mbc->nodes == 0) {
		fprintf(stderr, "need at least 1 node or reference node data\n");
		return -1;
	}

	switch (rot & MBC_ROT_MASK) {
	case MBC_ROT_NONE:
		break;

	case MBC_ROT_MAT:
		MBC_F_SET_ROT_MAT(mbc);
		break;

	case MBC_ROT_THETA:
		MBC_F_SET_ROT_THETA(mbc);
		break;

	case MBC_ROT_EULER_123:
		MBC_F_SET_ROT_EULER_123(mbc);
		break;

	default:
		fprintf(stderr, "unknown orientation parametrization 0x%lx in flags\n", (unsigned long)rot);
		return -1;
	}

	if (accels) {
		MBC_F_SET_ACCELS(mbc);
	}

	if (labels) {
		MBC_F_SET_LABELS(mbc);
	}

	if (mbc_rigid_init(&mbc->mbcr, refnode, labels, &rot, accels)) {
		return -1;
	}
	MBC_F_SET(mbc, (rot & MBC_REF_NODE_ROT_MASK));
	
	mbc->n_ptr = NULL;
	mbc->n_k_labels = NULL;
	mbc->n_k_x = NULL;
	mbc->n_k_theta = NULL;
	mbc->n_k_r = NULL;
	mbc->n_k_euler_123 = NULL;
	mbc->n_k_xp = NULL;
	mbc->n_k_omega = NULL;
	mbc->n_k_xpp = NULL;
	mbc->n_k_omegap = NULL;
	mbc->n_d_labels = NULL;
	mbc->n_d_f = NULL;
	mbc->n_d_m = NULL;

	if (mbc->nodes > 0) {
		char *ptr;

		mbc->k_size = (3 + 3); /* x, xp */

		switch (MBC_F_ROT(mbc)) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			mbc->k_size += 9 + 3;
			break;

		case MBC_ROT_THETA:
		case MBC_ROT_EULER_123:
			mbc->k_size += 3 + 3;
			break;
		}

		if (MBC_F_ACCELS(mbc)) {
			mbc->k_size += 3;
			if (MBC_F_ROT(mbc) != MBC_ROT_NONE) {
				mbc->k_size += 3;
			}
		}

		mbc->k_size *= mbc->nodes*sizeof(double);

		mbc->d_size = mbc->nodes*3*sizeof(double);
		if (MBC_F_ROT(mbc) != MBC_ROT_NONE) {
			mbc->d_size += mbc->nodes*3*sizeof(double);
		}

		if (MBC_F_LABELS(mbc)) {
			mbc->k_size += mbc->nodes*sizeof(uint32_t);
			mbc->d_size += mbc->nodes*sizeof(uint32_t);
		}

		mbc->n_ptr = (void *)malloc(MBC_N_SIZE(mbc));
		if (mbc->n_ptr == NULL) {
			fprintf(stderr, "nodal data malloc failed\n");
			return -1;
		}

		ptr = mbc->n_ptr;

		if (MBC_F_LABELS(mbc)) {
			mbc->n_k_labels = (uint32_t *)ptr;
			ptr += (nodes + nodes%2)*sizeof(uint32_t);
		}

		mbc->n_k_x = (double *)ptr;
		ptr += 3*sizeof(double)*nodes;

		switch (MBC_F_ROT(mbc)) {
		case MBC_ROT_NONE:
			break;

		case MBC_ROT_MAT:
			mbc->n_k_r = (double *)ptr;
			ptr += 9*sizeof(double)*nodes;
			break;

		case MBC_ROT_THETA:
			mbc->n_k_theta = (double *)ptr;
			ptr += 3*sizeof(double)*nodes;
			break;

		case MBC_ROT_EULER_123:
			mbc->n_k_euler_123 = (double *)ptr;
			ptr += 3*sizeof(double)*nodes;
			break;
		}

		mbc->n_k_xp = (double *)ptr;
		ptr += 3*sizeof(double)*nodes;

		if (MBC_F_ROT(mbc) != MBC_ROT_NONE) {
			mbc->n_k_omega = (double *)ptr;
			ptr += 3*sizeof(double)*nodes;
		}

		if (MBC_F_ACCELS(mbc)) {
			mbc->n_k_xpp = (double *)ptr;
			ptr += 3*sizeof(double)*nodes;

			if (MBC_F_ROT(mbc) != MBC_ROT_NONE) {
				mbc->n_k_omegap = (double *)ptr;
				ptr += 3*sizeof(double)*nodes;
			}
		}

		if (MBC_F_LABELS(mbc)) {
			mbc->n_d_labels = (uint32_t *)ptr;
			ptr += (nodes + nodes%2)*sizeof(uint32_t);
		}

		mbc->n_d_f = (double *)ptr;
		ptr += 3*sizeof(double)*nodes;

		if (MBC_F_ROT(mbc) != MBC_ROT_NONE) {
			mbc->n_d_m = (double *)ptr;
			ptr += 3*sizeof(double)*nodes;
		}
	}

	return 0;
}

/* negotiate nodal data
 *
 * mbc must be a pointer to a valid mbc_nodal_t structure
 *
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or nodes must be > 0
 *
 * the socket must be initialized and connected
 * sends a negotiation request to the master
 */
int
mbc_nodal_negotiate_request(mbc_nodal_t *mbc)
{
	int rc;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint32_t) + sizeof(uint32_t)];

	if (!MBC_F_REF_NODE(mbc) && mbc->nodes == 0) {
		fprintf(stderr, "need at least 1 node or reference node data\n");
		return -1;
	}

	if (!(mbc->mbc.sock_flags & MBC_SF_VALID)) {
		fprintf(stderr, "socket is not valid\n");
		return -1;
	}

	mbc->mbc.cmd = ES_NEGOTIATION;
	mbc_put_cmd((mbc_t *)mbc);

	uint32_ptr = (uint32_t *)&buf[0];
	uint32_ptr[0] = MBC_F(mbc);
	uint32_ptr[1] = mbc->nodes;

	rc = send(mbc->mbc.sock, (const void *)buf, sizeof(buf),
		mbc->mbc.send_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "send negotiate request failed (%ld)\n", (long)rc);
		return -1;
	}

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	switch (mbc->mbc.cmd) {
	case ES_ABORT:
		fprintf(stdout, "got ABORT from peer\n");
		return -1;

	case ES_OK:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%lu from peer\n", (unsigned long)mbc->mbc.cmd);
		return -1;
	}

	return 0;
}

int
mbc_nodal_negotiate_response(mbc_nodal_t *mbc)
{
	int rc;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint32_t) + sizeof(uint32_t)];
	unsigned uNodal;
	unsigned uRef;
	unsigned uLabels;
	unsigned uAccels;

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	switch (mbc->mbc.cmd) {
	case ES_NEGOTIATION:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%lu from peer\n", (unsigned long)mbc->mbc.cmd);
		return -1;
	}

	rc = recv(mbc->mbc.sock, (void *)buf, sizeof(buf), mbc->mbc.recv_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "recv negotiate request failed\n");
		return -1;
	}

	rc = 0;

	uint32_ptr = (uint32_t *)&buf[0];
	uNodal = (uint32_ptr[0] & MBC_MODAL_NODAL_MASK);
	uRef = (uint32_ptr[0] & MBC_REF_NODE);
	uLabels = (uint32_ptr[0] & MBC_LABELS);
	uAccels = (uint32_ptr[0] & MBC_ACCELS);

	if (uNodal != MBC_NODAL) {
		rc++;
	}

	if ((uRef && !MBC_F_REF_NODE(mbc)) || (!uRef && MBC_F_REF_NODE(mbc))) {
		rc++;
	}

	if (uLabels != MBC_F_LABELS(mbc)) {
		rc++;
	}

	if (uAccels != MBC_F_ACCELS(mbc)) {
		rc++;
	}

	if (uint32_ptr[1] != mbc->nodes) {
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
	if (mbc->n_ptr) {
		free(mbc->n_ptr);
		mbc->n_k_labels = NULL;
		mbc->n_k_x = NULL;
		mbc->n_k_theta = NULL;
		mbc->n_k_r = NULL;
		mbc->n_k_euler_123 = NULL;
		mbc->n_k_xp = NULL;
		mbc->n_k_omega = NULL;
		mbc->n_k_xpp = NULL;
		mbc->n_k_omegap = NULL;
		mbc->n_d_labels = NULL;
		mbc->n_d_f = NULL;
		mbc->n_d_m = NULL;
		mbc->n_ptr = NULL;
	}

	return mbc_destroy((mbc_t *)mbc);
}



/*
 * modal stuff
 */

/* get modal motion from peer
 *
 * if MBC_F_REF_NODE(mbc), access reference node motion
 * using macros MBC_X, MBC_R, MBC_V, MBC_W
 * if mbc->modes > 0, access modal motion using macros MBC_Q, MBC_QP
 */
int
mbc_modal_get_motion(mbc_modal_t *mbc)
{
	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	if (mbc->mbc.cmd == ES_ABORT) {
		fprintf(stdout, "got ABORT from peer\n");
		return -1;
	}

	if (mbc->mbc.cmd != ES_GOTO_NEXT_STEP) {
		if (MBC_F_REF_NODE(mbc)) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_R_KINEMATICS(mbc),
				MBC_R_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc == -1) {
				int save_errno = errno;
				const char *msg;

				msg = strerror(save_errno);
				fprintf(stderr, "recv(%lu) reference node failed (%ld: %s)\n",
					(unsigned long)MBC_R_KINEMATICS_SIZE(mbc),
					(long)save_errno, msg);
				return -1;

			} else if ((unsigned long)rc != MBC_R_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%lu) reference node failed (%ld)\n",
					(unsigned long)MBC_R_KINEMATICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}

		if (mbc->modes > 0) {
			ssize_t rc;

			rc = recv(mbc->mbc.sock, (void *)MBC_M_KINEMATICS(mbc),
				MBC_M_KINEMATICS_SIZE(mbc), mbc->mbc.recv_flags);
			if (rc == -1) {
				int save_errno = errno;
				const char *msg;

				msg = strerror(save_errno);
				fprintf(stderr, "recv(%lu) q, qP failed (%ld: %s)\n",
					(unsigned long)MBC_M_KINEMATICS_SIZE(mbc),
					(long)save_errno, msg);
				return -1;

			} else if ((unsigned long)rc != MBC_M_KINEMATICS_SIZE(mbc)) {
				fprintf(stderr, "recv(%lu) q, qP failed (%ld)\n",
					(unsigned long)MBC_M_KINEMATICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}
	}

	return 0;
}

/* put forces to peer
 *
 * if MBC_F_REF_NODE(mbc), force and moment must be set in storage pointed to
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
		/* reference node */
		if (MBC_F_REF_NODE(mbc)) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_R_DYNAMICS(mbc),
				MBC_R_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc == -1) {
				int save_errno = errno;
				const char *msg;

				msg = strerror(save_errno);
				fprintf(stderr, "send(%lu) reference node failed (%ld: %s)\n",
					(unsigned long)MBC_R_DYNAMICS_SIZE(mbc),
					(long)save_errno, msg);
				return -1;

			} else if ((unsigned long)rc != MBC_R_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%lu) reference node failed (%ld)\n",
					(unsigned long)MBC_R_DYNAMICS_SIZE(mbc), (long)rc);
				return -1;
			}
		}

		/* modal */
		if (mbc->modes > 0) {
			ssize_t	rc;

			rc = send(mbc->mbc.sock, (const void *)MBC_M_DYNAMICS(mbc),
				MBC_M_DYNAMICS_SIZE(mbc), mbc->mbc.send_flags);
			if (rc == -1) {
				int save_errno = errno;
				const char *msg;

				msg = strerror(save_errno);
				fprintf(stderr, "send(%lu) modes failed (%ld: %s)\n",
					(unsigned long)MBC_M_DYNAMICS_SIZE(mbc),
					(long)save_errno, msg);
				return -1;

			} else if ((unsigned long)rc != MBC_M_DYNAMICS_SIZE(mbc)) {
				fprintf(stderr, "send(%lu) modes failed (%ld)\n",
					(unsigned long)MBC_M_DYNAMICS_SIZE(mbc), (long)rc);
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
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or modes must be > 0
 *
 * if modes > 0, mallocs memory that needs to be freed calling
 * mbc_modal_destroy()
 */
int
mbc_modal_init(mbc_modal_t *mbc, int refnode, unsigned modes)
{
	unsigned rot = MBC_ROT_MAT;

	MBC_F_SET(mbc, MBC_MODAL);
	mbc->modes = modes;

	if (refnode) {
		MBC_F_SET_REF_NODE(mbc);
	}

	if (!MBC_F_REF_NODE(mbc) && modes == 0) {
		fprintf(stderr, "need at least 1 mode or reference node data\n");
		return -1;
	}

	/* FIXME: rotation configurable? */
	if (mbc_rigid_init(&mbc->mbcr, refnode, 0, &rot, 0)) {
		return -1;
	}
	MBC_F_SET(mbc, (rot & MBC_REF_NODE_ROT_MASK));

	if (mbc->modes > 0) {
		mbc->m = (double *)malloc(MBC_M_SIZE(mbc));
	}

	return 0;
}

/* negotiate modal data
 *
 * mbc must be a pointer to a valid mbc_modal_t structure
 *
 * at least reference node motion must be defined (MBC_F_REF_NODE(mbc)),
 * or modes must be > 0
 *
 * the socket must be initialized and connected
 * sends a negotiation request to the master
 */
int
mbc_modal_negotiate_request(mbc_modal_t *mbc)
{
	int rc;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint32_t) + sizeof(uint32_t)];

	if (!MBC_F_REF_NODE(mbc) && mbc->modes == 0) {
		fprintf(stderr, "need at least 1 mode or reference node data\n");
		return -1;
	}

	if (!(mbc->mbc.sock_flags & MBC_SF_VALID)) {
		fprintf(stderr, "socket is not valid\n");
		return -1;
	}

	mbc->mbc.cmd = ES_NEGOTIATION;
	mbc_put_cmd((mbc_t *)mbc);

	uint32_ptr = (uint32_t *)&buf[0];
	uint32_ptr[0] = (uint32_t)MBC_F(mbc);
	uint32_ptr[1] = mbc->modes;

	rc = send(mbc->mbc.sock, (const void *)buf, sizeof(buf),
		mbc->mbc.send_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "send negotiate request failed (%ld)\n", (long)rc);
		return -1;
	}

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	switch (mbc->mbc.cmd) {
	case ES_ABORT:
		fprintf(stdout, "got ABORT from peer\n");
		return -1;

	case ES_OK:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%lu from peer\n", (unsigned long)mbc->mbc.cmd);
		return -1;
	}

	return 0;
}

int
mbc_modal_negotiate_response(mbc_modal_t *mbc)
{
	int rc;
	uint32_t *uint32_ptr;
	char buf[sizeof(uint32_t) + sizeof(uint32_t)];
	unsigned uModal;
	unsigned uRef;

	if (mbc_get_cmd((mbc_t*)mbc)) {
		return -1;
	}

	if (mbc->mbc.verbose) {
		fprintf(stdout, "cmd from peer: %lu (%s)\n",
			(unsigned long)mbc->mbc.cmd, mbc_cmd2str(mbc->mbc.cmd));
	}

	switch (mbc->mbc.cmd) {
	case ES_NEGOTIATION:
		break;

	default:
		fprintf(stdout, "unexpected cmd=%lu from peer\n", (unsigned long)mbc->mbc.cmd);
		return -1;
	}

	rc = recv(mbc->mbc.sock, (void *)buf, sizeof(buf), mbc->mbc.recv_flags);
	if (rc != sizeof(buf)) {
		fprintf(stderr, "recv negotiate request failed\n");
		return -1;
	}

	rc = 0;

	uint32_ptr = (uint32_t *)&buf[0];
	uModal = uint32_ptr[0] & MBC_MODAL_NODAL_MASK;
	uRef = uint32_ptr[0] & MBC_REF_NODE;

	if (uModal != MBC_MODAL) {
		rc++;
	}

	if ((uRef && !MBC_F_REF_NODE(mbc)) || (!uRef && MBC_F_REF_NODE(mbc))) {
		rc++;
	}

	if (uint32_ptr[1] != mbc->modes) {
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

	return mbc_destroy((mbc_t *)mbc);
}

#endif /* USE_SOCKET */
