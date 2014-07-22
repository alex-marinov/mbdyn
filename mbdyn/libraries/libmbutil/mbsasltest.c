/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/*
 * The protocol is very simple.  It is client-initiated (of course),
 * by sending
 *
 * C: M
 * S: M<length><mechanism list>
 * C: S<length><mechanism chosen><length><additional data>
 * 
 * S: 	C<length><additional data>
 * C: 	C<length><additional data>
 *
 * S: O | F
 *
 * where M means "Methods", C means "Continuation" (with data), 
 * O means "OK" and F means "Fail"
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_SASL2) && defined(HAVE_THREADS) && (HAVE_SEMAPHORE_H)

#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <signal.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/poll.h>
#include <sys/socket.h>
#include "ac/getopt.h"

#include <sasl/sasl.h>
#include "mbsasl.h"

#include <pthread.h>
#include <semaphore.h>

static sem_t	sem;

static void *
server_f(void *arg)
{
	int			sock = ((int *)((void **)arg)[0])[0];
	struct mbdyn_sasl_t	*mbdyn_sasl = ((void **)arg)[1];

	/* server operations */
	mbdyn_sasl_init(mbdyn_sasl);

	if (mbdyn_sasl_auth(sock, NULL, mbdyn_sasl) == SASL_OK) {
		printf("[server] OK\n");
	} else {
		printf("[server] FAIL\n");
	}

	sem_post(&sem);
	
	return NULL;
}

int
main(int argc, char *argv[])
{
	int			sock, sockp[2];
	int			rc;
	void			*arg[2];
	pthread_t		th;
	struct mbdyn_sasl_t	mbdyn_sasl = MBDYN_SASL_INIT,
				mbdyn_sasl_client = MBDYN_SASL_INIT,
				mbdyn_sasl_server = MBDYN_SASL_INIT;

	while (1) {
		int opt = getopt(argc, argv, "s:" /* MBDYN_SASL_OPTIONS */ );

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 's':
			if (optarg[1] != '=' || mbdyn_sasl_parse_args(optarg[0], &optarg[2], &mbdyn_sasl)) {
				printf("UNKNOWN PARAMETER '%c'\n", opt);
				return 1;
			}
			break;

		default:
			printf("usage: %s [-s{ailmrsuw}=<value>]\n"
					"\ta=<authz>    client: authorization identity (optional)\n"
					"\tf=<flag>[=<value>]\n"
					"\ti=<remoteip> remote ip\n"
					"\tl=<localip>  local ip\n"
					"\tm=<method>   (list of) acceptable method(s)\n"
					"\tr=<realm>    client: user realm;\n"
					"\t             server: server realm\n"
#if 0
					"\ts={server|client}    use SASL to negotiate auth \n"
					"\t                     as server or client\n"
#endif
					"\tu=<user>     client: user identity\n"
					"\tw=<cred>     client: user credential\n", argv[0]);
			exit(EXIT_SUCCESS);
		}
	}

	/* validate server data */
	mbdyn_sasl_server.use_sasl = MBDYN_SASL_SERVER;
	mbdyn_sasl_server.sasl_mech = mbdyn_sasl.sasl_mech;
	mbdyn_sasl_server.sasl_realm = mbdyn_sasl.sasl_realm;

	if (mbdyn_sasl_validate(&mbdyn_sasl_server) != SASL_OK) {
		fprintf(stderr, "[server] SASL DATA DID NOT VALIDATE\n");
		return 1;
	}

	/* validate client data */
	mbdyn_sasl_client = mbdyn_sasl;
	mbdyn_sasl_client.use_sasl = MBDYN_SASL_CLIENT;

	if (mbdyn_sasl_validate(&mbdyn_sasl_client) != SASL_OK) {
		fprintf(stderr, "[client] SASL DATA DID NOT VALIDATE\n");
		return 1;
	}

	/* socketpair */
	rc = socketpair(PF_LOCAL, SOCK_STREAM, 0, sockp);
	if (rc != 0) {
		fprintf(stderr, "[client] socketpair() failed\n");
		exit(EXIT_FAILURE);
	}
	sock = sockp[0];

	/* server thread */
	sem_init(&sem, 0, 0);
	arg[0] = &sockp[1];
	arg[1] = &mbdyn_sasl_server;
	rc = pthread_create(&th, NULL, server_f, (void *)arg);
	if (rc != 0) {
		fprintf(stderr, "[client] pthread_create() failed\n");
		exit(EXIT_FAILURE);
	}
	pthread_detach(th);

	/* client operations */
	mbdyn_sasl_init(&mbdyn_sasl_client);

	if (mbdyn_sasl_auth(sock, NULL, &mbdyn_sasl_client) == SASL_OK) {
		printf("[client] OK\n");
	} else {
		printf("[client] FAIL\n");
	}

	/* wait for server */
	write(sock, "Q", 1);
	sem_wait(&sem);

	/* close the sasl session */
	mbdyn_sasl_fini();

	return 0;
}

#else /* ! defined(HAVE_SASL2) && defined(HAVE_THREADS) && (HAVE_SEMAPHORE_H) */

int
main(void)
{
	printf("need sasl2, pthreads and semaphores\n");
	exit(EXIT_FAILURE);
}

#endif /* ! defined(HAVE_SASL2) && defined(HAVE_THREADS) && (HAVE_SEMAPHORE_H) */


