/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/*
 * NOTE: this is all we use from MBDyn.  It was intentionally designed
 * to be configuration-independent.
 */
#include "mbc.h"

/* test tool specific stuff */
static volatile sig_atomic_t keep_going = 1;

static void
sh(int signum)
{
	keep_going = 0;
	signal(signum, SIG_DFL);
}

void
usage(void)
{
	fprintf(stderr,
		"usage: testsocket [options]\n"
		"\t-c [random:]<c>\t\tnumber of iterations\n"
		"\t-H <url>\tURL (local://path | inet://host:port)\n"
		"\t-M <modes>\tmodes number\n"
		"\t-r\t\tuse rigid body data\n"
		"\t-s <sleeptime>\t\tsleep time between tries\n"
		"\t-v\t\tverbose\n"
		"\t-x\t\tdata_and_next\n");
	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	int sleeptime = 1;
	int iters = 1;
	int iters_random = 0;
	unsigned steps;

	char *path = NULL;
	char *host = NULL;
	unsigned short int port = -1;

	mbc_modal_t	mbcx = { { 0 } };
	mbc_modal_t	*mbc = &mbcx;

	while (1) {
		int opt = getopt(argc, argv, "c:H:M:rs:vx");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'c':
			if (strncasecmp(optarg, "random:", sizeof("random:") -1) == 0) {
				iters_random = 1;
				iters = atoi(&optarg[sizeof("random:") -1]);

			} else {
				iters = atoi(optarg);
				printf("iterations: %d\n", iters);
			}
			if (iters < 1) {
				fprintf(stderr, "testedge: "
					"invalid sleep time %s\n",
					optarg);
				usage();
			}
			break;

		case 'H':
			if (strncasecmp(optarg, "inet://", sizeof("inet://") - 1) == 0) {
				char *ptr, *next;
				long l;

				host = optarg + sizeof("inet://") - 1;
				ptr = strchr(host, ':');
				if (ptr == NULL) {
					usage();
				}
				*ptr = '\0';
				ptr++;
				l = strtol(ptr, &next, 10);
				if (next == ptr || next[0] != '\0') {
					usage();
				}
				if (l <= 0) {
					usage();
				}
				port = (unsigned short)l;

			} else if (strncasecmp(optarg, "local://", sizeof("local://") - 1) == 0) {
				path = optarg + sizeof("local://") - 1;
				if (path[0] != '/') {
					usage();
				}

			} else {
				usage();
			}
				
			break;

		case 'M':
			mbc->modes = atoi(optarg);
			if (mbc->modes <= 0) {
				fprintf(stderr, "testedge: "
					"invalid mode number %s\n",
					optarg);
				usage();
			}
			break;

		case 'r':
			mbc->rigid = 1;
			break;

		case 's':
			sleeptime = atoi(optarg);
			if (sleeptime < 0) {
				fprintf(stderr, "testedge: "
					"invalid iters %s\n",
					optarg);
				usage();
			}
			break;

		case 'v':
			mbc->mbc.verbose = 1;
			break;

		case 'x':
			mbc->mbc.data_and_next = 1;
			break;

		default:
			usage();
		}
	}

	if (path) {
		if (mbc_unix_init((mbc_t *)mbc, path)) {
			exit(EXIT_FAILURE);
		}

	} else if (host) {
		if (mbc_inet_init((mbc_t *)mbc, host, port)) {
			exit(EXIT_FAILURE);
		}

	} else {
		usage();
	}

	if (mbc_modal_init(mbc, mbc->modes)) {
		exit(EXIT_FAILURE);
	}

	signal(SIGTERM, sh);
	signal(SIGINT, sh);

	for (steps = 0; keep_going > 0; steps++) {
		int iter;
		int niters;

		if (iters_random) {
			niters = rand() % iters + 1;
			printf("    iterations within this iter: %d\n", niters);

		} else {
			niters = iters;
		}

		for (iter = 0; iter < niters; iter++) {
			if (mbc_modal_get_motion(mbc)) {
				goto done;
			}

			if (mbc->rigid && mbc->mbc.verbose) {
				double *x = MBC_X(mbc);
				double *R = MBC_R(mbc);
				double *v = MBC_V(mbc);
				double *w = MBC_W(mbc);

				fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n", x[0], x[1], x[2]);
				fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n", R[0], R[3], R[6]);
				fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n", R[1], R[4], R[7]);
				fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}};\n", R[2], R[5], R[8]);
				fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n", v[0], v[1], v[2]);
				fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n", w[0], w[1], w[2]);
			}

			if (mbc->modes > 0 && mbc->mbc.verbose) {
				double *q = MBC_Q(mbc);
				double *qp = MBC_QP(mbc);
				int m;

				for (m = 0; m < mbc->modes; m++) {
					fprintf(stdout, "mode #%d: %+16.8e %+16.8e\n", m, q[m], qp[m]);
				}
			}

			if (sleeptime) {
				sleep(sleeptime);
			}

			/* set forces */
			if (mbc->rigid) {
				double *f = MBC_F(mbc);
				double *m = MBC_M(mbc);

				f[0] = 1.;
				f[1] = 2.;
				f[2] = 3.;

				m[0] = 4.;
				m[1] = 5.;
				m[2] = 6.;
			}

			if (mbc->modes > 0) {
				double *p = MBC_P(mbc);
				int m;

				for (m = 0; m < mbc->modes; m++) {
					p[m] = (double)(m + 1);
				}
			}

			if (mbc_modal_put_forces(mbc, (iter == niters - 1))) {
				goto done;
			}
		}
	}

done:;
	mbc_modal_destroy(mbc);
	mbc_destroy(&mbc->mbc);

	return 0;
}
