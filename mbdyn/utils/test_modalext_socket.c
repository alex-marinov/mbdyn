/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
#include <stdint.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

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
		"\t-c [random:]<c>\tnumber of iterations\n"
		"\t-f {fx,fy,fz,mx,my,mz} reference node force/moment\n"
		"\t-H <url>\tURL (local://path | inet://host:port)\n"
		"\t-M <modes>\tmodes number\n"
		"\t-p {p1,...,pM}\tmodal forces (need -M first)\n"
		"\t-r\t\tuse reference node data\n"
		"\t-s <sleeptime>\tsleep time between tries\n"
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

	int refnode = 0;
	unsigned modes = 0;

	mbc_modal_t	mbcx = { { 0 } };
	mbc_modal_t	*mbc = &mbcx;

	double fx[6], *f0 = NULL;
	double *p0 = NULL;

	while (1) {
		int opt = getopt(argc, argv, "c:f:H:M:p:rs:vx");

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
				fprintf(stderr, "test_modalext_socket: "
					"invalid iterations %s\n",
					optarg);
				usage();
			}
			break;

		case 'f': {
			char *s;
			int i;

			if (f0 != NULL) {
				fprintf(stderr, "test_modalext_socket: "
					"-f already provided\n");
				usage();
			}

			f0 = fx;

			s = optarg;
			for (i = 0; i < 6; i++) {
				char *next;

				f0[i] = strtod(s, &next);
				if (next == s) {
					fprintf(stderr, "test_modalext_socket: "
						"unable to parse f[%d]\n", i);
					usage();
				}

				if (i < 5) {
					if (next[0] != ',') {
						fprintf(stderr, "test_modalext_socket: "
							"unable to parse past f[%d]\n", i);
						usage();
					}

					s = &next[1];

				} else {
					if (next[0] != '\0') {
						fprintf(stderr, "test_modalext_socket: "
							"extra cruft past f[%d]\n", i);
						usage();
					}
				}
			}
			} break;

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

		case 'M': {
			int imodes;

			if (p0 != NULL) {
				fprintf(stderr, "test_modalext_socket: "
					"-M cannot follow -p\n");
				usage();
			}

			imodes = atoi(optarg);
			if (imodes <= 0) {
				fprintf(stderr, "test_modalext_socket: "
					"invalid mode number %s\n",
					optarg);
				usage();
			}
			modes = (unsigned)imodes;
			} break;

		case 'p': {
			char *s;
			int i;

			if (p0 != NULL) {
				fprintf(stderr, "test_modalext_socket: "
					"-p already provided\n");
				usage();
			}

			if (modes == 0) {
				fprintf(stderr, "test_modalext_socket: "
					"-p requires -M\n");
				usage();
			}

			p0 = (double *)calloc(sizeof(double), modes);
			if (p0 == NULL) {
				fprintf(stderr, "test_modalext_socket: "
					"malloc for modal force values failed\n");
				exit(EXIT_FAILURE);
			}

			s = optarg;
			for (i = 0; i < modes; i++) {
				char *next;

				p0[i] = strtod(s, &next);
				if (next == s) {
					fprintf(stderr, "test_modalext_socket: "
						"unable to parse p[%d]\n", i);
					usage();
				}

				if (i < modes - 1) {
					if (next[0] != ',') {
						fprintf(stderr, "test_modalext_socket: "
							"unable to parse past p[%d]\n", i);
						usage();
					}

					s = &next[1];

				} else {
					if (next[0] != '\0') {
						fprintf(stderr, "test_modalext_socket: "
							"extra cruft past p[%d]\n", i);
						usage();
					}
				}
			}
			} break;

		case 'r':
			refnode = 1;
			break;

		case 's':
			sleeptime = atoi(optarg);
			if (sleeptime < 0) {
				fprintf(stderr, "test_modalext_socket: "
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

	if (mbc_modal_init(mbc, refnode, modes)) {
		exit(EXIT_FAILURE);
	}

	if (mbc_modal_negotiate_request(mbc)) {
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

			if (refnode && mbc->mbc.verbose) {
				double *x = MBC_R_X(mbc);
				double *r;
				double *v = MBC_R_XP(mbc);
				double *w = MBC_R_OMEGA(mbc);

				fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n", x[0], x[1], x[2]);
				switch (MBC_F_ROT(mbc)) {
				case MBC_ROT_THETA:
					r = MBC_R_THETA(mbc);
					fprintf(stdout, "t={%+16.8e,%+16.8e,%+16.8e};\n", r[0], r[1], r[2]);
					break;

				case MBC_ROT_MAT:
					r = MBC_R_R(mbc);
					fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n", r[0], r[3], r[6]);
					fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n", r[1], r[4], r[7]);
					fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}};\n", r[2], r[5], r[8]);
					break;

				case MBC_ROT_EULER_123:
					r = MBC_R_EULER_123(mbc);
					fprintf(stdout, "e={%+16.8e,%+16.8e,%+16.8e};\n", r[0], r[1], r[2]);
					break;
				}
				fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n", v[0], v[1], v[2]);
				fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n", w[0], w[1], w[2]);
			}

			if (modes > 0 && mbc->mbc.verbose) {
				double *q = MBC_Q(mbc);
				double *qp = MBC_QP(mbc);
				int m;

				for (m = 0; m < modes; m++) {
					fprintf(stdout, "mode #%d: %+16.8e %+16.8e\n", m, q[m], qp[m]);
				}
			}

			if (sleeptime) {
				sleep(sleeptime);
			}

			/* set forces */
			if (refnode) {
				double *f = MBC_R_F(mbc);
				double *m = MBC_R_M(mbc);

				if (f0 != NULL) {
					f[0] = f0[0];
					f[1] = f0[1];
					f[2] = f0[2];

					m[0] = f0[3];
					m[1] = f0[4];
					m[2] = f0[5];

				} else {
					f[0] = 1.;
					f[1] = 2.;
					f[2] = 3.;

					m[0] = 4.;
					m[1] = 5.;
					m[2] = 6.;
				}
			}

			if (modes > 0) {
				double *p = MBC_P(mbc);
				int m;

				if (p0) {
					for (m = 0; m < modes; m++) {
						p[m] = p0[m];
					}

				} else {
					for (m = 0; m < modes; m++) {
						p[m] = (double)(m + 1);
					}
				}
			}

			if (mbc_modal_put_forces(mbc, (iter == niters - 1))) {
				goto done;
			}
		}
	}

done:;
	mbc_modal_destroy(mbc);

	if (p0) {
		free(p0);
	}

	return 0;
}
