/* $Header$ */
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
		"\t-a\t\tuse accelerations\n"
		"\t-c [random:]<c>\tnumber of iterations\n"
		"\t-f {fx,fy,fz,mx,my,mz} rigid body force/moment\n"
		"\t-H <url>\tURL (local://path | inet://host:port)\n"
		"\t-N <nodes>\tnodes number\n"
		"\t-p {f0x,f0y,f0z,m0x,m0y,m0z,...}\tnodal forces (need -N first)\n"
		"\t-r\t\tuse rigid body data\n"
		"\t-R {mat|theta|euler123}\torientation format\n"
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

	int accelerations = 0;
	unsigned rot = MBC_ROT_MAT;

	char *path = NULL;
	char *host = NULL;
	unsigned short int port = -1;

	mbc_nodal_t	mbcx = { { 0 } };
	mbc_nodal_t	*mbc = &mbcx;

	double fx[6], *f0 = NULL;
	double *p0 = NULL;

	while (1) {
		int opt = getopt(argc, argv, "ac:f:H:N:p:rR:s:vx");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'a':
			accelerations = 1;
			break;

		case 'c':
			if (strncasecmp(optarg, "random:", sizeof("random:") -1) == 0) {
				iters_random = 1;
				iters = atoi(&optarg[sizeof("random:") -1]);

			} else {
				iters = atoi(optarg);
				printf("iterations: %d\n", iters);
			}
			if (iters < 1) {
				fprintf(stderr, "test_strext_socket: "
					"invalid iterations %s\n",
					optarg);
				usage();
			}
			break;

		case 'f': {
			char *s;
			int i;

			if (f0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-f already provided\n");
				usage();
			}

			f0 = fx;

			s = optarg;
			for (i = 0; i < 6; i++) {
				char *next;
				const char fm[] = "fm";
				const char xyz[] = "xyz";

				f0[i] = strtod(s, &next);
				if (next == s) {
					fprintf(stderr, "test_strext_socket: "
						"unable to parse %c%c\n",
						fm[i/3], xyz[i%3]);
					usage();
				}

				if (i < 5) {
					if (next[0] != ',') {
						fprintf(stderr, "test_strext_socket: "
							"unable to parse %c%c\n",
							fm[i/3], xyz[i%3]);
						usage();
					}

					s = &next[1];

				} else {
					if (next[0] != '\0') {
						fprintf(stderr, "test_strext_socket: "
							"extra cruft past %c%c\n",
							fm[i/3], xyz[i%3]);
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

		case 'N':
			if (p0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-N cannot follow -p\n");
				usage();
			}

			mbc->nodes = atoi(optarg);
			if (mbc->nodes <= 0) {
				fprintf(stderr, "test_strext_socket: "
					"invalid node number %s\n",
					optarg);
				usage();
			}
			break;

		case 'p': {
			char *s;
			int i;

			if (p0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-p already provided\n");
				usage();
			}

			if (mbc->nodes <= 0) {
				fprintf(stderr, "test_strext_socket: "
					"-p requires -N\n");
				usage();
			}

			p0 = (double *)calloc(sizeof(double), 6*mbc->nodes);
			if (p0 == NULL) {
				fprintf(stderr, "test_strext_socket: "
					"malloc for nodal force values failed\n");
				exit(EXIT_FAILURE);
			}

			s = optarg;
			for (i = 0; i < 6*mbc->nodes; i++) {
				char *next;
				const char fm[] = "fm";
				const char xyz[] = "xyz";

				p0[i] = strtod(s, &next);
				if (next == s) {
					fprintf(stderr, "test_strext_socket: "
						"unable to parse %c%d%c\n",
						fm[(i/3)%2], i/6, xyz[i%3]);
					usage();
				}

				if (i < mbc->nodes - 1) {
					if (next[0] != ',') {
						fprintf(stderr, "test_strext_socket: "
							"unable to parse %c%d%c\n",
							fm[(i/3)%2], i/6, xyz[i%3]);
						usage();
					}

					s = &next[1];

				} else {
					if (next[0] != '\0') {
						fprintf(stderr, "test_strext_socket: "
							"extra cruft past %c%d%c\n",
							fm[(i/3)%2], i/6, xyz[i%3]);
						usage();
					}
				}
			}
			} break;

		case 'r':
			mbc->rigid = 1;
			break;

		case 'R':
			if (strcasecmp(optarg, "mat") == 0) {
				rot = MBC_ROT_MAT;

			} else if (strcasecmp(optarg, "theta") == 0) {
				rot = MBC_ROT_THETA;

			} else if (strcasecmp(optarg, "euler123") == 0) {
				rot = MBC_ROT_EULER_123;

			} else {
				fprintf(stderr, "test_strext_socket: "
					"unknown orientation format \"%s\"\n",
					optarg);
				usage();
			}
			break;

		case 's':
			sleeptime = atoi(optarg);
			if (sleeptime < 0) {
				fprintf(stderr, "test_strext_socket: "
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

	if (accelerations) {
		rot |= MBC_ACCELS;
	}
	if (mbc_nodal_init(mbc, mbc->nodes, rot)) {
		exit(EXIT_FAILURE);
	}

	if (mbc_nodal_negotiate_request(mbc)) {
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
			if (mbc_nodal_get_motion(mbc)) {
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

			if (mbc->nodes > 0 && mbc->mbc.verbose) {
				double *n_x = MBC_N_X(mbc);
				double *n_theta = MBC_N_THETA(mbc);
				double *n_xp = MBC_N_XP(mbc);
				double *n_omega = MBC_N_OMEGA(mbc);
				int n;

				for (n = 0; n < mbc->nodes; n++) {
					fprintf(stdout, "node #%d:\n"
							"    x=     %+16.8e %+16.8e %+16.8e\n"
							"    theta= %+16.8e %+16.8e %+16.8e\n"
							"    xp=    %+16.8e %+16.8e %+16.8e\n"
							"    omega= %+16.8e %+16.8e %+16.8e\n",
							n,
							n_x[3*n], n_x[3*n + 1], n_x[3*n + 2],
							n_theta[3*n], n_theta[3*n + 1], n_theta[3*n + 2],
							n_xp[3*n], n_xp[3*n + 1], n_xp[3*n + 2],
							n_omega[3*n], n_omega[3*n + 1], n_omega[3*n + 2]);
				}
			}

			if (sleeptime) {
				sleep(sleeptime);
			}

			/* set forces */
			if (mbc->rigid) {
				double *f = MBC_F(mbc);
				double *m = MBC_M(mbc);

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

			if (mbc->nodes > 0) {
				double *n_f = MBC_N_F(mbc);
				double *n_m = MBC_N_M(mbc);
				int n;

				if (p0) {
					for (n = 0; n < mbc->nodes; n++) {
						n_f[3*n] = p0[6*n];
						n_f[3*n + 1] = p0[6*n + 1];
						n_f[3*n + 2] = p0[6*n + 2];

						n_m[3*n] = p0[6*n + 3];
						n_m[3*n + 1] = p0[6*n + 4];
						n_m[3*n + 2] = p0[6*n + 5];
					}

				} else {
					for (n = 0; n < 3*mbc->nodes; n++) {
						n_f[n] = (double)(n + 1);
						n_m[n] = (double)(n + 1);
					}
				}
			}

			if (mbc_nodal_put_forces(mbc, (iter == niters - 1))) {
				goto done;
			}
		}
	}

done:;
	mbc_nodal_destroy(mbc);
	mbc_destroy(&mbc->mbc);

	if (p0) {
		free(p0);
	}

	return 0;
}
