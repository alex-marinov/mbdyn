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
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

#include "test_strext_socket_lib.h"

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

static void
usage(void)
{
	fprintf(stderr,
		"usage: testsocket [options]\n"
		"\t-a\t\tuse accelerations\n"
		"\t-c [random:]<c>\tnumber of iterations\n"
		"\t-f {fx,fy,fz,mx,my,mz} reference node force/moment\n"
		"\t-H <url>\tURL (local://path | inet://host:port)\n"
		"\t-l\t\tlabels\n"
		"\t-i <filename>\tinput file\n"
		"\t-n\t\tonly forces, no moments\n"
		"\t-N <nodes>\tnodes number\n"
		"\t-o <filename>\t output file\n"
		"\t-p {f0x,f0y,f0z,m0x,m0y,m0z,...}\tnodal forces (need -N first)\n"
		"\t-r\t\tuse reference node data\n"
		"\t-R {mat|theta|euler123}\torientation format\n"
		"\t-s <sleeptime>\tsleep time between tries\n"
		"\t-t <timeout>\thow long to wait for connection\n"
		"\t-v\t\tverbose\n"
		"\t-x\t\tdata_and_next\n");
	exit(EXIT_FAILURE);
}

static int sleeptime = 1;
static int iters = 1;
static int iters_random = 0;
static unsigned steps;

static int nomoments = 0;
static int refnode = 0;
static int nodes = 0;
static int labels = 0;
static int accelerations = 0;
static unsigned rot = MBC_ROT_MAT;

static char *path = NULL;
static char *host = NULL;
static unsigned short int port = -1;

static mbc_nodal_t	mbcx = { { 0 } };
static mbc_nodal_t	*mbc = &mbcx;

static double fx[6], *f0 = NULL;
static double *p0 = NULL;

static int inpfile = 0;
static int outfile = 0;
static FILE *inputfile  = NULL;
static FILE *outputfile = NULL;

void
test_init(int argc, char *argv[])
{
	while (1) {
		int opt = getopt(argc, argv, "ac:f:H:i:lnN:o:p:rR:s:t:vx");

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
			if (inpfile) {
				fprintf(stderr, "test_strext_socket: "
					"-i already provided\n");
				usage();
			}

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
#ifdef _WIN32
				fprintf(stderr, "test_strext_socket: "
					"local sockets are not supported in Windows\n");
				usage();
#else
				path = optarg + sizeof("local://") - 1;
				if (path[0] != '/') {
					usage();
				}
#endif /* _WIN32 */
			} else {
				usage();
			}
			break;

		case 'i' : {

			int size = 6;

			if (f0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-i error -f already provided\n");
				usage();
			}

			if (p0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-i error -p already provided\n");
				usage();
			}

			if (inputfile != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-i already provided\n");
				usage();
			}

			if (optarg == NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-i empty argument\n");
				usage();
			}

			inputfile = fopen(optarg, "r");
			if (inputfile == NULL) {
				fprintf(stderr,  "test_strext_socket: "
					"-i unable to open input file %s\n", optarg);
				usage();
			}

			if (nodes <= 0) {
				fprintf(stderr, "test_strext_socket: "
					"-o requires -N\n");
				usage();
			}

			if (nomoments) {
				size = 3;
		  	}

			f0 = fx;

			p0 = (double *)calloc(sizeof(double), size*nodes);
			if (p0 == NULL) {
				fprintf(stderr, "test_strext_socket: "
					"malloc for nodal force values failed\n");
				exit(EXIT_FAILURE);
			}


			inpfile = 1;
			} break;

		case 'l':
			labels = 1;
			break;

		case 'n':
			if (p0 != NULL) {
				fprintf(stderr, "-n must occur before -p\n");
				usage();
			}
			nomoments = 1;
			break;

		case 'N':
			if (p0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-N cannot follow -p\n");
				usage();
			}

			nodes = atoi(optarg);
			if (nodes <= 0) {
				fprintf(stderr, "test_strext_socket: "
					"invalid node number %s\n",
					optarg);
				usage();
			}
			break;

		case 'o':

			if (optarg == NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-o empty argument\n");
				usage();
			}

			outputfile = fopen(optarg, "w");
			if (outputfile == NULL) {
				fprintf(stderr, "unable to open output file %s\n", optarg);
				usage();
			}

			outfile = 1;
			break;

		case 'p': {
			char *s;
			int i, size = 6;

			if (p0 != NULL) {
				fprintf(stderr, "test_strext_socket: "
					"-p already provided\n");
				usage();
			}

			if (nodes <= 0) {
				fprintf(stderr, "test_strext_socket: "
					"-p requires -N\n");
				usage();
			}

			if (nomoments) {
				size = 3;
			}


			p0 = (double *)calloc(sizeof(double), size*nodes);
			if (p0 == NULL) {
				fprintf(stderr, "test_strext_socket: "
					"malloc for nodal force values failed\n");
				exit(EXIT_FAILURE);
			}

			s = optarg;
			for (i = 0; i < size*nodes; i++) {
				char *next;
				char fm[sizeof("fm")] = "fm";
				const char xyz[] = "xyz";

				if (nomoments) {
					fm[1] = 'f';
				}

				p0[i] = strtod(s, &next);
				if (next == s) {
					fprintf(stderr, "test_strext_socket: "
						"unable to parse %c%d%c\n",
						fm[(i/3)%2], i/size, xyz[i%3]);
					usage();
				}

				if (i < size*nodes - 1) {
					if (next[0] != ',') {
						fprintf(stderr, "test_strext_socket: "
							"unable to parse %c%d%c\n",
							fm[(i/3)%2], i/size, xyz[i%3]);
						usage();
					}

					s = &next[1];

				} else {
					if (next[0] != '\0') {
						fprintf(stderr, "test_strext_socket: "
							"extra cruft past %c%d%c\n",
							fm[(i/3)%2], i/size, xyz[i%3]);
						usage();
					}
				}
			}
			} break;

		case 'r':
			refnode = 1;
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

		case 't':
			if (strcasecmp(optarg, "forever") == 0) {
				mbc->mbc.timeout = -1;
			} else {
				mbc->mbc.timeout = atoi(optarg);
			}
			break;

		case 'v':
			mbc->mbc.verbose = 1;
			fprintf(stdout, "test_strext_socket: "
				"Verbose output option selected\n ");
			break;

		case 'x':
			mbc->mbc.data_and_next = 1;
			break;

		default:
			usage();
		}
	}

	if (nomoments) {
		rot = MBC_U_ROT_2_REF_NODE_ROT(rot);
	}

	if (mbc->mbc.verbose == 1){
        fprintf(stderr, "test_strext_socket: "
                    "Finished parsing input\n ");
	}

	if (path) {
#ifdef _WIN32
		fprintf(stderr, "test_strext_socket: "
				"Windows does not support local sockets, use inet sockets instead.\n ");
		exit(EXIT_FAILURE);
#else
		/* initialize UNIX socket (path) */
		if (mbc_unix_init((mbc_t *)mbc, path)) {
			exit(EXIT_FAILURE);
		}
#endif /* _WIN32 */

	} else if (host) {
		if (mbc->mbc.verbose == 1){
			fprintf(stdout, "test_strext_socket: "
				"Initialising inet socket (%s:%d)\n ", host, port);
		}

		/* initialize INET socket (host, port) */
		int retval = mbc_inet_init((mbc_t *)mbc, host, port);

		if (retval) {
		    fprintf(stderr, "test_strext_socket: "
				"mbc_inet_init call failed (return value %d)\n ", retval);
			exit(EXIT_FAILURE);
		}

	} else {
		usage();
	}

	/* initialize data structure:
	 */
	if (mbc_nodal_init(mbc, refnode, nodes, labels, rot, accelerations)) {
		fprintf(stderr, "test_strext_socket: "
			"mbc_nodal_init call failed\n");
		exit(EXIT_FAILURE);
	}

	/* "negotiate" configuration with MBDyn
	 * errors out if configurations are inconsistent */
	if (mbc_nodal_negotiate_request(mbc)) {
		exit(EXIT_FAILURE);
	}

	signal(SIGTERM, sh);
	signal(SIGINT, sh);
}

/*
 * specific to test_strext_socket
 */
void
test_run(void)
{
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
			/* receives motion when available
			 * errors out in case of problems
			 */
			if (mbc_nodal_get_motion(mbc)) {
				goto done;
			}

			if (outfile)  {
				fprintf(outputfile, "STEP %u ITERATION %d\n", steps, iter);
				fprintf(outputfile, "SIMULATION TIME %+16.8e\n", MBC_N_TIME(mbc)[0]);
			}

			if (refnode) {
				double *x = MBC_R_X(mbc);
				double *R;
				double *v = MBC_R_XP(mbc);
				double *w = MBC_R_OMEGA(mbc);
				if (outfile) {

					if (labels) {
						fprintf(outputfile, "REF (%u)\n", MBC_R_K_LABEL(mbc));
					} else {
						fprintf(outputfile, "REF \n");
					}

					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n", x[0], x[1], x[2]);

					switch (MBC_F_ROT(mbc)) {
					default:
						R =  MBC_R_R(mbc);
						fprintf(outputfile, "R %+16.8e %+16.8e %+16.8e", R[0], R[3], R[6]);
						fprintf(outputfile, "  %+16.8e %+16.8e %+16.8e", R[1], R[4], R[7]);
						fprintf(outputfile, "  %+16.8e %+16.8e %+16.8e\n", R[2], R[5], R[8]);
						break;
					case MBC_ROT_THETA:
						R =  MBC_R_THETA(mbc);
						fprintf(outputfile, "T %+16.8e %+16.8e %+16.8e\n",
								R[0], R[1], R[2]);
						break;

					case MBC_ROT_EULER_123:
						R =  MBC_R_EULER_123(mbc);
						fprintf(outputfile, "E %+16.8e %+16.8e %+16.8e\n",
								R[0], R[1], R[2]);
						break;
					}

					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n", v[0], v[1], v[2]);
					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n\n", w[0], w[1], w[2]);
				}
				else if (mbc->mbc.verbose) {

					if (labels) {
						fprintf(stdout, "reference node (%u):\n", MBC_R_K_LABEL(mbc));
					} else {
						fprintf(stdout, "reference node:\n");
					}
					fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n", x[0], x[1], x[2]);
					switch (MBC_F_ROT(mbc)) {
					default:
						R =  MBC_R_R(mbc);
						fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n", R[0], R[3], R[6]);
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n", R[1], R[4], R[7]);
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}}\n", R[2], R[5], R[8]);
						break;
					case MBC_ROT_THETA:
						R =  MBC_R_THETA(mbc);
						fprintf(stdout, " theta={%+16.8e,%+16.8e,%+16.8e\n}",
								R[0], R[1], R[2]);
						break;

					case MBC_ROT_EULER_123:
						R =  MBC_R_EULER_123(mbc);
						fprintf(stdout, "euler123={%+16.8e,%+16.8e,%+16.8e}\n",
								R[0], R[1], R[2]);
						break;
					}
					fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n", v[0], v[1], v[2]);
					fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n", w[0], w[1], w[2]);

				}
			}
			if (mbc->nodes > 0) {
				uint32_t *n_labels = MBC_N_K_LABELS(mbc);
				double *n_x = MBC_N_X(mbc);
				double *n_r;
				double *n_xp = MBC_N_XP(mbc);
				double *n_omega = MBC_N_OMEGA(mbc);
				unsigned n;

				if (outfile) {
					fprintf(outputfile, "POS %u\n", mbc->nodes);
					for (n = 0; n < mbc->nodes; n++) {
						if (labels) {
							fprintf(outputfile,"%d ", n_labels[n]);
						}
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							 n_x[3*n], n_x[3*n + 1], n_x[3*n + 2]);
					}
					if (nomoments == 0) {
						fprintf(outputfile, "ROT %u\n", mbc->nodes);
						for (n = 0; n < mbc->nodes; n++) {
							if (labels) {
								fprintf(outputfile, "%d ", n_labels[n]);
							}
							switch (MBC_F_ROT(mbc)) {
							default:
								n_r =  MBC_N_R(mbc);
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e"
										" %+16.8e %+16.8e %+16.8e"
										" %+16.8e %+16.8e %+16.8e\n",
									n_r[9*n], n_r[9*n + 3], n_r[9*n + 6],
									n_r[9*n + 1], n_r[9*n + 4], n_r[9*n + 7],
									n_r[9*n + 2], n_r[9*n + 5], n_r[9*n + 8]);
								break;

							case MBC_ROT_THETA:
								n_r =  MBC_N_THETA(mbc);
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									n_r[3*n], n_r[3*n + 1], n_r[3*n + 2]);
								break;

							case MBC_ROT_EULER_123:
								n_r =  MBC_N_EULER_123(mbc);
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									n_r[3*n], n_r[3*n + 1], n_r[3*n + 2]);
								break;
							}
						}
					}
					fprintf(outputfile, "VEL %u\n", mbc->nodes);
					for (n = 0; n < mbc->nodes; n++) {
						if (labels) {
							fprintf(outputfile, "%d ", n_labels[n]);
						}
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							 n_xp[3*n], n_xp[3*n + 1], n_xp[3*n + 2]);
					}
					if (nomoments == 0) {
						fprintf(outputfile, "W %u\n", mbc->nodes);
						for (n = 0; n < mbc->nodes; n++) {
							if (labels) {
								fprintf(outputfile, "%d ", n_labels[n]);
							}
							fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
								n_omega[3*n], n_omega[3*n + 1], n_omega[3*n + 2]);
						}
					}
				}
				else if (mbc->mbc.verbose) {
					for (n = 0; n < mbc->nodes; n++) {
						if (labels) {
							fprintf(stdout, "node #%d (%u):\n", n, n_labels[n]);
						} else {
							fprintf(stdout, "node #%d:\n", n);
						}
						fprintf(stdout, "    x=     %+16.8e %+16.8e %+16.8e\n",
							n_x[3*n], n_x[3*n + 1], n_x[3*n + 2]);
						if (nomoments == 0) {
							switch (MBC_F_ROT(mbc)) {
							default:
								n_r =  MBC_N_R(mbc);
								fprintf(stdout, "    R=     %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n",
									n_r[9*n], n_r[9*n + 3], n_r[9*n + 6],
									n_r[9*n + 1], n_r[9*n + 4], n_r[9*n + 7],
									n_r[9*n + 2], n_r[9*n + 5], n_r[9*n + 8]);
								break;

							case MBC_ROT_THETA:
								n_r =  MBC_N_THETA(mbc);
								fprintf(stdout, "    theta= %+16.8e %+16.8e %+16.8e\n",
									n_r[3*n], n_r[3*n + 1], n_r[3*n + 2]);
								break;

							case MBC_ROT_EULER_123:
								n_r =  MBC_N_EULER_123(mbc);
								fprintf(stdout, " euler123= %+16.8e %+16.8e %+16.8e\n",
									n_r[3*n], n_r[3*n + 1], n_r[3*n + 2]);
								break;
							}
						}
						fprintf(stdout, "    xp=    %+16.8e %+16.8e %+16.8e\n",
							n_xp[3*n], n_xp[3*n + 1], n_xp[3*n + 2]);
						if (nomoments == 0) {
							fprintf(stdout, "    omega= %+16.8e %+16.8e %+16.8e\n",
								n_omega[3*n], n_omega[3*n + 1], n_omega[3*n + 2]);
						}
					}
				}
			}

			if (sleeptime) {
				sleep(sleeptime);
			}

			/* set forces */
			if (inpfile && (iter == 0) && !feof(inputfile)) {
				unsigned i;
				unsigned n;
				int size = 6;
				if (fscanf(inputfile, "Step %u\n", &i) != 1) {
					fprintf(stderr, "Step: %u. Error while reading step"
						" number from input file\n", steps);
					exit(EXIT_FAILURE);
				}
				if (i != steps) {
					fprintf(stderr, "Error wrong step number from input file,"
						" is %u and shoul be %u\n", i, steps);
					exit(EXIT_FAILURE);
				}
				if (refnode) {
					if (fscanf(inputfile, "REF %lg %lg %lg %lg %lg %lg\n",
						&f0[0], &f0[1], &f0[2], &f0[3], &f0[4], &f0[5]) != 6) {
						fprintf(stderr, "Step: %u. Error while reading Reference Node"
							" forces from input file\n", steps);
						exit(EXIT_FAILURE);
					}
				}

				if (nomoments) {
					size = 3;
				}
				for (n = 0; n < mbc->nodes; n++) {
					if (nomoments == 0) {
						if (fscanf(inputfile, "%lg %lg %lg %lg %lg %lg\n",
							&p0[size*n], &p0[size*n +1], &p0[size*n + 2],
							&p0[size*n + 3], &p0[size*n +4], &p0[size*n + 5]) != 6) {
							fprintf(stderr, "Step: %u. Error while reading Force & Moments"
 								" for Node %u from input file\n", steps, n);
							exit(EXIT_FAILURE);
						}
					} else {
						if (fscanf(inputfile, "%lg %lg %lg\n",
							&p0[size*n], &p0[size*n + 1], &p0[size*n + 2]) != 3) {
							fprintf(stderr, "Step: %u. Error while reading Forces for Node %u"
								" from input file\n", steps, n);
							exit(EXIT_FAILURE);
						}
					}
				}
			}
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

			if (mbc->nodes > 0) {
				double *n_f = MBC_N_F(mbc);
				double *n_m = MBC_N_M(mbc);
				unsigned n;

				if (labels) {
					uint32_t *k_l = MBC_N_K_LABELS(mbc);
					uint32_t *d_l = MBC_N_D_LABELS(mbc);

					for (n = 0; n < mbc->nodes; n++) {
						d_l[n] = k_l[n];
					}
				}

				if (p0) {
					int size = 6;
					if (nomoments) {
						size = 3;
					}
					for (n = 0; n < mbc->nodes; n++) {
						n_f[3*n] = p0[size*n];
						n_f[3*n + 1] = p0[size*n + 1];
						n_f[3*n + 2] = p0[size*n + 2];

						if (nomoments == 0) {
							n_m[3*n] = p0[size*n + 3];
							n_m[3*n + 1] = p0[size*n + 4];
							n_m[3*n + 2] = p0[size*n + 5];
						}
					}

				} else {
					for (n = 0; n < 3*mbc->nodes; n++) {
						n_f[n] = (double)(n + 1);
						if (nomoments == 0) {
							n_m[n] = (double)(n + 1);
						}
					}
				}
			}

			/* sends forces
			 * second argument == 1 indicates convergence;
			 * otherwise MBDyn will send another solution
			 * and keep iterating */
			if (mbc_nodal_put_forces(mbc, (iter == niters - 1))) {
				goto done;
			}
		}
	}

done:;
	/* destroy data structure and close socket */
	mbc_nodal_destroy(mbc);

	if (p0) {
		free(p0);
	}
}

/*
 * specific to test_strext_socket_f
 */
void
tdata_(int32_t *REFNODE, int32_t *NODES, int32_t *ROT, int32_t *ITERS, int32_t *VERB,
	int32_t *RC_P)
{
	switch (MBC_F_ROT(mbc)) {
	case MBC_ROT_MAT:
		*ROT = 0;
		break;

	case MBC_ROT_THETA:
		*ROT = 1;
		break;

	case MBC_ROT_EULER_123:
		*ROT = 2;
		break;

	default:
		*RC_P = 1;
		return;
	}

	if (MBC_F_LABELS(mbc)) {
		*RC_P = 1;
		return;
	}

	if (MBC_F_ACCELS(mbc)) {
		*RC_P = 1;
		return;
	}

	*REFNODE = MBC_F_REF_NODE(mbc);
	*NODES = mbc->nodes;
	*VERB = mbc->mbc.verbose;

	*ITERS = iters;

	return;
}

void
tforce_(float *RF, float *RM, float *NF, float *NM)
{
	/* set forces */
	if (refnode) {
		if (f0 != NULL) {
			RF[0] = f0[0];
			RF[1] = f0[1];
			RF[2] = f0[2];

			RM[0] = f0[3];
			RM[1] = f0[4];
			RM[2] = f0[5];

		} else {
			RF[0] = 1.;
			RF[1] = 2.;
			RF[2] = 3.;

			RM[0] = 4.;
			RM[1] = 5.;
			RM[2] = 6.;
		}
	}

	if (mbc->nodes > 0) {
		unsigned n;

		if (p0) {
			for (n = 0; n < mbc->nodes; n++) {
				NF[3*n] = p0[6*n];
				NF[3*n + 1] = p0[6*n + 1];
				NF[3*n + 2] = p0[6*n + 2];

				NM[3*n] = p0[6*n + 3];
				NM[3*n + 1] = p0[6*n + 4];
				NM[3*n + 2] = p0[6*n + 5];
			}

		} else {
			for (n = 0; n < 3*mbc->nodes; n++) {
				NF[n] = (double)(n + 1);
				NM[n] = (double)(n + 1);
			}
		}
	}

	return;
}

void
tsend_(float *RF, float *RM, float *NF, float *NM,
        int32_t *CONVERGED_P, int32_t *RC_P)
{
	/* set forces */
	if (MBC_F_REF_NODE(mbc)) {
		double *f = MBC_R_F(mbc);
		double *m = MBC_R_M(mbc);

		f[0] = RF[0];
		f[1] = RF[1];
		f[2] = RF[2];

		m[0] = RM[0];
		m[1] = RM[1];
		m[2] = RM[2];
	}

	if (mbc->nodes > 0) {
		double *f = MBC_N_F(mbc);
		double *m = MBC_N_M(mbc);
		unsigned node;

		for (node = 0; node < mbc->nodes; node++) {
			f[3*node] = NF[3*node];
			f[3*node + 1] = NF[3*node + 1];
			f[3*node + 2] = NF[3*node + 2];

			m[3*node] = NM[3*node];
			m[3*node + 1] = NM[3*node + 1];
			m[3*node + 2] = NM[3*node + 2];
		}
	}

	/* NOTE: the flag indicates convergence when not 0 */
	if (mbc_nodal_put_forces(mbc, *CONVERGED_P)) {
		*RC_P = 1;
		return;
	}

	*RC_P = 0;
	return;
}

void
trecv_(float *RX, float *RR, float *RXP, float *ROMEGA,
	float *NX, float *NR, float *NXP, float *NOMEGA, int32_t *RC_P)
{
	if (mbc_nodal_get_motion(mbc)) {
		*RC_P = 1;
		return;
	}

	if (MBC_F_REF_NODE(mbc)) {
		double *x = MBC_R_X(mbc);
		double *r;
		double *v = MBC_R_XP(mbc);
		double *w = MBC_R_OMEGA(mbc);

		RX[0] = x[0];
		RX[1] = x[1];
		RX[2] = x[2];

		switch (MBC_F_ROT(mbc)) {
		case MBC_ROT_MAT:
			r = MBC_R_R(mbc);
			RR[0] = r[0];
			RR[1] = r[1];
			RR[2] = r[2];
			RR[3] = r[3];
			RR[4] = r[4];
			RR[5] = r[5];
			RR[6] = r[6];
			RR[7] = r[7];
			RR[8] = r[8];
			break;

		case MBC_ROT_THETA:
			r = MBC_R_THETA(mbc);
			RR[0] = r[0];
			RR[1] = r[1];
			RR[2] = r[2];
			break;

		case MBC_ROT_EULER_123:
			r = MBC_R_EULER_123(mbc);
			RR[0] = r[0];
			RR[1] = r[1];
			RR[2] = r[2];
			break;
		}

		RXP[0] = v[0];
		RXP[1] = v[1];
		RXP[2] = v[2];

		ROMEGA[0] = w[0];
		ROMEGA[1] = w[1];
		ROMEGA[2] = w[2];
	}

	if (mbc->nodes > 0) {
		double *x = MBC_N_X(mbc);
		double *r;
		double *v = MBC_N_XP(mbc);
		double *w = MBC_N_OMEGA(mbc);
		unsigned node;

		for (node = 0; node < mbc->nodes; node++) {
			NX[3*node] = x[3*node];
			NX[3*node + 1] = x[3*node + 1];
			NX[3*node + 2] = x[3*node + 2];

			switch (MBC_F_ROT(mbc)) {
			case MBC_ROT_MAT:
				r = MBC_N_R(mbc);
				NR[9*node + 0] = r[9*node + 0];
				NR[9*node + 1] = r[9*node + 1];
				NR[9*node + 2] = r[9*node + 2];
				NR[9*node + 3] = r[9*node + 3];
				NR[9*node + 4] = r[9*node + 4];
				NR[9*node + 5] = r[9*node + 5];
				NR[9*node + 6] = r[9*node + 6];
				NR[9*node + 7] = r[9*node + 7];
				NR[9*node + 8] = r[9*node + 8];
				break;

			case MBC_ROT_THETA:
				r = MBC_N_THETA(mbc);
				NR[3*node + 0] = r[3*node + 0];
				NR[3*node + 1] = r[3*node + 1];
				NR[3*node + 2] = r[3*node + 2];
				break;

			case MBC_ROT_EULER_123:
				r = MBC_N_EULER_123(mbc);
				NR[3*node + 0] = r[3*node + 0];
				NR[3*node + 1] = r[3*node + 1];
				NR[3*node + 2] = r[3*node + 2];
				break;
			}

			NXP[3*node] = v[3*node];
			NXP[3*node + 1] = v[3*node + 1];
			NXP[3*node + 2] = v[3*node + 2];

			NOMEGA[3*node] = w[3*node];
			NOMEGA[3*node + 1] = w[3*node + 1];
			NOMEGA[3*node + 2] = w[3*node + 2];
		}
	}

	*RC_P = 0;
	return;
}

