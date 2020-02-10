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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <signal.h>
#include <unistd.h>

/*
 * NOTE: this is all we use from MBDyn.  It was intentionally designed
 * to be configuration-independent.
 */
#include "mbcxxshared.h"

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
		"\t-H <region name> shared memory region name\n"
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
static MBCType refnoderot = MBC_ROT_NONE;
static int nodes = 0;
static int labels = 0;
static int accelerations = 0;
static MBCType rot = MBC_ROT_MAT;

static MBCSharedMemNodal *mbc;

static char *path;

static double fx[6], *f0 = NULL;
static double *p0 = NULL;

static int timeout;
static int verbose;
static int data_and_next;

static int inpfile = 0;
static int outfile = 0;
static FILE *inputfile  = NULL;
static FILE *outputfile = NULL;

/*
 * specific to test_strext_sharedmemory
 */
void test_run(void)
{
    if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: in test_run"  << std::endl; }

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
            if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: about to call mbc->GetMotion"  << std::endl; }
			if (mbc->GetMotion()) {
                if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: mbc->GetMotion was not zero"  << std::endl; }
				goto done;
			}

			if (outfile)  {
				fprintf(outputfile, "STEP %u ITERATION %d\n", steps, iter);
			}

			if (refnode) {
				if (outfile) {
					if (labels) {
						fprintf(outputfile, "REF (%u)\n",
							mbc->GetRefNodeKinematicsLabel());
					} else {
						fprintf(outputfile, "REF \n");
					}

					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
						mbc->GetRefNodeX(1), mbc->GetRefNodeX(2), mbc->GetRefNodeX(3));

					switch (mbc->GetRefNodeRot()) {
					default:
						fprintf(outputfile, "R %+16.8e %+16.8e %+16.8e",
							mbc->GetRefNodeR(1, 1), mbc->GetRefNodeR(1, 2), mbc->GetRefNodeR(1, 3));
						fprintf(outputfile, " %+16.8e %+16.8e %+16.8e",
							mbc->GetRefNodeR(2, 1), mbc->GetRefNodeR(2, 2), mbc->GetRefNodeR(2, 3));
						fprintf(outputfile, " %+16.8e %+16.8e %+16.8e\n",
							mbc->GetRefNodeR(3, 1), mbc->GetRefNodeR(3, 2), mbc->GetRefNodeR(3, 3));
						break;
					case MBC_ROT_THETA:
						fprintf(outputfile, "T %+16.8e %+16.8e %+16.8e\n",
							mbc->GetRefNodeTheta(1), mbc->GetRefNodeTheta(2), mbc->GetRefNodeTheta(3));
						break;

					case MBC_ROT_EULER_123:
						fprintf(outputfile, "E %+16.8e %+16.8e %+16.8e\n",
							mbc->GetRefNodeEuler123(1), mbc->GetRefNodeEuler123(2), mbc->GetRefNodeEuler123(3));
						break;
					}

					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
						mbc->GetRefNodeXP(1), mbc->GetRefNodeXP(2), mbc->GetRefNodeXP(3));
					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n\n",
						mbc->GetRefNodeOmega(1), mbc->GetRefNodeOmega(2), mbc->GetRefNodeOmega(3));

				} else if (mbc->bVerbose()) {

					if (mbc->bLabels()) {
						fprintf(stdout, "reference node (%u):\n", mbc->GetRefNodeKinematicsLabel());
					} else {
						fprintf(stdout, "reference node:\n");
					}
					fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->GetRefNodeX(1), mbc->GetRefNodeX(2), mbc->GetRefNodeX(3));
					switch (mbc->GetRefNodeRot()) {
					default:
						fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n",
							mbc->GetRefNodeR(1, 1), mbc->GetRefNodeR(1, 2), mbc->GetRefNodeR(1, 3));
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n",
							mbc->GetRefNodeR(2, 1), mbc->GetRefNodeR(2, 2), mbc->GetRefNodeR(2, 3));
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}}\n",
							mbc->GetRefNodeR(3, 1), mbc->GetRefNodeR(3, 2), mbc->GetRefNodeR(3, 3));
						break;
					case MBC_ROT_THETA:
						fprintf(stdout, " theta={%+16.8e,%+16.8e,%+16.8e\n}",
							mbc->GetRefNodeTheta(1), mbc->GetRefNodeTheta(2), mbc->GetRefNodeTheta(3));
						break;

					case MBC_ROT_EULER_123:
						fprintf(stdout, "euler123={%+16.8e,%+16.8e,%+16.8e}\n",
							mbc->GetRefNodeEuler123(1), mbc->GetRefNodeEuler123(2), mbc->GetRefNodeEuler123(3));
						break;
					}
					fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->GetRefNodeXP(1), mbc->GetRefNodeXP(2), mbc->GetRefNodeXP(3));
					fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->GetRefNodeOmega(1), mbc->GetRefNodeOmega(2), mbc->GetRefNodeOmega(3));
				}
			}
			if (mbc->GetNodes() > 0) {
				if (outfile) {
					fprintf(outputfile, "POS %u\n", mbc->GetNodes());
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(outputfile,"%d ", mbc->GetKinematicsLabel(n));
						}
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							 mbc->GetX(n, 1), mbc->GetX(n, 2), mbc->GetX(n, 3));
					}
					if (nomoments == 0) {
						fprintf(outputfile, "ROT %u\n", mbc->GetNodes());
						for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
							if (labels) {
								fprintf(outputfile, "%d ", mbc->GetKinematicsLabel(n));
							}
							switch (mbc->GetRot()) {
							default:
								fprintf(outputfile,
									"%+16.8e %+16.8e %+16.8e "
									"%+16.8e %+16.8e %+16.8e "
									"%+16.8e %+16.8e %+16.8e\n",
									mbc->GetR(n, 1, 1), mbc->GetR(n, 1, 2), mbc->GetR(n, 1, 3),
									mbc->GetR(n, 2, 1), mbc->GetR(n, 2, 2), mbc->GetR(n, 2, 3),
									mbc->GetR(n, 3, 1), mbc->GetR(n, 3, 2), mbc->GetR(n, 3, 3));
								break;

							case MBC_ROT_THETA:
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									mbc->GetTheta(n, 1), mbc->GetTheta(n, 2), mbc->GetTheta(n, 3));
								break;

							case MBC_ROT_EULER_123:
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									mbc->GetEuler123(n, 1), mbc->GetEuler123(n, 2), mbc->GetEuler123(n, 3));
								break;
							}
						}
					}
					fprintf(outputfile, "VEL %u\n", mbc->GetNodes());
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(outputfile, "%d ", mbc->GetKinematicsLabel(n));
						}
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							mbc->GetXP(n, 1), mbc->GetXP(n, 2), mbc->GetXP(n, 3));
					}
					if (nomoments == 0) {
						fprintf(outputfile, "W %u\n", mbc->GetNodes());
						for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
							if (labels) {
								fprintf(outputfile, "%d ", mbc->GetKinematicsLabel(n));
							}
							fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
								mbc->GetOmega(n, 1), mbc->GetOmega(n, 2), mbc->GetOmega(n, 3));
						}
					}

				} else if (mbc->bVerbose()) {
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(stdout, "node #%d (%u):\n", n - 1, mbc->GetKinematicsLabel(n));
						} else {
							fprintf(stdout, "node #%d:\n", n - 1);
						}
						fprintf(stdout, "    x=     %+16.8e %+16.8e %+16.8e\n",
							mbc->GetX(n, 1), mbc->GetX(n, 2), mbc->GetX(n, 3));
						if (nomoments == 0) {
							switch (mbc->GetRot()) {
							default:
								fprintf(stdout, "    R=     %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n",
									mbc->GetR(n, 1, 1), mbc->GetR(n, 1, 2), mbc->GetR(n, 1, 3),
									mbc->GetR(n, 2, 1), mbc->GetR(n, 2, 2), mbc->GetR(n, 2, 3),
									mbc->GetR(n, 3, 1), mbc->GetR(n, 3, 2), mbc->GetR(n, 3, 3));
								break;

							case MBC_ROT_THETA:
								fprintf(stdout, "    theta= %+16.8e %+16.8e %+16.8e\n",
									mbc->GetTheta(n, 1), mbc->GetTheta(n, 2), mbc->GetTheta(n, 3));
								break;

							case MBC_ROT_EULER_123:
								fprintf(stdout, " euler123= %+16.8e %+16.8e %+16.8e\n",
									mbc->GetEuler123(n, 1), mbc->GetEuler123(n, 2), mbc->GetEuler123(n, 3));
								break;
							}
						}
						fprintf(stdout, "    xp=    %+16.8e %+16.8e %+16.8e\n",
							mbc->GetXP(n, 1), mbc->GetXP(n, 2), mbc->GetXP(n, 3));
						if (nomoments == 0) {
							fprintf(stdout, "    omega= %+16.8e %+16.8e %+16.8e\n",
								mbc->GetOmega(n, 1), mbc->GetOmega(n, 2), mbc->GetOmega(n, 3));
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
				for (unsigned n = 0; n < mbc->GetNodes(); n++) {
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
				if (f0 != NULL) {
					mbc->SetRefNodeF(1, f0[0]);
					mbc->SetRefNodeF(2, f0[1]);
					mbc->SetRefNodeF(3, f0[2]);

					mbc->SetRefNodeM(1, f0[3]);
					mbc->SetRefNodeM(2, f0[4]);
					mbc->SetRefNodeM(3, f0[5]);

				} else {
					mbc->SetRefNodeF(1, 1.);
					mbc->SetRefNodeF(2, 2.);
					mbc->SetRefNodeF(3, 3.);

					mbc->SetRefNodeM(1, 4.);
					mbc->SetRefNodeM(1, 5.);
					mbc->SetRefNodeM(1, 6.);
				}
			}

			if (mbc->GetNodes() > 0) {
//				if (labels) {
//					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
//						mbc->GetDynamicsLabel(n) = mbc->GetKinematicsLabel(n);
//					}
//				}

				if (p0) {
					int size = 6;
					if (nomoments) {
						size = 3;
					}
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						mbc->SetF(n, 1, p0[size*(n - 1)]);
						mbc->SetF(n, 2, p0[size*(n - 1) + 1]);
						mbc->SetF(n, 3, p0[size*(n - 1) + 2]);

						if (nomoments == 0) {
							mbc->SetM(n, 1, p0[size*(n - 1) + 3]);
							mbc->SetM(n, 2, p0[size*(n - 1) + 4]);
							mbc->SetM(n, 3, p0[size*(n - 1) + 5]);
						}
					}

				} else {
					for (unsigned n = 1; n <= 3*mbc->GetNodes(); n++) {
						mbc->SetF((n - 1)/3 + 1, (n - 1)%3 + 1, (double)n);
						if (nomoments == 0) {
							mbc->SetM((n - 1)/3 + 1, (n - 1)%3 + 1, (double)n);
						}
					}
				}
			}

			/* sends forces
			 * second argument == 1 indicates convergence;
			 * otherwise MBDyn will send another solution
			 * and keep iterating */
            if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: about to call mbc->PutForces"  << std::endl; }
			if (mbc->PutForces((iter == niters - 1))) {
                if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: mbc->PutForces was not equal to niters - 1"  << std::endl; }
				goto done;
			}
		}
	}

done:;
	/* destroy class */
	if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: destroying mbc"  << std::endl; }
	delete mbc;

	if (p0) {
		free(p0);
	}
}


int main(int argc, char *argv[])
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
			if (optarg) {
				path = optarg;
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
				timeout = -1;
			} else {
				timeout = atoi(optarg);
			}
			break;

		case 'v':
			verbose = 1;
			break;

		case 'x':
			data_and_next = 1;
			break;

		default:
			usage();
		}
	}

	if (refnode && refnoderot == MBC_ROT_NONE) {
		refnoderot = rot;
	}

	if (nomoments) {
		rot = MBC_ROT_NONE;
	}

	/* initialize data structure:
	 */
	mbc = new MBCSharedMemNodal;
	std::cout << "Initializing MBCSharedMemNodal with path: " << path << std::endl;
	if (mbc->Initialize(refnoderot, nodes, labels, rot, accelerations, std::string(path))) {
		fprintf(stderr, "test_strext_socket: "
			"MBCNodal::Initialize() failed\n");
		usage();
	}

	if (verbose)
    {
        mbc->SetVerbose(true);
    }

	if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: after mbc->Initialize"  << std::endl; }


    mbc->Init();

    if (mbc->bVerbose()) { std::cout << "test_strext_sharedmemory_lib_cxx: after mbc->Init"  << std::endl; }

	/* "negotiate" configuration with MBDyn
	 * errors out if configurations are inconsistent */
	if (mbc->Negotiate()) {
		exit(EXIT_FAILURE);
	} else
    {
        test_run();
	}

	signal(SIGTERM, sh);
	signal(SIGINT, sh);
}

