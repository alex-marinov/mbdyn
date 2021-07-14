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
#include "mbcxx.h"

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
static MBCBase::Rot refnoderot = MBCBase::NONE;
static int nodes = 0;
static int labels = 0;
static int accelerations = 0;
static MBCBase::Rot rot = MBCBase::MAT;

static MBCNodal *mbc;

static char *path;
static char *host;
static unsigned short int port = -1;

static double fx[6], *f0 = NULL;
static double *p0 = NULL;

static int timeout;
static int verbose;
static int data_and_next;

static int inpfile = 0;
static int outfile = 0;
static FILE *inputfile  = NULL;
static FILE *outputfile = NULL;

extern "C" void
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
				rot = MBCBase::MAT;

			} else if (strcasecmp(optarg, "theta") == 0) {
				rot = MBCBase::THETA;

			} else if (strcasecmp(optarg, "euler123") == 0) {
				rot = MBCBase::EULER_123;

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

	if (refnode && refnoderot == MBCBase::NONE) {
		refnoderot = rot;
	}

	if (nomoments) {
		rot = MBCBase::NONE;
	}

	/* initialize data structure:
	 */
	mbc = new MBCNodal;
	if (mbc->Initialize(refnoderot, nodes, labels, rot, accelerations)) {
		fprintf(stderr, "test_strext_socket: "
			"MBCNodal::Initialize() failed\n");
		usage();
	}

	if (path) {
#ifdef _WIN32
		fprintf(stderr, "test_strext_socket: "
				"Windows does not support local sockets\n ");
		exit(EXIT_FAILURE);
#else
		/* initialize UNIX socket (path) */
		if (mbc->Init(path)) {
			exit(EXIT_FAILURE);
		}
#endif /* _WIN32 */
	} else if (host) {
		/* initialize INET socket (host, port) */
		if (mbc->Init(host, port)) {
			exit(EXIT_FAILURE);
		}

	} else {
		usage();
	}

	/* "negotiate" configuration with MBDyn
	 * errors out if configurations are inconsistent */
	if (mbc->Negotiate()) {
		exit(EXIT_FAILURE);
	}

	signal(SIGTERM, sh);
	signal(SIGINT, sh);
}

/*
 * specific to test_strext_socket
 */
extern "C" void
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
			if (mbc->GetMotion()) {
				goto done;
			}

			if (outfile)  {
				fprintf(outputfile, "STEP %u ITERATION %d\n", steps, iter);
			}

			if (refnode) {
				if (outfile) {
					if (labels) {
						fprintf(outputfile, "REF (%u)\n",
							mbc->KinematicsLabel());
					} else {
						fprintf(outputfile, "REF \n");
					}

					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
						mbc->X(1), mbc->X(2), mbc->X(3));
		
					switch (mbc->GetRefNodeRot()) {
					default:
						fprintf(outputfile, "R %+16.8e %+16.8e %+16.8e",
							mbc->R(1, 1), mbc->R(1, 2), mbc->R(1, 3));
						fprintf(outputfile, " %+16.8e %+16.8e %+16.8e",
							mbc->R(2, 1), mbc->R(2, 2), mbc->R(2, 3));
						fprintf(outputfile, " %+16.8e %+16.8e %+16.8e\n",
							mbc->R(3, 1), mbc->R(3, 2), mbc->R(3, 3));
						break;
					case MBCBase::THETA:
						fprintf(outputfile, "T %+16.8e %+16.8e %+16.8e\n",
							mbc->Theta(1), mbc->Theta(2), mbc->Theta(3));
						break;

					case MBCBase::EULER_123:
						fprintf(outputfile, "E %+16.8e %+16.8e %+16.8e\n",
							mbc->Euler123(1), mbc->Euler123(2), mbc->Euler123(3));
						break;
					}
	
					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
						mbc->XP(1), mbc->XP(2), mbc->XP(3));
					fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n\n",
						mbc->Omega(1), mbc->Omega(2), mbc->Omega(3));

				} else if (mbc->bVerbose()) {

					if (mbc->bLabels()) {
						fprintf(stdout, "reference node (%u):\n", mbc->KinematicsLabel());
					} else {
						fprintf(stdout, "reference node:\n");
					}
					fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->X(1), mbc->X(2), mbc->X(3));
					switch (mbc->GetRefNodeRot()) {
					default:
						fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n",
							mbc->R(1, 1), mbc->R(1, 2), mbc->R(1, 3));
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n",
							mbc->R(2, 1), mbc->R(2, 2), mbc->R(2, 3));
						fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}}\n",
							mbc->R(3, 1), mbc->R(3, 2), mbc->R(3, 3));
						break;
					case MBCBase::THETA:
						fprintf(stdout, " theta={%+16.8e,%+16.8e,%+16.8e\n}",
							mbc->Theta(1), mbc->Theta(2), mbc->Theta(3));
						break;

					case MBCBase::EULER_123:
						fprintf(stdout, "euler123={%+16.8e,%+16.8e,%+16.8e}\n",
							mbc->Euler123(1), mbc->Euler123(2), mbc->Euler123(3));
						break;
					}
					fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->XP(1), mbc->XP(2), mbc->XP(3));
					fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n",
						mbc->Omega(1), mbc->Omega(2), mbc->Omega(3));
				} 
			}
			if (mbc->GetNodes() > 0) {
				if (outfile) {
					fprintf(outputfile, "POS %u\n", mbc->GetNodes());
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(outputfile,"%d ", mbc->KinematicsLabel(n));
						} 
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							 mbc->X(n, 1), mbc->X(n, 2), mbc->X(n, 3));
					}
					if (nomoments == 0) {
						fprintf(outputfile, "ROT %u\n", mbc->GetNodes());
						for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
							if (labels) {
								fprintf(outputfile, "%d ", mbc->KinematicsLabel(n));
							}	
							switch (mbc->GetRot()) {
							default:
								fprintf(outputfile,
									"%+16.8e %+16.8e %+16.8e "
									"%+16.8e %+16.8e %+16.8e "
									"%+16.8e %+16.8e %+16.8e\n",
									mbc->R(n, 1, 1), mbc->R(n, 1, 2), mbc->R(n, 1, 3),
									mbc->R(n, 2, 1), mbc->R(n, 2, 2), mbc->R(n, 2, 3),
									mbc->R(n, 3, 1), mbc->R(n, 3, 2), mbc->R(n, 3, 3));
								break;

							case MBCBase::THETA:
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									mbc->Theta(n, 1), mbc->Theta(n, 2), mbc->Theta(n, 3));
								break;

							case MBCBase::EULER_123:
								fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
									mbc->Euler123(n, 1), mbc->Euler123(n, 2), mbc->Euler123(n, 3));
								break;
							}
						}	
					}
					fprintf(outputfile, "VEL %u\n", mbc->GetNodes());
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(outputfile, "%d ", mbc->KinematicsLabel(n));
						}
						fprintf(outputfile,"%+16.8e %+16.8e %+16.8e\n",
							mbc->XP(n, 1), mbc->XP(n, 2), mbc->XP(n, 3));
					}
					if (nomoments == 0) {
						fprintf(outputfile, "W %u\n", mbc->GetNodes());
						for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
							if (labels) {
								fprintf(outputfile, "%d ", mbc->KinematicsLabel(n));
							}	
							fprintf(outputfile, "%+16.8e %+16.8e %+16.8e\n",
								mbc->Omega(n, 1), mbc->Omega(n, 2), mbc->Omega(n, 3));
						}
					}

				} else if (mbc->bVerbose()) {
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						if (labels) {
							fprintf(stdout, "node #%d (%u):\n", n - 1, mbc->KinematicsLabel(n));
						} else {
							fprintf(stdout, "node #%d:\n", n - 1);
						}
						fprintf(stdout, "    x=     %+16.8e %+16.8e %+16.8e\n",
							mbc->X(n, 1), mbc->X(n, 2), mbc->X(n, 3));
						if (nomoments == 0) {
							switch (mbc->GetRot()) {
							default:
								fprintf(stdout, "    R=     %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n"
										"           %+16.8e %+16.8e %+16.8e\n",
									mbc->R(n, 1, 1), mbc->R(n, 1, 2), mbc->R(n, 1, 3),
									mbc->R(n, 2, 1), mbc->R(n, 2, 2), mbc->R(n, 2, 3),
									mbc->R(n, 3, 1), mbc->R(n, 3, 2), mbc->R(n, 3, 3));
								break;

							case MBCBase::THETA:
								fprintf(stdout, "    theta= %+16.8e %+16.8e %+16.8e\n",
									mbc->Theta(n, 1), mbc->Theta(n, 2), mbc->Theta(n, 3));
								break;

							case MBCBase::EULER_123:
								fprintf(stdout, " euler123= %+16.8e %+16.8e %+16.8e\n",
									mbc->Euler123(n, 1), mbc->Euler123(n, 2), mbc->Euler123(n, 3));
								break;
							}
						}
						fprintf(stdout, "    xp=    %+16.8e %+16.8e %+16.8e\n",
							mbc->XP(n, 1), mbc->XP(n, 2), mbc->XP(n, 3));
						if (nomoments == 0) {
							fprintf(stdout, "    omega= %+16.8e %+16.8e %+16.8e\n",
								mbc->Omega(n, 1), mbc->Omega(n, 2), mbc->Omega(n, 3));
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
					mbc->F(1) = f0[0];
					mbc->F(2) = f0[1];
					mbc->F(3) = f0[2];

					mbc->M(1) = f0[3];
					mbc->M(2) = f0[4];
					mbc->M(3) = f0[5];

				} else {
					mbc->F(1) = 1.;
					mbc->F(2) = 2.;
					mbc->F(3) = 3.;

					mbc->M(1) = 4.;
					mbc->M(2) = 5.;
					mbc->M(3) = 6.;
				}
			}

			if (mbc->GetNodes() > 0) {
				if (labels) {
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
					}
				}

				if (p0) {
					int size = 6;
					if (nomoments) {
						size = 3;
					}
					for (unsigned n = 1; n <= mbc->GetNodes(); n++) {
						mbc->F(n, 1) = p0[size*(n - 1)];
						mbc->F(n, 2) = p0[size*(n - 1) + 1];
						mbc->F(n, 3) = p0[size*(n - 1) + 2];

						if (nomoments == 0) {
							mbc->M(n, 1) = p0[size*(n - 1) + 3];
							mbc->M(n, 2) = p0[size*(n - 1) + 4];
							mbc->M(n, 3) = p0[size*(n - 1) + 5];
						}
					}

				} else {
					for (unsigned n = 1; n <= 3*mbc->GetNodes(); n++) {
						mbc->F((n - 1)/3 + 1, (n - 1)%3 + 1) = (double)n;
						if (nomoments == 0) {
							mbc->M((n - 1)/3 + 1, (n - 1)%3 + 1) = (double)n;
						}
					}
				}
			}

			/* sends forces
			 * second argument == 1 indicates convergence;
			 * otherwise MBDyn will send another solution
			 * and keep iterating */
			if (mbc->PutForces((iter == niters - 1))) {
				goto done;
			}
		}
	}

done:;
	/* destroy data structure and close socket */
	delete mbc;

	if (p0) {
		free(p0);
	}
}

/*
 * specific to test_strext_socket_f
 */
extern "C" void
tdata_(int32_t *REFNODE, int32_t *NODES, int32_t *ROT, int32_t *ITERS, int32_t *VERB,
	int32_t *RC_P)
{
	switch (mbc->GetRot()) {
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

	if (mbc->bLabels()) {
		*RC_P = 1;
		return;
	}

	if (mbc->bAccelerations()) {
		*RC_P = 1;
		return;
	}

	*REFNODE = mbc->bRefNode();
	*NODES = mbc->GetNodes();
	*VERB = mbc->bVerbose();

	*ITERS = iters;

	return;
}

extern "C" void
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

	if (mbc->GetNodes() > 0) {
		if (p0) {
			for (unsigned n = 0; n < mbc->GetNodes(); n++) {
				NF[3*n] = p0[6*n];
				NF[3*n + 1] = p0[6*n + 1];
				NF[3*n + 2] = p0[6*n + 2];

				NM[3*n] = p0[6*n + 3];
				NM[3*n + 1] = p0[6*n + 4];
				NM[3*n + 2] = p0[6*n + 5];
			}

		} else {
			for (unsigned n = 0; n < 3*mbc->GetNodes(); n++) {
				NF[n] = (double)(n + 1);
				NM[n] = (double)(n + 1);
			}
		}
	}

	return;
}

extern "C" void
tsend_(float *RF, float *RM, float *NF, float *NM,
        int32_t *CONVERGED_P, int32_t *RC_P)
{
	/* set forces */
	if (mbc->bRefNode()) {
		mbc->F(1) = RF[0];
		mbc->F(2) = RF[1];
		mbc->F(3) = RF[2];

		mbc->M(1) = RM[0];
		mbc->M(2) = RM[1];
		mbc->M(3) = RM[2];
	}

	if (mbc->GetNodes() > 0) {
		for (unsigned node = 1; node <= mbc->GetNodes(); node++) {
			mbc->F(node, 1) = NF[3*(node - 1)];
			mbc->F(node, 2) = NF[3*(node - 1) + 1];
			mbc->F(node, 3) = NF[3*(node - 1) + 2];

			mbc->M(node, 1) = NM[3*(node - 1)];
			mbc->M(node, 2) = NM[3*(node - 1) + 1];
			mbc->M(node, 3) = NM[3*(node - 1) + 2];
		}
	}

	/* NOTE: the flag indicates convergence when not 0 */
	if (mbc->PutForces(*CONVERGED_P)) {
		*RC_P = 1;
		return;
	}

	*RC_P = 0;
	return;
}

extern "C" void
trecv_(float *RX, float *RR, float *RXP, float *ROMEGA,
	float *NX, float *NR, float *NXP, float *NOMEGA, int32_t *RC_P)
{
	if (mbc->GetMotion()) {
		*RC_P = 1;
		return;
	}

	if (mbc->bRefNode()) {
		RX[0] = mbc->X(1);
		RX[1] = mbc->X(2);
		RX[2] = mbc->X(3);

		switch (mbc->GetRot()) {
		case MBCBase::MAT:
			RR[0] = mbc->R(1, 1);
			RR[1] = mbc->R(2, 1);
			RR[2] = mbc->R(3, 1);
			RR[3] = mbc->R(1, 2);
			RR[4] = mbc->R(2, 2);
			RR[5] = mbc->R(3, 2);
			RR[6] = mbc->R(1, 3);
			RR[7] = mbc->R(2, 3);
			RR[8] = mbc->R(3, 3);
			break;

		case MBCBase::THETA:
			RR[0] = mbc->Theta(1);
			RR[1] = mbc->Theta(2);
			RR[2] = mbc->Theta(3);
			break;

		case MBCBase::EULER_123:
			RR[0] = mbc->Euler123(1);
			RR[1] = mbc->Euler123(2);
			RR[2] = mbc->Euler123(3);
			break;

		default:
			*RC_P = 1;
			return;
		}

		RXP[0] = mbc->XP(1);
		RXP[1] = mbc->XP(2);
		RXP[2] = mbc->XP(3);

		ROMEGA[0] = mbc->Omega(1);
		ROMEGA[1] = mbc->Omega(2);
		ROMEGA[2] = mbc->Omega(3);
	}

	if (mbc->GetNodes() > 0) {
		for (unsigned node = 1; node <= mbc->GetNodes(); node++) {
			NX[3*(node - 1)] = mbc->X(node, 1);
			NX[3*(node - 1) + 1] = mbc->X(node, 2);
			NX[3*(node - 1) + 2] = mbc->X(node, 3);

			switch (mbc->GetRot()) {
			case MBCBase::MAT:
				NR[9*(node - 1)] = mbc->R(node, 1, 1);
				NR[9*(node - 1) + 1] = mbc->R(node, 2, 1);
				NR[9*(node - 1) + 2] = mbc->R(node, 3, 1);
				NR[9*(node - 1) + 3] = mbc->R(node, 1, 2);
				NR[9*(node - 1) + 4] = mbc->R(node, 2, 2);
				NR[9*(node - 1) + 5] = mbc->R(node, 3, 2);
				NR[9*(node - 1) + 6] = mbc->R(node, 1, 3);
				NR[9*(node - 1) + 7] = mbc->R(node, 2, 3);
				NR[9*(node - 1) + 8] = mbc->R(node, 3, 3);
				break;

			case MBCBase::THETA:
				NR[3*(node - 1)] = mbc->Theta(1);
				NR[3*(node - 1) + 1] = mbc->Theta(2);
				NR[3*(node - 1) + 2] = mbc->Theta(3);
				break;

			case MBCBase::EULER_123:
				NR[3*(node - 1)] = mbc->Euler123(1);
				NR[3*(node - 1) + 1] = mbc->Euler123(2);
				NR[3*(node - 1) + 2] = mbc->Euler123(3);
				break;

			default:
				*RC_P = 1;
				return;
			}

			NXP[3*(node - 1)] = mbc->XP(node, 1);
			NXP[3*(node - 1) + 1] = mbc->XP(node, 2);
			NXP[3*(node - 1) + 2] = mbc->XP(node, 3);

			NOMEGA[3*(node - 1)] = mbc->Omega(node, 1);
			NOMEGA[3*(node - 1) + 1] = mbc->Omega(node, 2);
			NOMEGA[3*(node - 1) + 2] = mbc->Omega(node, 3);
		}
	}

	*RC_P = 0;
	return;
}

