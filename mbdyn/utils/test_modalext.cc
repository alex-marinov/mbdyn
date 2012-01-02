/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <unistd.h>
#include <signal.h>


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
		"\t-H <url>\tURL (file://path)\n"
		"\t-M <modes>\tmodes number\n"
		"\t-r\t\tuse reference node data\n"
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

	std::string fname_in, fname_out, tmpfname_out;
	unsigned modes = 0;
	unsigned refnode = 0;
	char verbose = 0;
	char data_and_next = 0;

	double x[3], R[3][3], v[3], w[3], f[3], m[3];
	double *q = 0, *qp = 0, *p = 0;

	while (1) {
		int opt = getopt(argc, argv, "c:i:M:o:rs:vx");

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
					"invalid iterations %s\n",
					optarg);
				usage();
			}
			break;

		case 'i':
			fname_in = optarg;
			break;

		case 'M':
			modes = atoi(optarg);
			if (modes <= 0) {
				fprintf(stderr, "testedge: "
					"invalid mode number %s\n",
					optarg);
				usage();
			}
			break;

		case 'o':
			fname_out = optarg;
			break;

		case 'r':
			refnode = 1;
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
			verbose = 1;
			break;

		case 'x':
			data_and_next = 1;
			break;

		default:
			usage();
		}
	}

	if (fname_in.empty() || fname_out.empty()) {
		usage();
	}

	tmpfname_out = fname_out + ".tmp";

	if (modes > 0) {
		q = new double[3*modes];
		qp = &q[modes];
		p = &qp[modes];
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
			std::ifstream fin;

retry_1:;
			fin.open(fname_in.c_str());
			if (!fin) {
				if (!keep_going) {
					printf("interrupted\n");
					goto done;
				}

				usleep(1000000);
				goto retry_1;
			}

			/* get motion */
			if (refnode) {
				fin
					>> x[0] >> x[1] >> x[2]
					>> R[0][0] >> R[0][1] >> R[0][2]
					>> R[1][0] >> R[1][1] >> R[1][2]
					>> R[2][0] >> R[2][1] >> R[2][2]
					>> v[0] >> v[1] >> v[2]
					>> w[0] >> w[1] >> w[2];
				if (!fin) {
					goto done;
				}

				fprintf(stdout, "x={%+16.8e,%+16.8e,%+16.8e}\n", x[0], x[1], x[2]);
				fprintf(stdout, "R={{%+16.8e,%+16.8e,%+16.8e};\n", R[0][0], R[0][1], R[0][2]);
				fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e};\n", R[1][0], R[1][1], R[1][2]);
				fprintf(stdout, "   {%+16.8e,%+16.8e,%+16.8e}};\n", R[2][0], R[2][1], R[2][2]);
				fprintf(stdout, "v={%+16.8e,%+16.8e,%+16.8e}\n", v[0], v[1], v[2]);
				fprintf(stdout, "w={%+16.8e,%+16.8e,%+16.8e}\n", w[0], w[1], w[2]);
			}

			if (modes) {
				for (unsigned m = 0; m < modes; m++) {
					fin >> q[m] >> qp[m];
					if (!fin) {
						goto done;
					}
					fprintf(stdout, "mode #%d: %+16.8e %+16.8e\n", m, q[m], qp[m]);
				}
			}

			fin.close();
			unlink(fname_in.c_str());

			if (sleeptime) {
				usleep(1000000*sleeptime);
			}

			/* set forces */
			std::ofstream fout(tmpfname_out.c_str());

			if (refnode) {
				f[0] = 1.;
				f[1] = 2.;
				f[2] = 3.;

				m[0] = 4.;
				m[1] = 5.;
				m[2] = 6.;

				fout
					<< f[0] << ' ' << f[1] << ' ' << f[2] << std::endl
					<< m[0] << ' ' << m[1] << ' ' << m[2] << std::endl;
			}

			if (modes) {
				for (unsigned m = 0; m < modes; m++) {
					p[m] = (double)(m + 1);

					fout << p[m] << std::endl;
				}
			}

			rename(tmpfname_out.c_str(), fname_out.c_str());
		}
	}

done:;

	if (modes && q) {
		delete[] q;
	}

	return 0;
}
