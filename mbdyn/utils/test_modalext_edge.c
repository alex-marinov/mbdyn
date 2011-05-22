/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
#include <string.h>
#include <strings.h>
#include <signal.h>
#include <errno.h>

volatile sig_atomic_t keep_going = 1;

int do_rename;
int sleeptime = 1;
int verbose;
enum {
	RM,
	MR
} order = RM;

static void
sh(int signum)
{
	keep_going = 0;
	signal(signum, SIG_DFL);
}

static const char *cmd2str(int cmd)
{
	switch (cmd) {
	case 0:
		return "EDGE is initializing; MBDyn waits";

	case 1:
		return "EDGE is busy; MBDyn waits";

	case 2:
		return "EDGE waits (is ready to read kinematics); MBDyn iterates";

	case 3:
		return "EDGE is computing; MBDyn waits before reading forces";

	case 4:
		return "EDGE converged; MBDyn advances one step";

	case 5:
		return "EDGE wants to end simulation";

	default:
		return "unknown";
	}
}

static int
check_flag(const char *flag, int sleeptime)
{
	int rc;

	while (1) {
		char buf[BUFSIZ];
		FILE *f;
		char c = ' ';

		f = fopen(flag, "r");
		if (f == NULL && errno == ENOENT) {
			fprintf(stderr, "test_modalext_edge: file \"%s\" missing\n", flag);
			return 1;
		}

		fgets(buf, sizeof(buf), f);
		if (strcmp(buf, "UPDATE,N,0,0,1\n") != 0) {
			size_t len = strlen(buf);
			buf[len - 1] = '\0';
			fprintf(stderr, "test_modalext_edge: expecting \"UPDATE,N,0,0,1\", got \"%s\" from file \"%s\"\n", buf, flag);
			fclose(f);
			return -1;
		}
		fgets(buf, sizeof(buf), f);
		if (strcmp(buf, "FLAG,I,1,1,0\n") != 0) {
			size_t len = strlen(buf);
			buf[len - 1] = '\0';
			fprintf(stderr, "test_modalext_edge: expecting \"FLAG,I,1,1,0\", got \"%s\" from file \"%s\"\n", buf, flag);
			fclose(f);
			return -1;
		}
		rc = fread((void *)&c, 1, 1, f);
		fclose(f);
		if (rc == 1) {
			fprintf(stderr, "test_modalext_edge: got %c (%s) from file \"%s\"\n", c, cmd2str(c - '0'), flag);

			switch (c) {
			case '0':
			case '1':
			case '3':
				return 0;

			case '5':
				return 1;

			default:
				break;
			}
		}

		if (sleeptime) {
			fprintf(stderr, "test_modalext_edge: sleeping %d s\n", sleeptime);
			usleep(1000000*sleeptime);
		}
	}

	return 0;
}

static int
put_flag(const char *flag, int cmd)
{
	FILE *f;
	char ftmpname[] = "mbedgeXXXXXX";

#ifdef HAVE_MKSTEMP
	if (do_rename) {
		int filedes = mkstemp(ftmpname);
		f = fdopen(filedes, "w");

	} else
#endif // HAVE_MKSTEMP
	{
		f = fopen(flag, "w");
	}

	if (f == NULL) {
		int save_errno = errno;
		fprintf(stderr, "unable to open flag file \"%s\" for writing (%d: %s)\n",
			flag, save_errno, strerror(save_errno));
		exit(EXIT_FAILURE);
	}

	fprintf(f, "UPDATE,N,0,0,1\n");
	fprintf(f, "FLAG,I,1,1,0\n");
	fprintf(f, "%d", cmd);
	fclose(f);

	if (do_rename) {
retry:;
		if (rename(ftmpname, flag) == -1) {
			switch (errno) {
			case EBUSY:
				usleep(1000000*sleeptime);
				goto retry;

			default: {
				int save_errno = errno;
				fprintf(stderr, "unable to rename flag file \"%s\" (errno=%d: %s)\n",
					flag, save_errno, strerror(save_errno));
				exit(EXIT_FAILURE);
				}
			}
		}
	}

	return 0;
}

static int
put_rdata(const char *rdata, double fm[6])
{
	FILE *f;
	int i;

	f = fopen(rdata, "w");
	if (f == NULL) {
		int save_errno = errno;

		fprintf(stderr, "unable to open rigid data file \"%s\" (%d: %s)\n",
			rdata, save_errno, strerror(save_errno));
		exit(EXIT_FAILURE);
	}

	fprintf(f,
		"* rigid-body forces and moments\n"
		"body_forces,R,1,6,0\n");
	for (i = 0; i < 6; i++) {
		fprintf(f, "%e ", fm[i]);
	}
	fputc('\n', f);
	fclose(f);

	return 0;
}

static int
put_mdata(const char *mdata, int modes, double *fg)
{
	FILE *f;
	int i;

	f = fopen(mdata, "w");
	if (f == NULL) {
		int save_errno = errno;

		fprintf(stderr, "unable to open modal data file \"%s\" (%d: %s)\n",
			mdata, save_errno, strerror(save_errno));
		exit(EXIT_FAILURE);
	}

	fprintf(f,
		"* modal forces\n"
		"modal_force_flow,R,%d,1,0\n",
		modes);
	for (i = 0; i < modes; i++) {
		fprintf(f, "%e ", fg[i]);
	}
	fputc('\n', f);
	fclose(f);

	return 0;
}

static int
do_rigid0(const char *rflag, const char *rdata, double *fm)
{
	if (rflag != NULL) {
		FILE *f = NULL;
		int i;

		f = fopen(rflag, "r");
		if (f == NULL) {
			int save_errno = errno;
			if (save_errno == ENOENT) {
				put_flag(rflag, 0);

			} else {
				fprintf(stderr, "unable to open rigid flag file \"%s\" (%d: %s)\n",
					rflag, save_errno, strerror(save_errno));
				exit(EXIT_FAILURE);
			}

		} else {
			fclose(f);
		}

		for (i = 0; i < 6; i++) {
			fm[i] = 0.1*(i + 1);
		}

		put_rdata(rdata, fm);
	}

	return 0;
}

static int
do_modal0(const char *mflag, const char *mdata, int modes, double **fgp)
{
	if (mflag != NULL) {
		FILE *f = NULL;
		int i;

		f = fopen(mflag, "r");
		if (f == NULL) {
			int save_errno = errno;
			if (save_errno == ENOENT) {
				put_flag(mflag, 0);

			} else {
				fprintf(stderr, "unable to open modal flag file \"%s\" (%d: %s)\n",
					mflag, save_errno, strerror(save_errno));
				exit(EXIT_FAILURE);
			}

		} else {
			fclose(f);
		}

		*fgp = (double *)malloc(sizeof(double)*modes);
		for (i = 0; i < modes; i++) {
			(*fgp)[i] = ((double)i)/10.0;
		}

		put_mdata(mdata, modes, *fgp);
	}

	return 0;
}

int
do_rigid(const char *rflag, const char *rdata,
	int niters, int *iterp, int cmd,
	double *fm)
{
	/* rigid */
	if (rflag != NULL) {
		if (check_flag(rflag, sleeptime)) {
			*iterp = niters;
			keep_going = 0;
			return 0;
		}

		if (verbose) {
			char buf[BUFSIZ];
			FILE *f;

			f = fopen(rdata, "r");
			if (f == NULL) {
				int save_errno = errno;

				fprintf(stderr, "unable to open rigid data file \"%s\" (%d: %s)\n",
					rdata, save_errno, strerror(save_errno));
				exit(EXIT_FAILURE);
			}

			while (fgets(buf, sizeof(buf), f) != NULL) {
				fprintf(stderr, ">> rdata:%s", buf);
			}

			fclose(f);
		}

		put_rdata(rdata, fm);
		put_flag(rflag, cmd);
	}

	return 0;
}

int
do_modal(const char *mflag, const char *mdata,
	int niters, int *iterp, int cmd,
	int modes, double *fg)
{
	/* modal */
	if (mflag != NULL) {
		if (check_flag(mflag, sleeptime)) {
			*iterp = niters;
			keep_going = 0;
			return 0;
		}

		if (verbose) {
			char buf[BUFSIZ];
			FILE *f;

			f = fopen(mdata, "r");
			if (f == NULL) {
				int save_errno = errno;

				fprintf(stderr, "unable to open modal data file \"%s\" (%d: %s)\n",
					mdata, save_errno, strerror(save_errno));
				exit(EXIT_FAILURE);
			}

			while (fgets(buf, sizeof(buf), f) != NULL) {
				fprintf(stderr, ">> mdata:%s", buf);
			}

			fclose(f);
		}

		put_mdata(mdata, modes, fg);
		put_flag(mflag, cmd);
	}

	return 0;
}

void
usage(void)
{
	fprintf(stderr,
		"usage: test_modalext_edge [options]\n"
		"\t-c [random:]<c>\t\tnumber of iterations\n"
		"\t-m [flag|data]=<file>\tmodal file names (set both)\n"
		"\t-M <modes>\t\tmodes number\n"
		"\t-n\t\t\tuse \"rename\" when writing flag files\n"
		"\t-o {rm|mr}\tprocess rigid,modal (rm) or modal,rigid (mr)\n"
		"\t-r [flag|data]=<file>\trigid-body file names (set both)\n"
		"\t-s <sleeptime>\t\tsleep time between tries\n"
		"\t-v\t\t\tverbose\n" );
	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	char *rflag = NULL;
	char *rdata = NULL;
	char *mflag = NULL;
	char *mdata = NULL;
	int iters = 1;
	int iters_random = 0;
	unsigned steps;
	int modes = 5;
	double fm[6];
	double *fg = NULL;

	while (1) {
		int opt = getopt(argc, argv, "c:m:M:no:r:s:v");

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
				fprintf(stderr, "iterations: %d\n", iters);
			}
			if (iters < 1) {
				fprintf(stderr, "test_modalext_edge: "
					"invalid sleep time %s\n",
					optarg);
				usage();
			}
			break;

		case 'm':
			if (strncasecmp(optarg, "flag=", sizeof("flag=") - 1) == 0) {
				mflag = &optarg[sizeof("flag=") - 1];

			} else if (strncasecmp(optarg, "data=", sizeof("data=") - 1) == 0) {
				mdata = &optarg[sizeof("data=") - 1];

			} else {
				fprintf(stderr, "test_modalext_edge: "
					"unknown modal file \"%s\"\n",
					optarg);
				usage();
			}
			break;

		case 'M':
			modes = atoi(optarg);
			if (modes <= 0) {
				fprintf(stderr, "test_modalext_edge: "
					"invalid mode number %s\n",
					optarg);
				usage();
			}
			break;

		case 'n':
#ifdef HAVE_MKSTEMP
			do_rename++;
#else // ! HAVE_MKSTEMP
			fprintf(stderr, "test_modalext_edge: "
				"'-n' meaningless\n");
#endif // ! HAVE_MKSTEMP
			break;

		case 'o':
			if (strcmp(optarg, "rm") == 0) {
				order = RM;

			} else if (strcmp(optarg, "mr") == 0) {
				order = MR;

			} else {
				fprintf(stderr, "test_modalext_edge: "
					"invalid order \"%s\"\n",
					optarg);
				usage();
			}
			break;

		case 'r':
			if (strncasecmp(optarg, "flag=", sizeof("flag=") - 1) == 0) {
				rflag = &optarg[sizeof("flag=") - 1];

			} else if (strncasecmp(optarg, "data=", sizeof("data=") - 1) == 0) {
				rdata = &optarg[sizeof("data=") - 1];

			} else {
				fprintf(stderr, "test_modalext_edge: "
					"unknown rigid file \"%s\"\n",
					optarg);
				usage();
			}
			break;

		case 's':
			sleeptime = atoi(optarg);
			if (sleeptime < 0) {
				fprintf(stderr, "test_modalext_edge: "
					"invalid iters %s\n",
					optarg);
				usage();
			}
			break;

		case 'v':
			verbose++;
			break;

		default:
			usage();
		}
	}

	if (mflag == NULL && mdata != NULL) {
		fprintf(stderr, "test_modalext_edge: "
			"need modal flag file "
			"along with modal data file \"%s\"\n",
			mdata);
		usage();
	}

	if (mflag != NULL && mdata == NULL) {
		fprintf(stderr, "test_modalext_edge: "
			"need modal data file "
			"along with modal flag file \"%s\"\n",
			mflag);
		usage();
	}

	if (rflag == NULL && rdata != NULL) {
		fprintf(stderr, "test_modalext_edge: "
			"need rigid flag file "
			"along with rigid data file \"%s\"\n",
			rdata);
		usage();
	}

	if (rflag != NULL && rdata == NULL) {
		fprintf(stderr, "test_modalext_edge: "
			"need rigid data file "
			"along with rigid flag file \"%s\"\n",
			rflag);
		usage();
	}

	if (mflag == NULL && rflag == NULL) {
		fprintf(stderr, "test_modalext_edge: "
			"need at least rigid or modal files\n");
		usage();
	}

	signal(SIGTERM, sh);
	signal(SIGINT, sh);

	switch (order) {
	case RM:
		do_rigid0(rflag, rdata, fm);
		do_modal0(mflag, mdata, modes, &fg);
		break;

	case MR:
		do_modal0(mflag, mdata, modes, &fg);
		do_rigid0(rflag, rdata, fm);
		break;
	}

	for (steps = 0; keep_going > 0; steps++) {
		int iter;
		int niters;

		if (iters_random) {
			niters = rand() % iters + 1;
			fprintf(stderr, "    iterations within this iter: %d\n", niters);

		} else {
			niters = iters;
		}

		for (iter = 0; iter < niters; iter++) {
			int cmd = 2;

			if (iter == niters - 1) {
				cmd = 4;
			}

			switch (order) {
			case RM:
				do_rigid(rflag, rdata, niters, &iter, cmd, fm);
				do_modal(mflag, mdata, niters, &iter, cmd, modes, fg);
				break;

			case MR:
				do_modal(mflag, mdata, niters, &iter, cmd, modes, fg);
				do_rigid(rflag, rdata, niters, &iter, cmd, fm);
				break;
			}

		}
	}

	if (rflag != NULL) {
		put_flag(rflag, 5);
	}

	if (mflag != NULL) {
		put_flag(mflag, 5);
	}

	return 0;
}
