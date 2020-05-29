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
/* Portions Copyright 2008 Marco Morandini */

/* Driver del programma */

#include "mbconfig.h" 		/* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <fstream>

#include "ac/getopt.h"

extern "C" {
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif /* HAVE_SYS_TIMES_H */
}

#if 0
#define _GNU_SOURCE 1
#include <fenv.h>
static void __attribute__ ((constructor))
trapfpe ()
{
	/* Enable some exceptions.  At startup all exceptions are masked.  */

	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

#ifdef USE_MPI
#include "mbcomm.h"
MPI::Intracomm MBDynComm = MPI::COMM_SELF;
#ifdef USE_EXTERNAL
#include <list>
std::list<MPI::Intercomm>  InterfaceComms;
#include <external.h>
#include <vector>

#define MB_EXIT(what, err) \
	do { \
		if (mbp.using_mpi) { \
			External::SendClose();	\
			if ((err) != EXIT_SUCCESS) { \
				MBDynComm.Abort((err)); \
			} \
			MPI::Finalize(); \
		} \
		what((err)); \
	} while (0)

#else // !USE_EXTERNAL

#define MB_EXIT(what, err) \
	do { \
		if (mbp.using_mpi) { \
			if ((err) != EXIT_SUCCESS) { \
				MBDynComm.Abort((err)); \
			} \
			MPI::Finalize(); \
		} \
		what((err)); \
	} while (0)
#endif // !USE_EXTERNAL
#else // !USE_MPI

#define MB_EXIT(what, err) \
	what((err))

#endif // !USE_MPI

#ifdef USE_RTAI
#include "mbrtai_utils.h"
void *rtmbdyn_rtai_task = NULL;
#endif /* USE_RTAI */

#include <sys/stat.h>
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif // !PATH_MAX

#ifdef HAVE_FENV_H
#define _GNU_SOURCE 1
#include <fenv.h>
#endif // HAVE_FENV_H

const char sDefaultOutputFileName[] = "MBDyn";

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "solver.h"
#include "invsolver.h"
#include "modules.h"
#include "legalese.h"

#include "cleanup.h"

enum InputFormat {
	MBDYN,
	LASTFORMAT
};

enum InputSource {
	MBFILE_UNKNOWN,
	MBFILE_STDIN,
	MBFILE_OPT,
	MBFILE_ARGS
};

struct mbdyn_proc_t {
	int iSleepTime;
	bool bShowSymbolTable;
	std::ifstream FileStreamIn;
	std::istream* pIn;
	std::string sInputFileName;
	std::string sOutputFileName;
	bool bException;
	bool bRedefine;
	bool bTable;
	Table* pT;
	MathParser* pMP;
	InputFormat CurrInputFormat;
	InputSource CurrInputSource;
	unsigned int nThreads;
	bool using_mpi;
#ifdef USE_MPI
	int MyRank;
	char *ProcessorName;
	int WorldSize;
	int parallel_fSilent;
	int parallel_fPedantic;
#ifdef HAVE_GETPID
	pid_t pid;
#else // ! HAVE_GETPID
	int pid;
#endif // ! HAVE_GETPID
#endif // USE_MPI
};

/* Note: DEBUG_* codes are declared in "mbdyn.h" */
#ifdef DEBUG
const debug_array da[] = {
	{ "input",		MYDEBUG_INPUT			},
	{ "sol",		MYDEBUG_SOL			},
	{ "assembly",		MYDEBUG_ASSEMBLY		},
	{ "derivatives",	MYDEBUG_DERIVATIVES		},
	{ "fsteps",		MYDEBUG_FSTEPS			},
	{ "mem",		MYDEBUG_MEM			},
	{ "fnames",		MYDEBUG_FNAMES			},
#if defined(__GNUC__)
	{ "prettyfn",		MYDEBUG_FNAMES|MYDEBUG_PRETTYFN	},
#endif /* __GNUC__ */
	{ "mpi",		MYDEBUG_MPI			},
	{ "pred",		MYDEBUG_PRED			},
	{ "residual",		MYDEBUG_RESIDUAL		},
	{ "jacobian",		MYDEBUG_JAC			},
	{ "init",		MYDEBUG_INIT			},
	{ "output",		MYDEBUG_OUTPUT			},

	{ NULL,			MYDEBUG_NONE			}
};
#endif /* DEBUG */

static void
mbdyn_version(void)
{
	silent_cout(std::endl
		<< "MBDyn - MultiBody Dynamics " << VERSION << std::endl
		<< "configured on " << __DATE__ << " at " << __TIME__
		<< std::endl);
}

#ifdef HAVE_GETOPT
static void
mbdyn_usage(const char *sShortOpts)
{
	mbdyn_version();
	silent_cout(std::endl
		<< "mbdyn is a multibody simulation program." << std::endl
		<< std::endl
		<< "usage: mbdyn [" << sShortOpts << "] [input-file list] " << std::endl
		<< std::endl
		<< "  -f, --input-file {file}   reads from 'file' instead of stdin" << std::endl
		<< "  -o, --output-file {file}  writes to '{file}.<ext>'" << std::endl
		<< "                             instead of '{input-file}.<ext>'" << std::endl
		<< "                             (or 'Mbdyn.<ext>' when input from stdin)" << std::endl
		<< std::endl);
#ifdef DEBUG
	silent_cout(
		   "  -d, --debug {level[:level[...]]}  when using the debug version of the code," << std::endl
		<< "                            enables debug levels; available:" << std::endl
		<< "                                none" << std::endl);
	for (int i = 0; da[i].s != NULL; i++) {
		silent_cout(
		   "                                " << da[i].s << std::endl);
	}
	silent_cout(
		   "                                any" << std::endl);
#endif /* DEBUG */
	silent_cout(
		   "  -e, --exceptions          don't trap exceptions to ease debugging" << std::endl
		<< "  -E, --fp-mask[=...]       enable some floating point checks" << std::endl
		<< "  -h, --help                prints this message" << std::endl
		<< "  -H, --show-table          print symbol table and exit" << std::endl
		<< "  -l, --license             prints the licensing terms" << std::endl
		/*
		<< "  -N, --threads             number of threads (need multithread support)" << std::endl
		 */
		);
#ifdef USE_MPI
	silent_cout(
		   "  -p, --parallel            required when run in parallel (invoked by mpirun)" << std::endl);
#endif /* USE_MPI */
	silent_cout(
		   "  -P, --pedantic            pedantic warning messages (more p's: noisier)" << std::endl
		<< "  -r, --redefine" << std::endl
		<< "  -R, --no-redefine         redefine/don't redefine symbols in table" << std::endl
		<< "  -s, --silent              runs quietly (more s': quieter)" << std::endl
		<< "  -S, --sleep [rank={n}] {time}  sleep {time} (on process {n}) before" << std::endl
		<< "                            starting the analysis" << std::endl
		<< "  -t, --same-table" << std::endl
		<< "  -T, --no-same-table       use/don't use same symbol table for multiple runs" << std::endl
		<< "  -v, --version             show version and exit" << std::endl
		<< "  -w, --warranty            prints the warranty conditions" << std::endl
		<< "  -W, --working-dir {dir}   sets the working directory" << std::endl
		<< std::endl
		<< "Usually mbdyn reads the input from stdin and writes messages on stdout; a log" << std::endl
		<< "is put in '{file}.out', and data output is sent to various '{file}.<ext>'" << std::endl
		<< "files ('Mbdyn.<ext>' if input from stdin)" << std::endl
		<< std::endl
		<< std::endl);
}

static void
mbdyn_welcome(void)
{
	mbdyn_version();
	silent_cout(std::endl
		<< "Copyright 1996-2017 (C) Paolo Mantegazza "
			"and Pierangelo Masarati," << std::endl
		<< "Dipartimento di Ingegneria Aerospaziale "
			"<http://www.aero.polimi.it/>" << std::endl
		<< "Politecnico di Milano                   "
			"<http://www.polimi.it/>" << std::endl
		<< std::endl
		<< "MBDyn is free software, covered by the"
			" GNU General Public License," << std::endl
		<< "and you are welcome to change it and/or "
			"distribute copies of it" << std::endl
		<< "under certain conditions.  Use 'mbdyn --license' "
			"to see the conditions." << std::endl
		<< "There is absolutely no warranty for"
			" MBDyn.  Use \"mbdyn --warranty\"" << std::endl
		<< "for details." << std::endl
		<< std::endl);
}

/* Dati di getopt */
static char sShortOpts[] = "d:eE::f:hHlN:o:pPrRsS:tTvwW:";

#ifdef HAVE_GETOPT_LONG
static struct option LongOpts[] = {
	{ "debug",          required_argument, NULL,           int('d') },
	{ "exceptions",     no_argument,       NULL,           int('e') },
	{ "fp-mask",        optional_argument, NULL,           int('E') },
	{ "input-file",     required_argument, NULL,           int('f') },
	{ "help",           no_argument,       NULL,           int('h') },
	{ "show-table",     no_argument,       NULL,           int('H') },
	{ "license",        no_argument,       NULL,           int('l') },
	{ "threads",	    required_argument, NULL,	       int('N') },
	{ "output-file",    required_argument, NULL,           int('o') },
	{ "parallel",	    no_argument,       NULL,           int('p') },
	{ "pedantic",	    no_argument,       NULL,           int('P') },
	{ "redefine",       no_argument,       NULL,           int('r') },
	{ "no-redefine",    no_argument,       NULL,           int('R') },
	{ "silent",         no_argument,       NULL,           int('s') },
	{ "sleep",          required_argument, NULL,           int('S') },
	{ "same-table",     no_argument,       NULL,           int('t') },
	{ "no-same-table",  no_argument,       NULL,           int('T') },
	{ "version",        no_argument,       NULL,           int('v') },
	{ "warranty",       no_argument,       NULL,           int('w') },
	{ "working-dir",    required_argument, NULL,           int('W') },

	{ NULL,             0,                 NULL,           0        }
};
#endif /* HAVE_GETOPT_LONG */
#endif /* HAVE_GETOPT */

const std::string sDefaultInputFileName("MBDyn");

extern void GetEnviron(MathParser&);

static Solver* RunMBDyn(MBDynParser&, const std::string&, const std::string&, unsigned int, bool, bool);

#ifdef USE_MPI
static int
parse_parallel_args(mbdyn_proc_t& mbp, int argc, char *argv[])
{
	for (int i = 1; i < argc; i++) {
		if (!mbp.using_mpi && strncmp(argv[i], "-p", 2) == 0) {
			mbp.using_mpi = true;
			continue;
		}

		/* intercept silence/pedantic flags */
		if (strncmp(argv[i], "-s", 2) == 0) {
			mbp.parallel_fSilent++;

			for (unsigned j = 2; argv[i][j] == 's'; j++) {
				mbp.parallel_fSilent++;
			}

		} else if (strncmp(argv[i], "-P", 2) == 0) {
			mbp.parallel_fPedantic++;

			for (unsigned j = 2; argv[i][j] == 'P'; j++) {
				mbp.parallel_fPedantic++;
			}
		}
	}

	return 0;
}
#endif /* USE_MPI */


void
mbdyn_parse_arguments(mbdyn_proc_t& mbp, int argc, char *argv[], int& currarg)
{
#ifdef HAVE_GETOPT
	/* Dati acquisibili da linea di comando */
	int iIndexPtr = 0;

	mbp.iSleepTime = -1;

	/* Parsing della linea di comando */
	opterr = 0;
	while (true) {
#ifdef HAVE_GETOPT_LONG
		int iCurrOpt = getopt_long(argc, argv, sShortOpts,
			LongOpts, &iIndexPtr);
#else /* !HAVE_GETOPT_LONG */
		int iCurrOpt = getopt(argc, argv, sShortOpts);
#endif /* !HAVE_GETOPT_LONG */

		if (iCurrOpt == EOF) {
			break;
		}

		switch (iCurrOpt) {
		case int('d'):
#ifdef DEBUG
			if (get_debug_options(optarg, da)) {
				silent_cerr("Unable to interpret debug"
					" option argument;"
					" using default" << std::endl);
				::debug_level = DEFAULT_DEBUG_LEVEL;
				/* throw ErrGeneric(MBDYN_EXCEPT_ARGS); */
			}
#else /* !DEBUG */
			silent_cerr("Compile with '-DDEBUG'"
				" to use debug features" << std::endl);
#endif /* !DEBUG */
			break;

		case int('e'):
			mbp.bException = true;
			break;

		case int('E'): {
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
			int excepts = FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW;
			if (optarg != 0) {
				excepts = 0;
				char *start = optarg, *next;
				do {
					next = strchr(start, ',');
					size_t len;
					if (next != NULL) {
						len = next - start;
						next++;
					} else {
						len = strlen(start);
					}

					if (strncasecmp(start, "invalid", len) == 0) {
						excepts |= FE_INVALID;
					} else if (strncasecmp(start, "divbyzero", len) == 0) {
						excepts |= FE_DIVBYZERO;
					} else if (strncasecmp(start, "overflow", len) == 0) {
						excepts |= FE_OVERFLOW;
					} else {
						silent_cerr("option 'E': unrecognized arg \"" << std::string(start, len) << "\"" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					start = next;
				} while (next != 0);

				pedantic_cout("Option 'E': exceptions=" << excepts << std::endl);
			}
			/* Enable some exceptions.  At startup all exceptions are masked.  */
			feenableexcept(excepts);
#else // ! HAVE_FENV_H
			silent_cerr("Option 'E' unsupported; ignored" << std::endl);
#endif // ! HAVE_FENV_H
			} break;

		case int('f'):
			mbp.CurrInputFormat = MBDYN;
			mbp.CurrInputSource = MBFILE_OPT;
 			mbp.sInputFileName = optarg;
#ifdef _WIN32
			// open the file in non translated mode in order not to break seek operations
			mbp.FileStreamIn.open(mbp.sInputFileName.c_str(), std::ios::binary);
#else
			mbp.FileStreamIn.open(mbp.sInputFileName.c_str());
#endif
			if (!mbp.FileStreamIn) {
				int save_errno = errno;
				silent_cerr(std::endl
					<< "Unable to open file "
					"\"" << mbp.sInputFileName << "\" (" << save_errno << ": " << strerror(save_errno) << ");"
					" aborting..." << std::endl);
#ifdef USE_MPI
				if (mbp.using_mpi) {
					silent_cerr(" on " << mbp.ProcessorName);
				}
#endif /* USE_MPI */
				silent_cerr(";" << std::endl
					<< "aborting..."
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			mbp.pIn = dynamic_cast<std::istream *>(&mbp.FileStreamIn);
			break;

		case int('h'):
			mbdyn_usage(sShortOpts);
			throw NoErr(MBDYN_EXCEPT_ARGS);

		case int('H'):
			mbp.bShowSymbolTable = true;
			break;

		case int('l'):
			mbdyn_welcome();
			mbdyn_license();
			throw NoErr(MBDYN_EXCEPT_ARGS);

		case int('N'):
			if (strcmp(optarg, "auto") == 0) {
				mbp.nThreads = 0;

			} else {
				char *next;
				long n = strtoul(optarg, &next, 10);
				if (next[0] != '\0' || next == optarg) {
					silent_cerr("Unable to parse threads number, option \"-N " << optarg << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (n < 1 || n >= std::numeric_limits<unsigned>::max()) {
					silent_cerr("Invalid number of threads, option \"-N " << n << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				mbp.nThreads = unsigned(n);
			}
			break;

		case int('o'):
#ifdef HAVE_GETCWD
			if (!is_abs_path(optarg)) {
				char cwd[PATH_MAX];
				if (getcwd(cwd, sizeof(cwd)) == 0) {
					silent_cerr("Unable to set output file: getcwd failed" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				mbp.sOutputFileName = std::string(cwd) + DIR_SEP + optarg;

			} else
#endif // HAVE_GETCWD
			{
				mbp.sOutputFileName = optarg;
			}
			break;

		case int('p'):
#ifdef USE_MPI
			ASSERT(mbp.using_mpi);
#else /* !USE_MPI */
			silent_cerr("switch '-p' is meaningless "
				"without MPI" << std::endl);
#endif /* !USE_MPI */
			break;

		case int('P'):
#ifdef USE_MPI
			if (mbp.parallel_fPedantic > 0) {
				mbp.parallel_fPedantic--;
			} else
#endif /* USE_MPI */
			{
				::fPedantic++;
			}
			break;

		case int('r'):
			mbp.bRedefine = true;
			break;

		case int('R'):
			mbp.bRedefine = false;
			break;

		case int('s'):
#ifdef USE_MPI
			if (mbp.parallel_fSilent > 0) {
				mbp.parallel_fSilent--;
			} else
#endif /* USE_MPI */
			{
				::fSilent++;
			}
			break;

		case int('S'):
			if (optarg) {
				char	*s = optarg;

				if (strncasecmp(s, "rank=", STRLENOF("rank=")) == 0) {
#ifdef USE_MPI
					char	*next;
					long	r;

					s += STRLENOF("rank=");
					errno = 0;
					r = strtol(s, &next, 10);
					int save_errno = errno;
					if (next[0] != '\0') {
						if (next[0] != ',') {
							silent_cerr("error in argument -S " << optarg << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						s = &next[1];

					} else if (save_errno == ERANGE) {
						silent_cerr("rank=" << s << ": overflow" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (mbp.using_mpi && r != mbp.MyRank) {
						break;
					}
#else /* ! USE_MPI */
					silent_cerr("option -S " << optarg << " valid only when --with-mpi" << std::endl);
#endif /* ! USE_MPI */
				}

				if (s[0]) {
					char	*next;

					errno = 0;
					mbp.iSleepTime = strtol(s, &next, 10);
					int save_errno = errno;
					if (next[0] != '\0') {
						silent_cerr("error in argument -S " << optarg << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);

					} else if (save_errno == ERANGE) {
						silent_cerr("argument of -S " << s << " overflows" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

			} else {
				/* default: 10 s */
				mbp.iSleepTime = 10;
			}

			break;

		case int('t'):
			mbp.bTable = true;
			break;

		case int('T'):
			mbp.bTable = false;
			break;

		case int('v'):
			mbdyn_version();
			throw NoErr(MBDYN_EXCEPT_ARGS);

		case int('w'):
			mbdyn_welcome();
			mbdyn_warranty();
			throw NoErr(MBDYN_EXCEPT_ARGS);

		case int('W'):
#ifdef HAVE_CHDIR
			if (chdir(optarg)) {
				int save_errno = errno;
				silent_cerr("Error in chdir(\"" << optarg << "\") "
					"(" << save_errno << ": " << strerror(save_errno) << ")"
					<< std::endl);
				throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
			}
#else /* !HAVE_CHDIR */
			silent_cerr("chdir() not available"
				<< std::endl);
#endif /* !HAVE_CHDIR */
			break;

		case int('?'):
			silent_cerr("Unknown option -" << char(optopt) << std::endl);
			mbdyn_usage(sShortOpts);
			throw NoErr(MBDYN_EXCEPT_ARGS);

		default:
			silent_cerr(std::endl
				<< "Unrecoverable error; aborting..."
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	} // while (true) end of endless loop

	/*
	 * primo argomento utile (potenziale nome di file di ingresso)
	 */
	currarg = optind;
#endif /* HAVE_GETOPT */
}

static int
mbdyn_prepare_files(const std::string& sInputFileName, std::string& sOutputFileName)
{
	std::string sNormalizedInputFileName;
	const char *dirsep = std::strrchr(sInputFileName.c_str(), DIR_SEP);
	const char *dot = std::strrchr(sInputFileName.c_str(), '.');
	if (dot != 0 && (dirsep == 0 || dot > dirsep)) {
		sNormalizedInputFileName = std::string(sInputFileName, 0, dot - sInputFileName.c_str());

	} else {
		sNormalizedInputFileName = sInputFileName;
	}

	// fix output file name
	if (sOutputFileName.empty()) {
		if (!sInputFileName.empty()) {
			sOutputFileName = sNormalizedInputFileName;

		} else {
			sOutputFileName = ::sDefaultOutputFileName;
		}

	} else {
		struct stat	s;

		if (stat(sOutputFileName.c_str(), &s) != 0) {
			int	save_errno = errno;

			/* if does not exist, check path */
			if (save_errno != ENOENT) {
				char	*errmsg = strerror(save_errno);

				silent_cerr("stat(" << sOutputFileName << ") failed "
					"(" << save_errno << ": " << errmsg << ")" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			// note: use a reverse iterator find?
			const char *path = std::strrchr(sOutputFileName.c_str(), DIR_SEP);
			if (path != NULL) {
				std::string sOutputFilePath(sOutputFileName, 0, path - sOutputFileName.c_str());
				// sOutputFilePath.erase(path - sOutputFileName.c_str());

				if (stat(sOutputFilePath.c_str(), &s) != 0) {
					save_errno = errno;
					char	*errmsg = strerror(save_errno);

					silent_cerr("stat(" << sOutputFileName << ") failed because "
							"stat(" << sOutputFilePath << ") failed "
						"(" << save_errno << ": " << errmsg << ")" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);

				} else if (!S_ISDIR(s.st_mode)) {
					silent_cerr("path to file \"" << sOutputFileName << "\" is invalid ("
							"\"" << sOutputFilePath << "\" is not a dir)" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

		} else if (S_ISDIR(s.st_mode)) {
			std::string tmpIn;

			if (!sInputFileName.empty()) {
				if (dirsep != 0) {
					tmpIn = dirsep;

				} else {
					tmpIn = sNormalizedInputFileName;
				}

			} else {
				tmpIn = ::sDefaultOutputFileName;
			}

			if (sOutputFileName[sOutputFileName.size() - 1] != DIR_SEP && tmpIn[0] != DIR_SEP) {
				sOutputFileName += DIR_SEP;
			}
			sOutputFileName += tmpIn;
		}
	}

#ifdef HAVE_GETCWD
	if (!is_abs_path(sOutputFileName.c_str())) {
   		char cwd[PATH_MAX];
   		if (getcwd(cwd, sizeof(cwd)) == 0) {
			silent_cerr("Error in getcwd()" << std::endl);
      			throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
   		}
		sOutputFileName = std::string(cwd) + DIR_SEP + sOutputFileName;
	}
#endif // HAVE_GETCWD

	// chdir to input file's folder
	if (dirsep != 0) {
		std::string sInputDir(sInputFileName.c_str(), dirsep - sInputFileName.c_str());
#ifdef HAVE_CHDIR
		if (chdir(sInputDir.c_str())) {
			int save_errno = errno;
			silent_cerr("Error in chdir(" << sInputDir << ") (" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
			throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
		}
#else // !HAVE_CHDIR
		silent_cerr("warning: chdir(2) not available; chdir(" << sInputDir.c_str() << ") not performed" << std::endl);
#endif // !HAVE_CHDIR
	}

	return 0;
}

static int
mbdyn_program(mbdyn_proc_t& mbp, int argc, char *argv[], int& currarg)
{
#ifndef HAVE_SYS_TIMES_H
	const clock_t start = clock();
#endif
	mbdyn_welcome();
#ifdef USE_MPI
	if (mbp.using_mpi) {
		silent_cerr("PID=" << mbp.pid << " Process " << mbp.MyRank
			<< " (" << mbp.MyRank + 1 << " of " << mbp.WorldSize
			<< ") is alive on " << mbp.ProcessorName
			<< std::endl);
	}
#endif /* USE_MPI */

	/* Mostra la tabella dei simboli ed esce */
	if (mbp.bShowSymbolTable) {
#ifdef USE_MPI
		if (mbp.MyRank == 0)
#endif /* USE_MPI */
		{
			Table t(true);
			MathParser mp(t);
			GetEnviron(mp);
			silent_cout("default symbol table:"
				<< std::endl << mp.GetSymbolTable() << std::endl);
		}

		throw NoErr(MBDYN_EXCEPT_ARGS);
	}

#ifdef USE_SLEEP
	if (mbp.iSleepTime > -1) {
		silent_cerr("sleeping " << mbp.iSleepTime << "s" << std::endl);
		sleep(mbp.iSleepTime);
	}
#endif // USE_SLEEP

	/* risolve l'input */
	if (mbp.CurrInputSource == MBFILE_UNKNOWN) {
		if (currarg < argc) {
			mbp.CurrInputSource = MBFILE_ARGS;

		} else {
			/*
			 * se non e' un argomento prende
			 * lo standard input
			 */
			mbp.CurrInputSource = MBFILE_STDIN;
			mbp.CurrInputFormat = MBDYN;
			ASSERT(mbp.pIn == NULL);
			mbp.pIn = dynamic_cast<std::istream *>(&std::cin);
		}
	}

	/* Gestione di input/output */
	int last = 0;
	while (last == 0) {
		if (mbp.CurrInputSource == MBFILE_STDIN) {
			silent_cout("reading from stdin" << std::endl);
			last = 1;

		} else if (mbp.CurrInputSource == MBFILE_OPT) {
			silent_cout("reading from file \""
				<< mbp.sInputFileName
				<< "\"" << std::endl);
			last = 1;

		} else if (mbp.CurrInputSource == MBFILE_ARGS) {
			mbp.sInputFileName = argv[currarg];
			silent_cout("reading from file \""
				<< mbp.sInputFileName
				<< "\"" << std::endl);

			/* incrementa il numero di argomento */
			if (++currarg == argc) {
				last = 1;
			}

			{
				mbp.CurrInputFormat = MBDYN;
#ifdef _WIN32
				// open the file in non translated mode in order not to break seek operations
				mbp.FileStreamIn.open(mbp.sInputFileName.c_str(), std::ios::binary);
#else
				mbp.FileStreamIn.open(mbp.sInputFileName.c_str());
#endif
				if (!mbp.FileStreamIn) {
					int save_errno = errno;
					silent_cerr(std::endl
						<< "Unable to open file "
						"\"" << mbp.sInputFileName << "\" (" << save_errno << ": " << strerror(save_errno) << ");"
						" aborting..." << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
			mbp.pIn = &mbp.FileStreamIn;
		}

		Solver* pSolv = NULL;
		switch (mbp.CurrInputFormat) {
		case MBDYN: {
			if (mbp.pT == NULL) {
				SAFENEWWITHCONSTRUCTOR(mbp.pT, Table, Table(true));
			}
			if (mbp.pMP == NULL) {
				SAFENEWWITHCONSTRUCTOR(mbp.pMP,
					MathParser,
					MathParser(*mbp.pT, mbp.bRedefine));

					/* legge l'environment */
				GetEnviron(*mbp.pMP);
			}

#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
			char buf[PATH_MAX];
			if (getcwd(buf, sizeof(buf)) == 0) {
				int save_errno = errno;
				silent_cerr("getcwd() failed (" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
				throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
			}
			std::string sOrigCWD(buf);
#endif // HAVE_GETCWD && HAVE_CHDIR

			std::string sOutputFileName = mbp.sOutputFileName;
			mbdyn_prepare_files(mbp.sInputFileName, sOutputFileName);

			/* stream in ingresso */
			InputStream In(*mbp.pIn);
			MBDynParser HP(*mbp.pMP, In,
				mbp.sInputFileName == sDefaultInputFileName ? "initial file" : mbp.sInputFileName.c_str());

			pSolv = RunMBDyn(HP, mbp.sInputFileName,
				sOutputFileName,
				mbp.nThreads, mbp.using_mpi, mbp.bException);
			if (mbp.FileStreamIn.is_open()) {
				mbp.FileStreamIn.close();
			}
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
			if (chdir(sOrigCWD.c_str())) {
				int save_errno = errno;
				silent_cerr("chdir(" << sOrigCWD.c_str() << ") failed (" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
				throw ErrFileSystem(MBDYN_EXCEPT_ARGS);
			}
#endif // HAVE_GETCWD && HAVE_CHDIR
			break;
		}

		default:
			silent_cerr("You shouldn't be here!" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		} // switch (CurrInputFormat)

		clock_t ct = 0;

		if (pSolv != NULL) {
			ct += pSolv->GetCPUTime();
			SAFEDELETE(pSolv);
		}

		if (!mbp.bTable || currarg == argc) {
			if (mbp.pMP != NULL) {
				SAFEDELETE(mbp.pMP);
				mbp.pMP = NULL;
			}
			if (mbp.pT != NULL) {
				SAFEDELETE(mbp.pT);
				mbp.pT = NULL;
			}
		}

#ifdef HAVE_SYS_TIMES_H
		time_t tSecs = 0;
		time_t tMils = 0;

		/* Tempo di CPU impiegato */
		struct tms tmsbuf;
		times(&tmsbuf);

		ct += tmsbuf.tms_utime + tmsbuf.tms_cutime
			+ tmsbuf.tms_stime + tmsbuf.tms_cstime;

		long clk_tck = sysconf(_SC_CLK_TCK);
		tSecs = ct/clk_tck;
		tMils = ((ct*1000)/clk_tck)%1000;

		char timebuf[sizeof(" (1000000000h 59m 59s 999ms)")];
		char *s = timebuf;
		size_t l = sizeof(timebuf), n;

		n = snprintf(s, l, "%ld.%03ld", tSecs, tMils);

		silent_cout(std::endl << "The simulation required "
			<< timebuf << " seconds of CPU time");

		if (tSecs > 60) {
			n = snprintf(s, l, " (" /* ) */ );
			s += n;
			l -= n;

			if (tSecs > 3600) {
				n = snprintf(s, l, "%ldh ", tSecs/3600);
				if (n >= l) {
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				s += n;
				l -= n;
			}

			n = snprintf(s, l,
					/* ( */ "%02ldm %02lds %03ldms)",
					(tSecs%3600)/60,
					(tSecs%3600)%60, tMils);

			s += n;
			l -= n;

			silent_cout(timebuf);
		}
#ifdef USE_MPI
		if (mbp.using_mpi) {
			silent_cout(" on " << mbp.ProcessorName
				<< " (" << mbp.MyRank + 1
				<< " of " << mbp.WorldSize << ")");
		}
#endif /* USE_MPI */
		silent_cout(std::endl);

#else /* ! HAVE_SYS_TIMES_H */
		silent_cout(std::endl << "Elapsed time "
			<< double(clock() - start)/CLOCKS_PER_SEC << " seconds\n");
#endif /* ! HAVE_SYS_TIMES_H */
	} // while (last == 0)
	throw NoErr(MBDYN_EXCEPT_ARGS);
}

int
main(int argc, char* argv[])
{
	int	rc = EXIT_SUCCESS;

	mbdyn_proc_t mbp;
       	mbp.bException = false;
       	mbp.bRedefine = false;
       	mbp.bTable = false;
	mbp.pT = 0;
	mbp.pMP = 0;
       	mbp.bShowSymbolTable = false;
	mbp.pIn = NULL;
       	mbp.sInputFileName = sDefaultInputFileName;
	mbp.nThreads = 0;
	mbp.iSleepTime = -1;
	mbp.CurrInputFormat = MBDYN;
	mbp.CurrInputSource = MBFILE_UNKNOWN;

	atexit(mbdyn_cleanup_destroy);

#ifdef USE_SOCKET
#ifdef _WIN32
    // winsock requires initialisation
    WSADATA wsa_data;
    if (WSAStartup(MAKEWORD(2,2),&wsa_data) != 0)
    {
        silent_cerr("Winsock initialisation failed. Error Code: "
                    << WSAGetLastError() << std::endl);
        return -1;
    }
    silent_cout("Winsock initialised successfully" << std::endl);
#endif /* _WIN32 */
#endif /* USE_SOCKET */

#ifdef USE_MPI
	char	ProcessorName_[1024] = "localhost";

	mbp.ProcessorName = ProcessorName_;
	mbp.WorldSize = 1;
	mbp.MyRank = 0;
	mbp.parallel_fSilent = 0,
	mbp.parallel_fPedantic = 0;

#ifdef HAVE_GETPID
	mbp.pid = getpid();
#endif // HAVE_GETPID

    	/*
	 * FIXME: this is a hack to detect whether mbdyn has been
	 * invoked thru mpirun (which means we need to initialize
	 * MPI) or not (which means we don't); need to check how
	 * portable it is ...
	 *
	 * the check is on the first two chars because "most" of
	 * the mpirun/MPI extra args start with -p<something>
	 */
	parse_parallel_args(mbp, argc, argv);

	::fSilent = mbp.parallel_fSilent;
	::fPedantic = mbp.parallel_fPedantic;

	if (mbp.using_mpi) {
		MPI::Init(argc, argv);
		mbp.WorldSize = MPI::COMM_WORLD.Get_size();
		mbp.MyRank = MPI::COMM_WORLD.Get_rank();

		if (mbp.MyRank > 0) {
			/*
			 * need a second take because MPI::Init()
			 * restores the inital args, so if Get_rank() > 0
			 * the eventual -s/-P flags have been restored
			 */
			parse_parallel_args(mbp, argc, argv);

			::fSilent = mbp.parallel_fSilent;
			::fPedantic = mbp.parallel_fPedantic;
		}

		/*
		 * all these temporaries are to avoid complains from
		 * the compiler (MPI's API is really messy ):
		 */
		int ProcessorNameLength = sizeof(ProcessorName_);
		MPI::Get_processor_name(mbp.ProcessorName, ProcessorNameLength);

		silent_cerr("using MPI (explicitly required by '-p*' switch)"
			<< std::endl);
	}
#endif /* USE_MPI */

    	/* primo argomento valido (potenziale nome di file di ingresso) */
    	int currarg = 0;

   	if (argc > 0) {
        	currarg = 1;
    	}

    	try {
		mbdyn_parse_arguments(mbp, argc, argv, currarg);
	} catch (NoErr& e) {
		silent_cout("MBDyn terminated normally" << std::endl);
		rc = EXIT_SUCCESS;
		MB_EXIT(return, rc);
	} catch (ErrInterrupted& e) {
		silent_cout("MBDyn was interrupted" << std::endl);
		rc = 2;
		MB_EXIT(exit, rc);
    	} catch (...) {
                silent_cerr("An error occurred during the execution of MBDyn" << std::endl);
                rc = EXIT_FAILURE;
                MB_EXIT(exit, rc);
        }

	if (mbp.bException) {
		mbdyn_program(mbp, argc, argv, currarg);

	} else {
	    	/* The program is a big try block */
	    	try {
			mbdyn_program(mbp, argc, argv, currarg);
		} catch (NoErr& e) {
			silent_cout("MBDyn terminated normally" << std::endl);
			rc = EXIT_SUCCESS;
		} catch (ErrInterrupted& e) {
			silent_cout("MBDyn was interrupted" << std::endl);
			rc = 2;
		} catch (std::ios::failure& err) {
			silent_cerr("An IO error occurred during the execution of MBDyn (" << err.what() << ");"
				" aborting... " << std::endl);
			rc = EXIT_FAILURE;
		} catch (std::bad_alloc& err) {
			silent_cerr("An allocation error occurred during the execution of MBDyn (" << err.what() << ");\n"
				"This means that - somewhere deep into the program - there was an allocation request\\"
				"for an amount of memory larger than what is currently available in your machine.\n\n\\"
				"Check your input file\n\n"
				"As a last resort, and to debug the issue, you can try to re-launch MBDyn \n"
				"within a debugger and with the run-time option -e .\n\n\\"
				"Aborting... " << std::endl);
			rc = EXIT_FAILURE;
#ifdef DEBUG
		} catch (const std::exception& e) {
			silent_cerr("An error occurred during the execution of MBDyn (" << e.what() << ");"
				" aborting..." << std::endl);
			rc = EXIT_FAILURE;
#endif // DEBUG
		} catch (...) {
			silent_cerr("An error occurred during the execution of MBDyn;"
				" aborting... " << std::endl);
			rc = EXIT_FAILURE;
			MB_EXIT(exit, rc);
		}
 	}

	if (mbp.pMP != 0) {
		SAFEDELETE(mbp.pMP);
		mbp.pMP = 0;
	}

	if (mbp.pT != 0) {
		SAFEDELETE(mbp.pT);
		mbp.pT = 0;
	}

#ifdef USE_RTAI
	if (::rtmbdyn_rtai_task) {
		(void)rtmbdyn_rt_task_delete(&::rtmbdyn_rtai_task);
		::rtmbdyn_rtai_task = NULL;
	}
#endif /* USE_RTAI */


#ifdef USE_SOCKET
#ifdef _WIN32
    // winsock requires cleanup
    WSACleanup();
#endif /* _WIN32 */
#endif /* USE_SOCKET */

	mbdyn_cleanup();

    	MB_EXIT(return, rc);
} // main() end


static Solver*
RunMBDyn(MBDynParser& HP,
	 const std::string& sInputFileName,
	 const std::string& sOutputFileName,
	 unsigned int nThreads,
	 bool using_mpi,
	 bool bException)
{
    	DEBUGCOUTFNAME("RunMBDyn");

    	Solver* pSolv(0);

    	/* NOTE: the flag "bParallel" states whether the analysis
	 * must be run in parallel or not; the flag "using_mpi"
	 * can be true because the "-p" switch was detected;
	 * it can be switched on by the "parallel" keyword in the
	 * "data" block, but the analysis can still be scalar if
	 * only one machine is available */
    	bool bParallel(false);

    	/* parole chiave */
    	const char* sKeyWords[] = {
        	"begin",
		"end",
        	"data",

		/* problem */
		"problem",
        	"integrator",			/* deprecated */

			/* problem types */
			"initial" "value",
	 		"multistep",		/* deprecated */
			"parallel" "initial" "value",
			"inverse" "dynamics",

		"parallel",			/* deprecated */

		NULL
    	};

    	/* enum delle parole chiave */
    	enum KeyWords {
        	UNKNOWN = -1,

        	BEGIN = 0,
		END,
        	DATA,

		PROBLEM,
        	INTEGRATOR,			/* deprecated */

			/* problem types */
			INITIAL_VALUE,
        		MULTISTEP,		/* deprecated */
			PARALLEL_INITIAL_VALUE,
			INVERSE_DYNAMICS,

		PARALLEL,			/* deprecated */

        	LASTKEYWORD
    	};

    	/* tabella delle parole chiave */
    	KeyTable K(HP, sKeyWords);

    	/* legge i dati della simulazione */
	/* NOTE: if the first GetDescription() fails because
	 * of an end-of-file, don't treat it as an error */
	KeyWords cd;
	try {
		cd = KeyWords(HP.GetDescription());
	} catch (EndOfFile& e) {
		throw NoErr(MBDYN_EXCEPT_ARGS);
	}
	/* looking for "begin"... */
	if (cd != BEGIN) {
        	silent_cerr(std::endl
	    		<< "Error: <begin> expected at line "
	    		<< HP.GetLineData() << "; aborting..." << std::endl);
        	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    	}

	/* looking for "data"... */
    	if (KeyWords(HP.GetWord()) != DATA) {
        	silent_cerr(std::endl
	    		<< "Error: <begin: data;> expected at line "
	    		<< HP.GetLineData() << "; aborting..." << std::endl);
        	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    	}

    	KeyWords CurrInt = INITIAL_VALUE;

    	/* Ciclo infinito */
    	while (true) {
        	switch (KeyWords(HP.GetDescription())) {
		case INTEGRATOR:
			pedantic_cout("statement \"integrator\" is deprecated; "
				"use \"problem\" instead." << std::endl);
        	case PROBLEM:
            		switch (KeyWords(HP.GetWord())) {
            		case MULTISTEP:
				pedantic_cout("deprecated \"multistep\" problem; "
					"use \"initial value\" instead" << std::endl);
			case INITIAL_VALUE:
	        		CurrInt = INITIAL_VALUE;
	        		break;

	        	case PARALLEL_INITIAL_VALUE:
	        		CurrInt = PARALLEL_INITIAL_VALUE;
#ifdef USE_MPI
				/* NOTE: use "parallel" in "data" block
				 * for models that should always be solved
				 * in parallel; otherwise this directive
				 * is superseded by the "-p" command-line
				 * switch */
				using_mpi = true;
	    			break;
#else /* !USE_MPI */
            			silent_cerr("compile with -DUSE_MPI to enable "
					"parallel solution" << std::endl);
	    			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_MPI */

			case INVERSE_DYNAMICS:
	        		CurrInt = INVERSE_DYNAMICS;
				break;

            		default:
	        		silent_cerr(std::endl
		    			<< "Unknown problem at line "
	            			<< HP.GetLineData()
					<< "; aborting..." << std::endl);
	        		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            		}
            		break;

        	case END:
	    		if (KeyWords(HP.GetWord()) != DATA) {
	        		silent_cerr(std::endl
		    			<< "Error: <end: data;> expected"
					" at line " << HP.GetLineData()
					<< "; aborting..." << std::endl);
	        		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	    		}
	    		goto endofcycle;

        	default:
	    		silent_cerr(std::endl
	        		<< "Unknown description at line "
	        		<< HP.GetLineData()
				<< "; aborting..." << std::endl);
	    		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        	}
    	}

   	/* Uscita dal ciclo infinito */
endofcycle:

#ifdef USE_MPI
	/* using_mpi is detected from command line switch "-p"
	 * or set after "parallel" config statement */
	if (using_mpi) {
		unsigned int size = MPI::COMM_WORLD.Get_size();
		unsigned int Bcount = 0;

#ifdef USE_EXTERNAL
		int* pBuff = NULL;
		SAFENEWARR(pBuff, int, size+1);
		pBuff[0] = 0;
		MPI::COMM_WORLD.Allgather(pBuff, 1, MPI::INT,
				pBuff+1, 1, MPI::INT);
		std::vector<unsigned int> iInterfaceProc(5);
		unsigned int Icount = 0;
		for (unsigned int i = 1; i <= size; i++) {
			if (pBuff[i] == pBuff[0]) {
				Bcount++;
			}
			if (pBuff[i] > 9) {
				if (Icount <= 5) {
					iInterfaceProc[Icount++] = i-1;

				} else {
					iInterfaceProc.push_back(i-1);
				}
			}
		}
		if (Bcount == size) {
			/* l'unica cosa che c'e' e' MBDyn */
			MBDynComm = MPI::COMM_WORLD.Dup();
		} else {
			MBDynComm = MPI::COMM_WORLD.Split(pBuff[0], 1);
		}
		if (Icount != 0) {
			for (unsigned int ii = 0; ii < Icount; ii++) {
				InterfaceComms.push_back(MBDynComm.Create_intercomm(0, MPI::COMM_WORLD, iInterfaceProc[ii], INTERF_COMM_LABEL));
			}
		}
		SAFEDELETEARR(pBuff);
#else /* ! USE_EXTERNAL */
		MBDynComm = MPI::COMM_WORLD.Dup();
		Bcount = size;
#endif /* ! USE_EXTERNAL */
		if (MBDynComm.Get_rank()) {
			silent_cout("MBDyn will run on " << Bcount
					<< " processors starting from "
					"processor " /* ??? */
					<< std::endl);
		}
       		bParallel = (MBDynComm.Get_size() != 1);

		if (bParallel) {
			/*
			 * E' necessario poter determinare in questa routine
			 * quale e' il master in modo da far calcolare la soluzione
			 * solo su di esso
			 */
			MyRank = MBDynComm.Get_rank();
			/* chiama il gestore dei dati generali della simulazione */

			/*
			 * I file di output vengono stampati localmente
			 * da ogni processo aggiungendo al termine
			 * dell'OutputFileName il rank del processo
			 */
			ASSERT(MPI::COMM_WORLD.Get_size() > 1);
			int iRankLength = 1 + (int)log10(MPI::COMM_WORLD.Get_size() - 1);

			char buf[STRLENOF(".") + iRankLength + 1];
			snprintf(buf, sizeof(buf), ".%.*d", iRankLength, MyRank);
			sOutputFileName += buf;
		}
	}
#endif /* USE_MPI */

    	switch (CurrInt) {
    	case INITIAL_VALUE:
    	case PARALLEL_INITIAL_VALUE:
        	SAFENEWWITHCONSTRUCTOR(pSolv,
				Solver,
				Solver(HP, sInputFileName,
					sOutputFileName, nThreads, bParallel));
		break;

    	case INVERSE_DYNAMICS:
        	SAFENEWWITHCONSTRUCTOR(pSolv,
				InverseSolver,
				InverseSolver(HP, sInputFileName,
					sOutputFileName, nThreads, bParallel));
		break;

	default:
        	silent_cerr("Unknown integrator; aborting..." << std::endl);
        	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

    	/* Runs the simulation */
	if (bException) {
		pSolv->Run();

	} else {
		try {
    			/* Runs the simulation */
	    		pSolv->Run();
		} catch (...) {
			if (pSolv) {
				SAFEDELETE(pSolv);
				pSolv = 0;
			}
			throw;
		}
	}

    	return pSolv;
} // RunMBDyn
