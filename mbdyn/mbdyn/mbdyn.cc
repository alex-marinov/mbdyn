/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* Driver del programma */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h> 		/* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/fstream>
#include <ac/getopt.h>

extern "C" {
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif /* HAVE_SYS_TIMES_H */
}

#ifdef USE_MPI 
#include <mbcomm.h>
MPI::Intracomm MBDynComm = MPI::COMM_SELF;
#ifdef USE_EXTERNAL
#include<list>
std::list<MPI::Intercomm>  InterfaceComms;
#include<external.h>
#include<vector>

#define MB_EXIT(err) \
	do { \
		if (using_mpi) { \
			External::SendClose();	\
			if ((err) != EXIT_SUCCESS) { \
				MBDynComm.Abort((err)); \
			} \
			MPI::Finalize(); \
		} \
		exit((err)); \
	} while (0)

#else /* USE_EXTERNAL */

#define MB_EXIT(err) \
	do { \
		if (using_mpi) { \
			if ((err) != EXIT_SUCCESS) { \
				MBDynComm.Abort((err)); \
			} \
			MPI::Finalize(); \
		} \
		exit((err)); \
	} while (0)
#endif /* USE_EXTERNAL */
#else /* !USE_MPI */

#define MB_EXIT(err) \
	exit((err))

#endif /* !USE_MPI */

#ifdef USE_RTAI
#include <mbrtai_utils.h>
void *mbdyn_rtai_task = NULL;
#endif /* USE_RTAI */

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

#include <solver.h>

/* Note: DEBUG_* codes are declared in "mbdyn.h" */
#ifdef DEBUG
const debug_array da[] = {
    { "input",          MYDEBUG_INPUT           },
    { "sol",            MYDEBUG_SOL             },
    { "assembly",       MYDEBUG_ASSEMBLY        },
    { "derivatives",    MYDEBUG_DERIVATIVES     },
    { "fsteps",         MYDEBUG_FSTEPS          },
    { "mem",            MYDEBUG_MEM             },
    { "fnames",         MYDEBUG_FNAMES          },
#if defined(__GNUC__)
    { "prettyfn",       MYDEBUG_FNAMES|MYDEBUG_PRETTYFN },
#endif /* __GNUC__ */
    { "mpi",            MYDEBUG_MPI             },
    { "pred",           MYDEBUG_PRED            },
    { "residual",       MYDEBUG_RESIDUAL        },
    { "jacobian",       MYDEBUG_JAC             },
    { "init",           MYDEBUG_INIT            },
    { "output",         MYDEBUG_OUTPUT          },
     
    { NULL,             MYDEBUG_NONE            }
};
#endif /* DEBUG */

#ifdef HAVE_GETOPT
static void
mbdyn_usage(const char *sShortOpts)
{
    silent_cout(std::endl
	<< "MBDyn - MultiBody Dynamics " << VERSION << std::endl
	<< "compiled on " << __DATE__ << " at " << __TIME__ << std::endl 
	<< std::endl
	<< "mbdyn is a multibody simulation program." << std::endl
	<< std::endl
	<< "usage: mbdyn [" << sShortOpts << "] [input-file list] " << std::endl 
	<< std::endl
	<< "  -f, --input-file {file}  :"
	" reads from 'file' instead of stdin" << std::endl
	<< "  -o, --output-file {file} : writes to '{file}.xxx'" << std::endl
	<< "                             instead of '{input-file}.xxx'" << std::endl
	<< "                            "
	" (or 'Mbdyn.xxx' if input from stdin)" << std::endl
	<< "  -W, --working-dir {dir}  :"
	" sets the working directory" << std::endl
	<< "  -m, --mail {address}     :"
	" mails to {address} at job completion" << std::endl
	<< "  -n, --nice [level]       :"
	" change the execution priority of the process" << std::endl);
#ifdef DEBUG
    silent_cout("  -d, --debug {level[:level[...]]} :"
        " when using the debug version of the code," << std::endl
        << "                            "
	" enables debug levels; available:" << std::endl
        << "                                 none" << std::endl);
    for (int i = 0; da[i].s != NULL; i++) {
        silent_cout("                                 " << da[i].s << std::endl);
    }      
    silent_cout("                                 any" << std::endl);
#endif /* DEBUG */
    silent_cout("  -t, --same-table" << std::endl
        << "  -T, --no-same-table      :"
        " use/don't use same symbol table for multiple runs" << std::endl
        << "  -r, --redefine" << std::endl
        << "  -R, --no-redefine        :"
        " redefine/don't redefine symbols in table" << std::endl
	<< "  -H, --show-table         :"
	" print symbol table and exit" << std::endl
        << "  -s, --silent             :"
        " runs quietly" << std::endl
        << "  -P, --pedantic           :"
        " pedantic warning messages" << std::endl);
#ifdef USE_MPI
    silent_cout("  -p, --parallel           :"
        " required when run in parallel (invoked by mpirun)" << std::endl);
#endif /* USE_MPI */
    silent_cout("  -h, --help               :"
        " prints this message" << std::endl
        << "  -l, --license            :"
        " prints the licensing terms" << std::endl
        << "  -w, --warranty           :"
        " prints the warranty conditions" << std::endl
        << std::endl
        << "Usually mbdyn reads the input from stdin"
        " and writes messages on stdout; a log" << std::endl
        << "is put in '{file}.out', and data output"
        " is sent to various '{file}.xxx' files" << std::endl
        << "('Mbdyn.xxx' if input from stdin)" << std::endl
        << std::endl
	<< std::endl);
}

static void
mbdyn_welcome(void)
{
       	silent_cout(std::endl
		<< "MBDyn - Multi-Body Dynamics " << VERSION 
		<< std::endl
		<< "compiled on " << __DATE__
			<< " at " << __TIME__ << std::endl
		<< std::endl
		<< "Copyright 1997-2004 (C) Paolo Mantegazza "
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
static char sShortOpts[] = "a:d:f:hHlm:n:N::o:pPrRsS:tTwW:";

#ifdef HAVE_GETOPT_LONG
static struct option LongOpts[] = {
	{ "adams-file",     required_argument, NULL,           int('a') },
	{ "debug",          required_argument, NULL,           int('d') },
	{ "input-file",     required_argument, NULL,           int('f') },
	{ "help",           no_argument,       NULL,           int('h') },
	{ "show-table",     no_argument,       NULL,           int('H') },
	{ "license",        no_argument,       NULL,           int('l') },
	{ "mail",           required_argument, NULL,           int('m') },
	{ "nice",           optional_argument, NULL,           int('n') },
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
	{ "warranty",       no_argument,       NULL,           int('w') },
	{ "working-dir",    required_argument, NULL,           int('W') },
	
	{ NULL,             0,                 NULL,           0        }
};
#endif /* HAVE_GETOPT_LONG */
#endif /* HAVE_GETOPT */


extern void GetEnviron(MathParser&);
const char* sDefaultInputFileName = "MBDyn";



Solver* RunMBDyn(MBDynParser&, const char* const, const char* const, bool);
#ifdef MBDYN_X_MAIL_MESSAGE
static void SendMessage(const char* const, const char* const, time_t, time_t);
#endif /* MBDYN_X_MAIL_MESSAGE */

#ifdef USE_MPI
static int
parse_args(int argc, char *argv[], bool &using_mpi,
		int &parallel_fSilent, int &parallel_fPedantic)
{
	for (int i = 1; i < argc; i++) {
		if (!using_mpi && strncmp(argv[i], "-p", 2) == 0) {
			using_mpi = true;
			continue;
		}
		
		/* intercept silence/pedantic flags */
		if (strncmp(argv[i], "-s", 2) == 0) {
			parallel_fSilent++;

			for (unsigned j = 2; argv[i][j] == 's'; j++) {
				parallel_fSilent++;
			}

		} else if (strncmp(argv[i], "-P", 2) == 0) {
			parallel_fPedantic++;

			for (unsigned j = 2; argv[i][j] == 'P'; j++) {
				parallel_fPedantic++;
			}
		}
	}

	return 0;
}
#endif /* USE_MPI */

int
main(int argc, char* argv[])
{
	int	rc = EXIT_SUCCESS;

	bool	using_mpi = false;

#ifdef USE_MPI
	int	WorldSize = 1;
	int	MyRank = 0;
	char	ProcessorName_[1024] = "localhost",
		*ProcessorName = ProcessorName_;
	int	parallel_fSilent = 0,
		parallel_fPedantic = 0;

    	/* 
	 * FIXME: this is a hack to detect whether mbdyn has been
	 * invoked thru mpirun (which means we need to initialize
	 * MPI) or not (which means we don't); need to check how 
	 * portable it is ...
	 *
	 * the check is on the first two chars because "most" of
	 * the mpirun/MPI extra args start with -p<something>
	 */
	parse_args(argc, argv, using_mpi, parallel_fSilent, parallel_fPedantic);

	::fSilent = parallel_fSilent;
	::fPedantic = parallel_fPedantic;

	if (using_mpi) {
		MPI::Init(argc, argv);	   
		WorldSize = MPI::COMM_WORLD.Get_size();
		MyRank = MPI::COMM_WORLD.Get_rank();

		if (MyRank > 0) {
			/*
			 * need a second take because MPI::Init()
			 * restores the inital args, so if Get_rank() > 0
			 * the eventual -s/-P flags have been restored
			 */
			parse_args(argc, argv, using_mpi,
					parallel_fSilent, parallel_fPedantic);

			::fSilent = parallel_fSilent;
			::fPedantic = parallel_fPedantic;
		}
		
		/*
		 * all these temporaries are to avoid complains from
		 * the compiler (MPI's API is really messy ):
		 */
		int ProcessorNameLength = sizeof(ProcessorName_);
		MPI::Get_processor_name(ProcessorName, ProcessorNameLength); 

		silent_cerr("using MPI (explicitly required by '-p*' switch)"
			<< std::endl);
	}
#endif /* USE_MPI */
   
    	/* primo argomento valido (potenziale nome di file di ingresso) */
    	int currarg = 0;
    	if (argc > 0) {
        	currarg = 1;
    	}
   
    	/* The program is a big try block */
    	try {
        	enum InputFormat {
	    		MBDYN,
	    		ADAMS,
	    		LASTFORMAT
        	} CurrInputFormat = MBDYN;
      
        	/* Stream di ingresso dati */
		std::istream* pIn = NULL;
		std::ifstream FileStreamIn;
        	char* sInputFileName = (char *)sDefaultInputFileName;
        	char* sOutputFileName = NULL;
      
        	enum InputSource {
	    		FILE_UNKNOWN,
	    		FILE_STDIN,
	    		FILE_OPT,
	    		FILE_ARGS	   
        	} CurrInputSource = FILE_UNKNOWN;
      
        	/* Gestione dei parametri da linea di comando */
        	int fRedefine = 0;
        	int fTable = 0;
      
        	/* Mostra la tabella dei simboli ed esce */
        	int fShowSymbolTable = 0;
      
#ifdef HAVE_GETOPT
        	/* Dati acquisibili da linea di comando */
#ifdef MBDYN_X_MAIL_MESSAGE
        	char* sMailToAddress = NULL;
#endif /* MBDYN_X_MAIL_MESSAGE */
      
        	int iIndexPtr = 0;

#ifdef HAVE_NICE
		int niceIncr = 0;
#endif /* HAVE_NICE */

		int iSleepTime = -1;

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
	    		case int('m'):
#ifdef MBDYN_X_MAIL_MESSAGE
	        		sMailToAddress = optarg;
#else /* ! MBDYN_X_MAIL_MESSAGE */
				silent_cerr("warning: option -m has been "
					"disabled because of potential "
					"vulnerabilities" << std::endl);
#endif /* ! MBDYN_X_MAIL_MESSAGE */
	        		break;

			case int('N'):
				break;

#ifdef HAVE_NICE
			case int('n'):
				if (optarg != 0) {
#ifdef HAVE_STRTOL
					char *eptr = NULL;
					niceIncr = strtol(optarg, &eptr, 10);
					if (eptr != NULL && eptr[0] != '\0') {
		    				silent_cerr("Unable to "
							"parse nice level <" 
							<< optarg << ">"
							<< std::endl);
		    				throw ErrGeneric();
					}
#else /* !HAVE_STRTOL */
					niceIncr = atoi(optarg);
#endif /* !HAVE_STRTOL */
				} else {
					niceIncr = 10;
				}
				break;
#endif /* HAVE_NICE */
	    
	    		case int('f'):
	        		CurrInputFormat = MBDYN;
	        		CurrInputSource = FILE_OPT;
	        		sInputFileName = optarg;
	        		FileStreamIn.open(sInputFileName);
	        		if (!FileStreamIn) {
		    			silent_cerr(std::endl 
		        			<< "Unable to open file \""
						<< sInputFileName << "\"");
#ifdef USE_MPI
					if (using_mpi) {
						silent_cerr(" on "
							<< ProcessorName);
					}
#endif /* USE_MPI */
					silent_cerr(";" << std::endl 
						<< "aborting ..."
						<< std::endl);
		    			throw ErrGeneric();
	        		}
	        		pIn = (std::istream*)&FileStreamIn;
	        		break;
	    
	    		case int('a'):
#ifdef USE_ADAMS_PP
	        		CurrInputFormat = ADAMS;
	        		CurrInputSource = FILE_OPT;
	        		sInputFileName = optarg;
	     
	        		silent_cerr("ADAMS input not implemented yet,"
		    			" cannot open file '"
					<< sInputFileName << "'" << std::endl);
				throw ErrGeneric();
	        		break;
#else /* !USE_ADAMS_PP */
	        		silent_cerr("Illegal option -a" << std::endl);
	        		throw ErrGeneric();
				break;
#endif /* !USE_ADAMS_PP */
	    
	    		case int('o'):
	        		sOutputFileName = optarg;
	        		break;

	    		case int('d'):
#ifdef DEBUG
	        		if (get_debug_options(optarg, da)) {
		    			silent_cerr("Unable to interpret debug"
						" option argument;"
						" using default" << std::endl);
		    			::debug_level = DEFAULT_DEBUG_LEVEL;
		    			/* throw ErrGeneric(); */
	        		}
#else /* !DEBUG */
	        		silent_cerr("Compile with '-DDEBUG'"
					" to use debug features" << std::endl);
#endif /* !DEBUG */
	        		break;
	       
	    		case int('t'):
	        		fTable = 1;
	        		break;
	    
	    		case int('T'):
	        		fTable = 0;
	        		break;	  

			case int('p'):
#ifdef USE_MPI
				ASSERT(using_mpi);
#else /* !USE_MPI */
				silent_cerr("switch '-p' is meaningless "
						"without MPI" << std::endl);
#endif /* !USE_MPI */
				break;
	    
	    		case int('P'):
#ifdef USE_MPI
				if (parallel_fPedantic> 0) {
					parallel_fPedantic--;
				} else 
#endif /* USE_MPI */
				{
	        			::fPedantic++;
				}
	        		break;

	    		case int('r'):
	        		fRedefine = 1;
	        		break;
	    
	    		case int('R'):
	        		fRedefine = 0;
	        		break;
	    
	    		case int('s'):
#ifdef USE_MPI
				if (parallel_fSilent > 0) {
					parallel_fSilent--;
				} else 
#endif /* USE_MPI */
				{
	        			::fSilent++;
				}
	        		break;

			case int('S'):
				if (optarg) {
					char	*s = optarg;

					if (strncasecmp(s, "rank=", sizeof("rank=") - 1) == 0) {
#ifdef USE_MPI
						char	*next;
						long	r;
						
						s += sizeof("rank=") - 1;
						r = strtol(s, &next, 10);
						if (next[0] != '\0') {
							if (next[0] != ',') {
								silent_cerr("error in argument -S " << optarg << std::endl);
								throw ErrGeneric();
							}
							s = &next[1];
						}

						if (using_mpi && r != MyRank) {
							break;
						}
#else /* ! USE_MPI */
						silent_cerr("option -S " << optarg << " valid only when --with-mpi" << std::endl);
#endif /* ! USE_MPI */
					}

					if (s[0]) {
						char	*next;

						iSleepTime = strtol(s, &next, 10);
						if (next[0] != '\0') {
							silent_cerr("error in argument -S " << optarg << std::endl);
							throw ErrGeneric();
						}
					}

				} else {
					/* default: 10 s */
					iSleepTime = 10;
				}

				break;

	    		case int('l'):
				mbdyn_welcome();
	        		mbdyn_license();
	        		throw NoErr();
	    
	    		case int('w'):
				mbdyn_welcome();
				mbdyn_warranty();
	        		throw NoErr();

			case int('W'):
#ifdef HAVE_CHDIR
      				if (chdir(optarg)) {
					silent_cerr("Error in chdir(\""
						<< optarg << "\")"
						<< std::endl);
	 				throw ErrFileSystem();
      				}
#else /* !HAVE_CHDIR */
				silent_cerr("chdir() not available"
					<< std::endl);
#endif /* !HAVE_CHDIR */
				break;
	    
	    		case int('?'):
	        		silent_cerr("Unknown option -"
					<< char(optopt) << std::endl);

	    		case int('h'):
		    		mbdyn_usage(sShortOpts);
	        		throw NoErr();
	    
	    		case int('H'):
	        		fShowSymbolTable++;
	        		break;
	    
	    		default:
	        		silent_cerr(std::endl 
	            			<< "Unrecoverable error; aborting ..."
					<< std::endl);
	        		throw ErrGeneric();
	    		}
        	}

		if (iSleepTime > -1) {
			sleep(iSleepTime);
		}
		
        	/*
		 * primo argomento utile (potenziale nome di file di ingresso)
		 */
        	currarg = optind;
#endif /* HAVE_GETOPT */
		mbdyn_welcome();
#ifdef USE_MPI
		if (using_mpi) {
        		silent_cerr("Process " << MyRank 
	    			<< " (" << MyRank + 1 << " of " << WorldSize
            			<< ") is alive on " << ProcessorName
				<< std::endl);
		}
#endif /* USE_MPI */
      
      		/* Mostra la tabella dei simboli ed esce */
        	if (fShowSymbolTable > 0) {
#ifdef USE_MPI
	    		if (MyRank == 0) {
#endif /* USE_MPI */
	        		Table t(31, 1);
				MathParser mp(t);
				GetEnviron(mp);
	        		silent_cout("default symbol table:"
					<< std::endl << mp.GetSymbolTable() << std::endl);
#ifdef USE_MPI
	    		}
#endif /* USE_MPI */
	 
	    		throw NoErr();
        	}

#ifdef HAVE_NICE
		if (niceIncr != 0) {
			silent_cout("setting nice(" << niceIncr << ")" 
					<< std::endl);
			if (nice(niceIncr)) {
				silent_cerr("nice(" << niceIncr 
					<< ") failed; ignored" << std::endl);
			}
		}
#endif /* HAVE_NICE */

        	/* risolve l'input */
        	if (CurrInputSource == FILE_UNKNOWN) {
	    		if (currarg < argc) {
	        		CurrInputSource = FILE_ARGS;

	    		} else {
	        		/*
				 * se non e' un argomento prende
				 * lo standard input
				 */
	        		CurrInputSource = FILE_STDIN;
	        		CurrInputFormat = MBDYN;
	        		ASSERT(pIn == NULL);
	        		pIn = (std::istream*)&std::cin;
	    		}
        	}
      
        	/* Gestione di input/output */      
        	Table* pT = NULL;
        	MathParser* pMP = NULL;
      
        	int last = 0;
        	while (last == 0) {
	    		if (CurrInputSource == FILE_STDIN) {
				silent_cout("reading from stdin" << std::endl);
				last = 1;

			} else if (CurrInputSource == FILE_OPT) {
				silent_cout("reading from file \"" 
						<< sInputFileName 
						<< "\"" << std::endl);
	        		last = 1;

	    		} else if (CurrInputSource == FILE_ARGS) {
	        		sInputFileName = argv[currarg];
				silent_cout("reading from file \""
						<< sInputFileName
						<< "\"" << std::endl);
				
	        		/* incrementa il numero di argomento */
	        		if (++currarg == argc) {
					last = 1; 
				}
	    
#ifdef USE_ADAMS_PP
	        		/* ADAMS extension */
	        		char* p = strrchr(sInputFileName, int('.'));
	        		if (p != NULL 
		    		    && strlen(p+1) == 3 
				    && !strncasecmp(p+1, "adm", 3)) {
	            			CurrInputFormat = ADAMS;
	       
	            			silent_cerr("ADAMS input "
							"not implemented yet, "
							"cannot open file '"
							<< sInputFileName 
							<< "'" << std::endl);
	            			throw ErrGeneric();
	        		} else {
#endif /* USE_ADAMS_PP */
	            			CurrInputFormat = MBDYN;
	       
	            			FileStreamIn.open(sInputFileName);
	            			if (!FileStreamIn) {
		        			silent_cerr(std::endl 
			    				<< "Unable to open"
							" file \""
							<< sInputFileName 
			    				<< "\"; aborting ..."
							<< std::endl);
		        			throw ErrGeneric();
	            			}
#ifdef USE_ADAMS_PP
	        		}
#endif /* USE_ADAMS_PP */
	        		pIn = &FileStreamIn;
	    		}
	 	 
	    		Solver* pSolv = NULL;
	    		switch (CurrInputFormat) {
	    		case MBDYN: {	     
	        		if (pT == NULL) {
		    			SAFENEWWITHCONSTRUCTOR(pT,
							       Table,
							       Table(31, 1));
	        		}
	        		if (pMP == NULL) {
		    			SAFENEWWITHCONSTRUCTOR(pMP, 
		                           		       MathParser, 
					   			MathParser(*pT,
								    fRedefine));
		
		    			/* legge l'environment */
		    			GetEnviron(*pMP);
	        		} 
		
	        		/* stream in ingresso */
	        		InputStream In(*pIn);
	        		MBDynParser HP(*pMP, In, 
						sInputFileName == sDefaultInputFileName ? "initial file" : sInputFileName);
	     
	        		pSolv = RunMBDyn(HP, sInputFileName, 
						 sOutputFileName, using_mpi);
				if (FileStreamIn.is_open()) {
	        			FileStreamIn.close();
				}
	        		break;
	    		}
	       
	    		case ADAMS:
	        		silent_cerr("ADAMS input not implemented yet!"
					<< std::endl);
	        		throw ErrNotImplementedYet();
	    
	    		default:
	        		silent_cerr("You shouldn't be here!"
						<< std::endl);
	        		throw ErrGeneric();
	    		}
	    
			clock_t ct = 0;

	    		if (pSolv != NULL) {
				ct += pSolv->GetCPUTime();
	        		SAFEDELETE(pSolv);
	    		}
	 
	    		if (fTable == 0 || currarg == argc) {
	        		if (pMP != NULL) {
	            			SAFEDELETE(pMP);
	            			pMP = NULL;
	        		}
	        		if (pT != NULL) {
	            			SAFEDELETE(pT);
	            			pT = NULL;
	        		}
	    		}
	 
	    		time_t tSecs = 0;
	    		time_t tMils = 0;
#ifdef HAVE_SYS_TIMES_H	 
	    		/* Tempo di CPU impiegato */
	    		struct tms tmsbuf;
	    		times(&tmsbuf);

	    		ct += tmsbuf.tms_utime + tmsbuf.tms_cutime
				+ tmsbuf.tms_stime + tmsbuf.tms_cstime;

			long clk_tck = sysconf(_SC_CLK_TCK);
	    		tSecs = ct/clk_tck;
	    		tMils = ((ct*1000)/clk_tck)%1000;

			char timebuf[] = " (1000000000h 59m 59s 999ms)";
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
						throw ErrGeneric();
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
			if (using_mpi) {
	    			silent_cout(" on " << ProcessorName
					<< " (" << MyRank + 1 
					<< " of " << WorldSize << ")");
			}
#endif /* USE_MPI */
	    		silent_cout(std::endl);
#endif /* HAVE_SYS_TIMES_H */

#ifdef MBDYN_X_MAIL_MESSAGE
#ifdef HAVE_GETOPT
	    		/* E-mail all'utente */
	    		if (sMailToAddress) {
	        		SendMessage(sInputFileName, sMailToAddress,
					    tSecs, tMils);
				sMailToAddress = NULL;
	    		}
#endif /* HAVE_GETOPT */
#endif /* MBDYN_X_MAIL_MESSAGE */
        	}

        	throw NoErr();

    	} catch (NoErr) {     
        	silent_cout("MBDyn terminated normally" << std::endl);
        	rc = EXIT_SUCCESS;

	} catch (ErrInterrupted) {
        	silent_cout("MBDyn was interrupted" << std::endl);
        	rc = 2;

    	} catch (...) {
        	silent_cerr("An error occurred during the execution of MBDyn;"
	    		" aborting ... " << std::endl);
        	rc = EXIT_FAILURE;
    	}

#ifdef USE_RTAI
	if (mbdyn_rtai_task) {
		(void)mbdyn_rt_task_delete(&mbdyn_rtai_task);
		mbdyn_rtai_task = NULL;
	}
#endif /* USE_RTAI */
   
    	MB_EXIT(rc);
}


Solver* 
RunMBDyn(MBDynParser& HP, 
	 const char* const sInputFileName,
	 const char* const sOutputFileName,
	 bool using_mpi)
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
        	"integrator",
        	"multistep",
        	"rungekutta",
        	"parallel",
        	"schur",
		NULL
    	};

    	/* enum delle parole chiave */
    	enum KeyWords {
        	UNKNOWN = -1,
		
        	BEGIN = 0,
		END,
        	DATA,
        	INTEGRATOR,
        	MULTISTEP,
        	RUNGEKUTTA,
        	PARALLEL,
        	SSCHUR,
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
	} catch (EndOfFile) {
		throw NoErr();
	}
	/* looking for "begin"... */	
	if (cd != BEGIN) {
        	silent_cerr(std::endl 
	    		<< "Error: <begin> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << std::endl);
        	throw ErrGeneric();
    	}

	/* looking for "data"... */	
    	if (KeyWords(HP.GetWord()) != DATA) {
        	silent_cerr(std::endl 
	    		<< "Error: <begin: data;> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << std::endl);
        	throw ErrGeneric();
    	}
   
    	KeyWords CurrInt = MULTISTEP;
   
    	/* Ciclo infinito */
    	while (true) {	
        	switch (KeyWords(HP.GetDescription())) {
        	case INTEGRATOR:
            		switch (KeyWords(HP.GetWord())) {
            		case RUNGEKUTTA:
	        		CurrInt = RUNGEKUTTA;
	        		break;
		
            		case MULTISTEP:
	        		CurrInt = MULTISTEP;
	        		break;
		
            		case SSCHUR:
				silent_cerr("warning: \"schur\" solver "
					"is deprecated;" << std::endl);
#ifdef USE_MPI
				silent_cerr("use \"parallel\" with "
					"\"multistep\" solver instead"
					<< std::endl);
	        		CurrInt = MULTISTEP;
				using_mpi = true;
#else /* !USE_MPI */
				silent_cerr("compile with -DUSE_MPI "
					"to enable parallel solution" 
					<< std::endl);
				throw ErrGeneric();
#endif /* !USE_MPI */
	        		break;
		
            		default:
	        		silent_cerr(std::endl 
		    			<< "Unknown integrator at line " 
	            			<< HP.GetLineData()
					<< "; aborting ..." << std::endl);
	        		throw ErrGeneric();
            		}
            		break;    

        	case PARALLEL:
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
	    		throw ErrGeneric();
#endif /* !USE_MPI */

        	case END:
	    		if (KeyWords(HP.GetWord()) != DATA) {
	        		silent_cerr(std::endl 
		    			<< "Error: <end: data;> expected"
					" at line " << HP.GetLineData()
					<< "; aborting ..." << std::endl);
	        		throw ErrGeneric();
	    		}
	    		goto endofcycle;        
	 
        	default:
	    		silent_cerr(std::endl 
	        		<< "Unknown description at line " 
	        		<< HP.GetLineData()
				<< "; aborting ..." << std::endl);
	    		throw ErrGeneric();      
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
	}
#endif /* USE_MPI */

    	switch (CurrInt) {
    	case MULTISTEP: 
        	SAFENEWWITHCONSTRUCTOR(pSolv,
				Solver,
				Solver(HP, sInputFileName, 
					sOutputFileName, bParallel));
        	break;
 
    	case RUNGEKUTTA:
		/* FIXME: this is rather obsolete, since our 
		 * integration scheme incorporates implicit 
		 * Runge-Kutta and multistep/single-lag methods
		 * in one scheme; the "thirdorder" method
		 * is actually an implicit Runge-Kutta with
		 * tunable algorithmic dissipation */
        	silent_cerr("Sorry, implicit Runge-Kutta isn't supported yet;"
	    		<< std::endl << "aborting ..." << std::endl);
        	throw ErrNotImplementedYet();

	default:
        	silent_cerr("Unknown integrator; aborting ..." << std::endl);
        	throw ErrGeneric();   
    	}

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
    
    	return pSolv;
}

#ifdef MBDYN_X_MAIL_MESSAGE
/*
 * Plenty of potential vulnerabilities; never enable
 */
static void 
SendMessage(const char* const sInputFileName,
	    const char* const sMailToAddress,
	    time_t tSecs,
	    time_t tMils)
{
    	DEBUGCOUTFNAME("SendMessage");
   
    	/* Scrive il messaggio in un file temporaneo */
	std::ofstream Msg("mbdyn.msg");
    	Msg << "MBDyn terminated job ";
    	if (sInputFileName != NULL) {
        	Msg << "'" << sInputFileName << "' ";
    	}
#ifdef HAVE_SYS_TIMES_H
    	Msg << "in " << tSecs << '.' << tMils << " seconds of CPU time";
#endif /* HAVE_SYS_TIMES_H */
    	Msg << '.' << std::endl;
    	Msg.close();
   
    	/* Crea la linea di comando */
    	char* sCmd = NULL;
    	SAFENEWARR(sCmd, char, (29+strlen(sMailToAddress)+11+1));
    	char* s = sCmd;
    	strcpy(s, "/bin/mail -s 'mbdyn message' ");
    	s = sCmd+strlen(sCmd);
    	strcpy(s, sMailToAddress);
    	s = sCmd+strlen(sCmd);
    	strcpy(s, " <mbdyn.msg");
    	system(sCmd);
   
    	/* Manda il messagio e cancella il file temporaneo */
    	system("rm mbdyn.msg");
    	SAFEDELETEARR(sCmd);
}
#endif /* MBDYN_X_MAIL_MESSAGE */

