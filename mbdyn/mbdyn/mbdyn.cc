/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
#include <mpi++.h>

MPI::Intracomm MBDynComm = MPI::COMM_SELF;

#define MB_EXIT(err) \
	do { \
		if (using_mpi) { \
			MPI::Finalize(); \
		} \
		exit((err)); \
	} while (0)

#else /* !USE_MPI */

#define MB_EXIT(err) \
	exit((err))

#endif /* !USE_MPI */

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

#include <multistp.h>

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
mbdyn_usage(std::ostream& out, const char *sShortOpts)
{
    std::cout 
        << std::endl
	<< "MBDyn - Multi-Body Dynamics " << VERSION << std::endl
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
	" change the execution priority of the process" << std::endl;
#ifdef DEBUG
    out
        << "  -d, --debug {level[:level[...]]} :"
        " when using the debug version of the code," << std::endl
        << "                            "
	" enables debug levels; available:" << std::endl
        << "                                 none" << std::endl;
    for (int i = 0; da[i].s != NULL; i++) {
        out << "                                 " << da[i].s << std::endl;
    }      
    out 
        << "                                 any" << std::endl;
#endif /* DEBUG */
    out    
        << "  -t, --same-table" << std::endl
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
        " pedantic warning messages" << std::endl
#ifdef USE_MPI
        << "  -p, --parallel           :"
        " required when run in parallel (invoked by mpirun)" << std::endl
#endif /* USE_MPI */
        << "  -h, --help               :"
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
        << std::endl;
    std::cout << std::endl;   
}

/* Dati di getopt */
static char sShortOpts[] = "a:d:f:hHlm:n::o:pPrRstTwW:";
enum MyOptions {
	MAIL = 0,
	INPUT_FILE,
	ADAMS_FILE,
	OUTPUT_FILE,
	DEBUG_LEVEL,
	ONE_TABLE,
	MULTIPLE_TABLE,
	SILENT,
	HELP,
	
	LASTOPTION
} /* MyOptions */ ;

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
	{ "output-file",    required_argument, NULL,           int('o') },
	{ "parallel",	    no_argument,       NULL,           int('p') },
	{ "pedantic",	    no_argument,       NULL,           int('P') },
	{ "redefine",       no_argument,       NULL,           int('r') },
	{ "no-redefine",    no_argument,       NULL,           int('R') },
	{ "silent",         no_argument,       NULL,           int('s') },
	{ "same-table",     no_argument,       NULL,           int('t') },
	{ "no-same-table",  no_argument,       NULL,           int('T') },
	{ "warranty",       no_argument,       NULL,           int('w') },
	{ "working-dir",    required_argument, NULL,           int('W') },
	
	{ NULL,             0,                 NULL,           0        }
};
#endif /* HAVE_GETOPT_LONG */
#endif /* HAVE_GETOPT */


extern void GetEnviron(MathParser&);

/* flag di silent run (no output su stdout) */
int fSilent = 0;
int fPedantic = 0;
const char* sDefaultInputFileName = "MBDyn";



Integrator* RunMBDyn(MBDynParser&, const char* const, const char* const);
#ifdef MBDYN_X_MAIL_MESSAGE
static void SendMessage(const char* const, const char* const, time_t, time_t);
#endif /* MBDYN_X_MAIL_MESSAGE */


int
main(int argc, char* argv[])
{
#ifdef USE_MPI
	int using_mpi = 0;
	int WorldSize = 1;
	int myrank = 0;
	char ProcessorName_[1024] = "localhost", *ProcessorName = ProcessorName_;

    	/* 
	 * FIXME: this is a hack to detect whether mbdyn has been
	 * invoked thru mpirun (which means we need to initialize
	 * MPI) or not (which means we don't); need to check how 
	 * portable it is ...
	 */
	for (char **s = argv; s[0]; s++) {
		if (strncmp(s[0], "-p", 2) == 0) {
			MPI::Init(argc, argv);	   
			WorldSize = MPI::COMM_WORLD.Get_size();
			myrank = MPI::COMM_WORLD.Get_rank();
		
			/*
			 * all these temporaries are to avoid complains from
			 * the compiler (MPI's API is really messy ):
			 */
			int ProcessorNameLength = sizeof(ProcessorName_);
			MPI::Get_processor_name(ProcessorName, 
					ProcessorNameLength); 
			using_mpi = 1;
			break;
		}
	}

	if ( using_mpi == 1 ) {
		std::cerr << "using MPI (required by '-p' switch)" << std::endl;
	}
#endif /* USE_MPI */
   
    	/* primo argomento valido (potenziale nome di file di ingresso) */
    	int currarg = 0;
    	if (argc > 0) {
        	currarg = 1;
    	}
   
    	/* The program is a big try block */
#ifdef USE_EXCEPTIONS
    	try {
#endif /* USE_EXCEPTIONS */
      
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

        	/* Parsing della linea di comando */
        	opterr = 0;
        	while (1) {
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
#else /* !MBDYN_X_MAIL_MESSAGE */
				std::cerr << "warning: option -m has been "
					"disabled because of potential "
					"vulnerabilities" << std::endl;
#endif /* MBDYN_X_MAIL_MESSAGE */
	        		break;

#ifdef HAVE_NICE
			case int('n'):
				if (optarg != 0) {
#ifdef HAVE_STRTOL
					char *eptr = NULL;
					niceIncr = strtol(optarg, &eptr, 10);
					if (eptr != NULL && eptr[0] != '\0') {
		    				std::cerr << "Unable to "
							"parse nice level <" 
							<< optarg << ">"
							<< std::endl;
		    				THROW(ErrGeneric());
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
		    			std::cerr 
		        			<< std::endl 
		        			<< "Unable to open file '"
						<< sInputFileName;
#ifdef USE_MPI
					if (using_mpi) {
						std::cerr << "' on " << ProcessorName;
					}
#endif /* USE_MPI */
					std::cerr
						<< ";" << std::endl 
						<< "aborting ..." << std::endl;
		    			THROW(ErrGeneric());
	        		}
	        		pIn = (std::istream*)&FileStreamIn;
	        		break;
	    
	    		case int('a'):
#ifdef USE_ADAMS_PP
	        		CurrInputFormat = ADAMS;
	        		CurrInputSource = FILE_OPT;
	        		sInputFileName = optarg;
	     
	        		std::cerr << "ADAMS input not implemented yet,"
		    			" cannot open file '"
					<< sInputFileName << "'" << std::endl;
				THROW(ErrGeneric());
	        		break;
#else /* !USE_ADAMS_PP */
	        		std::cerr << "Illegal option -a" << std::endl;
	        		THROW(ErrGeneric());
				break;
#endif /* !USE_ADAMS_PP */
	    
	    		case int('o'):
	        		sOutputFileName = optarg;
	        		break;

	    		case int('d'):
#ifdef DEBUG
	        		if (get_debug_options(optarg, da)) {
		    			std::cerr << "Unable to interpret debug"
						" option argument;"
						" using default" << std::endl;
		    			::debug_level = DEFAULT_DEBUG_LEVEL;
		    			/* THROW(ErrGeneric()); */
	        		}
#else /* !DEBUG */
	        		std::cerr << "Compile with '-DDEBUG'"
					" to use debug features" << std::endl;
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
				ASSERT(using_mpi == 1);
#else /* !USE_MPI */
				std::cerr << "switch '-p' is meaningless without MPI" << std::endl;
#endif /* !USE_MPI */
				break;
	    
	    		case int('P'):
	        		::fPedantic++;
	        		break;

	    		case int('r'):
	        		fRedefine = 1;
	        		break;
	    
	    		case int('R'):
	        		fRedefine = 0;
	        		break;
	    
	    		case int('s'):
	        		::fSilent++;
	        		break;

	    		case int('l'):
	        		mbdyn_license(std::cout);
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
	        		MB_EXIT(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
	    
	    		case int('w'):
				mbdyn_warranty(std::cout);
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
	        		MB_EXIT(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */

			case int('W'):
#ifdef HAVE_CHDIR
      				if (chdir(optarg)) {
					std::cerr << "Error in chdir(\""
						<< optarg << "\")"
						<< std::endl;
	 				THROW(ErrFileSystem());
      				}
#else /* !HAVE_CHDIR */
				std::cerr << "chdir() not available"
					<< std::endl;
#endif /* !HAVE_CHDIR */
				break;
	    
	    		case int('?'):
	        		std::cerr << "Unknown option -"
					<< char(optopt) << std::endl;

	    		case int('h'):
#ifdef USE_MPI
	        		if (myrank == 0) {
#endif /* USE_MPI */
		    			mbdyn_usage(std::cout, sShortOpts);
#ifdef USE_MPI
	        		}
#endif /* USE_MPI */
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
	        		MB_EXIT(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
	    
	    		case int('H'):
	        		fShowSymbolTable++;
	        		break;
	    
	    		default:
	        		std::cerr << std::endl 
	            			<< "Unrecoverable error; aborting ..."
					<< std::endl;
	        		THROW(ErrGeneric());
	    		}
        	}
		
        	/*
		 * primo argomento utile (potenziale nome di file di ingresso)
		 */
        	currarg = optind;
#endif /* HAVE_GETOPT */

        	silent_cout(std::endl
                    	    << "MBDyn - Multi-Body Dynamics " << VERSION 
		    	    << std::endl
		    	    << "compiled on " << __DATE__ << " at " << __TIME__ 
		    	    << std::endl
			    << std::endl
		    	    << "Copyright 1997-2002 (C) Paolo Mantegazza and Pierangelo Masarati,"
			    << std::endl
		    	    << "Dipartimento di Ingegneria Aerospaziale,"
		    	    " Politecnico di Milano." << std::endl
 			    << std::endl
		    	    << "MBDyn is free software, covered by the"
		    	    " GNU General Public License, and you are" << std::endl
		    	    << "welcome to change it and/or distribute"
		    	    " copies of it under certain conditions." << std::endl
		    	    << "Use 'mbdyn --license' to see the conditions."
			    << std::endl
		    	    << "There is absolutely no warranty for"
		    	    " MBDyn.  Use \"mbdyn --warranty\" for details."
			    << std::endl 
		    	    << std::endl);
#ifdef USE_MPI
		if (using_mpi) {
        		std::cerr << "Process " << myrank 
	    			<< " (" << myrank+1 << " of " << WorldSize
            			<< ") is alive on " << ProcessorName << std::endl;
		}
#endif /* USE_MPI */
      
      		/* Mostra la tabella dei simboli ed esce */
        	if (fShowSymbolTable > 0) {
#ifdef USE_MPI
	    		if (myrank == 0) {
#endif /* USE_MPI */
	        		Table t(31, 1);
	        		std::cout << "default symbol table:"
					<< std::endl << t << std::endl;
#ifdef USE_MPI
	    		}
#endif /* USE_MPI */
	 
#ifdef USE_EXCEPTIONS
	    		throw NoErr();
#else /* !USE_EXCEPTIONS */
	    		MB_EXIT(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
        	}

#ifdef HAVE_NICE
		if (niceIncr != 0) {
			silent_cout("setting nice(" << niceIncr << ")" 
					<< std::endl);
			if (nice(niceIncr)) {
				std::cerr << "nice(" << niceIncr 
					<< ") failed; ignored" << std::endl;
			}
		}
#endif /* HAVE_NICE */

        	/* risolve l'input */
        	if (CurrInputSource == FILE_UNKNOWN) {
	    		if (argv[currarg] != NULL) {
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
				silent_cout("reading from file '" 
						<< sInputFileName 
						<< "'" << std::endl);
	        		last = 1;

	    		} else if (CurrInputSource == FILE_ARGS) {
	        		sInputFileName = argv[currarg];
				silent_cout("reading from file '"
						<< sInputFileName
						<< "'" << std::endl);
				
	        		/* incrementa il numero di argomento */
	        		currarg++; 
	        		if (argv[currarg] == NULL) {
	            			last = 1;
	        		}
	    
#ifdef USE_ADAMS_PP
	        		/* ADAMS extension */
	        		char* p = strrchr(sInputFileName, int('.'));
	        		if (p != NULL 
		    		    && strlen(p+1) == 3 
				    && !strncasecmp(p+1, "adm", 3)) {
	            			CurrInputFormat = ADAMS;
	       
	            			std::cerr << "ADAMS input not implemented"
						" yet, cannot open file '"
						<< sInputFileName << "'"
						<< std::endl;
	            			THROW(ErrGeneric());
	        		} else {
#endif /* USE_ADAMS_PP */
	            			CurrInputFormat = MBDYN;
	       
	            			FileStreamIn.open(sInputFileName);
	            			if (!FileStreamIn) {
		        			std::cerr << std::endl 
			    				<< "Unable to open"
							" file '"
							<< sInputFileName 
			    				<< "'; aborting ..."
							<< std::endl;
		        			THROW(ErrGeneric());
	            			}
#ifdef USE_ADAMS_PP
	        		}
#endif /* USE_ADAMS_PP */
	        		pIn = &FileStreamIn;
	    		}
	 	 
	    		Integrator* pIntg = NULL;
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
		
	        		/* parser del blocco di controllo */
	        		KeyTable K(0, NULL);
	     
	        		/* stream in ingresso */
	        		InputStream In(*pIn);
	        		MBDynParser HP(*pMP, K, In, 
						sInputFileName == sDefaultInputFileName ? "initial file" : sInputFileName);
	     
	        		pIntg = RunMBDyn(HP, sInputFileName, 
						 sOutputFileName);
				if (FileStreamIn.is_open()) {
	        			FileStreamIn.close();
				}
	        		break;
	    		}
	       
	    		case ADAMS:
	        		std::cerr << "ADAMS input not implemented yet!"
					<< std::endl;
	        		THROW(ErrNotImplementedYet());
	    
	    		default:
	        		std::cerr << "You shouldn't be here!" << std::endl;
	        		THROW(ErrGeneric());
	    		}
	    
	    		if (pIntg != NULL) {
	        		SAFEDELETE(pIntg);
	    		}
	 
	    		if (fTable == 0 || argv[currarg] == NULL) {
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
	    		struct tms buf;
	    		times(&buf);
	    		clock_t ct = buf.tms_utime+buf.tms_cutime
				+buf.tms_stime+buf.tms_cstime;
			long clk_tck = sysconf(_SC_CLK_TCK);
	    		tSecs = ct/clk_tck;
	    		tMils = ((ct*1000)/clk_tck)%1000;
	    		std::cout << std::endl << "The simulation required " 
	        		<< tSecs << '.' << tMils 
	        		<< " seconds of CPU time";
			if (tSecs > 60) {
				std::cout << " (";
				if (tSecs > 3600) {
					std::cout << tSecs/3600 << "h ";
				}
				std::cout << (tSecs%3600)/60 << "m "
					<< (tSecs%3600)%60 << "s "
					<< tMils << "ms)";
			}
#ifdef USE_MPI
			if (using_mpi) {
	    			std::cout << " on " << ProcessorName;
			}
#endif /* USE_MPI */
	    		std::cout << std::endl;
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

#ifdef USE_EXCEPTIONS
        	throw NoErr();
    	} catch (NoErr) {     
        	silent_cout("MBDyn terminated normally" << std::endl);
        	MB_EXIT(EXIT_SUCCESS);
    	} catch (...) {
        	std::cerr << "An error occurred during the execution of MBDyn;"
	    		" aborting ... " << std::endl;
        	MB_EXIT(EXIT_FAILURE);
    	}
#endif /* USE_EXCEPTIONS */
   
    	MB_EXIT(EXIT_SUCCESS);
}


Integrator* 
RunMBDyn(MBDynParser& HP, 
	 const char* const sInputFileName,
	 const char* const sOutputFileName)
{
    	DEBUGCOUTFNAME("RunMBDyn");
   
    	Integrator* pIntg = NULL;

#ifdef USE_MPI
    	/* flag di parallelo */
    	flag fParallel(0);
#endif /* USE_MPI */
	flag fIterative(0);
	
    	/* parole chiave */
    	const char* sKeyWords[] = { 
        	"begin",
		"end",
        	"data",
        	"integrator",
        	"multistep",
        	"rungekutta",
        	"parallel",
		"iterative",
        	"schur"
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
		ITERATIVE,
        	SSCHUR,
        	LASTKEYWORD
    	};
   
    	/* tabella delle parole chiave */
    	KeyTable K((int)LASTKEYWORD, sKeyWords);
   
    	/* Attacca la tabella al parser */
    	HP.PutKeyTable(K);
   
    	/* legge i dati della simulazione */
    	if (KeyWords(HP.GetDescription()) != BEGIN) {
        	std::cerr << std::endl 
	    		<< "Error: <begin> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << std::endl;
        	THROW(ErrGeneric());
    	}
   
    	if (KeyWords(HP.GetWord()) != DATA) {
        	std::cerr << std::endl 
	    		<< "Error: <begin: data;> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << std::endl;
        	THROW(ErrGeneric());
    	}
   
    	KeyWords CurrInt = MULTISTEP;
   
    	/* Ciclo infinito */
    	while (1) {	
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
				std::cerr << "warning: \"schur\" solver "
					"is deprecated;" << std::endl;
#ifdef USE_MPI
				std::cerr << "use \"parallel\" with "
					"\"multistep\" solver instead"
					<< std::endl;
	        		CurrInt = MULTISTEP;
				fParallel = 1;
#else /* !USE_MPI */
				std::cerr << "compile with -DUSE_MPI "
					"to enable parallel solution" 
					<< std::endl;
				THROW(ErrGeneric());
#endif /* !USE_MPI */
	        		break;
		
            		default:
	        		std::cerr << std::endl 
		    			<< "Unknown integrator at line " 
	            			<< HP.GetLineData()
					<< "; aborting ..." << std::endl;
	        		THROW(ErrGeneric());
            		}
            		break;    

        	case PARALLEL:
#ifdef USE_MPI
			if (HP.fIsArg()) {
				int size = MPI::COMM_WORLD.Get_size();
				int NMbProc = HP.GetInt(); /* numero di porcessi assegnati ad mbdyn */
				if (NMbProc > size) {
					std::cerr << "WARNING: the requested number of processors is" 
						<< "higher than the comm_world size.\n  The new processor number is "  
						<< size << std::endl;
					NMbProc = size;
				}
				if (NMbProc == size) {
					/* c'e' solo MBDyn */
					MBDynComm = MPI::COMM_WORLD.Dup();
				} else {
					/* creo il gruppo che va da 1 a NMbProc */
					MPI::Group WorldGroup = MPI::COMM_WORLD.Get_group();
					int ran[3] = {1, NMbProc, 1};
					MPI::Group MBDynGroup = WorldGroup.Range_incl(1, &ran);
					MBDynComm = MPI::COMM_WORLD.Create(MBDynGroup);
					MBDynGroup.Free();
					WorldGroup.Free();
				}
			} else {
			 	MBDynComm = MPI::COMM_WORLD.Dup();
			}
            		fParallel = 1;
	    		break;
#else /* !USE_MPI */
            		std::cerr << "complile with -DUSE_MPI to enable "
				"parallel solution" << std::endl;
	    		THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
	    		break;
#endif /* USE_EXCEPTIONS */
#endif /* !USE_MPI */
		case ITERATIVE:
			fIterative = 1;
			break;

        	case END:
	    		if (KeyWords(HP.GetWord()) != DATA) {
	        		std::cerr << std::endl 
		    			<< "Error: <end: data;> expected"
					" at line " << HP.GetLineData()
					<< "; aborting ..." << std::endl;
	        		THROW(ErrGeneric());
	    		}
	    		goto endofcycle;        
	 
        	default:
	    		std::cerr << std::endl 
	        		<< "Unknown description at line " 
	        		<< HP.GetLineData()
				<< "; aborting ..." << std::endl;
	    		THROW(ErrGeneric());      
        	}
    	}
   
   	/* Uscita dal ciclo infinito */
endofcycle:   
   
    	switch (CurrInt) {
    	case MULTISTEP: 
#ifdef USE_MPI
		if (fParallel) {
        		SAFENEWWITHCONSTRUCTOR(pIntg,
					MultiStepIntegrator,
					MultiStepIntegrator(HP, sInputFileName, 
						sOutputFileName, fParallel));
		} else {
#endif /* USE_MPI */
			if (fIterative){
				SAFENEWWITHCONSTRUCTOR(pIntg,
					MultiStepIntegrator,
					MultiStepIntegrator(HP, sInputFileName,
						sOutputFileName, 0, fIterative));
			} else {
				SAFENEWWITHCONSTRUCTOR(pIntg,
					MultiStepIntegrator,
					MultiStepIntegrator(HP, sInputFileName,
						sOutputFileName));
			}
#ifdef USE_MPI
		}
#endif /* USE_MPI */
        	break;
      
    	case RUNGEKUTTA:
        	std::cerr << "Sorry, implicit Runge-Kutta isn't supported yet;"
	    		<< std::endl << "aborting ..." << std::endl;
        	THROW(ErrNotImplementedYet());

	default:
        	std::cerr << "Unknown integrator; aborting ..." << std::endl;
        	THROW(ErrGeneric());   
    	}
   
    	/* Runs the simulation */
    	pIntg->Run();
    
    	return pIntg;
}

#ifdef MBDYN_X_MAIL_MESSAGE
/*
 * Plenty of potential vulnerabilities
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
    	Msg << "Mbdyn terminated job ";
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

