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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

#ifdef DEBUG_MEMMANAGER
clMemMan MBDynMM("mbdyn");
#endif /* DEBUG_MEMMANAGER */

extern "C" {
#ifdef HAVE_TIMES_H
#include <sys/times.h>
#endif /* HAVE_TIMES_H */

/* Per il parsing della linea di comando */
#include <unistd.h>
#include <fstream.h>
#ifdef HAVE_GETOPT
#include <getopt.h>
#endif /* HAVE_GETOPT */ 
}

#ifdef USE_MPI 
#include <mpi++.h>
#include <schur.h>
#endif /* USE_MPI */

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
    { "init",           MYDEBUG_INIT            },
    { "output",         MYDEBUG_OUTPUT          },
     
    { NULL,             MYDEBUG_NONE            }
};
#endif /* DEBUG */

#ifdef HAVE_GETOPT
static void
mbdyn_usage( ostream& out, const char *sShortOpts )
{
    cout 
        << endl
	<< "MBDyn - Multi-Body Dynamics " << VERSION << endl
	<< "compiled on " << __DATE__ << " at " << __TIME__ << endl 
	<< endl
	<< "mbdyn is a multibody simulation program." << endl
	<< endl
	<< "usage: mbdyn [" << sShortOpts << "] [input-file list] " << endl 
	<< endl
	<< "  -f, --input-file {file}  :"
	" reads from 'file' instead of stdin" << endl
	<< "  -o, --output-file {file} : writes to '{file}.xxx'" << endl
	<< "                             instead of '{input-file}.xxx'" << endl
	<< "                            "
	" (or 'Mbdyn.xxx' if input from stdin)" << endl
	<< "  -m, --mail {address}     :"
	" mails to {address} at job completion" << endl;
#ifdef DEBUG
    out
        << "  -d, --debug {level[:level[...]]} :"
        " when using the debug version of the code," << endl
        << "                            "
	" enables debug levels; available:" << endl
        << "                                 none" << endl;
    for (int i = 0; da[i].s != NULL; i++) {
        out << "                                 " << da[i].s << endl;
    }      
    out 
        << "                                 any" << endl;
#endif /* DEBUG */
    out    
        << "  -t, --same-table" << endl
        << "  -T, --no-same-table      :"
        " use/don't use same symbol table for multiple runs" << endl
        << "  -r, --redefine" << endl
        << "  -R, --no-redefine        :"
        " redefine/don't redefine symbols in table" << endl
        << "  -s, --silent             :"
        " runs quietly" << endl
        << "  -h, --help               :"
        " prints this message" << endl
        << "  -l, --license            :"
        " prints the licensing terms" << endl
        << "  -w, --warranty           :"
        " prints the warranty conditions" << endl
        << endl
        << "Usually mbdyn reads the input from stdin"
        " and writes messages on stdout; a log" << endl
        << "is put in '{file}.out', and data output"
        " is put in various '{file}.xxx' files" << endl
        << "('Mbdyn.xxx' if input from stdin)" << endl
        << endl;
    cout << endl;   
}

/* Dati di getopt */
static char sShortOpts[] = "a:d:f:hHlm:o:rRstTw";
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
	{ "output-file",    required_argument, NULL,           int('o') },
	{ "redefine",       no_argument,       NULL,           int('r') },
	{ "no-redefine",    no_argument,       NULL,           int('R') },
	{ "silent",         no_argument,       NULL,           int('s') },
	{ "same-table",     no_argument,       NULL,           int('t') },
	{ "no-same-table",  no_argument,       NULL,           int('T') },
	{ "warranty",       no_argument,       NULL,           int('w') },
	
	{ NULL,             0,                 NULL,           0        }
};
#endif /* HAVE_GETOPT_LONG */
#endif /* HAVE_GETOPT */


extern void GetEnviron(MathParser&);

/* flag di silent run (no output su stdout) */
int fSilent = 0;


Integrator* RunMBDyn(MBDynParser&, const char* const, const char* const);
void SendMessage(const char* const, const char* const, time_t, time_t);


int
main(int argc, char* argv[])
{   
#ifdef USE_MPI
    	/* Inizializza i processi */
    	MPI::Init(argc, argv);	   
    	int WorldSize = MPI::COMM_WORLD.Get_size();
    	int myrank = MPI::COMM_WORLD.Get_rank();
    	int ProcessorNameLenght = 1024;
    	char* ProcessorName = NULL;
    	SAFENEWARR(ProcessorName, char, ProcessorNameLenght, MBDynMM);
    	MPI::Get_processor_name(ProcessorName, ProcessorNameLenght); 
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
        	istream* pIn = NULL;
        	ifstream FileStreamIn;
        	char* sInputFileName = NULL;
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
        	int iMailToBeSent = 0;
        	char* sMailToAddress = NULL;
      
        	int iIndexPtr = 0;

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
	        		iMailToBeSent = 1;
	        		SAFESTRDUP(sMailToAddress, optarg, MBDynMM);
	        		break;
	    
	    		case int('f'):
	        		CurrInputFormat = MBDYN;
	        		CurrInputSource = FILE_OPT;
				if (sInputFileName != NULL) {
					SAFEDELETEARR(sInputFileName, MBDynMM);
					sInputFileName = NULL;
				}
	        		SAFESTRDUP(sInputFileName, optarg, MBDynMM);
	        		FileStreamIn.open(sInputFileName);
	        		if (!FileStreamIn) {
		    			cerr 
		        			<< endl 
		        			<< "Unable to open file '"
						<< sInputFileName
#ifdef USE_MPI
		        			<< " on " << ProcessorName
#endif /* USE_MPI */
						<< ";" << endl 
						<< "aborting ..." << endl;
		    			THROW(ErrGeneric());
	        		}
	        		pIn = (istream*)&FileStreamIn;
	        		break;
	    
	    		case int('a'):
#ifdef USE_ADAMS_PP
	        		CurrInputFormat = ADAMS;
	        		CurrInputSource = FILE_OPT;
	        		SAFESTRDUP(sInputFileName, optarg, MBDynMM);
	     
	        		cerr << "ADAMS input not implemented yet,"
		    			" cannot open file '"
					<< sInputFileName << "'" << endl;
				THROW(ErrGeneric());
	        		break;
#else /* !USE_ADAMS_PP */
	        		cerr << "Illegal option -a" << endl;
	        		THROW(ErrGeneric());
				break;
#endif /* !USE_ADAMS_PP */
	    
	    		case int('o'):
	        		if (sOutputFileName != NULL) {
	            			SAFEDELETEARR(sOutputFileName, MBDynMM);
	            			sOutputFileName = NULL;
	        		}
	        		SAFESTRDUP(sOutputFileName, optarg, MBDynMM);
	        		break;

	    		case int('d'):
#ifdef DEBUG
	        		if (get_debug_options(optarg, da)) {
		    			cerr << "Unable to interpret debug"
						" option argument;"
						" using default" << endl;
		    			::debug_level = DEFAULT_DEBUG_LEVEL;
		    			/* THROW(ErrGeneric()); */
	        		}
#else /* !DEBUG */
	        		cerr << "Compile with '-DDEBUG'"
					" to use debug features" << endl;
#endif /* !DEBUG */
	        		break;
	       
	    		case int('t'):
	        		fTable = 1;
	        		break;
	    
	    		case int('T'):
	        		fTable = 0;
	        		break;	  
	    
	    		case int('r'):
	        		fRedefine = 1;
	        		break;
	    
	    		case int('R'):
	        		fRedefine = 0;
	        		break;
	    
	    		case int('s'):
	        		::fSilent = 1;
	        		break;

	    		case int('l'):
	        		cout << "license not available yet;"
					" see GPL" << endl;
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
#ifdef USE_MPI
	        		MPI::Finalize();
#endif /* USE_MPI */
	        		exit(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
	    
	    		case int('w'):
	        		cout << "warranty not available yet;"
	            			" see warranty coming with GPL"
					<< endl;
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
#ifdef USE_MPI
	        		MPI::Finalize();
#endif /* USE_MPI */
	        		exit(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
	    
	    		case int('?'):
	        		cerr << "Unknown option -"
					<< char(optopt) << endl;
	    		case int('h'):
	     
#ifdef USE_MPI
	        		if (myrank == 0) {
#endif /* USE_MPI */
		    			mbdyn_usage( cout, sShortOpts );
#ifdef USE_MPI
	        		}
#endif /* USE_MPI */
#ifdef USE_EXCEPTIONS
	        		throw NoErr();
#else /* !USE_EXCEPTIONS */
#ifdef USE_MPI
	        		MPI::Finalize();
#endif /* USE_MPI */
	        		exit(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
	    
	    		case int('H'):
	        		fShowSymbolTable++;
	        		break;
	    
	    		default:
	        		cerr << endl 
	            			<< "Unrecoverable error; aborting ..."
					<< endl;
	        		THROW(ErrGeneric());
	    		}
        	}
      
        	/*
		 * primo argomento utile (potenziale nome di file di ingresso)
		 */
        	currarg = optind;
#endif /* HAVE_GETOPT */

        	silent_cout(endl
                    	    << "MBDyn - Multi-Body Dynamics " << VERSION 
		    	    << endl
		    	    << "compiled on " << __DATE__ << " at " << __TIME__ 
		    	    << endl
			    << endl
		    	    << "Copyright 1997-2000 Pierangelo Masarati,"
			    << endl
		    	    << "Dipartimento di Ingegneria Aerospaziale,"
		    	    " Politecnico di Milano." << endl
		    	    << "MBDyn is free software, covered by the"
		    	    " GNU General Public License, and you are" << endl
		    	    << "welcome to change it and/or distribute"
		    	    " copies of it under certain conditions." << endl
		    	    << "Use 'mbdyn --license' to see the conditions."
			    << endl
		    	    << "There is absolutely no warranty for"
		    	    " MBDyn.  Use \"mbdyn --warranty\" for details."
			    << endl 
		    	    << endl);
#ifdef USE_MPI
        	cerr << "Process " << myrank 
	    		<< " (" << myrank+1 << " of " << WorldSize
            		<< ") is alive on " << ProcessorName << endl;
#endif /* USE_MPI */
      
      		/* Mostra la tabella dei simboli ed esce */
        	if (fShowSymbolTable > 0) {
#ifdef USE_MPI
	    		if (myrank == 0) {
#endif /* USE_MPI */
	        		Table t(31, 1);
	        		cout << "default symbol table:"
					<< endl << t << endl;
#ifdef USE_MPI
	    		}
#endif /* USE_MPI */
	 
#ifdef USE_EXCEPTIONS
	    		throw NoErr();
#else /* !USE_EXCEPTIONS */
#ifdef USE_MPI
	    		MPI::Finalize();
#endif /* USE_MPI */
	    		exit(EXIT_SUCCESS);
#endif /* !USE_EXCEPTIONS */
        	}
            
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
	        		pIn = (istream*)&cin;
	    		}
        	}
      
        	/* Gestione di input/output */      
        	Table* pT = NULL;
        	MathParser* pMP = NULL;
      
        	int last = 0;
        	while (last == 0) {
	    		if (CurrInputSource == FILE_STDIN
			    || CurrInputSource == FILE_OPT) {
	        		last = 1;
	    		} else if (CurrInputSource == FILE_ARGS) {
	        		sInputFileName = argv[currarg];
	        		DEBUGLCOUT(MYDEBUG_INPUT, 
		           		   "input file: "
					   << argv[currarg] << endl);
	    
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
	       
	            			cerr << "ADAMS input not implemented"
						" yet, cannot open file '"
						<< sInputFileName << "'"
						<< endl;
	            			THROW(ErrGeneric());
	        		} else {
#endif /* USE_ADAMS_PP */
	            			CurrInputFormat = MBDYN;
	       
	            			FileStreamIn.open(sInputFileName);
	            			if (!FileStreamIn) {
		        			cerr << endl 
			    				<< "Unable to open"
							" file '"
							<< sInputFileName 
			    				<< "'; aborting ..."
							<< endl;
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
							       Table(31, 1),
							       MBDynMM);
	        		}
	        		if (pMP == NULL) {
		    			SAFENEWWITHCONSTRUCTOR(pMP, 
		                           		       MathParser, 
					   			MathParser(*pT,
								    fRedefine), 
					   			MBDynMM);
		
		    			/* legge l'environment */
		    			GetEnviron(*pMP);
	        		} 
		
	        		/* parser del blocco di controllo */
	        		KeyTable K(0, NULL);
	     
	        		/* stream in ingresso */
	        		InputStream In(*pIn);
	        		MBDynParser HP(*pMP, K, In);
	     
	        		pIntg = RunMBDyn(HP, sInputFileName, 
						 sOutputFileName);
	        		FileStreamIn.close();
	        		break;
	    		}
	       
	    		case ADAMS:
	        		cerr << "ADAMS input not implemented yet!"
					<< endl;
	        		THROW(ErrNotImplementedYet());
	    
	    		default:
	        		cerr << "You shouldn't be here!" << endl;
	        		THROW(ErrGeneric());
	    		}
	    
	    		if (pIntg != NULL) {
	        		SAFEDELETE(pIntg, MBDynMM);
	    		}
	 
	    		if (fTable == 0 || argv[currarg] == NULL) {
	        		if (pMP != NULL) {
	            			SAFEDELETE(pMP, MBDynMM);
	            			pMP = NULL;
	        		}
	        		if (pT != NULL) {
	            			SAFEDELETE(pT, MBDynMM);
	            			pT = NULL;
	        		}
	    		}
	 
	    		time_t tSecs = 0;
	    		time_t tCents = 0;
#ifdef HAVE_TIMES_H	 
	    		/* Tempo di CPU impiegato */
	    		struct tms buf;
	    		times(&buf);
	    		clock_t ct = buf.tms_utime+buf.tms_cutime
				+buf.tms_stime+buf.tms_cstime;
	    		tSecs = ct/CLK_TCK;
	    		tCents = ((ct*100)/CLK_TCK)%100;
	    		cout << endl << "The simulation required " 
	        		<< tSecs << '.' << tCents 
	        		<< " seconds of CPU time";
#ifdef USE_MPI
	    		cout << " on " << ProcessorName;
#endif /* USE_MPI */
	    		cout << endl;
#endif /* HAVE_TIMES_H */
	 
#ifdef HAVE_GETOPT
	    		/* E-mail all'utente */
	    		if (iMailToBeSent) {
	        		SendMessage(sInputFileName, sMailToAddress,
					    tSecs, tCents);
	    		}
#endif /* HAVE_GETOPT */	    
        	}

#ifdef HAVE_GETOPT
		if (sMailToAddress) {
			SAFEDELETEARR(sMailToAddress, MBDynMM);
		}
#endif /* HAVE_GETOPT */
		if (sInputFileName) {
			SAFEDELETEARR(sInputFileName, MBDynMM);
		}
		if (sOutputFileName) {
			SAFEDELETEARR(sOutputFileName, MBDynMM);
		}

#ifdef USE_EXCEPTIONS
        	throw NoErr();
    	}
   
    	catch (NoErr) {     
        	silent_cout("MBDyn terminated normally" << endl);
#ifdef USE_MPI 
        	MPI::Finalize();
#endif /* USE_MPI */
      
        	exit(EXIT_SUCCESS);
    	}
    	catch (...) {
        	cerr << "An error occurred during the execution of MBDyn;"
	    		" aborting ... " << endl;
#ifdef USE_MPI
        	MPI::Finalize();
#endif /* USE_MPI */
        	exit(EXIT_FAILURE);
    	}
#endif /* USE_EXCEPTIONS */
   
#ifdef USE_MPI 
    	MPI::Finalize();
#endif /* USE_MPI */
    	return EXIT_SUCCESS;
}


Integrator* 
RunMBDyn(MBDynParser& HP, 
	 const char* const sInputFileName,
	 const char* const sOutputFileName)
{
    	DEBUGCOUTFNAME("RunMBDyn");
   
    	Integrator* pIntg = NULL;

    	/* flag di parallelo */
    	flag fParallel(0);

    	/* parole chiave */
    	const char* sKeyWords[] = { 
        	"begin",
		"end",
        	"data",
        	"integrator",
        	"multistep",
        	"rungekutta",
        	"parallel",
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
        	SSCHUR,
        	LASTKEYWORD
    	};
   
    	/* tabella delle parole chiave */
    	KeyTable K((int)LASTKEYWORD, sKeyWords);
   
    	/* Attacca la tabella al parser */
    	HP.PutKeyTable(K);
   
    	/* legge i dati della simulazione */
    	if (KeyWords(HP.GetDescription()) != BEGIN) {
        	cerr << endl 
	    		<< "Error: <begin> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << endl;
        	THROW(ErrGeneric());
    	}
   
    	if (KeyWords(HP.GetWord()) != DATA) {
        	cerr << endl 
	    		<< "Error: <begin: data;> expected at line " 
	    		<< HP.GetLineData() << "; aborting ..." << endl;
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
	        		CurrInt = SSCHUR;
	        		break;
		
            		default:
	        		cerr << endl 
		    			<< "Unknown integrator at line " 
	            			<< HP.GetLineData()
					<< "; aborting ..." << endl;
	        		THROW(ErrGeneric());
            		}
            		break;    

        	case PARALLEL:
#ifdef USE_MPI
            		fParallel = 1;
	    		break;
#else /* !USE_MPI */
            		cerr << "complile with -DUSE_MPI to have parallel"
				<< endl;
	    		THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
	    		break;
#endif /* USE_EXCEPTIONS */
#endif /* !USE_MPI */

        	case END:
	    		if (KeyWords(HP.GetWord()) != DATA) {
	        		cerr << endl 
		    			<< "Error: <end: data;> expected"
					" at line " << HP.GetLineData()
					<< "; aborting ..." << endl;
	        		THROW(ErrGeneric());
	    		}
	    		goto endofcycle;        
	 
        	default:
	    		cerr << endl 
	        		<< "Unknown description at line " 
	        		<< HP.GetLineData()
				<< "; aborting ..." << endl;
	    		THROW(ErrGeneric());      
        	}
    	}
   
   	/* Uscita dal ciclo infinito */
endofcycle:   
   
    	switch (CurrInt) {
    	case MULTISTEP:
		if (fParallel == 1) {
	    		cerr << "Sorry, multistep method cannot be parallel;"
	        		" aborting ..." << endl;
	    		THROW(ErrGeneric());
        	}
        	SAFENEWWITHCONSTRUCTOR(pIntg,
			               MultiStepIntegrator,
				       MultiStepIntegrator(HP, 
				       			   sInputFileName, 
							   sOutputFileName),
			       	       MBDynMM);
        	break;
      
    	case RUNGEKUTTA:
        	cerr << "Sorry, implicit Runge-Kutta isn't supported yet;"
	    		<< endl << "aborting ..." << endl;
        	THROW(ErrNotImplementedYet());

    	case SSCHUR:
        	if (!fParallel) {
	    		cerr << "Sorry Schur method is supported only"
	        		" for parallel computations; aborting ..."
				<< endl;
	    		THROW(ErrNotImplementedYet());
        	}
#ifdef USE_MPI
        	if (MPI::COMM_WORLD.Get_size() == 1) {
            		cerr << "Schur method is inefficient"
	        		" if used with only one processor;" << endl 
				<< "multistep method will be used instead"
				<< endl;

            		SAFENEWWITHCONSTRUCTOR(pIntg,
                                   	       MultiStepIntegrator,
					       MultiStepIntegrator(HP,
					       		sInputFileName,
							sOutputFileName),
                                   	       MBDynMM);
	    		break;
        	}
        	SAFENEWWITHCONSTRUCTOR(pIntg, 
			       	       SchurMultiStepIntegrator, 
				       SchurMultiStepIntegrator(HP, 
				       			sInputFileName, 
							sOutputFileName,
				       fParallel),
			       MBDynMM);
        	break;
#else /* !USE_MPI */
        	cerr << "compile with -DUSE_MPI to have parallel" << endl;
        	THROW(ErrGeneric());
#endif /* !USE_MPI */
    
    	default:
        	cerr << "Unknown integrator; aborting ..." << endl;
        	THROW(ErrGeneric());   
    	}
   
    	/* Runs the simulation */
    	pIntg->Run();
    
    	return pIntg;
}

void 
SendMessage(const char* const sInputFileName,
	    const char* const sMailToAddress,
	    time_t tSecs,
	    time_t tCents)
{
    	DEBUGCOUTFNAME("SendMessage");
   
    	/* Scrive il messaggio in un file temporaneo */
    	ofstream Msg("mbdyn.msg");
    	Msg << "Mbdyn terminated job ";
    	if (sInputFileName != NULL) {
        	Msg << "'" << sInputFileName << "' ";
    	}
#ifdef HAVE_TIMES_H
    	Msg << "in " << tSecs << '.' << tCents << " seconds of CPU time";
#endif /* HAVE_TIMES_H */
    	Msg << '.' << endl;
    	Msg.close();
   
    	/* Crea la linea di comando */
    	char* sCmd = NULL;
    	SAFENEWARR(sCmd, char, (29+strlen(sMailToAddress)+11+1), MBDynMM);
    	char* s = sCmd;
    	strcpy(s, "/bin/mail -s 'mbdyn message' ");
    	s = sCmd+strlen(sCmd);
    	strcpy(s, sMailToAddress);
    	s = sCmd+strlen(sCmd);
    	strcpy(s, " <mbdyn.msg");
    	system(sCmd);
   
    	/* Manda il messagio e cancella il file temporaneo */
    	system("rm mbdyn.msg");
    	SAFEDELETEARR(sCmd, MBDynMM);
}

