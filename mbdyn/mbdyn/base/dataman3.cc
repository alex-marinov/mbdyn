/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* Continua il DataManager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <sys/stat.h>
#include <unistd.h>

#if defined(HAVE_LOADABLE) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* HAVE_LOADABLE && HAVE_LTDL_H */

#include <dataman.h>
#include <dataman_.h>

#include <drive.h>
#include <presnode.h>
#include <readclaw.h>
#include <j2p.h>
#include <sah.h>

class NotAllowed {};

/* Legge i dati di controllo */

void DataManager::ReadControl(MBDynParser& HP, 
			      const char* sInputFileName, 
			      const char* sOutputFileName)
{
   DEBUGCOUTFNAME("DataManager::ReadControl");

   /* Nome del file di output */
   char* sOutName = NULL;
   if (sOutputFileName != NULL) {
      SAFESTRDUP(sOutName, sOutputFileName);
   }
   
   /* parole chiave del blocco di controllo */
   const char* sKeyWords[] = { 
      "end",
      "controldata",
      
      psReadControlNodes[Node::STRUCTURAL],      
      psReadControlNodes[Node::ELECTRIC],
      psReadControlNodes[Node::ABSTRACT],      
      psReadControlNodes[Node::PARAMETER],      
      psReadControlNodes[Node::HYDRAULIC],      
      
      psReadControlElems[Elem::AUTOMATICSTRUCTURAL],
      psReadControlElems[Elem::GRAVITY],	
      psReadControlElems[Elem::BODY],
      psReadControlElems[Elem::JOINT],
      psReadControlElems[Elem::BEAM],
      psReadControlElems[Elem::PLATE],      
      psReadControlElems[Elem::AIRPROPERTIES],
      psReadControlElems[Elem::ROTOR],
      psReadControlElems[Elem::AEROMODAL],
      psReadControlElems[Elem::AERODYNAMIC],      
      psReadControlElems[Elem::FORCE],      
      psReadControlElems[Elem::GENEL],
      psReadControlElems[Elem::ELECTRICBULK],
      psReadControlElems[Elem::ELECTRIC],      
      psReadControlElems[Elem::HYDRAULIC],      
      psReadControlElems[Elem::BULK],
      psReadControlElems[Elem::LOADABLE],
      psReadControlElems[Elem::RTAI_OUTPUT],

      psReadControlDrivers[Drive::FILEDRIVE],

      "loadable" "path",
      
      "skip" "initial" "joint" "assembly",
      "use",
      "in" "assembly",
      "initial" "stiffness",
      "omega" "rotates",
      "no",
      "yes",	
      "initial" "tolerance",
      "max" "initial" "iterations",

      "print",
      "dof" "stats",
      
      "title",
      "make" "restart" "file",
      "output" "file" "name",
      "output" "precision",
      "output" "frequency",
      "output" "results",
      "default" "output",
      "all",	
      "none",
      "reference" "frames",

      "default" "scale",
      
      NULL
   };
   
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
      END = 0,
      CONTROLDATA,
      STRUCTURALNODES,
      ELECTRICNODES,
      ABSTRACTNODES,
      PARAMETERNODES,
      HYDRAULICNODES,       
      
      AUTOMATICSTRUCTURAL,
      GRAVITY,
      RIGIDBODIES,
      JOINTS,
      BEAMS,
      PLATES,
      AIRPROPERTIES,
      ROTORS,
      AEROMODALS,
      AERODYNAMICELEMENTS,
      FORCES,
      GENELS,
      ELECTRICBULKELEMENTS,
      ELECTRICELEMENTS,
      HYDRAULICELEMENTS,
      BULKELEMENTS,
      LOADABLEELEMENTS,
      RTAIOUTPUTELEMENTS,

      FILEDRIVERS,

      LOADABLEPATH,
      
      SKIPINITIALJOINTASSEMBLY,
      USE,
      INASSEMBLY,
      INITIALSTIFFNESS,
      OMEGAROTATES,
      NO,
      YES,
      INITIALTOLERANCE,
      MAXINITIALITERATIONS,

      PRINT,
      DOFSTATS,
      
      TITLE,
      MAKERESTARTFILE,
      OUTPUTFILENAME,
      OUTPUTPRECISION,
      OUTPUTFREQUENCY,

      OUTPUTRESULTS,
      DEFAULTOUTPUT,
      ALL,
      NONE,
      REFERENCEFRAMES,

      DEFAULTSCALE,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);

   KeyWords CurrDesc;
   while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
      switch (CurrDesc) {
	 
	 /******** Nodes *********/

	 /* Numero di nodi strutturali attesi */
       case STRUCTURALNODES: {
#if defined(USE_STRUCT_NODES)	 
	  int iDmy = HP.GetInt();
	  NodeData[Node::STRUCTURAL].iNum = iDmy;
	  DofData[DofOwner::STRUCTURALNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Structural nodes: " << iDmy << std::endl);
	  break;
#else /* USE_STRUCT_NODES */
	  std::cerr << "you're not allowed to use structural nodes" 
		  << std::endl;
	  THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
       }
	 
	 /* Numero di nodi elettrici attesi */
       case ELECTRICNODES: {
#if defined(USE_ELECTRIC_NODES)	  
	  int iDmy = HP.GetInt();
	  NodeData[Node::ELECTRIC].iNum = iDmy;
	  DofData[DofOwner::ELECTRICNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Electric nodes: " << iDmy << std::endl);
	  break;
#else /* USE_ELECTRIC_NODES */
	  std::cerr << "you're not allowed to use electric nodes" << std::endl;
	  THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
       }	     
	 
	 /* Numero di nodi astratti attesi */
       case ABSTRACTNODES: {
#if defined(USE_ELECTRIC_NODES)	  
	  int iDmy = HP.GetInt();
	  NodeData[Node::ABSTRACT].iNum = iDmy;
	  DofData[DofOwner::ABSTRACTNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Abstract nodes: " << iDmy << std::endl);
	  break;
#else /* USE_ELECTRIC_NODES */
	  std::cerr << "you're not allowed to use abstract nodes" << std::endl;
	  THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
       }
	 
	 /* Numero di nodi astratti attesi */
       case PARAMETERNODES: {		  
	  int iDmy = HP.GetInt();
	  NodeData[Node::PARAMETER].iNum = iDmy;	     	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Parameter nodes: " << iDmy << std::endl);
	  break;
       }	     
	 
	 /* Numero di nodi idraulici attesi */
       case HYDRAULICNODES: {
#if defined(USE_HYDRAULIC_NODES)
	  int iDmy = HP.GetInt();
	  NodeData[Node::HYDRAULIC].iNum = iDmy;
	  DofData[DofOwner::HYDRAULICNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic nodes: " << iDmy << std::endl);
	  break;
#else /* defined(USE_HYDRAULIC_NODES) */
	  std::cerr << "you're not allowed to use hydraulic nodes" << std::endl;
	  THROW(ErrGeneric());
#endif /* defined(USE_HYDRAULIC_NODES) */
       }	     
	 
	 
	 /******** Elements *********/

	 /* Numero di corpi rigidi attesi */
#if defined(USE_STRUCT_NODES)	  
       case AUTOMATICSTRUCTURAL: {
	  int iDmy = HP.GetInt();
#ifdef DEBUG
#if 0
	  ElemData[Elem::AUTOMATICSTRUCTURAL].iNum = iDmy;
#endif /* 0 */
	  DEBUGLCOUT(MYDEBUG_INPUT, "Automatic structural elements expected: " 
		    << iDmy << std::endl);
#else
          iDmy = 0;
#endif
	  break;
       }
	 
	 /* Accelerazione di gravita' */
       case GRAVITY: {
	  if(ElemData[Elem::GRAVITY].iNum > 0) {
	     std::cerr
	       << "warning: gravity acceleration already defined;" << std::endl
	       << "only one definition will be considered" << std::endl;
	  }
	  ElemData[Elem::GRAVITY].iNum = 1;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Gravity acceleration expected in elements data" << std::endl);
	  break;
       }
	 
	 /* Numero di corpi rigidi attesi */
       case RIGIDBODIES: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BODY].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Rigid bodies: " << iDmy << std::endl);
	  break;
       }	     
	 
	 /* Numero di vincoli attesi */
       case JOINTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::JOINT].iNum = iDmy;	     
	  DofData[DofOwner::JOINT].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Joints: " << iDmy << std::endl);
	  if (iDmy > 0 ) {		       
	     fInitialJointAssemblyToBeDone = flag(1);
	  }
	  break;
       }	     
	 
	 /* Numero di travi attese */
       case BEAMS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BEAM].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Beams: " << iDmy << std::endl);
	  if (iDmy > 0 ) {		       
	     fInitialJointAssemblyToBeDone = flag(1);
	  }
	  break;
       }	     

#if defined(USE_AERODYNAMIC_ELEMS)
	 /* Elementi aerodinamici: proprieta' dell'aria */
       case AIRPROPERTIES: {
	  if (ElemData[Elem::AIRPROPERTIES].iNum > 0) {
	     std::cerr
	       << "warning: air properties already defined;" << std::endl
	       << "only one definition will be considered" << std::endl;
	  }
	  ElemData[Elem::AIRPROPERTIES].iNum = 1;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Air properties expected in elements data" << std::endl);
	  break;
       }
	 
	 /* Elementi aerodinamici: rotori */
       case ROTORS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::ROTOR].iNum = iDmy;	     
	  DofData[DofOwner::ROTOR].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Rotors: " << iDmy << std::endl);
	  break;
       }	     	     
	 
       case AEROMODALS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::AEROMODAL].iNum = iDmy;	     
	  DofData[DofOwner::AEROMODAL].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Aeromodals: " << iDmy << std::endl);
	  break;
       }	     	     
	 
	 /* Elementi aerodinamici: vari elementi aerodinamici senza dof */
       case AERODYNAMICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::AERODYNAMIC].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements: " << iDmy << std::endl);
	  break;
       }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

	 /* Numero di forze e coppie attese */
       case FORCES: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::FORCE].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Forces: " << iDmy << std::endl);
	  break;
       }	     
	 
#if defined(USE_ELECTRIC_NODES)
	 /* Numero di vincoli attesi */
       case GENELS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::GENEL].iNum = iDmy;	     
	  DofData[DofOwner::GENEL].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Genels: " << iDmy << std::endl);
	  break;
       }	     
	 
	 /* Numero di elementi elettrici attesi */
       case ELECTRICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::ELECTRIC].iNum = iDmy;	     
	  DofData[DofOwner::ELECTRIC].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Electric elements: " << iDmy << std::endl);
	  break;
       }	     
#endif /* USE_ELECTRIC_NODES */
	 
	 /* Numero di elementi idraulici attesi */
       case HYDRAULICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::HYDRAULIC].iNum = iDmy;	     
	  DofData[DofOwner::HYDRAULIC].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic elements: " << iDmy << std::endl);
	  break;
       }	     
	 
	 /* Numero di elementi elettrici attesi */
       case BULKELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BULK].iNum = iDmy;	     	    
	  DEBUGLCOUT(MYDEBUG_INPUT, "Bulk elements: " << iDmy << std::endl);
	  break;
       }	     
	 
#if defined(HAVE_LOADABLE)
	 /* Numero di elementi caricabili attesi */
       case LOADABLEELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::LOADABLE].iNum = iDmy;	     	    
	  DofData[DofOwner::LOADABLE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Loadable elements: " << iDmy << std::endl);
	  break;
       }

#if defined(HAVE_LTDL_H)
       case LOADABLEPATH: {
          int mode = 0;

	  if (loadableElemInitialized == false) {
	     if (lt_dlinit()) {
	  	std::cerr << "unable to initialize loadable elements" << std::endl;
      		THROW(ErrGeneric());
	     }
	     loadableElemInitialized = true;
	  }

	  if (HP.IsKeyWord("set")) {
	     mode = 0;
	  } else if (HP.IsKeyWord("add")) {
	     mode = 1;
	  }

	  const char *s = HP.GetFileName();
	  if (s == NULL) {
	     std::cerr << "missing path in \"loadable path\" statement "
		     "at line " << HP.GetLineData() << std::endl;
	     THROW(ErrGeneric());
	  }

	  if (mode == 0) {
	     if (lt_dlsetsearchpath(s) != 0) {
	        std::cerr << "unable to set path \"" << s 
			<< "in \"loadable path\" statement at line "
			<< HP.GetLineData() << std::endl;
	        THROW(ErrGeneric());
	     }
	  } else {
	     if (lt_dladdsearchdir(s) != 0) {
	        std::cerr << "unable to add path \"" << s 
			<< "\" in \"loadable path\" statement at line "
			<< HP.GetLineData() << std::endl;
	        THROW(ErrGeneric());
	     }
	  }

	  break;
       }
#endif /* HAVE_LTDL_H */
#endif /* defined(HAVE_LOADABLE) */

       case RTAIOUTPUTELEMENTS: {
#ifdef USE_RTAI
	  int iDmy = HP.GetInt();
	  ElemData[Elem::RTAI_OUTPUT].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "RTAI output elements: " << iDmy
		<< std::endl);
#else /* ! USE_RTAI */
          std::cerr << "cannot use RTAI output elements when not configured --with-rtai" << std::endl;
#endif /* ! USE_RTAI */
	  break;
       }

	 /* Numero di drivers attesi */
       case FILEDRIVERS: {
	  int iDmy = HP.GetInt();
	  DriveData[Drive::FILEDRIVE].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "File drivers: " << iDmy << std::endl);
	  break;
       }	     
	 
	    
	 /********* Miscellaneous *********/

#if defined(USE_STRUCT_NODES)	    
	 /* Spegne il flag di assemblaggio iniziale; 
	  * di default viene eseguito solo se sono definiti vincoli */
       case SKIPINITIALJOINTASSEMBLY: {
	  fSkipInitialJointAssembly = flag(1);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Skipping initial joint assembly" << std::endl);
	  break;
       }	     
	 
	 /* Uso di diversi tipi di elementi nell'assemblaggio iniziale */
       case USE: {
	  while (1) {		       
	     switch (KeyWords(HP.GetWord())) {
		/* Esce dal ciclo */
	      case INASSEMBLY: {
		 goto EndOfUse;
	      }
		
	      case RIGIDBODIES: {			    
		 ElemData[Elem::BODY].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Rigid bodies will be used in initial joint assembly" << std::endl);
		 break;
	      }
		
	      case GRAVITY: {			    
		 ElemData[Elem::GRAVITY].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Gravity will be used in initial joint assembly" << std::endl);
		 break;
	      }
		
	      case FORCES: {
		 ElemData[Elem::FORCE].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Forces will be used in initial joint assembly" << std::endl);
		 break;
	      }
		
		/* Lo lascio per backwards compatibility */
	      case BEAMS: {			    
#if 0
		 ElemData[Elem::BEAM].fToBeUsedInAssembly = flag(1);
#endif /* 0 */
		 DEBUGLCOUT(MYDEBUG_INPUT, "Beams are used in initial joint assembly by default" << std::endl);
		 break;
	      }			    

#if defined(USE_AERODYNAMIC_ELEMS)
	      case AERODYNAMICELEMENTS: {
		 ElemData[Elem::AERODYNAMIC].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements will be used in initial joint assembly" << std::endl);
		 
		 if (ElemData[Elem::AIRPROPERTIES].fToBeUsedInAssembly == flag(0)) {
		    ElemData[Elem::AIRPROPERTIES].fToBeUsedInAssembly = flag(1);
		 }
		 
		 break;
	      }
#endif /* USE_AERODYNAMIC_ELEMS */
		
#if defined(HAVE_LOADABLE)
	      case LOADABLEELEMENTS: {
		 ElemData[Elem::LOADABLE].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Loadable Elements will be used in initial joint assembly" << std::endl);
		 break;
	      }
#endif /* defined(HAVE_LOADABLE) */
		
		
		
		/* Elemento non autorizzato */
	      default: {
	         std::cerr << "Element type at line "
		   << HP.GetLineData() 
		     << " is not allowed; aborting ..." << std::endl;
		 
		 THROW(DataManager::ErrElemNotAllowedInAssembly());
	      }		       
		
		/* Errore */
	      case UNKNOWN: {
		 std::cerr << "Unknown element type at line "
		   << HP.GetLineData() << "; aborting ..." << std::endl;
		 
		 THROW(DataManager::ErrUnknownElem());
	      }		       			      		       
	     }
	  }		  
	  
	  EndOfUse:
	  break;		    
       }
	 
	 
	 /* Rigidezza delle molle fittizie usate 
	  * nell'assemblaggio iniziale; 
	  * se viene fornito un solo valore, viene usato sia per posizione
	  * che per velocita';
	  * se ne vengono forniti due, il primo vale per la posizione ed
	  * il secondo per la velocita' */
       case INITIALSTIFFNESS: {
	  dInitialPositionStiffness = 
	    HP.GetReal(dDefaultInitialStiffness);
	  
	  if (HP.fIsArg()) {		 
	     dInitialVelocityStiffness = 
	       HP.GetReal(dDefaultInitialStiffness);
	  } else {		 
	     dInitialVelocityStiffness = dInitialPositionStiffness;
	  }
	  
	  DEBUGLCOUT(MYDEBUG_INPUT, "Initial position stiffness: " 
		     << dInitialPositionStiffness << std::endl
		     << "Initial velocity stiffness: "
		     << dInitialVelocityStiffness << std::endl);
	  
	  break;
       }
	 
	 /* Omega solidale con il nodo o con il rif. globale */
       case OMEGAROTATES: {
	  if (HP.IsKeyWord("yes")) {
	     fOmegaRotates = flag(1);
	  } else if (HP.IsKeyWord("no")) {
	     fOmegaRotates = flag(0);
	  } else {
	     std::cerr << "Invalid option at line " 
		     << HP.GetLineData() << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  break;
       }
	 
	 /* Tolleranza nell'assemblaggio iniziale; viene calcolata come:
	  * sqrt(sum(res^2)/(1.+sum(sol^2)) */
       case INITIALTOLERANCE: {
	  dInitialAssemblyTol = 
	    HP.GetReal(dDefaultInitialAssemblyTol);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Initial assembly tolerance: " 
		     << dInitialAssemblyTol << std::endl);
	  break;
       }	     
	 
	 /* Numero massimo di iterazioni nell'assemblaggio iniziale;
	  * di default ne e' consentita solo una, indice di condizioni 
	  * iniziali corrette */
       case MAXINITIALITERATIONS: {
	  iMaxInitialIterations = 
	    HP.GetInt(iDefaultMaxInitialIterations);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Max initial iterations: " 
		     << iMaxInitialIterations << std::endl);
	  break;
       }
#endif /* USE_STRUCT_NODES */

       case PRINT:
	  if (HP.IsKeyWord("dof" "stats")) {
	     fPrintDofStats = 1;
	  }
	  break;
	 
	 /* Titolo */
       case TITLE: {
	  ASSERT(sSimulationTitle == NULL);
	  if (sSimulationTitle != NULL) {
	     SAFEDELETEARR(sSimulationTitle);
	  }
	  const char* sTmp(HP.GetStringWithDelims());
	  SAFESTRDUP(sSimulationTitle, sTmp);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Simulation title: \"" 
		     << sSimulationTitle << '"' << std::endl);
	  break;
       }	     
	 
	 /* Crea il file di restart */
       case MAKERESTARTFILE: {
	  DEBUGLCOUT(MYDEBUG_INPUT, "Restart file will be generated " << std::endl);
	  if (HP.fIsArg()) {
	     if (HP.IsKeyWord("iterations")) {
		RestartEvery = ITERATIONS;
		iRestartIterations = HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT, "every " << iRestartIterations << " iterations" << std::endl);
	     } else if (HP.IsKeyWord("time")) {
		RestartEvery = TIME;
		dRestartTime = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "every " << dRestartTime << " time units" << std::endl);
	     } else {
		std::cerr << "Error: unrecognized restart option at line "
		  << HP.GetLineData() << std::endl;

		THROW(DataManager::ErrGeneric());
	     }
	  } else {
	     RestartEvery = ATEND;
	  }
	  break;
       }
	 
       case OUTPUTFILENAME: {
	  if (sOutName != NULL) {
	     SAFEDELETEARR(sOutName);
	     sOutName = NULL;
	  }
	  const char* sTmp(HP.GetFileName());
	  if (sTmp == NULL) {
	     std::cerr << "Null file name at line "
		     << HP.GetLineData() << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  SAFESTRDUP(sOutName, sTmp);
	  break;
       }

       case OUTPUTPRECISION: {
          int iPrec = HP.GetInt();
	  OutHdl.SetPrecision(iPrec);
	  break;
       }
	 
       case OUTPUTFREQUENCY: {
          integer iFreq = HP.GetInt();
	  if (iFreq < 1) {
		  std::cerr << "Illegal output frequency " << iFreq
			  << " at line " << HP.GetLineData() << std::endl;
		  THROW(ErrGeneric());
	  }
	  iOutputFrequency = iFreq;
	  break;
       }
	 
       case OUTPUTRESULTS: {
	while (HP.fIsArg()) {
	
		/* require support for ADAMS/View .res output */
	  	if (HP.IsKeyWord("adams")) {
	 	 	ResMode |= RES_ADAMS;

	  		if (HP.fIsArg() && HP.IsKeyWord("model" "name")) {
	     			if (sAdamsModelName != NULL) {
					SAFEDELETEARR(sAdamsModelName);
					sAdamsModelName = NULL;
	     			}
	     
	     			const char *tmp = HP.GetStringWithDelims();
	     			SAFESTRDUP(sAdamsModelName, tmp);
	  		} else {
	     			SAFESTRDUP(sAdamsModelName, "mbdyn");
	  		}

		/* require support for MotionView output */
		} else if (HP.IsKeyWord("motion" "view")) {
			ResMode |= RES_MOTIONVIEW;

			/*
			 * add output info
			 */
		}
	}
	break;
       }
	 
       case DEFAULTOUTPUT: {
	  while (HP.fIsArg()) {
	     KeyWords CurrDefOut(KeyWords(HP.GetWord()));
	     switch (CurrDefOut) {
	      case ALL: {
		 for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		    ElemData[iCnt].fDefaultOut = flag(1);
		 }			 
		 for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		    NodeData[iCnt].fDefaultOut = flag(1);
		 }			 
		 break;
	      }
		
	      case NONE: {
		 for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		    ElemData[iCnt].fDefaultOut = flag(0);
		 }			 
		 for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		    NodeData[iCnt].fDefaultOut = flag(0);
		 }			 
		 break;
	      }

#if defined(USE_STRUCT_NODES)
	      case REFERENCEFRAMES:
		 fOutputFrames = 1;
		 break;
#endif /* USE_STRUCT_NODES */

#if defined(USE_STRUCT_NODES)
	      case STRUCTURALNODES: {
		 NodeData[Node::STRUCTURAL].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_STRUCT_NODES */
		
#if defined(USE_ELECTRIC_NODES)
	      case ELECTRICNODES: {			 
		 NodeData[Node::ELECTRIC].fDefaultOut = flag(1);
		 break;
	      }
		
	      case ABSTRACTNODES: {			 
		 NodeData[Node::ABSTRACT].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
	      case HYDRAULICNODES: {
		 NodeData[Node::HYDRAULIC].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_HYDRAULIC_NODES */
		
#if defined(USE_STRUCT_NODES)
	      case GRAVITY: {			 
		 ElemData[Elem::GRAVITY].fDefaultOut = flag(1);
		 break;
	      }
		
	      case RIGIDBODIES: {			 
		 ElemData[Elem::BODY].fDefaultOut = flag(1);
		 break;
	      }
		
	      case JOINTS: {			 
		 ElemData[Elem::JOINT].fDefaultOut = flag(1);
		 break;
	      }
		
	      case BEAMS: {			 
		 ElemData[Elem::BEAM].fDefaultOut = flag(1);
		 break;
	      }
		
	      case PLATES: {
		 ElemData[Elem::PLATE].fDefaultOut = flag(1);
		 break;
	      }

#if defined(USE_AERODYNAMIC_ELEMS)
	      case AIRPROPERTIES: {			 
		 ElemData[Elem::AIRPROPERTIES].fDefaultOut = flag(1);
		 break;
	      }
		
	      case ROTORS: {			 
		 ElemData[Elem::ROTOR].fDefaultOut = flag(1);
		 break;
	      }
		
	      case AEROMODALS: {			 
		 ElemData[Elem::AEROMODAL].fDefaultOut = flag(1);
		 break;
	      }
	      case AERODYNAMICELEMENTS: {			 
		 ElemData[Elem::AERODYNAMIC].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

	      case FORCES: {			 
		 ElemData[Elem::FORCE].fDefaultOut = flag(1);
		 break;
	      }
		
#if defined(USE_ELECTRIC_NODES)
	      case GENELS: {			 
		 ElemData[Elem::GENEL].fDefaultOut = flag(1);
		 break;
	      }
		
	      case ELECTRICBULKELEMENTS: {			 
		 ElemData[Elem::ELECTRICBULK].fDefaultOut = flag(1);
		 break;
	      }
		
	      case ELECTRICELEMENTS: {			 
		 ElemData[Elem::ELECTRIC].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC)
	      case HYDRAULICELEMENTS: {
		 ElemData[Elem::HYDRAULIC].fDefaultOut = flag(1);
		 break;
	      }
#endif /* USE_HYDRAULIC */
#if defined(HAVE_LOADABLE)
	      case LOADABLEELEMENTS: {
		 ElemData[Elem::LOADABLE].fDefaultOut = flag(1);
		 break;
	      }
#endif /* defined(HAVE_LOADABLE) */
		
		
	      case UNKNOWN: {
		 std::cerr << "warning: unknown output case at line " 
		   << HP.GetLineData() << std::endl;
		 ASSERT(0);		 
		 break;
	      }
		
	      default: {
		 std::cerr << "case " << sKeyWords[CurrDesc] << " at line "
		   << HP.GetLineData() << " is not allowed" << std::endl;
		 ASSERT(0);		 
		 break;
	      }		    		    
	     }		 
	  }
	  
	  break;
       }
	 
	 
       case DEFAULTSCALE: {
#ifndef __HACK_SCALE_RES__
	  std::cerr << "warning: residual and solution scaling is disabled"
		  << std::endl;
#endif /* !__HACK_SCALE_RES__ */

	  while (HP.fIsArg()) {
	     KeyWords CurrDefOut(KeyWords(HP.GetWord()));
	     doublereal dScale = HP.GetReal(1.);

	     switch (CurrDefOut) {
	      case ALL: {
		 for (int iCnt = 0; iCnt < DofOwner::LASTDOFTYPE; iCnt++) {
		    DofData[iCnt].dDefScale= dScale;
		 }			 
	      }
		
#if defined(USE_STRUCT_NODES)
	      case STRUCTURALNODES: {
		 DofData[DofOwner::STRUCTURALNODE].dDefScale = dScale;
		 break;
	      }
#endif /* USE_STRUCT_NODES */
		
#if defined(USE_ELECTRIC_NODES)
	      case ELECTRICNODES: {
		 DofData[DofOwner::ELECTRICNODE].dDefScale = dScale;
		 break;
	      }
		
	      case ABSTRACTNODES: {
		 DofData[DofOwner::ABSTRACTNODE].dDefScale = dScale;
		 break;
	      }
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
	      case HYDRAULICNODES: {
		 DofData[DofOwner::HYDRAULICNODE].dDefScale = dScale;
		 break;
	      }
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_STRUCT_NODES)
	      case JOINTS: {
		 DofData[DofOwner::JOINT].dDefScale = dScale;
		 break;
	      }
		
#if defined(USE_AERODYNAMIC_ELEMS)
	      case ROTORS: {
		 DofData[DofOwner::ROTOR].dDefScale = dScale;
		 break;
	      }
		
	      case AEROMODALS: {
		 DofData[DofOwner::AEROMODAL].dDefScale = dScale;
		 break;
	      }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
	      case GENELS: {
		 DofData[DofOwner::GENEL].dDefScale = dScale;
		 break;
	      }
		
	      case ELECTRICBULKELEMENTS: {
		 DofData[DofOwner::ELECTRICBULK].dDefScale = dScale;
		 break;
	      }
		
	      case ELECTRICELEMENTS: {
		 DofData[DofOwner::ELECTRIC].dDefScale = dScale;
		 break;
	      }
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
	      case HYDRAULICELEMENTS: {
		 DofData[DofOwner::HYDRAULIC].dDefScale = dScale;
		 break;
	      }
#endif /* USE_HYDRAULIC */

#if defined(HAVE_LOADABLE)
	      case LOADABLEELEMENTS: {
		 DofData[DofOwner::LOADABLE].dDefScale = dScale;
		 break;
	      }
#endif /* defined(HAVE_LOADABLE) */
		
	      case UNKNOWN: {
		 std::cerr << "warning: unknown output case at line " 
		   << HP.GetLineData() << std::endl;
		 ASSERT(0);		 
		 break;
	      }
		
	      default: {
		 std::cerr << "case " << sKeyWords[CurrDesc] << " at line "
		   << HP.GetLineData() << " is not allowed" << std::endl;
		 ASSERT(0);		 
		 break;
	      }		    		    
	     }		 
	  }
	  
	  break;
       }
	 
	 
	 /* add more entries ... */
	 
       case UNKNOWN: { /*
			* If description is not in key table the parser
	                * returns UNKNONW, so "default" can be used to 
	                * intercept control cases that are not allowed.
			*/
	  DEBUGCERR("");
	  std::cerr << "unknown description at line " 
		  << HP.GetLineData() << std::endl;
	  ASSERT(0);
	  break;
       }
	 
       default: {
	  std::cerr << "case " << sKeyWords[CurrDesc] << " at line " 
	    << HP.GetLineData() << " is not allowed" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }	 
      }
   }
   
   if (KeyWords(HP.GetWord()) != CONTROLDATA) {	
      DEBUGCERR("");
      std::cerr << "<end: control data;> expected at line " 
	<< HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      std::cerr << "semicolon expected at line "
	      << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }   
   
   /* Inizializza l'OutputHandler */
   if (sOutName == NULL) {
      if (sInputFileName == NULL) {
	 OutHdl.Init(sDefaultOutputFileName, 0);
      } else {
	 OutHdl.Init(sInputFileName, 0);
      }
   } else {

             /* FIXME: ora in caso di solutore parallelo, se sOutName
	      * e' una directory il nome dei file viene incasinato */
	     
      struct stat s;
      if (stat(sOutName, &s) == 0) {
	 char *tmpout = NULL;
	 const char *tmpin = sInputFileName ? sInputFileName : sDefaultOutputFileName;
	 
	 if (S_ISDIR(s.st_mode)) {
	    int l = strlen(sOutName)+strlen(tmpin)+1;
	    SAFENEWARR(tmpout, char, l);
	    strcpy(tmpout, sOutName);
	    strcat(tmpout, tmpin);
	    SAFEDELETEARR(sOutName);
	    sOutName = tmpout;
	 }
      }
      
      OutHdl.Init(sOutName, 0);
      SAFEDELETEARR(sOutName);
   }

   OutHdl.Log() << "output frequency: " << iOutputFrequency << std::endl;

   DEBUGLCOUT(MYDEBUG_INPUT, "End of control data" << std::endl);
} /* End of DataManager::ReadControl() */


/* Legge se un item deve scrivere sull'output handler */

flag DataManager::fReadOutput(MBDynParser& HP, enum Elem::Type t)
{
   flag fDef = fGetDefaultOutputFlag(t);
   if (!HP.IsKeyWord("output")) {
      return fDef;
   }
   
   if (HP.IsKeyWord("no")) {    
      return flag(0);
   } else if (HP.IsKeyWord("yes")) {         
      return flag(1);
   } else if (HP.IsKeyWord("default")) {
      return fDef;
   } else {    
      std::cerr << "Unknown output flag for element \""
	<< psElemNames[t] << "\" at line " << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return fDef;
#endif /* USE_EXCEPTIONS */
} /* End of DataManager::fReadOutput */


flag DataManager::fReadOutput(MBDynParser& HP, enum Node::Type t)
{
   flag fDef = fGetDefaultOutputFlag(t);
   if (!HP.IsKeyWord("output")) {
      return fDef;
   }
   
   if (HP.IsKeyWord("no")) {    
      return flag(0);
   } else if (HP.IsKeyWord("yes")) {         
      return flag(1);
   } else if (HP.IsKeyWord("default")) {
      return fDef;
   } else {    
      std::cerr << "Unknown output flag for node \"" 
	<< psNodeNames[t] << "\" at line " << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return fDef;
#endif /* USE_EXCEPTIONS */
} /* End of DataManager::fReadOutput */


doublereal
DataManager::dReadScale(MBDynParser& HP, enum DofOwner::Type t)
{
	doublereal d = dGetDefaultScale(t);

	if (!HP.IsKeyWord("scale")) {
		return d;
	}

#ifndef __HACK_SCALE_RES__
	std::cerr << "warning: residual and solution scaling is disabled"
		<< std::endl;
#endif /* !__HACK_SCALE_RES__ */

	if (!HP.IsKeyWord("default")) {
		d = HP.GetReal(d);
	}

	return d;
}


/* legge i nodi e li costruisce */

int
DataManager::ReadScalarAlgebraicNode(MBDynParser& HP, 
		unsigned int uLabel, Node::Type type,
		doublereal& dX)
{
	/* verifica di esistenza del nodo */
	if (pFindNode(type, uLabel) != NULL) {
		std::cerr << "line " << HP.GetLineData() 
      			<< ": " << psNodeNames[type] << "(" << uLabel
      			<< ") already defined" << std::endl;
		THROW(DataManager::ErrGeneric());
	}		  

	if (HP.fIsArg()) {
		/* eat keyword "value" */
		if (!HP.IsKeyWord("value")) {
			pedantic_cerr(psNodeNames[type] << "(" << uLabel 
     				<< "): initial value specified without "
     				"\"value\" keyword (deprecated)" << std::endl);
		}
		dX = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x = " << dX 
				<< " is supplied" << std::endl);
		return 1;
	}

	return 0;
}

int
DataManager::ReadScalarDifferentialNode(MBDynParser& HP, 
		unsigned int uLabel, Node::Type type,
		doublereal& dX, doublereal& dXP)
{
	if (ReadScalarAlgebraicNode(HP, uLabel, type, dX) == 1) {
		if (HP.IsKeyWord("derivative")) {
			dX = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, 
					"Initial derivative value xp = " 
					<< dXP << " is supplied" << std::endl);
			return 1;
		}
	}
	
	return 0;	
}

void DataManager::ReadNodes(MBDynParser& HP)
{   
   DEBUGCOUTFNAME("DataManager::ReadNodes");
     
   /* parole chiave del blocco di controllo */
   const char* sKeyWords[] = { 
      "end",
	"nodes",

	psReadNodesNodes[Node::STRUCTURAL],
	psReadNodesNodes[Node::ELECTRIC],
	psReadNodesNodes[Node::ABSTRACT],
	psReadNodesNodes[Node::PARAMETER],
	psReadNodesNodes[Node::HYDRAULIC],

	"output"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	
	END = 0,
	NODES,
	
	STRUCTURAL,
	ELECTRIC,
	ABSTRACT,
	PARAMETER,
	HYDRAULIC,

	OUTPUT,
	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   
   /* struttura di servizio che conta i nodi tipo per tipo */
   int iNumTypes[Node::LASTNODETYPE];
   for (int i = 0; i < Node::LASTNODETYPE; i++) {      
      iNumTypes[i] = NodeData[i].iNum;
   }      
         
   int iMissingNodes = iTotNodes;
   DEBUGLCOUT(MYDEBUG_INPUT, "Expected nodes: " << iMissingNodes << std::endl);

   KeyWords CurrDesc;   
   while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
      
      if (CurrDesc == OUTPUT) {
	 DEBUGLCOUT(MYDEBUG_INPUT, "nodes to be output: ");
	 
	 Node::Type Typ;
	 switch (KeyWords(HP.GetWord())) {

	  case STRUCTURAL: {
#if defined(USE_STRUCT_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "structural" << std::endl);
	     Typ = Node::STRUCTURAL;
	     break;
#else /* USE_STRUCT_NODES */
	     std::cerr << "you're not allowed to use structural nodes"
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
	  }
	    
	  case ELECTRIC: {	
#if defined(USE_ELECTRIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "electric" << std::endl);
	     Typ = Node::ELECTRIC;	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     std::cerr << "you're not allowed to use electric nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	  case ABSTRACT: {	   
#if defined(USE_ELECTRIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "abstract" << std::endl);
	     Typ = Node::ABSTRACT;
	     break;
#else /* USE_ELECTRIC_NODES */
	     std::cerr << "you're not allowed to use abstract nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	  case PARAMETER: {	    
	     DEBUGLCOUT(MYDEBUG_INPUT, "parameter" << std::endl);
	     Typ = Node::PARAMETER;	     
	     break;
	  }

	  case HYDRAULIC: {
#if defined (USE_HYDRAULIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic" << std::endl);
	     Typ = Node::HYDRAULIC;	     
	     break;
#else /* defined (USE_HYDRAULIC_NODES) */
	     std::cerr << "you're not allowed to use hydraulic nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* defined (USE_HYDRAULIC_NODES) */
	  }
	    
	  default: {
	     std::cerr << "Error: unknown node type, cannot modify output" 
		     << std::endl;

	     THROW(DataManager::ErrUnknownNode());
	  }
	 }
	 
	 while (HP.fIsArg()) {
	    unsigned int uL = (unsigned int)HP.GetInt();	    	 
	    Node* pN = pFindNode(Typ, uL);
	    if (pN == NULL) {
	       std::cerr << "Error: " << psNodeNames[Typ] << "(" 
		       << uL << ") is not defined; output cannot be modified" 
		       << std::endl;
	    } else {
	       DEBUGLCOUT(MYDEBUG_INPUT, "node " << uL << std::endl);
	       pN->SetOutputFlag(flag(1));
	    }
	 }	 
      } else {      
	 /* puntatore al puntatore al nodo */
	 Node** ppN = NULL;      
	 
	 /* legge la label */
	 unsigned int uLabel = 0;
	 if (CurrDesc >= 0) {	 
	    uLabel = unsigned(HP.GetInt());
	 }
	 
	 /* Nome del nodo */
	 const char *sName = NULL;
	 if (HP.IsKeyWord("name")) {
	    const char *sTmp = HP.GetStringWithDelims();
	    SAFESTRDUP(sName, sTmp);
	 }
	 
	 /* in base al tipo, avviene l'allocazione */
	 switch (CurrDesc) {
	    
	    /* Struttura del blocco di lettura dei nodi:
	     * 
	     * - test sul numero di nodi letti rispetto a quelli dichiarati
	     *   nel blocco di controllo:
	     *       if(iNumTypes[Node::??]-- <= 0)
	     *         < gestione dell'errore >
	     * 
	     * - verifica di esistenza del nodo 
	     *       if(pFindNode(Node::??, uLabel) != NULL)
	     *         < gestione dell'errore >
	     * 
	     * - lettura dati specifici
	     * 
	     * - allocazione e costruzione del nodo:
	     *   - assegnazione del puntatore:
	     *       ppN = NodeData[Node::??].ppFirstNode+
	     *         iNumTypes[Node::??];
	     * 
	     *   - allocazione e costruzione:
	     *       SAFENEW((??Node*)*ppN, ??Node(uLabel));
	     * 
	     * - correzione del DofOwner relativo al nodo:
	     *       (DofData[DofOwner::??].pFirstDofOwner+
	     *         iNumTypes[Node::??])->iNumDofs = 
	     *         DofData[DofOwner::??].iSize;
	     * 
	     * - scrittura dei dati specifici dell'oggetto creato.
	     *   In alternativa i dati possono essere passato tutti 
	     *   attraverso il costruttore.
	     */
	    

	    /* Nodi strutturali */
	  case STRUCTURAL: {
#if defined(USE_STRUCT_NODES)	   
	     silent_cout("Reading structural node " << uLabel << std::endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::STRUCTURAL]-- <= 0) {
		DEBUGCERR("");
		std::cerr << "line " << HP.GetLineData() 
		  << ": structural node " << uLabel 
		  << " exceedes structural nodes number" << std::endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindStructNode(uLabel) != NULL) {
		DEBUGCERR("");
		std::cerr << "line " << HP.GetLineData() 
		  << ": structural node " << uLabel
		  << " already defined" << std::endl;

		THROW(DataManager::ErrGeneric());
	     }		  
	     
	     
	     /* lettura dei dati specifici */
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::STRUCTURAL].iNum
	       -iNumTypes[Node::STRUCTURAL]-1;
	     ppN = NodeData[Node::STRUCTURAL].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner+i;
	     
	     *ppN = ReadStructNode(this, HP, pDO, uLabel);
	     HP.PutKeyTable(K);
	     
	     break;
#else /* USE_STRUCT_NODES */
	     std::cerr << "you're not allowed to use structural nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
	  }
	    
	    /* nodi elettrici */
	  case ELECTRIC: {
#if defined(USE_ELECTRIC_NODES)	   
	     silent_cout("Reading electric node " << uLabel << std::endl);

	     /* 
	      * verifica che non siano gia' stati letti tutti 
	      * quelli previsti 
	      */
	     if (iNumTypes[Node::ELECTRIC]-- <= 0) {
		std::cerr << "line " << HP.GetLineData() 
      			<< ": " << psNodeNames[Node::ELECTRIC] << "(" << uLabel
      			<< ") exceeds declared number" << std::endl;
		THROW(DataManager::ErrGeneric());
	     }

	     /* Initial values */
	     doublereal dx(0.);
	     doublereal dxp(0.);

	     ReadScalarDifferentialNode(HP, uLabel, Node::ELECTRIC, dx, dxp);
	     doublereal dScale = dReadScale(HP, DofOwner::ELECTRICNODE);
	     flag fOut = fReadOutput(HP, Node::ELECTRIC);
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::ELECTRIC].iNum
	       -iNumTypes[Node::ELECTRIC]-1;
	     ppN = NodeData[Node::ELECTRIC].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::ELECTRICNODE].pFirstDofOwner+i;
	     pDO->SetScale(dScale);
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    ElectricNode,
				    ElectricNode(uLabel, pDO, dx, dxp, fOut));
	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     std::cerr << "you're not allowed to use electric nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	    /* nodi astratti */
	  case ABSTRACT: {
#if defined(USE_ELECTRIC_NODES)	   
	     silent_cout("Reading abstract node " << uLabel << std::endl);
	     
	     /* 
	      * verifica che non siano gia' stati letti tutti 
	      * quelli previsti 
	      */
	     if (iNumTypes[Node::ABSTRACT]-- <= 0) {
		std::cerr << "line " << HP.GetLineData() 
      			<< ": " << psNodeNames[Node::ABSTRACT] << "(" << uLabel
      			<< ") exceeds declared number" << std::endl;
		THROW(DataManager::ErrGeneric());
	     }

	     /* lettura dei dati specifici */
	     doublereal dx(0.);
	     doublereal dxp(0.);
	     ReadScalarDifferentialNode(HP, uLabel, Node::ABSTRACT, dx, dxp);
	     doublereal dScale = dReadScale(HP, DofOwner::ABSTRACTNODE);
	     flag fOut = fReadOutput(HP, Node::ABSTRACT);
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::ABSTRACT].iNum
	       -iNumTypes[Node::ABSTRACT]-1;
	     ppN = NodeData[Node::ABSTRACT].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::ABSTRACTNODE].pFirstDofOwner+i;
	     pDO->SetScale(dScale);
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    AbstractNode,
				    AbstractNode(uLabel, pDO, dx, dxp, fOut));
	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     std::cerr << "you're not allowed to use abstract nodes" 
		     << std::endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	    /* parametri */
	  case PARAMETER: {	     
	     silent_cout("Reading parameter " << uLabel << std::endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::PARAMETER]-- <= 0) {
		DEBUGCERR("");
		std::cerr << "line " << HP.GetLineData() 
		  << ": parameter " << uLabel
		  << " exceedes parameters number" << std::endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindNode(Node::PARAMETER, uLabel) != NULL) {
		DEBUGCERR("");
		std::cerr << "line " << HP.GetLineData() 
		  << ": parameter " << uLabel
		  << " already defined" << std::endl;
		
		THROW(DataManager::ErrGeneric());
	     }		  
	     
	     /* bound a elemento */
	     if (HP.IsKeyWord("element")) {
		DEBUGLCOUT(MYDEBUG_INPUT, "parameter node " << uLabel 
			   << "is linked to an element" << std::endl);
		flag fOut = fReadOutput(HP, Node::PARAMETER);
		
		/* allocazione e creazione */
		int i = NodeData[Node::PARAMETER].iNum
		  -iNumTypes[Node::PARAMETER]-1;
		ppN = NodeData[Node::PARAMETER].ppFirstNode+i;
		
		SAFENEWWITHCONSTRUCTOR(*ppN, 
				       Elem2Param,
				       Elem2Param(uLabel, 
						  &DummyDofOwner, 
						  fOut));

	     } else if (HP.IsKeyWord("sample" "and" "hold") ||
			     HP.IsKeyWord("sample'n'hold")) {

		DEBUGLCOUT(MYDEBUG_INPUT, "parameter node " << uLabel 
			   << "is a sample-and-hold" << std::endl);
		
		ScalarDof SD(ReadScalarDof(this, HP, 0));
		HP.PutKeyTable(K);

		DriveCaller *pDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller,
				TimeDriveCaller(&DrvHdl));

		doublereal dSP = HP.GetReal();
		if (dSP <= 0.) {
			std::cerr << "illegal sample period for SampleAndHold("
				<< uLabel << ") at line " << HP.GetLineData()
				<< std::endl;
			THROW(ErrGeneric());
		}
		
		flag fOut = fReadOutput(HP, Node::PARAMETER);
		
		/* allocazione e creazione */
		int i = NodeData[Node::PARAMETER].iNum
		  -iNumTypes[Node::PARAMETER]-1;
		ppN = NodeData[Node::PARAMETER].ppFirstNode+i;
		
		SAFENEWWITHCONSTRUCTOR(*ppN, 
				       SampleAndHold,
				       SampleAndHold(uLabel, 
						       &DummyDofOwner,
						       SD.pNode,
						       pDC,
						       dSP,
						       fOut));
		
	     /* strain gage */
	     } else if (HP.IsKeyWord("straingage")) {
		DEBUGLCOUT(MYDEBUG_INPUT, "parameter node " << uLabel 
			   << "is a strain gage" << std::endl);
		
		doublereal dY = HP.GetReal();
		doublereal dZ = HP.GetReal();
		
		flag fOut = fReadOutput(HP, Node::PARAMETER);
		
		/* allocazione e creazione */
		int i = NodeData[Node::PARAMETER].iNum
		  -iNumTypes[Node::PARAMETER]-1;
		ppN = NodeData[Node::PARAMETER].ppFirstNode+i;
		
		SAFENEWWITHCONSTRUCTOR(*ppN, 
				       StrainGageParam,
				       StrainGageParam(uLabel, 
						       &DummyDofOwner,
						       dY, dZ,
						       fOut));
		
		/* parametro generico */
	     } else {	     
			     
		/* lettura dei dati specifici */
		doublereal dX(0.);	    
		if (HP.fIsArg()) {
		   dX = HP.GetReal();
		   DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x = " << dX
			      << " is supplied" << std::endl);
		}	      
		
		flag fOut = fReadOutput(HP, Node::PARAMETER);
		
		/* allocazione e creazione */
		int i = NodeData[Node::PARAMETER].iNum
		  -iNumTypes[Node::PARAMETER]-1;
		ppN = NodeData[Node::PARAMETER].ppFirstNode+i;
		
		SAFENEWWITHCONSTRUCTOR(*ppN, 
				       ParameterNode,
				       ParameterNode(uLabel,
						     &DummyDofOwner, 
						     dX, fOut));
	     }
	     
	     break;
	  }

	    
	    
#if defined(USE_HYDRAULIC_NODES)
	    
	    /* nodi idraulici */
	  case HYDRAULIC: {
	     silent_cout("Reading hydraulic node " << uLabel << std::endl);
	     
	     /* 
	      * verifica che non siano gia' stati letti tutti 
	      * quelli previsti 
	      */
	     if (iNumTypes[Node::HYDRAULIC]-- <= 0) {
		std::cerr << "line " << HP.GetLineData() 
      			<< ": " << psNodeNames[Node::HYDRAULIC] << "(" << uLabel
      			<< ") exceeds declared number" << std::endl;
		THROW(DataManager::ErrGeneric());
	     }

	     /* lettura dei dati specifici */
	     doublereal dx(0.);
	     ReadScalarAlgebraicNode(HP, uLabel, Node::ABSTRACT, dx);
	     doublereal dScale = dReadScale(HP, DofOwner::HYDRAULICNODE);
	     flag fOut = fReadOutput(HP, Node::HYDRAULIC);
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::HYDRAULIC].iNum
	       -iNumTypes[Node::HYDRAULIC]-1;
	     ppN = NodeData[Node::HYDRAULIC].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::HYDRAULICNODE].pFirstDofOwner+i;
	     pDO->SetScale(dScale);
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    PressureNode,
				    PressureNode(uLabel, pDO, dx, fOut));
	     
	     break;
	  }
#endif /* USE_HYDRAULIC_NODES */
	    
	    
	    /* aggiungere eventuali nuovi tipi di nodo */
	    
	    
	    /* in caso di tipo errato */
	  case UNKNOWN: {
	     DEBUGCERR("");
	     std::cerr << "unknown node type at line " 
		     << HP.GetLineData() << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
	    
	  default: {
	     DEBUGCERR("");
	     std::cerr << "node type " << sKeyWords[CurrDesc] 
	       << " at line " << HP.GetLineData() << " is not allowed" 
	       << std::endl;

	     THROW(DataManager::ErrGeneric());
	  }
	 }
	 
	 /* verifica dell'allocazione - comune a tutti i casi */
	 if (*ppN == NULL) {
	    DEBUGCERR("");
	    std::cerr << "error in allocation, "
		    "item DataManager::*ppNodes, node " 
		    << uLabel << std::endl;

	    THROW(ErrMemory());
	 }

	 if (sName != NULL) {
 	    (*ppN)->PutName(sName);
	    SAFEDELETEARR(sName);
	 }
	 
	 /* Decrementa i nodi attesi */
	 iMissingNodes--;
      }
   }   
   
   if (KeyWords(HP.GetWord()) != NODES) {
      DEBUGCERR("");
      std::cerr << "<end: nodes;> expected at line "
	<< HP.GetLineData() << std::endl;

      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      std::cerr << "semicolon expected at line " 
	      << HP.GetLineData() << std::endl;

      THROW(DataManager::ErrGeneric());
   }   
   
   if (iMissingNodes > 0) {
      DEBUGCERR("");
      std::cerr << "warning: " << iMissingNodes
	<< " nodes are missing;" << std::endl;
      for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
	 if (iNumTypes[iCnt] > 0) {
	    std::cerr << "  " << iNumTypes[iCnt] 
	      << ' ' << psNodeNames[iCnt] << std::endl;
	 }	 
      }      

      THROW(DataManager::ErrMissingNodes());
   }      
   
   DEBUGLCOUT(MYDEBUG_INPUT, "End of nodes data" << std::endl);
} /* End of DataManager::ReadNodes() */


void DataManager::ReadDrivers(MBDynParser& HP)
{
   DEBUGCOUTFNAME("DataManager::ReadDrivers");
     
   /* parole chiave del blocco di controllo */
   const char* sKeyWords[] = { 
      "end",
	"drivers",
	"file"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	
	END = 0,
	DRIVERS,
	FILEDRIVE,
	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   
   /* strutture di conteggio dei drivers letti */
   int iNumTypes[Drive::LASTDRIVETYPE];
   for (int i = 0; 
	i < Drive::LASTDRIVETYPE; 
	iNumTypes[i] = DriveData[i++].iNum) { 
      NO_OP; 
   }   
   
   int iMissingDrivers = iTotDrive;
   DEBUGLCOUT(MYDEBUG_INPUT, "Expected drivers: " << iMissingDrivers 
		   << std::endl);

   KeyWords CurrDesc;   
   while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {		
      /* puntatore al puntatore al driver */
      Drive** ppD = NULL;
      
      /* legge la label */	
      unsigned int uLabel = 0;
      if (CurrDesc >= 0) {	 
	 uLabel = unsigned(HP.GetInt());
      }      

      /* in base al tipo, avviene l'allocazione */
      switch (CurrDesc) {
	 
	 /* drivers */
       case FILEDRIVE: {	    
	  silent_cout("Reading file driver " << uLabel << std::endl);
	  
	  if (iNumTypes[Drive::FILEDRIVE]-- <= 0) {
	     DEBUGCERR("");
	     std::cerr << "line " << HP.GetLineData() 
	       << ": driver " << uLabel
	       << " exceedes file drivers number" << std::endl;	     
	     THROW(DataManager::ErrGeneric());
	  }	  
	  
	  /* allocazione e creazione */
	  int i = DriveData[Drive::FILEDRIVE].iNum
	    -iNumTypes[Drive::FILEDRIVE]-1;
	  ppD = DriveData[Drive::FILEDRIVE].ppFirstDrive+i;
	  
	  *ppD = ReadFileDriver(this, HP, uLabel);
	  HP.PutKeyTable(K);
	  break;
       }
	 	
	 /* aggiungere qui i nuovi tipi */
	 
	 /* in caso di tipo sconosciuto */
       default: {
	  DEBUGCERR("");
	  std::cerr << "unknown drive type at line " 
		  << HP.GetLineData() << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
      }  
      
      /* verifica dell'allocazione */
      if (*ppD == NULL) {
	 DEBUGCERR("");
	 std::cerr << "error in allocation, "
		 "item DataManager::*ppDrive, element " 
		 << uLabel << std::endl;
	 THROW(ErrMemory());
      }
      
      /* decrementa il totale degli elementi mancanti */
      iMissingDrivers--;
   }
   
   if (KeyWords(HP.GetWord()) != DRIVERS) {
      DEBUGCERR("");
      std::cerr << "<end: drivers;> expected at line " 
	<< HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      std::cerr << "semicolon expected at line " 
	      << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }   
   
   if (iMissingDrivers > 0) {
      std::cerr << "warning, " << iMissingDrivers
	<< " drivers are missing" << std::endl;
      THROW(DataManager::ErrGeneric());
   }      
   
   DEBUGLCOUT(MYDEBUG_INPUT, "End of drivers data" << std::endl);
} /* End of ReadDrivers */


/* Legge un legame costitutivo monodimensionale */
ConstitutiveLaw1D* 
DataManager::ReadConstLaw1D(MBDynParser& HP, DefHingeType::Type& T)
{
   return ReadConstLaw(this, HP, &DrvHdl, T, (ConstitutiveLaw1D*)NULL);
}


/* Legge un legame costitutivo tridimensionale */
ConstitutiveLaw3D* 
DataManager::ReadConstLaw3D(MBDynParser& HP, DefHingeType::Type& T)
{  
   return ReadConstLaw(this, HP, &DrvHdl, T, (ConstitutiveLaw3D*)NULL);
}


/* Legge un legame costitutivo esadimensionale */
ConstitutiveLaw6D* 
DataManager::ReadConstLaw6D(MBDynParser& HP, DefHingeType::Type& T)
{
   return ReadConstLaw(this, HP, &DrvHdl, T, (ConstitutiveLaw6D*)NULL);
}

/* DataManager - end */

Node* DataManager::ReadNode(MBDynParser& HP, Node::Type type)
{
	unsigned int uNode = (unsigned int)HP.GetInt();
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << std::endl);
	
	/* verifica di esistenza del nodo */
	Node* pNode;
	if ((pNode = pFindNode(type, uNode)) == NULL) {
		std::cerr << ": " << psNodeNames[type] << " node " << uNode
			<< " not defined at line " 
			<< HP.GetLineData() << std::endl;
		THROW(DataManager::ErrGeneric());
	}

	return pNode;
}

Elem* DataManager::ReadElem(MBDynParser& HP, Elem::Type type)
{
	unsigned int uElem = (unsigned int)HP.GetInt();
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Element " << uElem << std::endl);
	
	/* verifica di esistenza dell'elemento */
	Elem* pElem;

	if ((pElem = (Elem*)pFindElem(type, uElem)) == NULL) {
		std::cerr << ": " << psElemNames[type] << uElem
			<< " not defined at line " 
			<< HP.GetLineData() << std::endl;
		THROW(DataManager::ErrGeneric());
	}
	
	return pElem;
}


int GetDofOrder(MBDynParser& HP, Node* pNode, int iIndex)
{
   /*
    * Ordine del grado di liberta' da considerare 
    * (stato, oppure derivata se esiste)
    */
   if (HP.IsKeyWord("differential")) {
      /* 
       * Questo check e' pleonastico, viene gia' fatto dal chiamante
       */
      if (pNode->SetDof(iIndex-1) != DofOrder::DIFFERENTIAL) {
	 std::cerr << psNodeNames[pNode->GetNodeType()] 
		 << "(" << pNode->GetLabel() 
		 << "): invalid order for index " << iIndex 
      		 << " variable at line " << HP.GetLineData()
	   << std::endl;
	 THROW(DataManager::ErrGeneric());
      }
      return 1;
   } 

   if (HP.IsKeyWord("algebraic")) {
      return 0;

   } /* else */

   std::cerr << psNodeNames[pNode->GetNodeType()] 
	   << "(" << pNode->GetLabel() 
	   << "): unknown or illegal order for index " << iIndex << std::endl
      	   << "(hint: you may need to specify "
	   "\"differential\" or \"algebraic\" when referencing " << std::endl
      	   << "a generic degree of freedom at line " 
      	   << HP.GetLineData() << ")" << std::endl;

   THROW(DataManager::ErrGeneric());
#ifndef USE_EXCEPTIONS
   return 0;
#endif /* USE_EXCEPTIONS */
}


ScalarDof ReadScalarDof(const DataManager* pDM, MBDynParser& HP, flag fOrder)
{
   /* tabella delle parole chiave */
   KeyTable KDof((int)Node::LASTNODETYPE, psReadNodesNodes);
   HP.PutKeyTable(KDof);	     	        

   /* Label del nodo */
   unsigned int uNode = (unsigned int)HP.GetInt();
   DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << std::endl);
   
   /* Tipo del nodo */
   Node::Type Type = Node::Type(HP.GetWord());
   if (Type == Node::UNKNOWN) {
	   std::cerr << "line " << HP.GetLineData() 
		   << ": unknown node type" << std::endl;
      THROW(ErrGeneric());
   }   
   DEBUGLCOUT(MYDEBUG_INPUT, "Node type: " << psNodeNames[Type] << std::endl);
   
   /* verifica di esistenza del nodo */
   Node* pNode;
   if ((pNode = pDM->pFindNode(Type, uNode)) == NULL) {
      std::cerr << psNodeNames[Type] << "(" << uNode << ") not defined" 
	" at line " << HP.GetLineData() << std::endl;	 
      THROW(DataManager::ErrGeneric());
   }

   /* si procura il numero di dof del nodo */
   unsigned int iMaxIndex = pNode->iGetNumDof();
   DEBUGLCOUT(MYDEBUG_INPUT, "max index: " << iMaxIndex << std::endl);
   
   /* se il nodo ha piu' di un dof, chiede quale dof si desidera */
   unsigned int iIndex = 1;
   if (iMaxIndex > 1) {
      iIndex = HP.GetInt();
      if (iIndex > iMaxIndex) {
	 std::cerr << "Illegal index " << iIndex << ", " 
	   << psNodeNames[Type] << "(" << uNode << ") has only " 
	   << iMaxIndex << " dofs at line " 
	   << HP.GetLineData() << std::endl;
	 THROW(ErrGeneric());
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "index: " << iIndex << std::endl);
   }
   
   /* se e' richiesto l'order del dof e se il dof e' differenziale ... */
   int iOrder = 0;
   if (fOrder && pNode->iGetNumDof() > 0 && pNode->SetDof(iIndex-1) == DofOrder::DIFFERENTIAL) {
      iOrder = GetDofOrder(HP, pNode, iIndex);
   }
   DEBUGLCOUT(MYDEBUG_INPUT, "order: " << iOrder << std::endl);

   /* se il nodo non e' scalare, alloca un Node2Scalar che lo wrappa */
   if (iMaxIndex > 1) {
       NodeDof nd(pNode->GetLabel(), iIndex-1, pNode);
       
       pNode = NULL;
       /* Chi dealloca questa memoria? ci vorrebbe l'handle */	 
       SAFENEWWITHCONSTRUCTOR(pNode, Node2Scalar, Node2Scalar(nd));
       pedantic_cerr(psNodeNames[Type] << "(" << uNode 
	       << "): warning, possibly allocating a NodeDof "
	       "that nobody will delete until handles will be used"
	       "at line " << HP.GetLineData() << std::endl);
   } 
   
   return ScalarDof((ScalarNode*)pNode, iOrder);
}
   
   
   
   
   
   
  

/* Legge una shape1D; 
 * NOTA: il proprietario del puntatore alla Shape la deve distruggere */

#ifdef DEBUG_MEMMANAGER
#undef DEBUG_MEMMANAGER
#endif	   

#if (defined(USE_STRUCT_NODES) && defined(USE_AERODYNAMIC_ELEMS))
Shape* ReadShape(MBDynParser& HP)
{
   DEBUGCOUTFNAME("ReadShape");
   
   const char* sKeyWords[] = {
      "const",
	"linear",
	"piecewise" "linear",
	"parabolic",
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	CONST = 0,
	LINEAR,
	PIECEWISELINEAR,
	PARABOLIC,
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);   
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord;
   if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
      CurrKeyWord = CONST;
   }
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      std::cout << "shape type: " << sKeyWords[CurrKeyWord] << std::endl;
   }   
#endif /* DEBUG */

   Shape* pS = NULL;
   
   switch (CurrKeyWord) {
      
      /* forma costante */
    case CONST: {
       /* lettura dei dati specifici */
       doublereal dConst = HP.GetReal();
       DEBUGLCOUT(MYDEBUG_INPUT, "Const value: " << dConst << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pS,
			      ConstShape1D,
			      ConstShape1D(dConst));
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* forma lineare */
    case LINEAR: {
       /* lettura dei dati specifici */
       doublereal da0;
       doublereal da1;
       if (HP.IsKeyWord("coefficients")) {
          da0 = HP.GetReal();
	  da1 = HP.GetReal();
       } else {
	  doublereal dm = HP.GetReal();
	  doublereal dp = HP.GetReal();
          da0 = (dp+dm)/2.;
          da1 = (dp-dm)/2;
       }

       DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: " 
		       << da0 << ", " << da1 << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pS,
			      LinearShape1D,
			      LinearShape1D(da0, da1));
       
       /* scrittura dei dati specifici */	     
       
       break;
    }

      /* forma lineare a tratti (costante al di fuori del dominio definito) */
    case PIECEWISELINEAR: {
       int np = HP.GetInt();
       if (np <= 0) {
	  std::cerr << "Illegal number of points " << np 
            << " for piecewise linear shape at line " 
	    << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }

       doublereal *px = NULL;
       doublereal *pv = NULL;

       SAFENEWARR(px, doublereal, np);
       SAFENEWARR(pv, doublereal, np);

       px[0] = HP.GetReal();
       if (px[0] < -1. || px[0] > 1.) {
	  std::cerr << "Illegal value " << px[0] 
	    << "for first point abscissa (must be -1. < x < 1.) "
	    "in piecewise linear shape at line " 
	    << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }
       pv[0] = HP.GetReal();

       for (int i = 1; i < np; i++) {
	  px[i] = HP.GetReal();
	  if (px[i] <= px[i-1] || px[i] > 1.) {
	     std::cerr << "Illegal value " << px[i]
	       << "for point " << i+1 << " abscissa (must be " << px[i-1]
	       << " < x < 1.) in piecewise linear shape at line "
	       << HP.GetLineData() << std::endl;
	     THROW(ErrGeneric());
	  }
	  pv[i] = HP.GetReal();
       }

       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pS,
		       PiecewiseLinearShape1D,
		       PiecewiseLinearShape1D(np, px, pv));
  
       break;
    }

      /* forma lineare */
    case PARABOLIC: {
       /* lettura dei dati specifici */
       doublereal dm = HP.GetReal();
       doublereal da0 = HP.GetReal();
       doublereal dp = HP.GetReal();
       doublereal da1 = (dp-dm)/2.;
       doublereal da2 = (dp+dm)/2.-da0;
       DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: " 
		       << da0 << ", " << da1 << ", " << da2 << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pS,
			      ParabolicShape1D,
			      ParabolicShape1D(da0, da1, da2));
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* Non c'e' default in quanto all'inizio il default e' stato messo
       * pari a CONST */
    default: {
       ASSERTMSG(0, "You shouldn't have reached this point");
       THROW(ErrGeneric());       
    }
   }
   
   ASSERT(pS != NULL);
   return pS;   
} /* ReadShape */

#endif /* STRUCT && AERODYNAMIC */

