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

/* Continua il DataManager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <sys/stat.h>
#include <unistd.h>

#include <dataman.h>
#include <dataman_.h>

#include <drive.h>
#include <presnode.h>
#include <readclaw.h>
#include <j2p.h>

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
      psReadControlElems[Elem::AERODYNAMIC],      
      psReadControlElems[Elem::FORCE],      
      psReadControlElems[Elem::GENEL],
      psReadControlElems[Elem::ELECTRICBULK],
      psReadControlElems[Elem::ELECTRIC],      
      psReadControlElems[Elem::HYDRAULIC],      
      psReadControlElems[Elem::BULK],
      psReadControlElems[Elem::LOADABLE],
      psReadControlDrivers[Drive::FILEDRIVE],       
      
      
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
      "adams" "res" "output",
      "default" "output",
      "all",	
      "none"
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
      AERODYNAMICELEMENTS,
      FORCES,
      GENELS,
      ELECTRICBULKELEMENTS,
      ELECTRICELEMENTS,
      HYDRAULICELEMENTS,
      BULKELEMENTS,
      LOADABLEELEMENTS,
      FILEDRIVERS,
      
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
      ADAMSRESOUTPUT,
      DEFAULTOUTPUT,
      ALL,
      NONE,
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
	  DEBUGLCOUT(MYDEBUG_INPUT, "Structural nodes: " << iDmy << endl);
	  break;
#else /* USE_STRUCT_NODES */
	  cerr << "you're not allowed to use structural nodes" << endl;
	  THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
       }
	 
	 /* Numero di nodi elettrici attesi */
       case ELECTRICNODES: {
#if defined(USE_ELECTRIC_NODES)	  
	  int iDmy = HP.GetInt();
	  NodeData[Node::ELECTRIC].iNum = iDmy;
	  DofData[DofOwner::ELECTRICNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Electric nodes: " << iDmy << endl);
	  break;
#else /* USE_ELECTRIC_NODES */
	  cerr << "you're not allowed to use electric nodes" << endl;
	  THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
       }	     
	 
	 /* Numero di nodi astratti attesi */
       case ABSTRACTNODES: {
#if defined(USE_ELECTRIC_NODES)	  
	  int iDmy = HP.GetInt();
	  NodeData[Node::ABSTRACT].iNum = iDmy;
	  DofData[DofOwner::ABSTRACTNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Abstract nodes: " << iDmy << endl);
	  break;
#else /* USE_ELECTRIC_NODES */
	  cerr << "you're not allowed to use abstract nodes" << endl;
	  THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
       }
	 
	 /* Numero di nodi astratti attesi */
       case PARAMETERNODES: {		  
	  int iDmy = HP.GetInt();
	  NodeData[Node::PARAMETER].iNum = iDmy;	     	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Parameter nodes: " << iDmy << endl);
	  break;
       }	     
	 
	 /* Numero di nodi idraulici attesi */
       case HYDRAULICNODES: {
#if defined(USE_HYDRAULIC_NODES)
	  int iDmy = HP.GetInt();
	  NodeData[Node::HYDRAULIC].iNum = iDmy;
	  DofData[DofOwner::HYDRAULICNODE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic nodes: " << iDmy << endl);
	  break;
#else /* defined(USE_HYDRAULIC_NODES) */
	  cerr << "you're not allowed to use hydraulic nodes" << endl;
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
		    << iDmy << endl);
#else
          iDmy = 0;
#endif
	  break;
       }
	 
	 /* Accelerazione di gravita' */
       case GRAVITY: {
	  if(ElemData[Elem::GRAVITY].iNum > 0) {
	     cerr
	       << "warning: gravity acceleration already defined;" << endl
	       << "only one definition will be considered" << endl;
	  }
	  ElemData[Elem::GRAVITY].iNum = 1;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Gravity acceleration expected in elements data" << endl);
	  break;
       }
	 
	 /* Numero di corpi rigidi attesi */
       case RIGIDBODIES: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BODY].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Rigid bodies: " << iDmy << endl);
	  break;
       }	     
	 
	 /* Numero di vincoli attesi */
       case JOINTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::JOINT].iNum = iDmy;	     
	  DofData[DofOwner::JOINT].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Joints: " << iDmy << endl);
	  if (iDmy > 0 ) {		       
	     fInitialJointAssemblyToBeMade = flag(1);
	  }
	  break;
       }	     
	 
	 /* Numero di travi attese */
       case BEAMS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BEAM].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Beams: " << iDmy << endl);
	  if (iDmy > 0 ) {		       
	     fInitialJointAssemblyToBeMade = flag(1);
	  }
	  break;
       }	     

#if defined(USE_AERODYNAMIC_ELEMS)
	 /* Elementi aerodinamici: proprieta' dell'aria */
       case AIRPROPERTIES: {
	  if (ElemData[Elem::AIRPROPERTIES].iNum > 0) {
	     cerr
	       << "warning: air properties already defined;" << endl
	       << "only one definition will be considered" << endl;
	  }
	  ElemData[Elem::AIRPROPERTIES].iNum = 1;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Air properties expected in elements data" << endl);
	  break;
       }
	 
	 /* Elementi aerodinamici: rotori */
       case ROTORS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::ROTOR].iNum = iDmy;	     
	  DofData[DofOwner::ROTOR].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Rotors: " << iDmy << endl);
	  break;
       }	     	     
	 
	 /* Elementi aerodinamici: vari elementi aerodinamici senza dof */
       case AERODYNAMICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::AERODYNAMIC].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements: " << iDmy << endl);
	  break;
       }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

	 /* Numero di forze e coppie attese */
       case FORCES: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::FORCE].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "Forces: " << iDmy << endl);
	  break;
       }	     
	 
#if defined(USE_ELECTRIC_NODES)
	 /* Numero di vincoli attesi */
       case GENELS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::GENEL].iNum = iDmy;	     
	  DofData[DofOwner::GENEL].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Genels: " << iDmy << endl);
	  break;
       }	     
	 
	 /* Numero di elementi elettrici attesi */
       case ELECTRICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::ELECTRIC].iNum = iDmy;	     
	  DofData[DofOwner::ELECTRIC].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Electric elements: " << iDmy << endl);
	  break;
       }	     
#endif /* USE_ELECTRIC_NODES */
	 
	 /* Numero di elementi idraulici attesi */
       case HYDRAULICELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::HYDRAULIC].iNum = iDmy;	     
	  DofData[DofOwner::HYDRAULIC].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic elements: " << iDmy << endl);
	  break;
       }	     
	 
	 /* Numero di elementi elettrici attesi */
       case BULKELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::BULK].iNum = iDmy;	     	    
	  DEBUGLCOUT(MYDEBUG_INPUT, "Bulk elements: " << iDmy << endl);
	  break;
       }	     
	 
#if defined(HAVE_LOADABLE)
	 /* Numero di elementi elettrici attesi */
       case LOADABLEELEMENTS: {
	  int iDmy = HP.GetInt();
	  ElemData[Elem::LOADABLE].iNum = iDmy;	     	    
	  DofData[DofOwner::LOADABLE].iNum = iDmy;
	  DEBUGLCOUT(MYDEBUG_INPUT, "Loadable elements: " << iDmy << endl);
	  break;
       }
#endif /* defined(HAVE_LOADABLE) */
	 	 
	 /* Numero di drivers attesi */
       case FILEDRIVERS: {
	  int iDmy = HP.GetInt();
	  DriveData[Drive::FILEDRIVE].iNum = iDmy;	     
	  DEBUGLCOUT(MYDEBUG_INPUT, "File drivers: " << iDmy << endl);
	  break;
       }	     
	 
	    
	 /********* Miscellaneous *********/

#if defined(USE_STRUCT_NODES)	    
	 /* Spegne il flag di assemblaggio iniziale; 
	  * di default viene eseguito solo se sono definiti vincoli */
       case SKIPINITIALJOINTASSEMBLY: {
	  fSkipInitialJointAssembly = flag(1);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Skipping initial joint assembly" << endl);
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
		 DEBUGLCOUT(MYDEBUG_INPUT, "Rigid bodies will be used in initial joint assembly" << endl);
		 break;
	      }
		
	      case GRAVITY: {			    
		 ElemData[Elem::GRAVITY].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Gravity will be used in initial joint assembly" << endl);
		 break;
	      }
		
	      case FORCES: {
		 ElemData[Elem::FORCE].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Forces will be used in initial joint assembly" << endl);
		 break;
	      }
		
		/* Lo lascio per backwards compatibility */
	      case BEAMS: {			    
#if 0
		 ElemData[Elem::BEAM].fToBeUsedInAssembly = flag(1);
#endif /* 0 */
		 DEBUGLCOUT(MYDEBUG_INPUT, "Beams are used in initial joint assembly by default" << endl);
		 break;
	      }			    

#if defined(USE_AERODYNAMIC_ELEMS)
	      case AERODYNAMICELEMENTS: {
		 ElemData[Elem::AERODYNAMIC].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements will be used in initial joint assembly" << endl);
		 
		 if (ElemData[Elem::AIRPROPERTIES].fToBeUsedInAssembly == flag(0)) {
		    ElemData[Elem::AIRPROPERTIES].fToBeUsedInAssembly = flag(1);
		 }
		 
		 break;
	      }
#endif /* USE_AERODYNAMIC_ELEMS */
		
#if defined(HAVE_LOADABLE)
	      case LOADABLEELEMENTS: {
		 ElemData[Elem::LOADABLE].fToBeUsedInAssembly = flag(1);
		 DEBUGLCOUT(MYDEBUG_INPUT, "Loadable Elements will be used in initial joint assembly" << endl);
		 break;
	      }
#endif /* defined(HAVE_LOADABLE) */
		
		
		
		/* Elemento non autorizzato */
	      default: {
		 cerr << "Element type at line "
		   << HP.GetLineData() 
		     << " is not allowed; aborting ..." << endl;
		 
		 THROW(DataManager::ErrElemNotAllowedInAssembly());
	      }		       
		
		/* Errore */
	      case UNKNOWN: {
		 cerr << "Unknown element type at line "
		   << HP.GetLineData() << "; aborting ..." << endl;
		 
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
		     << dInitialPositionStiffness << endl
		     << "Initial velocity stiffness: "
		     << dInitialVelocityStiffness << endl);
	  
	  break;
       }
	 
	 /* Omega solidale con il nodo o con il rif. globale */
       case OMEGAROTATES: {
	  if (HP.IsKeyWord("yes")) {
	     fOmegaRotates = flag(1);
	  } else if (HP.IsKeyWord("no")) {
	     fOmegaRotates = flag(0);
	  } else {
	     cerr << "Invalid option at line " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  break;
       }
	 
	 /* Tolleranza nell'assemblaggio iniziale; viene calcolata come:
	  * sqrt(sum(res^2)/(1.+sum(sol^2)) */
       case INITIALTOLERANCE: {
	  dInitialAssemblyToll = 
	    HP.GetReal(dDefaultInitialAssemblyToll);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Initial assembly tolerance: " 
		     << dInitialAssemblyToll << endl);
	  break;
       }	     
	 
	 /* Numero massimo di iterazioni nell'assemblaggio iniziale;
	  * di default ne e' consentita solo una, indice di condizioni 
	  * iniziali corrette */
       case MAXINITIALITERATIONS: {
	  iMaxInitialIterations = 
	    HP.GetInt(iDefaultMaxInitialIterations);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Max initial iterations: " 
		     << iMaxInitialIterations << endl);
	  break;
       }
#endif /* USE_STRUCT_NODES */

       case PRINT:
	  if (HP.IsKeyWord("dofstats")) {
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
		     << sSimulationTitle << '"' << endl);
	  break;
       }	     
	 
	 /* Crea il file di restart */
       case MAKERESTARTFILE: {
	  DEBUGLCOUT(MYDEBUG_INPUT, "Restart file will be generated " << endl);
	  if (HP.fIsArg()) {
	     if (HP.IsKeyWord("iterations")) {
		RestartEvery = ITERATIONS;
		iRestartIterations = HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT, "every " << iRestartIterations << " iterations" << endl);
	     } else if (HP.IsKeyWord("time")) {
		RestartEvery = TIME;
		dRestartTime = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "every " << dRestartTime << " time units" << endl);
	     } else {
		cerr << "Error: unrecognized restart option at line "
		  << HP.GetLineData() << endl;

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
	     cerr << "Null file name at line " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  SAFESTRDUP(sOutName, sTmp);
	  break;
       }
	 
       case ADAMSRESOUTPUT: {
	  fAdamsResOutput = 1;
	  if (HP.fIsArg()) {
	     if (sAdamsModelName != NULL) {
		SAFEDELETEARR(sAdamsModelName);
		sAdamsModelName = NULL;
	     }
	     
	     const char *tmp = HP.GetStringWithDelims();
	     SAFESTRDUP(sAdamsModelName, tmp);
	  } else {
	     SAFESTRDUP(sAdamsModelName, "mbdyn");
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

#if defined(HAVE_LOADABLE)
	      case LOADABLEELEMENTS: {
		 ElemData[Elem::LOADABLE].fDefaultOut = flag(1);
		 break;
	      }
#endif /* defined(HAVE_LOADABLE) */
		
		
	      case UNKNOWN: {
		 cerr << "warning: unknown output case at line " 
		   << HP.GetLineData() << endl;
		 ASSERT(0);		 
		 break;
	      }
		
	      default: {
		 cerr << "case " << sKeyWords[CurrDesc] << " at line "
		   << HP.GetLineData() << " is not allowed" << endl;
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
	  cerr << "unknown description at line " << HP.GetLineData() << endl;
	  ASSERT(0);
	  break;
       }
	 
       default: {
	  cerr << "case " << sKeyWords[CurrDesc] << " at line " 
	    << HP.GetLineData() << " is not allowed" << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
      }
   }
   
   if (KeyWords(HP.GetWord()) != CONTROLDATA) {	
      DEBUGCERR("");
      cerr << "<end: control data;> expected at line " 
	<< HP.GetLineData() << endl;
      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;
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

   DEBUGLCOUT(MYDEBUG_INPUT, "End of control data" << endl);
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
      cerr << "Unknown output flag for element \""
	<< psElemNames[t] << "\" at line " << HP.GetLineData() << endl;
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
      cerr << "Unknown output flag for node \"" 
	<< psNodeNames[t] << "\" at line " << HP.GetLineData() << endl;
      THROW(DataManager::ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return fDef;
#endif /* USE_EXCEPTIONS */
} /* End of DataManager::fReadOutput */


/* legge i nodi e li costruisce */

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
   DEBUGLCOUT(MYDEBUG_INPUT, "Expected nodes: " << iMissingNodes << endl);

   KeyWords CurrDesc;   
   while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
      
      if (CurrDesc == OUTPUT) {
	 DEBUGLCOUT(MYDEBUG_INPUT, "nodes to be output: ");
	 
	 Node::Type Typ;
	 switch (KeyWords(HP.GetWord())) {

	  case STRUCTURAL: {
#if defined(USE_STRUCT_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "structural" << endl);
	     Typ = Node::STRUCTURAL;
	     break;
#else /* USE_STRUCT_NODES */
	     cerr << "you're not allowed to use structural nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
	  }
	    
	  case ELECTRIC: {	
#if defined(USE_ELECTRIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "electric" << endl);
	     Typ = Node::ELECTRIC;	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     cerr << "you're not allowed to use electric nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	  case ABSTRACT: {	   
#if defined(USE_ELECTRIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "abstract" << endl);
	     Typ = Node::ABSTRACT;
	     break;
#else /* USE_ELECTRIC_NODES */
	     cerr << "you're not allowed to use abstract nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	  case PARAMETER: {	    
	     DEBUGLCOUT(MYDEBUG_INPUT, "parameter" << endl);
	     Typ = Node::PARAMETER;	     
	     break;
	  }

	  case HYDRAULIC: {
#if defined (USE_HYDRAULIC_NODES)
	     DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic" << endl);
	     Typ = Node::HYDRAULIC;	     
	     break;
#else /* defined (USE_HYDRAULIC_NODES) */
	     cerr << "you're not allowed to use hydraulic nodes" << endl;
	     THROW(ErrGeneric());
#endif /* defined (USE_HYDRAULIC_NODES) */
	  }
	    
	  default: {
	     cerr << "Error: unknown node type, cannot modify output" << endl;

	     THROW(DataManager::ErrUnknownNode());
	  }
	 }
	 
	 while (HP.fIsArg()) {
	    unsigned int uL = (unsigned int)HP.GetInt();	    	 
	    Node* pN = pFindNode(Typ, uL);
	    if (pN == NULL) {
	       cerr << "Error: node " << uL << ", type " << psNodeNames[Typ]
		 << "is not defined; output cannot be modified" << endl;
	    } else {
	       DEBUGLCOUT(MYDEBUG_INPUT, "node " << uL << endl);
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
	     silent_cout("Reading structural node " << uLabel << endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::STRUCTURAL]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": structural node " << uLabel 
		  << " exceedes structural nodes number" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindStructNode(uLabel) != NULL) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": structural node " << uLabel
		  << " already defined" << endl;

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
	     cerr << "you're not allowed to use structural nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
	  }
	    
	    /* nodi elettrici */
	  case ELECTRIC: {
#if defined(USE_ELECTRIC_NODES)	   
	     silent_cout("Reading electric node " << uLabel << endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::ELECTRIC]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": electric node " << uLabel
		  << " exceedes electric nodes number" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindElectricNode(uLabel) != NULL) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": electric node " << uLabel
		  << " already defined" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* Initial values */
	     doublereal dx = HP.GetReal();
	     doublereal dxp = HP.GetReal();
	     
	     flag fOut = fReadOutput(HP, Node::ELECTRIC);
	     
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::ELECTRIC].iNum
	       -iNumTypes[Node::ELECTRIC]-1;
	     ppN = NodeData[Node::ELECTRIC].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::ELECTRICNODE].pFirstDofOwner+i;
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    ElectricNode,
				    ElectricNode(uLabel, pDO, dx, dxp, fOut));
	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     cerr << "you're not allowed to use electric nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	    /* nodi astratti */
	  case ABSTRACT: {
#if defined(USE_ELECTRIC_NODES)	   
	     silent_cout("Reading abstract node " << uLabel << endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::ABSTRACT]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": abstract node " << uLabel
		  << " exceedes abstract nodes number" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindNode(Node::ABSTRACT, uLabel) != NULL) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": abstract node " << uLabel
		  << " already defined" << endl;

		THROW(DataManager::ErrGeneric());
	     }		  
	     
	     
	     /* lettura dei dati specifici */
	     doublereal dX(0.);
	     doublereal dXP(0.);
	     if (HP.fIsArg()) {
		dX = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x = " << dX
			   << " is supplied" << endl);
		if (HP.fIsArg()) {
		   dXP = HP.GetReal();
		   DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x' = " << dXP
			      << " is supplied" << endl);
		}		 
	     }		  
	     
	     flag fOut = fReadOutput(HP, Node::ABSTRACT);
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::ABSTRACT].iNum
	       -iNumTypes[Node::ABSTRACT]-1;
	     ppN = NodeData[Node::ABSTRACT].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::ABSTRACTNODE].pFirstDofOwner+i;
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    AbstractNode,
				    AbstractNode(uLabel, pDO, dX, dXP, fOut));
	     
	     break;
#else /* USE_ELECTRIC_NODES */
	     cerr << "you're not allowed to use abstract nodes" << endl;
	     THROW(ErrGeneric());
#endif /* USE_ELECTRIC_NODES */
	  }
	    
	    /* parametri */
	  case PARAMETER: {	     
	     silent_cout("Reading parameter " << uLabel << endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::PARAMETER]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": parameter " << uLabel
		  << " exceedes parameters number" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindNode(Node::PARAMETER, uLabel) != NULL) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": parameter " << uLabel
		  << " already defined" << endl;
		
		THROW(DataManager::ErrGeneric());
	     }		  
	     
	     /* bound a elemento */
	     if (HP.IsKeyWord("element")) {
		DEBUGLCOUT(MYDEBUG_INPUT, "parameter node " << uLabel 
			   << "is linked to an element" << endl);
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
		
	     /* strain gage */
	     } else if (HP.IsKeyWord("straingage")) {
		DEBUGLCOUT(MYDEBUG_INPUT, "parameter node " << uLabel 
			   << "is a strain gage" << endl);
		
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
			      << " is supplied" << endl);
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
	     silent_cout("Reading hydraulic node " << uLabel << endl);
	     
	     /* verifica che non siano gia' stati letti tutti 
	      * quelli previsti */
	     if (iNumTypes[Node::HYDRAULIC]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": hydraulic node " << uLabel
		  << " exceedes hydraulic nodes number" << endl;

		THROW(DataManager::ErrGeneric());
	     }
	     
	     /* verifica di esistenza del nodo */
	     if (pFindNode(Node::HYDRAULIC, uLabel) != NULL) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": hydraulic node " << uLabel
		  << " already defined" << endl;

		THROW(DataManager::ErrGeneric());
	     }		  
	     
	     
	     /* lettura dei dati specifici */
	     doublereal dX(0.);	  
	     if (HP.fIsArg()) {
		dX = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x = " << dX
			  << " is supplied" << endl);
	     }		  
	     
	     flag fOut = fReadOutput(HP, Node::HYDRAULIC);
	     
	     /* allocazione e creazione */
	     int i = NodeData[Node::HYDRAULIC].iNum
	       -iNumTypes[Node::HYDRAULIC]-1;
	     ppN = NodeData[Node::HYDRAULIC].ppFirstNode+i;
	     DofOwner* pDO = DofData[DofOwner::HYDRAULICNODE].pFirstDofOwner+i;
	     
	     SAFENEWWITHCONSTRUCTOR(*ppN, 
				    PressureNode,
				    PressureNode(uLabel, pDO, dX, fOut));
	     
	     break;
	  }
#endif /* USE_HYDRAULIC_NODES */
	    
	    
	    /* aggiungere eventuali nuovi tipi di nodo */
	    
	    
	    /* in caso di tipo errato */
	  case UNKNOWN: {
	     DEBUGCERR("");
	     cerr << "unknown node type at line " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	    
	  default: {
	     DEBUGCERR("");
	     cerr << "node type " << sKeyWords[CurrDesc] 
	       << " at line " << HP.GetLineData() << " is not allowed" 
	       << endl;

	     THROW(DataManager::ErrGeneric());
	  }
	 }
	 
	 /* verifica dell'allocazione - comune a tutti i casi */
	 if (*ppN == NULL) {
	    DEBUGCERR("");
	    cerr << "error in allocation, item DataManager::*ppNodes, node "
	      << uLabel << endl;

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
      cerr << "<end: nodes;> expected at line "
	<< HP.GetLineData() << endl;

      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;

      THROW(DataManager::ErrGeneric());
   }   
   
   if (iMissingNodes > 0) {
      DEBUGCERR("");
      cerr << "warning: " << iMissingNodes
	<< " nodes are missing;" << endl;
      for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
	 if (iNumTypes[iCnt] > 0) {
	    cerr << "  " << iNumTypes[iCnt] 
	      << ' ' << psNodeNames[iCnt] << endl;
	 }	 
      }      

      THROW(DataManager::ErrMissingNodes());
   }      
   
   DEBUGLCOUT(MYDEBUG_INPUT, "End of nodes data" << endl);
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
   DEBUGLCOUT(MYDEBUG_INPUT, "Expected drivers: " << iMissingDrivers << endl);

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
	  silent_cout("Reading file driver " << uLabel << endl);
	  
	  if (iNumTypes[Drive::FILEDRIVE]-- <= 0) {
	     DEBUGCERR("");
	     cerr << "line " << HP.GetLineData() 
	       << ": driver " << uLabel
	       << " exceedes file drivers number" << endl;	     
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
	  cerr << "unknown drive type at line " << HP.GetLineData() << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
      }  
      
      /* verifica dell'allocazione */
      if (*ppD == NULL) {
	 DEBUGCERR("");
	 cerr << "error in allocation, item DataManager::*ppDrive, element "
	   << uLabel << endl;
	 THROW(ErrMemory());
      }
      
      /* decrementa il totale degli elementi mancanti */
      iMissingDrivers--;
   }
   
   if (KeyWords(HP.GetWord()) != DRIVERS) {
      DEBUGCERR("");
      cerr << "<end: drivers;> expected at line " 
	<< HP.GetLineData() << endl;
      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;
      THROW(DataManager::ErrGeneric());
   }   
   
   if (iMissingDrivers > 0) {
      cerr << endl << "warning, " << iMissingDrivers
	<< " drivers are missing" << endl;
      THROW(DataManager::ErrGeneric());
   }      
   
   DEBUGLCOUT(MYDEBUG_INPUT, "End of drivers data" << endl);
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
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << endl);
	
	/* verifica di esistenza del nodo */
	Node* pNode;
	if ((pNode = pFindNode(type, uNode)) == NULL) {
		cerr << ": " << psNodeNames[type] << " node " << uNode
			<< " not defined at line " << HP.GetLineData() << endl;
		THROW(DataManager::ErrGeneric());
	}

	return pNode;
}

Elem* DataManager::ReadElem(MBDynParser& HP, Elem::Type type)
{
	unsigned int uElem = (unsigned int)HP.GetInt();
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Element " << uElem << endl);
	
	/* verifica di esistenza dell'elemento */
	Elem* pElem;

	if ((pElem = (Elem*)pFindElem(type, uElem)) == NULL) {
		cerr << ": " << psElemNames[type] << uElem
			<< " not defined at line " << HP.GetLineData() << endl;
		THROW(DataManager::ErrGeneric());
	}
	
	return pElem;
}


int GetDofOrder(MBDynParser& HP, Node* pNode, int iIndex)
{
   /* Ordine del grado di liberta' da considerare 
    * (stato, oppure derivata se esiste) */
   if (HP.IsKeyWord("differential")) {
      /* Questo check e' pleonastico, viene gia' fatto dal chiamante */
      if (pNode->SetDof(iIndex-1) != DofOrder::DIFFERENTIAL) {
	 cerr << endl << sDMClassName 
	   << " at line " << HP.GetLineData()
	     << ": invalid order for index " << iIndex 
	   << " variable" << endl;
	 THROW(DataManager::ErrGeneric());
      }
      return 1;
   } 
   if (HP.IsKeyWord("algebraic")) {
#ifdef DEBUG
      if ((pNode->SetDof(iIndex-1) != DofOrder::DIFFERENTIAL)
	  && (pNode->SetDof(iIndex-1) != DofOrder::ALGEBRAIC)) {
	 cerr << endl << sDMClassName 
	   << " at line " << HP.GetLineData()
	     << ": invalid order for index " << iIndex 
	   << " variable" << endl;
	 THROW(DataManager::ErrGeneric());
      }	     
#endif      
      return 0;
   } /* else */
   cerr << endl << sDMClassName
     << " at line " << HP.GetLineData()
       << ": unknown or illegal order for index " << iIndex
     << " variable" << endl;
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
   DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << endl);
   
   /* Tipo del nodo */
   Node::Type Type = Node::Type(HP.GetWord());
   if (Type == Node::UNKNOWN) {
      cerr << "line " << HP.GetLineData() << ": unknown node type" << endl;
      THROW(ErrGeneric());
   }   
   DEBUGLCOUT(MYDEBUG_INPUT, "Node type: " << psNodeNames[Type] << endl);
   
   /* verifica di esistenza del nodo */
   Node* pNode;
   if ((pNode = pDM->pFindNode(Type, uNode)) == NULL) {
      cerr << endl << sDMClassName
	<< " at line " << HP.GetLineData() 
	<< ": node " << uNode
	<< " not defined" << endl;	 
      THROW(DataManager::ErrGeneric());
   }

   /* si procura il numero di dof del nodo */
   unsigned int iMaxIndex = pNode->iGetNumDof();
   DEBUGLCOUT(MYDEBUG_INPUT, "max index: " << iMaxIndex << endl);
   
   /* se il nodo ha piu' di un dof, chiede quale dof si desidera */
   unsigned int iIndex = 1;
   if (iMaxIndex > 1) {
      iIndex = HP.GetInt();
      if (iIndex > iMaxIndex) {
	 cerr << "Illegal index " << iIndex << ", " 
	   << psNodeNames[Type] << "(" << uNode << ") has only " 
	   << iMaxIndex << " dofs" << endl;
	 THROW(ErrGeneric());
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "index: " << iIndex << endl);
   }
   
   /* se e' richiesto l'order del dof e se il dof e' differenziale ... */
   int iOrder = 0;
   if (fOrder && pNode->iGetNumDof() > 0 && pNode->SetDof(iIndex-1) == DofOrder::DIFFERENTIAL) {
      iOrder = GetDofOrder(HP, pNode, iIndex);
   }
   DEBUGLCOUT(MYDEBUG_INPUT, "order: " << iOrder << endl);

   /* se il nodo non e' scalare, alloca un Node2Scalar che lo wrappa */
   if (iMaxIndex > 1) {
       NodeDof nd(pNode->GetLabel(), iIndex-1, pNode);
       
       pNode = NULL;
       /* Chi dealloca questa memoria? ci vorrebbe l'handle */	 
       SAFENEWWITHCONSTRUCTOR(pNode, Node2Scalar, Node2Scalar(nd));
       cerr << sDMClassName
	 << ": warning, possibly allocating a NodeDof that nobody will delete until Handles will be used"
	 << endl;
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
   if(CurrKeyWord >= 0) {      
      cout << "shape type: " << sKeyWords[CurrKeyWord] << endl;
   }   
#endif /* DEBUG */

   Shape* pS = NULL;
   
   switch (CurrKeyWord) {
      
      /* forma costante */
    case CONST: {
       /* lettura dei dati specifici */
       doublereal dConst = HP.GetReal();
       DEBUGLCOUT(MYDEBUG_INPUT, "Const value: " << dConst << endl);
       
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
	  if (dm < 0 || dp) {
             cerr << "Illegal chord boundaries {" << dm << "," << dp
	       << "} for linear shape at line " << HP.GetLineData() << endl;
	     THROW(ErrGeneric());
	  }
          da0 = (dp+dm)/2.;
          da1 = (dp-dm)/2;
       }

       DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: " << da0 << ", " << da1 << endl);
       
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
	  cerr << "Illegal number of points " << np 
            << " for piecewise linear shape at line " 
	    << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }

       doublereal *px = NULL;
       doublereal *pv = NULL;

       SAFENEWARR(px, doublereal, np);
       SAFENEWARR(pv, doublereal, np);

       px[0] = HP.GetReal();
       if (px[0] < -1. || px[0] > 1.) {
	  cerr << "Illegal value " << px[0] 
	    << "for first point abscissa (must be -1. < x < 1.) "
	    "in piecewise linear shape at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }
       pv[0] = HP.GetReal();

       for (int i = 1; i < np; i++) {
	  px[i] = HP.GetReal();
	  if (px[i] <= px[i-1] || px[i] > 1.) {
             cerr << "Illegal value " << px[i]
	       << "for point " << i+1 << " abscissa (must be " << px[i-1]
	       << " < x < 1.) in piecewise linear shape at line "
	       << HP.GetLineData() << endl;
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
       if (dm < 0. || da0 < 0. || dp < 0.) {
             cerr << "Illegal chord boundaries {" 
	       << dm << "," << da0 << "," << dp
	       << "} for parabolic shape at line " << HP.GetLineData() << endl;
	     THROW(ErrGeneric());
       }
       doublereal da1 = (dp-dm)/2.;
       doublereal da2 = (dp+dm)/2.-da0;
       DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: " << da0 << ", " << da1 << ", " << da2 <<endl);
       
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

