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

/* Elementi elettrici */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_ELECTRIC_NODES

#include <elec.h>
#include <elecnode.h>
#include <drive.h>
#include <strnode.h>
#include <accelerometer.h>
#include <displacement.h>
#include <motor.h>
#include <dataman.h>
#include <discctrl.h>

/* Electric - begin */

Electric::Electric(unsigned int uL, Electric::Type T, 
		   const DofOwner* pDO, flag fOut)
: Elem(uL, Elem::ELECTRIC, fOut), 
ElemWithDofs(uL, Elem::ELECTRIC, pDO, fOut), 
ElecT(T) 
{
   NO_OP; 
}


Electric::~Electric(void) 
{
   NO_OP;
}

   
/* Contributo al file di restart 
 * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
 * relativa alla classe derivata */
std::ostream& Electric::Restart(std::ostream& out) const {
   return out << "  electric: " << GetLabel();
}


/* Tipo dell'elemento (usato solo per debug ecc.) */
Elem::Type Electric::GetElemType(void) const 
{
   return Elem::ELECTRIC; 
}

/* Electric - end */


/* Legge un forgetting factor */

ForgettingFactor* ReadFF(MBDynParser& HP, integer iNumOutputs)
{
   ForgettingFactor* pFF = NULL;
   
   if (HP.IsKeyWord("forgettingfactor")) {
      if (HP.IsKeyWord("const")) {
	 doublereal d = HP.GetReal();
	 
	 SAFENEWWITHCONSTRUCTOR(pFF,
				ConstForgettingFactor,
				ConstForgettingFactor(d));
	 
      } else if (HP.IsKeyWord("dynamic")) {
	 /* uso la 2^a versione */
	 integer n1 = HP.GetInt();
	 integer n2 = HP.GetInt();
	 doublereal dRho = HP.GetReal();
	 doublereal dFact = HP.GetReal();
	 doublereal dKRef = HP.GetReal();
	 doublereal dKLim = HP.GetReal();
	 
	 SAFENEWWITHCONSTRUCTOR(pFF,
				DynamicForgettingFactor2,
				DynamicForgettingFactor2(n1, n2, 
							 iNumOutputs,
							 dRho, dFact,
							 dKRef, dKLim));
	 
      } else {	      
	 std::cerr << "line " << HP.GetLineData() 
	   << ": unknown forgetting factor" << std::endl;
	 THROW(ErrGeneric());
      }	      
   } else {
      /* default */
      SAFENEWWITHCONSTRUCTOR(pFF,
			     ConstForgettingFactor,
			     ConstForgettingFactor(1.));
   }
   
   ASSERT(pFF != NULL);
   return pFF;
}


/* Legge un eccitatore persistente */

PersistentExcitation* ReadPX(DataManager* pDM, MBDynParser& HP, integer iNumInputs)
{
   PersistentExcitation* pPX = NULL;
   
   if (HP.IsKeyWord("excitation")) {
      if (iNumInputs == 1) {
	 DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	 SAFENEWWITHCONSTRUCTOR(pPX, ScalarPX, ScalarPX(pDC));
      } else {
	 DriveCaller** ppDC = NULL;
	 SAFENEWARR(ppDC, DriveCaller*, iNumInputs);
	 
	 for (integer i = iNumInputs; i-- > 0; ) {
	    ppDC[i] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());	    
	 }
	 SAFENEWWITHCONSTRUCTOR(pPX, VectorPX, VectorPX(iNumInputs, ppDC));
      }
   } else {
      /* Null excitation */
      SAFENEW(pPX, NullPX);
   }
   
   ASSERT(pPX != NULL);
   return pPX;
}

                   


/* Legge un elemento elettrico */
   
Elem* ReadElectric(DataManager* pDM,
		   MBDynParser& HP, 
		   const DofOwner* pDO, 
		   unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadElectric()");
   
   const char* sKeyWords[] = {
      "accelerometer",
      "displacement",
      "motor",
      "discretecontrol",	
      "identification",
          "const",
          "dynamic",
      "control",
      "adaptivecontrol"	                    
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	
      ACCELEROMETER = 0,
      DISPLACEMENT,

      MOTOR,

      DISCRETECONTROL,	
      IDENTIFICATION,
          CONST,
          DYNAMIC,
      CONTROL,
      ADAPTIVECONTROL,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   /* lettura del tipo di elemento elettrico */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      std::cout << "electric element type: " 
	<< sKeyWords[CurrKeyWord] << std::endl;
   }   
#endif   

   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {
      /*  */
      
    case ACCELEROMETER: {
#if defined(USE_STRUCT_NODES)

       int f = 0;
       if (HP.IsKeyWord("translational")) {
	  f = 1;
       } else if(HP.IsKeyWord("rotational")) {
	  f = 2;
       }
       
       if (f) {
	  
	  /* nodo strutturale collegato */
	  unsigned int uNode = (unsigned int)HP.GetInt();	     
	  DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
	  
	  /* verifica di esistenza del nodo strutturale */
	  StructNode* pStrNode;
	  if ((pStrNode = pDM->pFindStructNode(uNode)) == NULL) {
	     std::cerr << "line " << HP.GetLineData() 
	       << ": structural node " << uNode
	       << " not defined" << std::endl;	  
	     THROW(DataManager::ErrGeneric());
	  }		  
	  
	  /* nodo astratto collegato */
	  uNode = (unsigned int)HP.GetInt();	     
	  DEBUGCOUT("Linked to Abstract Node " << uNode << std::endl);
	  
	  /* verifica di esistenza del nodo astratto */
	  AbstractNode* pAbsNode;
	  if ((pAbsNode = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	     std::cerr << "line " << HP.GetLineData() 
	       << ": abstract node " << uNode
	       << " not defined" << std::endl;	  
	     THROW(DataManager::ErrGeneric());
	  }		  
	  
	  /* Direzione */
	  Vec3 Dir(HP.GetVecRel(ReferenceFrame(pStrNode)));
	  doublereal d = Dir.Dot();
	  if (d > 0.) {
	     Dir /= d;
	     DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);
	  } else {
	     std::cerr << "Warning, null direction in accelerometer "
	       << uLabel << std::endl;
	  }
	  
	  switch (f) {
	   case 1: {
	      /* offset */
	      Vec3 Tmpf(HP.GetPosRel(ReferenceFrame(pStrNode)));
	      flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
	      
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     TraslAccel,
				     TraslAccel(uLabel, pDO, 
						pStrNode, pAbsNode,
						Dir, Tmpf, fOut));
	      break;
	   } 
	   case 2: {
	      flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
	      
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     RotAccel,
				     RotAccel(uLabel, pDO, 
					      pStrNode, pAbsNode,
					      Dir, fOut));
	      break;
	   }
	   default: {
	      std::cerr << "you shouldn't be here!" << std::endl;
	      THROW(ErrGeneric());
	   }
	  }	  
	  
       } else {
	  
	  /* nodo strutturale collegato */
	  unsigned int uNode = (unsigned int)HP.GetInt();	     
	  DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
	  
	  /* verifica di esistenza del nodo strutturale */
	  StructNode* pStrNode;
	  if ((pStrNode = pDM->pFindStructNode(uNode)) == NULL) {
	     std::cerr << "line " << HP.GetLineData() 
	       << ": structural node " << uNode
	       << " not defined" << std::endl;	  
	     THROW(DataManager::ErrGeneric());
	  }		  
	  
	  /* nodo astratto collegato */
	  uNode = (unsigned int)HP.GetInt();	     
	  DEBUGCOUT("Linked to Abstract Node " << uNode << std::endl);
	  
	  /* verifica di esistenza del nodo astratto */
	  AbstractNode* pAbsNode;
	  if ((pAbsNode = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	     std::cerr << "line " << HP.GetLineData() 
	       << ": abstract node " << uNode
	       << " not defined" << std::endl;	  
	     THROW(DataManager::ErrGeneric());
	  }		  
	  
	  /* Direzione */
	  Vec3 Dir(HP.GetVecRel(ReferenceFrame(pStrNode)));
	  doublereal d = Dir.Dot();
	  if (d > 0.) {
	     Dir /= d;
	     DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);
	  } else {
	     std::cerr << "Warning, null direction in accelerometer "
	       << uLabel << std::endl;
	  }
	  
	  /* Parametri */
	  doublereal dOmega = HP.GetReal();
	  if (dOmega <= 0.) {		  
	     std::cerr << "Warning, illegal Omega in accelerometer " 
	       << uLabel << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  doublereal dTau = HP.GetReal();
	  if (dTau <= 0.) {		  
	     std::cerr << "Warning, illegal Tau in accelerometer " 
	       << uLabel << "; aborting ..." << std::endl;	  
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  doublereal dCsi = HP.GetReal();
	  if (dCsi <= 0. || dCsi > 1.) {		  
	     std::cerr << "Warning, illegal Csi in accelerometer " 
	       << uLabel << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  doublereal dKappa = HP.GetReal();
	  if (dKappa == 0.) {		  
	     std::cerr << "Warning, null Kappa in accelerometer " 
	       << uLabel << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }	     
	  
	  DEBUGCOUT("Omega: " << dOmega 
		    << ", Tau: " << dTau
		    << ", Csi: " << dCsi
		    << ", Kappa: " << dKappa << std::endl);
	  
	  flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
	  
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 Accelerometer,
				 Accelerometer(uLabel, pDO, pStrNode, pAbsNode,
					       Dir, dOmega, dTau, dCsi, dKappa, 
					       fOut));
       }
       
       break;
       
#else /* USE_STRUCT_NODES */
       std::cerr << "you're not allowed to use accelerometer elements" << std::endl;
       THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
    }
      
    case DISPLACEMENT: {
#if defined(USE_STRUCT_NODES)
	  
       /* nodo strutturale collegato 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo strutturale */
       StructNode* pStrNode1;
       if ((pStrNode1 = pDM->pFindStructNode(uNode)) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": structural node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* offset 1 */
       Vec3 Tmpf1(HP.GetPosRel(ReferenceFrame(pStrNode1)));
	  
       /* nodo strutturale collegato 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo strutturale */
       StructNode* pStrNode2;
       if ((pStrNode2 = pDM->pFindStructNode(uNode)) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": structural node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
	  
       /* offset 2 */
       Vec3 Tmpf2(HP.GetPosRel(ReferenceFrame(pStrNode2)));
       
       /* nodo astratto collegato */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Abstract Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo astratto */
       AbstractNode* pAbsNode;
       if ((pAbsNode = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      DispMeasure,
			      DispMeasure(uLabel, pDO, 
					  pStrNode1, pStrNode2, pAbsNode,
					  Tmpf1, Tmpf2, fOut));
       break;
       
#else /* USE_STRUCT_NODES */
       std::cerr << "you're not allowed to use displacement measure elements" << std::endl;
       THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
    }

    case MOTOR: {
#if defined(USE_STRUCT_NODES)
	  
       /* nodo strutturale collegato 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo strutturale */
       StructNode* pStrNode1;
       pStrNode1 = pDM->pFindStructNode(uNode);
       if (pStrNode1 == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": structural node " << uNode
	    << " not defined" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }
	
       /* nodo strutturale collegato 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Structural Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo strutturale */
       StructNode* pStrNode2;
       pStrNode2 = pDM->pFindStructNode(uNode);
       if (pStrNode2 == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": structural node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  

       /* direzione */
       Vec3 TmpDir(HP.GetVecRel(ReferenceFrame(pStrNode1)));
       if (TmpDir.Norm() < DBL_EPSILON) {
	       std::cerr << "motor direction is illegal at line "
		       << HP.GetLineData() << std::endl;
	       THROW(ErrGeneric());
       }

       /* nodo astratto collegato */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Abstract Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo astratto */
       ElectricNode* pVoltage1;
       pVoltage1 = (ElectricNode *)pDM->pFindNode(Node::ELECTRIC, uNode);
       if (pVoltage1 == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo astratto collegato */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Abstract Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo astratto */
       ElectricNode* pVoltage2;
       pVoltage2 = (ElectricNode *)pDM->pFindNode(Node::ELECTRIC, uNode);
       if (pVoltage2 == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	      << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }

       doublereal dG = HP.GetReal();
       doublereal dl = HP.GetReal();
       doublereal dr = HP.GetReal();
       
       flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, Motor, Motor(uLabel, pDO, 
	       		       pStrNode1, pStrNode2, pVoltage1, pVoltage2,
			       TmpDir, dG, dl, dr, fOut));
       break;
       
#else /* USE_STRUCT_NODES */
       std::cerr << "you're not allowed to use electric motors" << std::endl;
       THROW(ErrGeneric());
#endif /* USE_STRUCT_NODES */
    }

      /*  */
    case DISCRETECONTROL: {
       /* lettura dei dati specifici */
       
       /* Dati generali del controllore */
       integer iNumOutputs = HP.GetInt();
       integer iNumInputs = HP.GetInt();
       integer iOrderA = HP.GetInt();
       integer iOrderB = iOrderA;
       if (HP.IsKeyWord("fir")) {	 
	  iOrderB = HP.GetInt();
       }
              
       integer iNumIter = HP.GetInt();

       DEBUGCOUT("Discrete controller of order " << iOrderA);
       if (iOrderB != iOrderA) {
	  DEBUGCOUT(" (fir order " << iOrderB << ')' << std::endl);
       }
       DEBUGCOUT(": " << iNumOutputs << " output(s) and " 
		 << iNumInputs << " input(s)" << std::endl
		 << "Update every " << iNumIter << " iterations" << std::endl);
       
       /* Tipo di controllo */
       DiscreteControlProcess* pDCP = NULL;
       switch (HP.GetWord()) {
	case CONTROL: {
	   /* Add the file with control data */
	   const char* sControlFile(HP.GetFileName());
	   
	   DEBUGCOUT("Getting control matrices from file <"
		     << sControlFile << '>' << std::endl);
	   
	   std::ifstream iFIn(sControlFile);
	   if (!iFIn) {
	      std::cerr << "Error in opening control file <" 
		<< sControlFile << '>' << std::endl;	      
	      THROW(DataManager::ErrGeneric());
	   }
	   
	   /* Construction of controller */
	   SAFENEWWITHCONSTRUCTOR(pDCP, 
				  DiscreteControlARXProcess_Debug,
				  DiscreteControlARXProcess_Debug(iNumOutputs, 
								  iNumInputs,
								  iOrderA,
								  iOrderB,
								  iFIn));
	   
	   iFIn.close();	   
	   break;
	}
	  
	case IDENTIFICATION: {
	   flag f_ma = 0;
	   if (HP.IsKeyWord("arx")) {
	      f_ma = 0;
	   } else if (HP.IsKeyWord("armax")) {
	      f_ma = 1;
	   }

	   /* Forgetting factor */
	   ForgettingFactor* pFF = ReadFF(HP, iNumOutputs);
	   
	   /* Persistent excitation */
	   PersistentExcitation* pPX = ReadPX(pDM, HP, iNumInputs);
	   HP.PutKeyTable(K);
	   
	   char* s = NULL;
	   if (HP.IsKeyWord("file")) {
	      s = (char*)HP.GetFileName();
	   }
	   	   
	   /* Construction of controller */
	   SAFENEWWITHCONSTRUCTOR(pDCP,
				  DiscreteIdentProcess_Debug,
				  DiscreteIdentProcess_Debug(iNumOutputs,
							     iNumInputs,
							     iOrderA,
							     iOrderB,
							     pFF, pPX,
							     f_ma, s));
	      
	   break;
	}	   
	   
	case ADAPTIVECONTROL:  {
#ifdef USE_DBC
	   flag f_ma = 0;
	   doublereal dPeriodicFactor(0.);
	   
	   if (HP.IsKeyWord("arx")) {
	      DEBUGCOUT("ARX adaptive control" << std::endl);
	      f_ma = 0;
	   } else if (HP.IsKeyWord("armax")) {
	      DEBUGCOUT("ARMAX adaptive control" << std::endl);
	      f_ma = 1;
	   }

	   if (HP.IsKeyWord("periodic")) {
	      dPeriodicFactor = HP.GetReal();
           }
	   
	   GPCDesigner* pCD = NULL;
	   if (HP.IsKeyWord("gpc")) {
	      DEBUGCOUT("GPC adaptive control" << std::endl);

	      integer iPredS = HP.GetInt();
	      integer iContrS = HP.GetInt();	 
	      integer iPredH = HP.GetInt();
	      integer iContrH = 0;
	      
	      DEBUGCOUT("prediction advancing horizon: " << iPredS << std::endl
			<< "prediction receding horizon: " << iPredH << std::endl
			<< "control advancing horizon: " << iContrS << std::endl
			<< "control receding horizon: " << iContrH << std::endl);
	      
	      if (iPredS < 0) {
		 std::cerr << "Prediction advancing horizon (" << iPredS 
		   << ") must be positive" << std::endl;
		 THROW(ErrGeneric());
	      }
	      if (iPredH < 0) {
		 std::cerr << "Prediction receding horizon (" << iPredH
		   << ") must be positive" << std::endl;
		 THROW(ErrGeneric());
	      }
	      if (iPredH >= iPredS) {
		 std::cerr << "Prediction receding horizon (" << iPredH 
		   << ") must be smaller than prediction advancing horizon ("
		   << iPredS << ")" << std::endl;
		 THROW(ErrGeneric());
	      }
	      if (iContrS < 0) {
		 std::cerr << "Control advancing horizon (" << iContrS
		   << ") must be positive" << std::endl;
		 THROW(ErrGeneric());
	      }
	      
	      doublereal* pW = NULL;
	      doublereal* pR = NULL;
	      SAFENEWARR(pW, doublereal, iPredS-iPredH);
	      SAFENEWARR(pR, doublereal, iContrS-iContrH);
	      
	      if (HP.IsKeyWord("predictionweights")) {
		 DEBUGCOUT("prediction weights:" << std::endl);
		 for (integer i = iPredS-iPredH; i-- > 0; ) {
		    pW[i] = HP.GetReal();
		    DEBUGCOUT("W[" << i+1 << "] = " << pW[i] << std::endl);
		 }
	      } else {
		 for (integer i = 0; i < iPredS-iPredH; i++) {
		    pW[i] = 1.;
		 }
	      }
	      
	      if (HP.IsKeyWord("controlweights")) {
		 DEBUGCOUT("control weights:" << std::endl);
		 for (integer i = iContrS-iContrH; i-- > 0; ) {
		    pR[i] = HP.GetReal();
		    DEBUGCOUT("R[" << i+1 << "] = " << pR[i] << std::endl);
		 }
	      } else {
		 for (integer i = 0; i < iContrS-iContrH; i++) {
		    pR[i] = 1.;
		 }
	      }
	      
	      DEBUGCOUT("Weight Drive:" << std::endl);
	      DriveCaller* pLambda = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	      
	      SAFENEWWITHCONSTRUCTOR(pCD,
				     GPC,
				     GPC(iNumOutputs, iNumInputs,
					 iOrderA, iOrderB, 
					 iPredS, iContrS, 
					 iPredH, iContrH,
					 pW, pR, pLambda,
					 dPeriodicFactor, f_ma));
	      
	   } else if (HP.IsKeyWord("deadbeat")) {
	      DEBUGCOUT("DeadBeat adaptive control" << std::endl);
	      
	      int iPredS = HP.GetInt();
	      int iContrS = HP.GetInt();
	      SAFENEWWITHCONSTRUCTOR(pCD,
				     DeadBeat,
				     DeadBeat(iNumOutputs, iNumInputs,
					      iOrderA, iOrderB, 
					      iPredS, iContrS,
					      dPeriodicFactor, f_ma));
	   }
	   
	   /* Forgetting factor */
	   DEBUGCOUT("Forgetting Factor:" << std::endl);
	   ForgettingFactor* pFF = ReadFF(HP, iNumOutputs);
	   
	   /* Persistent excitation */
	   DEBUGCOUT("Persistent Excitation:" << std::endl);
	   PersistentExcitation* pPX = ReadPX(pDM, HP, iNumInputs);
	   HP.PutKeyTable(K);
	   
	   DriveCaller* pTrig = NULL;
	   if (HP.IsKeyWord("trigger")) {	      
	      DEBUGCOUT("Trigger:" << std::endl);
	      pTrig = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	   } else {
	      SAFENEWWITHCONSTRUCTOR(pTrig, 
				     OneDriveCaller,
				     OneDriveCaller(pDM->pGetDrvHdl()));
	   }
	   
	   /* desired output */
	   DriveCaller** pvDesiredOut = NULL;
	   if (HP.IsKeyWord("desiredoutput")) {
	      DEBUGCOUT("Desired output:" << std::endl);
	      SAFENEWARR(pvDesiredOut, DriveCaller*, iNumOutputs);
	      
	      for (integer i = 0; i < iNumOutputs; i++) {
		 DEBUGCOUT("output[" << i+1 << "]:" << std::endl);
		 pvDesiredOut[i] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	      }
	   }
	   	   
	   char* s = NULL;
	   if (HP.IsKeyWord("file")) {
	      s = (char*)HP.GetFileName();
	      DEBUGCOUT("Identified matrices will be output in file <" << s << '>' << std::endl);
	   }
	   	   
	   /* Construction of controller */
	   ASSERT(f_ma == 0 || f_ma == 1);
	   SAFENEWWITHCONSTRUCTOR(pDCP,
				  DAC_Process_Debug,
				  DAC_Process_Debug(iNumOutputs,
						    iNumInputs,
						    iOrderA,
						    iOrderB,
						    pFF,
						    pCD,
						    pPX,
						    pTrig,
						    pvDesiredOut,
						    s,
						    f_ma));
	   
	   break;
#else /* !USE_DBC */
	      std::cerr << "GPC/deadbeat control is not available" << std::endl;
	      THROW(ErrGeneric());
#endif /* !USE_DBC */
	}	   

	  
	  
	default: {
	   std::cerr << "Sorry, not implemented yed" << std::endl;	   
	   THROW(ErrNotImplementedYet());
	}
       }
       
       
       
       if (!HP.IsKeyWord("outputs")) {
	  std::cerr << "Error, outputs expected at line " 
	    << HP.GetLineData() << std::endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       ScalarDof* pOutputs = NULL;
       SAFENEWARR(pOutputs, ScalarDof, iNumOutputs);
       DriveCaller** ppOutScaleFact = NULL;
       SAFENEWARR(ppOutScaleFact, DriveCaller*, iNumOutputs);
       
       /* Allocazione nodi e connessioni */
       for (int i = 0; i < iNumOutputs; i++) {
	  pOutputs[i] = ReadScalarDof(pDM, HP, 1);
	  if (HP.IsKeyWord("scale")) {
	     ppOutScaleFact[i] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	  } else {
	     ppOutScaleFact[i] = NULL;
	     SAFENEWWITHCONSTRUCTOR(ppOutScaleFact[i], 
				    OneDriveCaller, 
				    OneDriveCaller(pDM->pGetDrvHdl()));
	  }
       }
                            
       if (!HP.IsKeyWord("inputs")) {
	  std::cerr << "Error, inputs expected at line "
	    << HP.GetLineData() << std::endl;
	  THROW(DataManager::ErrGeneric());
       }	   
       
       /* Same thing for input nodes */
       ScalarDof* pInputs = NULL;
       SAFENEWARR(pInputs, ScalarDof, iNumInputs);
       
       /* Allocazione nodi e connessioni */
       for (int i = 0; i < iNumInputs; i++) {
	  pInputs[i] = ReadScalarDof(pDM, HP, 1);
	  if (pInputs[i].pNode->GetNodeType() ==  Node::PARAMETER) {
	     std::cerr << "Sorry, parameters are not allowed as input nodes" 
	       << std::endl;	     
	     THROW(DataManager::ErrGeneric());	      
	  }	      
       }
       
       HP.PutKeyTable(K);

       flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      DiscreteControlElem,
			      DiscreteControlElem(uLabel, pDO,
						  iNumOutputs,
						  pOutputs,
						  ppOutScaleFact,
						  iNumInputs,
						  pInputs,
						  pDCP,
						  iNumIter,
						  fOut));
       
       break;
    }

      /* Aggiungere altri elementi elettrici */
      
    default: {
       std::cerr << "unknown electric element type in electric element " << uLabel 
	 << " at line " << HP.GetLineData() << std::endl;       
       THROW(DataManager::ErrGeneric());
    }	
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      std::cerr << "semicolon expected at line " << HP.GetLineData() << std::endl;     
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* ReadElectric() */

#endif /* USE_ELECTRIC_NODES */

