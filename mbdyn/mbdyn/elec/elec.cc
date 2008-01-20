/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

Electric::Electric(unsigned int uL,
		   const DofOwner* pDO, flag fOut)
: Elem(uL, fOut),
ElemWithDofs(uL, pDO, fOut)
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
	 silent_cerr("line " << HP.GetLineData()
	   << ": unknown forgetting factor" << std::endl);
	 throw ErrGeneric();
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
	 DriveCaller* pDC = HP.GetDriveCaller();
	 SAFENEWWITHCONSTRUCTOR(pPX, ScalarPX, ScalarPX(pDC));
      } else {
	 DriveCaller** ppDC = NULL;
	 SAFENEWARR(ppDC, DriveCaller*, iNumInputs);

	 for (integer i = iNumInputs; i-- > 0; ) {
	    ppDC[i] = HP.GetDriveCaller();
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
      CONTROL,
      ADAPTIVECONTROL,

      LASTKEYWORD
   };

   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);

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
       int f = 0;
       if (HP.IsKeyWord("translational")) {
	  f = 1;
       } else if(HP.IsKeyWord("rotational")) {
	  f = 2;
       }

       if (f) {
	  /* TODO: check if downgradable to ScalarNode */

	  /* nodo strutturale collegato */
	  StructNode* pStrNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

	  /* nodo astratto collegato */
	  ScalarDifferentialNode* pAbsNode = dynamic_cast<ScalarDifferentialNode *>(pDM->ReadNode(HP, Node::ABSTRACT));

	  /* Direzione */
	  Vec3 Dir(HP.GetVecRel(ReferenceFrame(pStrNode)));
	  doublereal d = Dir.Dot();
	  if (d > 0.) {
	     Dir /= d;
	     DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);
	  } else {
	     silent_cerr("Warning, null direction in accelerometer "
	       << uLabel << std::endl);
	     throw ErrGeneric();
	  }

	  /* offset */
	  Vec3 Tmpf(0.);
	  if (HP.IsKeyWord("position") || HP.IsKeyWord("offset") || f == 1) {
	     Tmpf = HP.GetPosRel(ReferenceFrame(pStrNode));
	  }
	  flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
	  switch (f) {
	   case 1:
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     TraslAccel,
				     TraslAccel(uLabel, pDO,
						pStrNode, pAbsNode,
						Dir, Tmpf, fOut));
	      break;

	   case 2:
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     RotAccel,
				     RotAccel(uLabel, pDO,
					      pStrNode, pAbsNode,
					      Dir, fOut));
	      break;

	   default:
	      silent_cerr("you shouldn't be here!" << std::endl);
	      throw ErrGeneric();
	  }

       } else {

	  /* nodo strutturale collegato */
	  StructNode* pStrNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

	  /* nodo astratto collegato */
	  ScalarDifferentialNode* pAbsNode = dynamic_cast<ScalarDifferentialNode *>(pDM->ReadNode(HP, Node::ABSTRACT));

	  /* Direzione */
	  Vec3 Dir(HP.GetVecRel(ReferenceFrame(pStrNode)));
	  doublereal d = Dir.Dot();
	  if (d > 0.) {
	     Dir /= d;
	     DEBUGCOUT("Direction: " << std::endl << Dir << std::endl);
	  } else {
	     silent_cerr("Warning, null direction in accelerometer "
	       << uLabel << std::endl);
	     throw ErrGeneric();
	  }

	  /* Parametri */
	  doublereal dOmega = HP.GetReal();
	  if (dOmega <= 0.) {
	     silent_cerr("Warning, illegal Omega in accelerometer "
	       << uLabel << std::endl);
	     throw DataManager::ErrGeneric();
	  }

	  doublereal dTau = HP.GetReal();
	  if (dTau <= 0.) {
	     silent_cerr("Warning, illegal Tau in accelerometer "
	       << uLabel << "; aborting ..." << std::endl);
	     throw DataManager::ErrGeneric();
	  }

	  doublereal dCsi = HP.GetReal();
	  if (dCsi <= 0. || dCsi > 1.) {
	     silent_cerr("Warning, illegal Csi in accelerometer "
	       << uLabel << std::endl);
	     throw DataManager::ErrGeneric();
	  }

	  doublereal dKappa = HP.GetReal();
	  if (dKappa == 0.) {
	     silent_cerr("Warning, null Kappa in accelerometer "
	       << uLabel << std::endl);
	     throw DataManager::ErrGeneric();
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
    }

    case DISPLACEMENT: {
       /* nodo strutturale collegato 1 */
       StructNode* pStrNode1 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

       /* offset 1 */
       Vec3 Tmpf1(HP.GetPosRel(ReferenceFrame(pStrNode1)));

       /* nodo strutturale collegato 2 */
       StructNode* pStrNode2 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

       /* offset 2 */
       Vec3 Tmpf2(HP.GetPosRel(ReferenceFrame(pStrNode2)));

       /* nodo astratto collegato */
       ScalarDifferentialNode* pAbsNode = dynamic_cast<ScalarDifferentialNode *>(pDM->ReadNode(HP, Node::ABSTRACT));

       flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

       SAFENEWWITHCONSTRUCTOR(pEl,
			      DispMeasure,
			      DispMeasure(uLabel, pDO,
					  pStrNode1, pStrNode2, pAbsNode,
					  Tmpf1, Tmpf2, fOut));
       break;
    }

    case MOTOR: {
       /* nodo strutturale collegato 1 */
       StructNode* pStrNode1 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

       /* nodo strutturale collegato 2 */
       StructNode* pStrNode2 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

       /* direzione */
       Vec3 TmpDir(HP.GetVecRel(ReferenceFrame(pStrNode1)));
       if (TmpDir.Norm() < DBL_EPSILON) {
	       silent_cerr("motor direction is illegal at line "
		       << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }

       /* nodo elettrico1 collegato */
       ElectricNode* pVoltage1 = dynamic_cast<ElectricNode *>(pDM->ReadNode(HP, Node::ELECTRIC));

       /* nodo elettrico2 collegato */
       ElectricNode* pVoltage2 = dynamic_cast<ElectricNode *>(pDM->ReadNode(HP, Node::ELECTRIC));

       doublereal dG = HP.GetReal();
       doublereal dl = HP.GetReal();
       doublereal dr = HP.GetReal();

       flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);

       SAFENEWWITHCONSTRUCTOR(pEl, Motor, Motor(uLabel, pDO,
	       		       pStrNode1, pStrNode2, pVoltage1, pVoltage2,
			       TmpDir, dG, dl, dr, fOut));
       break;
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
	      silent_cerr("Error in opening control file <"
		<< sControlFile << '>' << std::endl);
	      throw DataManager::ErrGeneric();
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
		 silent_cerr("Prediction advancing horizon (" << iPredS
		   << ") must be positive" << std::endl);
		 throw ErrGeneric();
	      }
	      if (iPredH < 0) {
		 silent_cerr("Prediction receding horizon (" << iPredH
		   << ") must be positive" << std::endl);
		 throw ErrGeneric();
	      }
	      if (iPredH >= iPredS) {
		 silent_cerr("Prediction receding horizon (" << iPredH
		   << ") must be smaller than prediction advancing horizon ("
		   << iPredS << ")" << std::endl);
		 throw ErrGeneric();
	      }
	      if (iContrS < 0) {
		 silent_cerr("Control advancing horizon (" << iContrS
		   << ") must be positive" << std::endl);
		 throw ErrGeneric();
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
	      DriveCaller* pLambda = HP.GetDriveCaller();

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

	   DriveCaller* pTrig = NULL;
	   if (HP.IsKeyWord("trigger")) {
	      DEBUGCOUT("Trigger:" << std::endl);
	      pTrig = HP.GetDriveCaller();
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
		 pvDesiredOut[i] = HP.GetDriveCaller();
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
	      silent_cerr("GPC/deadbeat control is not available" << std::endl);
	      throw ErrGeneric();
#endif /* !USE_DBC */
	}



	default: {
	   silent_cerr("Sorry, not implemented yed" << std::endl);
	   throw ErrNotImplementedYet();
	}
       }



       if (!HP.IsKeyWord("outputs")) {
	  silent_cerr("Error, outputs expected at line "
	    << HP.GetLineData() << std::endl);

	  throw DataManager::ErrGeneric();
       }

       ScalarDof* pOutputs = NULL;
       SAFENEWARRNOFILL(pOutputs, ScalarDof, iNumOutputs);
       DriveCaller** ppOutScaleFact = NULL;
       SAFENEWARR(ppOutScaleFact, DriveCaller*, iNumOutputs);

       /* Allocazione nodi e connessioni */
       for (int i = 0; i < iNumOutputs; i++) {
	  pOutputs[i] = ReadScalarDof(pDM, HP, 1);
	  if (HP.IsKeyWord("scale")) {
	     ppOutScaleFact[i] = HP.GetDriveCaller();
	  } else {
	     ppOutScaleFact[i] = NULL;
	     SAFENEW(ppOutScaleFact[i], OneDriveCaller);
	  }
       }

       if (!HP.IsKeyWord("inputs")) {
	  silent_cerr("Error, inputs expected at line "
	    << HP.GetLineData() << std::endl);
	  throw DataManager::ErrGeneric();
       }

       /* Same thing for input nodes */
       ScalarDof* pInputs = NULL;
       SAFENEWARRNOFILL(pInputs, ScalarDof, iNumInputs);

       /* Allocazione nodi e connessioni */
       for (int i = 0; i < iNumInputs; i++) {
	  pInputs[i] = ReadScalarDof(pDM, HP, 1);
	  if (pInputs[i].pNode->GetNodeType() ==  Node::PARAMETER) {
	     silent_cerr("Sorry, parameters are not allowed as input nodes"
	       << std::endl);
	     throw DataManager::ErrGeneric();
	  }
       }

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
       silent_cerr("unknown electric element type in electric element " << uLabel
	 << " at line " << HP.GetLineData() << std::endl);
       throw DataManager::ErrGeneric();
    }
   }

   /* Se non c'e' il punto e virgola finale */
   if (HP.IsArg()) {
      silent_cerr("semicolon expected at line " << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric();
   }

   return pEl;
} /* ReadElectric() */

