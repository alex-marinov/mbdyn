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

/* Genel */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <genel_.h>
#include <genfilt.h>
#include <swashpl.h>
#include <rottrim.h>
#include <dataman.h>

/* genel - begin */

Genel::Genel(unsigned int uL,
	     GenelType::Type T, 
	     const DofOwner* pDO, 
	     flag fOut)
: Elem(uL, ElemType::GENEL, fOut), 
ElemWithDofs(uL, ElemType::GENEL, pDO, fOut),
GenelT(T)
{
   NO_OP;
}


Genel::~Genel(void)
{
   NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& Genel::Restart(ostream& out) const
{
   return out << "  genel: " << GetLabel();
}

/* Genel - end */


/* Legge un genel */

Elem* ReadGenel(DataManager* pDM, 
		MBDynParser& HP,
		const DofOwner* pDO, 
		unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadGenel()");
   
   const char* sKeyWords[] = {
      "swashplate",
	"rotortrim",
	"clamp",
	"distance",
	"spring",
	"springsupport",
	"crossspringsupport",
	"mass",
	"scalarfilter",
	"statespaceSISO",
	"statespaceMIMO"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	SWASHPLATE = 0,
	ROTORTRIM,
	CLAMP,
	DISTANCE,
	SPRING,
	SPRINGSUPPORT,
	CROSSSPRINGSUPPORT,
	MASS,
	SCALARFILTER,
	STATESPACESISO,
	STATESPACEMIMO,
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   /* lettura del tipo di vincolo */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {
      cout << "genel type: " << sKeyWords[CurrKeyWord] << endl;
   }
#endif   
   
   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {      
      /* genel piatto oscillante */
    case SWASHPLATE: {
       /* nodo Collettivo */
       unsigned int uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pCollIn;
       if ((pCollIn = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       // DriveCaller* pColl = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // HP.PutKeyTable(K);
       
       flag fCollLimits(0);
       doublereal dCollMax(0.);
       doublereal dCollMin(0.);
       if (HP.IsKeyWord("limits")) {
	  fCollLimits = flag(1);
	  dCollMin = HP.GetReal();
	  dCollMax = HP.GetReal();
       }
	   
       
       /* nodo Longitudinale */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pLongIn;
       if ((pLongIn = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       // DriveCaller* pLong = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // HP.PutKeyTable(K);
       
       flag fForeAftLimits(0);
       doublereal dForeAftMax(0.);
       doublereal dForeAftMin(0.);
       if (HP.IsKeyWord("limits")) {
	  fForeAftLimits = flag(1);
	  dForeAftMin = HP.GetReal();
	  dForeAftMax = HP.GetReal();
       }	   
       
       
       /* nodo Laterale */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pLatIn;
       if ((pLatIn = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  	   	  	   
       
       // DriveCaller* pLat = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // HP.PutKeyTable(K);
       
       flag fLatLimits(0);
       doublereal dLatMax(0.);
       doublereal dLatMin(0.);
       if (HP.IsKeyWord("limits")) {
	  fLatLimits = flag(1);
	  dLatMin = HP.GetReal();
	  dLatMax = HP.GetReal();
       }
       
       
       /* nodo collegato 1 */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode1;
       if ((pNode1 = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
              
       /* nodo collegato 2 */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode2;
       if ((pNode2 = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData()
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       
       /* nodo collegato 3 */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode3;
       if ((pNode3 = (AbstractNode*)(pDM->pFindNode(NodeType::ABSTRACT, uNode))) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       doublereal dDynCoef = 0.;
       if (HP.fIsArg()) {
	  dDynCoef = HP.GetReal(dDynCoef);
       }
       
       doublereal dCyclFact = 1.;
       if (HP.fIsArg()) {
	  dCyclFact = HP.GetReal(dCyclFact);
       }
       
       doublereal dCollFact = 1.;
       if (HP.fIsArg()) {
	  dCollFact = HP.GetReal(dCollFact);
       }
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      SwashPlate,
			      SwashPlate(uLabel, pDO,
					 pCollIn, // pColl, 
					 pLongIn, // pLong,
					 pLatIn,  // pLat,
					 pNode1, pNode2, pNode3,
					 dDynCoef,
					 dCyclFact, 
					 dCollFact,
					 fCollLimits,
					 dCollMin,
					 dCollMax,
					 fForeAftLimits,
					 dForeAftMin,
					 dForeAftMax,
					 fLatLimits,
					 dLatMin,
					 dLatMax,
					 fOut), 
			      DMmm);
       break;
    }
      

      
    case ROTORTRIM: {
#if defined(USE_AERODYNAMIC_ELEMS)
       unsigned int uL = HP.GetInt();
       Rotor* pRot = (Rotor*)(((Elem*)pDM->pFindElem(ElemType::ROTOR, uL))->pGet());
       if (pRot == NULL) {
	  cerr << "line " << HP.GetLineData() << ": can't find rotor "
	    << uL << endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDifferentialNode* pvNodes[3];
       uL = HP.GetInt();
       pvNodes[0] = (AbstractNode*)pDM->pFindNode(NodeType::ABSTRACT, uL);
       if (pvNodes[0] == NULL) {
	  cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << endl;
	  THROW(ErrGeneric());
       }
       uL = HP.GetInt();
       pvNodes[1] = (AbstractNode*)pDM->pFindNode(NodeType::ABSTRACT, uL);
       if (pvNodes[1] == NULL) {
	  cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << endl;
	  THROW(ErrGeneric());
       }
       uL = HP.GetInt();
       pvNodes[2] = (AbstractNode*)pDM->pFindNode(NodeType::ABSTRACT, uL);
       if (pvNodes[2] == NULL) {
	  cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << endl;
	  THROW(ErrGeneric());
       }
       
       DEBUGCOUT("Rotor trim " << uLabel 
		 << " linked to rotor " << pRot->GetLabel() << endl
		 << "abstract nodes: " 
		 << pvNodes[0]->GetLabel() << ", "
		 << pvNodes[1]->GetLabel() << ", "
		 << pvNodes[2]->GetLabel() << endl);
       
       DriveCaller* pvDrives[3];
       pvDrives[0] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       pvDrives[1] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       pvDrives[2] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       doublereal dSigma = HP.GetReal();
       DEBUGCOUT("Sigma: " << dSigma << endl);

       doublereal dGamma = HP.GetReal();
       DEBUGCOUT("Gamma: " << dGamma << endl);
       
       doublereal dP = HP.GetReal();
       DEBUGCOUT("P: " << dP << endl);
       
       doublereal dTau0 = HP.GetReal();
       DEBUGCOUT("Tau0: " << dTau0 << endl);

       doublereal dTau1 = HP.GetReal();
       DEBUGCOUT("Tau1: " << dTau1 << endl);
       
       doublereal dKappa0 = HP.GetReal();
       DEBUGCOUT("Kappa0: " << dKappa0 << endl);

       doublereal dKappa1 = HP.GetReal();
       DEBUGCOUT("Kappa1: " << dKappa1 << endl);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      RotorTrim,
			      RotorTrim(uLabel, pDO, pRot,
					pvNodes[0], pvNodes[1], pvNodes[2],
					pvDrives[0], pvDrives[1], pvDrives[2],
					dSigma, dGamma, dP, 
					dTau0, dTau1, dKappa0, dKappa1,
					0 /* fOut */ ),
			      DMmm);

#else // defined(USE_AERODYNAMIC_ELEMS)
       cerr << "can't use a rotor trim element without rotors" << endl;
       THROW(ErrGeneric());       
#endif // defined(USE_AERODYNAMIC_ELEMS)
       
       break;
    }
      
    case CLAMP: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel clamp" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // DEBUGCOUT("Stiffness: " << dK << endl);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelClamp,
			      GenelClamp(uLabel, pDO, pDC, SD, fOut),
			      DMmm);
       
       break;
    }
      
    case DISTANCE: {
       ScalarDof SD1 = ReadScalarDof(pDM, HP, 1);
       if (SD1.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel distance" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  

       ScalarDof SD2 = ReadScalarDof(pDM, HP, 1);
       if (SD2.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel distance" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // DEBUGCOUT("Stiffness: " << dK << endl);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelDistance,
			      GenelDistance(uLabel, pDO, pDC, SD1, SD2, fOut),
			      DMmm);
       
       
       break;
    }
      
    case SPRING: {
       ScalarDof SD1 = ReadScalarDof(pDM, HP, 1);
       if (SD1.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel springs" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       ScalarDof SD2 = ReadScalarDof(pDM, HP, 1);
       if (SD2.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel springs" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       if (ConstLawType != DefHingeType::ELASTIC) {
	  cerr << "Error at line " << HP.GetLineData() 
	    << ": elastic constitutive laws only are allowed" << endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelSpring,
			      GenelSpring(uLabel, pDO, pCL, SD1, SD2, fOut),
			      DMmm);
       
       
       break;
    }
      
    case SPRINGSUPPORT: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel spring supports" << endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if ((SD.iOrder != 0) 
		  || (SD.pNode->SetDof(0) != DofOrder::DIFFERENTIAL)) {
	  cerr << "Sorry, a spring support must be linked to the algebraic value of a differential node" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	         
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);

       switch (ConstLawType) {
	case DefHingeType::ELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelSpringSupport,	
				  GenelSpringSupport(uLabel, pDO, pCL, 
						     SD, fOut),
				  DMmm);
	   break;
	}
	case DefHingeType::VISCOUS:
	case DefHingeType::VISCOELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelSpringDamperSupport,
				  GenelSpringDamperSupport(uLabel, pDO, pCL, 
							   SD, fOut),
				  DMmm);
	   break;
	}
	default: {
	   cerr << "You shouldn't be here!" << endl;
	   THROW(DataManager::ErrGeneric());
	}	  
       }
       
       break;
    }      
      
    case CROSSSPRINGSUPPORT: {
       ScalarDof SDRow = ReadScalarDof(pDM, HP, 1);
       if (SDRow.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel spring supports" << endl;
	  THROW(DataManager::ErrGeneric());	      
       } 
       
       ScalarDof SDCol = ReadScalarDof(pDM, HP, 1);
       if ((SDCol.iOrder != 0) 
	   || (SDCol.pNode->SetDof(0) != DofOrder::DIFFERENTIAL)) {
	  cerr << "Sorry, a spring support must be linked to the algebraic value of a differential node" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	         
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);

       switch (ConstLawType) {
	case DefHingeType::ELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelCrossSpringSupport,	
				  GenelCrossSpringSupport(uLabel, pDO, pCL, 
							  SDRow, SDCol, fOut),
				  DMmm);
	   break;
	}
	case DefHingeType::VISCOUS:
	case DefHingeType::VISCOELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelCrossSpringDamperSupport,
				  GenelCrossSpringDamperSupport(uLabel, pDO, pCL, 
								SDRow, SDCol, fOut),
				  DMmm);
	   break;
	}
	default: {
	   cerr << "You shouldn't be here!" << endl;
	   THROW(DataManager::ErrGeneric());
	}	  
       }
       
       break;
    }      
      
    case MASS: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel mass" << endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  cerr << "Sorry, only differential dofs are allowed for genel mass" << endl;
	  THROW(DataManager::ErrGeneric());
       }

       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelMass,
			      GenelMass(uLabel, pDO, pDC, SD, fOut),
			      DMmm);
       
       
       break;
    }

#if 0
    case SCALARFILTER: {
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  NodeType::PARAMETER
	   || SD_u.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel scalar filter" << endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD_y.pNode->SetDof(0) != DofOrder::DIFFERENTIAL
		  || SD_u.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  cerr << "Sorry, only differential dofs are allowed for genel scalar filter" << endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       unsigned int na = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << na << " order denominator" << endl);
       doublereal* pdP = NULL;
       SAFENEWARR(pdP, doublereal, na+1, DMmm);
       
       pdP[0] = 1.;
       
       if (na > 0) {	
	  for (unsigned int iCnt = 1; iCnt <= na; iCnt++) {
	     pdP[iCnt] = HP.GetReal();
	     DEBUGCOUT("a(" << iCnt << ") = " << pdP[iCnt] << endl);
	  }
       }
	     
       unsigned int nb = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << nb << " order numerator" << endl);
       doublereal* pdTau = NULL;
       SAFENEWARR(pdTau, doublereal, nb+1, DMmm);
             
       for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	  pdTau[iCnt] = HP.GetReal();
	  DEBUGCOUT("b(" << iCnt << ") = " << pdTau[iCnt] << endl);
       }
            
       if (HP.IsKeyWord("gain")) {
	  doublereal gain = HP.GetReal();
	  DEBUGCOUT("Gain is: " << gain << endl);
	  for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	     pdTau[iCnt] *= gain;
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelFilter,
			      GenelFilter(uLabel, pDO, SD_y, SD_u,
					  na, nb, pdP, pdTau, fOut),
			      DMmm);
              
       break;
    }
#endif // 0
      
    case SCALARFILTER: {
       /* output */
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel scalar filter output" << endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       /* input */
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);

       /* ordine del denominatore */
       unsigned int na = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << na << " order denominator" << endl);
       doublereal* pdP = NULL;
       SAFENEWARR(pdP, doublereal, na+1, DMmm);
       pdP[0] = 1.;
       
       /* i coefficienti vengono letti a partire da a1, perche' il polinomio 
	* e' assunto monico (a0 = 1) */
       if (na > 0) {	
	  for (unsigned int iCnt = 1; iCnt <= na; iCnt++) {
	     pdP[iCnt] = HP.GetReal();
	     DEBUGCOUT("a(" << iCnt << ") = " << pdP[iCnt] << endl);
	  }
       }
	     
       /* ordine del numeratore */
       unsigned int nb = HP.GetInt();
       if (nb > na) {
	  cerr << "illegal (non proper) transfer function for scalar filter " 
	    << uLabel << " at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }

       DEBUGCOUT("GenelFilter " << uLabel << " has a " << nb << " order numerator" << endl);
       doublereal* pdTau = NULL;
       SAFENEWARR(pdTau, doublereal, nb+1, DMmm);
             
       /* i coefficienti vengono letti a partire da b0; e' possibile inserire 
	* il polinomio in forma monica e modificare direttamente il guadagno */
       for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	  pdTau[iCnt] = HP.GetReal();
	  DEBUGCOUT("b(" << iCnt << ") = " << pdTau[iCnt] << endl);
       }
            
       /* guadagno (opzionale) */
       if (HP.IsKeyWord("gain")) {
	  doublereal gain = HP.GetReal();
	  DEBUGCOUT("Gain is: " << gain << endl);
	  for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	     pdTau[iCnt] *= gain;
	  }
       }
       
       flag fState(0);
       if (HP.IsKeyWord("state")) {
	  if (HP.IsKeyWord("steady")) {
	     fState = 1;
	  } else {
	     cerr << "unknown state option at line " << HP.GetLineData() << endl;
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelFilterEq,
			      GenelFilterEq(uLabel, pDO, SD_y, SD_u,
					    na, nb, pdP, pdTau, 
					    fState, fOut),
			      DMmm);
              
       break;
    }
      
    case STATESPACESISO: {
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  NodeType::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for genel state space SISO output" << endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD_y.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  cerr << "Sorry, only differential dofs are allowed for genel state space SISO output" << endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       unsigned int Order = HP.GetInt();       
       DEBUGCOUT("State Space SISO " << uLabel << " is of order " << Order << endl);
       
       if (!HP.IsKeyWord("matrixA")) {
	  cerr << "matrix A expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdA = NULL;
       SAFENEWARR(pdA, doublereal, Order*Order, DMmm);
       doublereal* pd = pdA;
       for (unsigned int i = 0; i < Order*Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixB")) {
	  cerr << "matrix B expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdB = NULL;
       SAFENEWARR(pdB, doublereal, Order, DMmm);
       pd = pdB;
       for (unsigned int i = 0; i < Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixC")) {
	  cerr << "matrix C expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdC = NULL;
       SAFENEWARR(pdC, doublereal, Order, DMmm);
       pd = pdC;
       for (unsigned int i = 0; i < Order; i++) {
	  *pd++ = HP.GetReal();
       }

       doublereal dD = 0.;
       if (HP.IsKeyWord("matrixD")) {
	  dD = HP.GetReal();
       }
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelStateSpaceSISO,
			      GenelStateSpaceSISO(uLabel, pDO, SD_y, SD_u,
						  Order, 
						  pdA, pdB, pdC, dD, fOut),
			      DMmm);
              
       break;
    }
      
      
    case STATESPACEMIMO: {
       int iNumOutputs = HP.GetInt();
       if (iNumOutputs <= 0) {
	  cerr << "illegal number of outputs for state space MIMO " << uLabel
	    << " at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDof* pvSD_y = NULL;
       SAFENEWARR(pvSD_y, ScalarDof, iNumOutputs, DMmm);
       for (int i = 0; i < iNumOutputs; i++) {	  
	  pvSD_y[i] = ReadScalarDof(pDM, HP, 1);
	  if (pvSD_y[i].pNode->GetNodeType() ==  NodeType::PARAMETER) {
	     cerr << "line " << HP.GetLineData()
	       << ": sorry, parameters are not allowed for genel state space MIMO output" << endl;
	     THROW(DataManager::ErrGeneric());	      
	  } else if (pvSD_y[i].pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	     cerr << "line " << HP.GetLineData() 
	       << ": sorry, only differential dofs are allowed for genel state space MIMO output" << endl;
	     THROW(DataManager::ErrGeneric());
	  }
       }
       
       int iNumInputs = HP.GetInt();
       if (iNumInputs <= 0) {
	  cerr << "illegal number of inputs for state space MIMO " << uLabel
	    << " at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDof* pvSD_u = NULL;
       SAFENEWARR(pvSD_u, ScalarDof, iNumInputs, DMmm);
       for (int i = 0; i < iNumInputs; i++) {	  
	  pvSD_u[i] = ReadScalarDof(pDM, HP, 1);
       }
       
       unsigned int Order = HP.GetInt();       
       DEBUGCOUT("State Space MIMO " << uLabel << " is of order " << Order << endl);
       
       if (!HP.IsKeyWord("matrixA")) {
	  cerr << "matrix A expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdA = NULL;
       SAFENEWARR(pdA, doublereal, Order*Order, DMmm);
       doublereal* pd = pdA;
       for (unsigned int i = 0; i < Order*Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixB")) {
	  cerr << "matrix B expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdB = NULL;
       SAFENEWARR(pdB, doublereal, Order*iNumInputs, DMmm);
       pd = pdB;
       for (unsigned int i = 0; i < Order*iNumInputs; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixC")) {
	  cerr << "matrix C expected at line " << HP.GetLineNumber() << endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdC = NULL;
       SAFENEWARR(pdC, doublereal, iNumOutputs*Order, DMmm);
       pd = pdC;
       for (unsigned int i = 0; i < iNumOutputs*Order; i++) {
	  *pd++ = HP.GetReal();
       }

       doublereal* pdD = NULL;
       if (HP.IsKeyWord("matrixD")) {
	  SAFENEWARR(pdD, doublereal, iNumOutputs*iNumInputs, DMmm);
	  pd = pdD;
	  for (int i = 0; i < iNumOutputs*iNumInputs; i++) {
	     *pd++ = HP.GetReal();
	  }
       }       
       
       flag fOut = pDM->fReadOutput(HP, ElemType::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelStateSpaceMIMO,
			      GenelStateSpaceMIMO(uLabel, pDO, 
						  iNumOutputs, pvSD_y, 
						  iNumInputs, pvSD_u,
						  Order, 
						  pdA, pdB, pdC, pdD, fOut),
			      DMmm);

       break;
    }
      
      
      
      
      
      
      /* Aggiungere altri genel */
      
    default: {
       cerr << "unknown genel type in genel " << uLabel
	 << " at line " << HP.GetLineData() << endl;
       
       THROW(DataManager::ErrGeneric());
    }	
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;      
      THROW(DataManager::ErrGeneric());
   }      
   
   return pEl;
} /* ReadGenel() */
