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

/* Genel */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_ELECTRIC_NODES

#include <genel_.h>
#include <genfilt.h>
#include <swashpl.h>
#ifdef USE_AERODYNAMIC_ELEMS
#include <rottrim.h>
#endif /* USE_AERODYNAMIC_ELEMS */
#include <dataman.h>

/* genel - begin */

Genel::Genel(unsigned int uL,
	     Genel::Type T, 
	     const DofOwner* pDO, 
	     flag fOut)
: Elem(uL, Elem::GENEL, fOut), 
ElemWithDofs(uL, Elem::GENEL, pDO, fOut),
GenelT(T)
{
   NO_OP;
}


Genel::~Genel(void)
{
   NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream& Genel::Restart(std::ostream& out) const
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
      std::cout << "genel type: " << sKeyWords[CurrKeyWord] << std::endl;
   }
#endif   
   
   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {      
      /* genel piatto oscillante */
    case SWASHPLATE: {
       /* nodo Collettivo */
       unsigned int uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pCollIn;
       if ((pCollIn = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
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
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pLongIn;
       if ((pLongIn = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
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
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pLatIn;
       if ((pLatIn = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
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
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode1;
       if ((pNode1 = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
              
       /* nodo collegato 2 */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode2;
       if ((pNode2 = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData()
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       
       /* nodo collegato 3 */
       uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << std::endl);
       
       /* verifica di esistenza del nodo */
       AbstractNode* pNode3;
       if ((pNode3 = (AbstractNode*)(pDM->pFindNode(Node::ABSTRACT, uNode))) == NULL) {
	  std::cerr << "line " << HP.GetLineData() 
	    << ": abstract node " << uNode
	    << " not defined" << std::endl;	  
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
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
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
					 fOut));
       break;
    }
      

      
    case ROTORTRIM: {
#ifdef USE_AERODYNAMIC_ELEMS
       unsigned int uL = HP.GetInt();
       Rotor* pRot = (Rotor*)(((Elem*)pDM->pFindElem(Elem::ROTOR, uL))->pGet());
       if (pRot == NULL) {
	  std::cerr << "line " << HP.GetLineData() << ": can't find rotor "
	    << uL << std::endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDifferentialNode* pvNodes[3];
       uL = HP.GetInt();
       pvNodes[0] = (AbstractNode*)pDM->pFindNode(Node::ABSTRACT, uL);
       if (pvNodes[0] == NULL) {
	  std::cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << std::endl;
	  THROW(ErrGeneric());
       }
       uL = HP.GetInt();
       pvNodes[1] = (AbstractNode*)pDM->pFindNode(Node::ABSTRACT, uL);
       if (pvNodes[1] == NULL) {
	  std::cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << std::endl;
	  THROW(ErrGeneric());
       }
       uL = HP.GetInt();
       pvNodes[2] = (AbstractNode*)pDM->pFindNode(Node::ABSTRACT, uL);
       if (pvNodes[2] == NULL) {
	  std::cerr << "line " << HP.GetLineData() << ": can't find abstract node "
	    << uL << std::endl;
	  THROW(ErrGeneric());
       }
       
       DEBUGCOUT("Rotor trim " << uLabel 
		 << " linked to rotor " << pRot->GetLabel() << std::endl
		 << "abstract nodes: " 
		 << pvNodes[0]->GetLabel() << ", "
		 << pvNodes[1]->GetLabel() << ", "
		 << pvNodes[2]->GetLabel() << std::endl);
       
       DriveCaller* pvDrives[3];
       pvDrives[0] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       pvDrives[1] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       pvDrives[2] = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       doublereal dSigma = HP.GetReal();
       DEBUGCOUT("Sigma: " << dSigma << std::endl);

       doublereal dGamma = HP.GetReal();
       DEBUGCOUT("Gamma: " << dGamma << std::endl);
       
       doublereal dP = HP.GetReal();
       DEBUGCOUT("P: " << dP << std::endl);
       
       doublereal dTau0 = HP.GetReal();
       DEBUGCOUT("Tau0: " << dTau0 << std::endl);

       doublereal dTau1 = HP.GetReal();
       DEBUGCOUT("Tau1: " << dTau1 << std::endl);
       
       doublereal dKappa0 = HP.GetReal();
       DEBUGCOUT("Kappa0: " << dKappa0 << std::endl);

       doublereal dKappa1 = HP.GetReal();
       DEBUGCOUT("Kappa1: " << dKappa1 << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      RotorTrim,
			      RotorTrim(uLabel, pDO, pRot,
					pvNodes[0], pvNodes[1], pvNodes[2],
					pvDrives[0], pvDrives[1], pvDrives[2],
					dSigma, dGamma, dP, 
					dTau0, dTau1, dKappa0, dKappa1,
					0 /* fOut */ ));

#else /* !USE_AERODYNAMIC_ELEMS */
       std::cerr << "can't use a rotor trim element without rotors" << std::endl;
       THROW(ErrGeneric());       
#endif /* !USE_AERODYNAMIC_ELEMS */
       
       break;
    }
      
    case CLAMP: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel clamp" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // DEBUGCOUT("Stiffness: " << dK << std::endl);
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelClamp,
			      GenelClamp(uLabel, pDO, pDC, SD, fOut));
       
       break;
    }
      
    case DISTANCE: {
       ScalarDof SD1 = ReadScalarDof(pDM, HP, 1);
       if (SD1.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel distance" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  

       ScalarDof SD2 = ReadScalarDof(pDM, HP, 1);
       if (SD2.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel distance" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       // DEBUGCOUT("Stiffness: " << dK << std::endl);
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelDistance,
			      GenelDistance(uLabel, pDO, pDC, SD1, SD2, fOut));
       
       
       break;
    }
      
    case SPRING: {
       ScalarDof SD1 = ReadScalarDof(pDM, HP, 1);
       if (SD1.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel springs" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       ScalarDof SD2 = ReadScalarDof(pDM, HP, 1);
       if (SD2.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel springs" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       if (ConstLawType != DefHingeType::ELASTIC) {
	  std::cerr << "Error at line " << HP.GetLineData() 
	    << ": elastic constitutive laws only are allowed" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelSpring,
			      GenelSpring(uLabel, pDO, pCL, SD1, SD2, fOut));
       
       
       break;
    }
      
    case SPRINGSUPPORT: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel spring supports" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if ((SD.iOrder != 0) 
		  || (SD.pNode->SetDof(0) != DofOrder::DIFFERENTIAL)) {
	  std::cerr << "Sorry, a spring support must be linked to the algebraic value of a differential node" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	         
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

       switch (ConstLawType) {
	case DefHingeType::ELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelSpringSupport,	
				  GenelSpringSupport(uLabel, pDO, pCL, 
						     SD, fOut));
	   break;
	}
	case DefHingeType::VISCOUS:
	case DefHingeType::VISCOELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelSpringDamperSupport,
				  GenelSpringDamperSupport(uLabel, pDO, pCL, 
							   SD, fOut));
	   break;
	}
	default: {
	   std::cerr << "You shouldn't be here!" << std::endl;
	   THROW(DataManager::ErrGeneric());
	}	  
       }
       
       break;
    }      
      
    case CROSSSPRINGSUPPORT: {
       ScalarDof SDRow = ReadScalarDof(pDM, HP, 1);
       if (SDRow.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel spring supports" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       } 
       
       ScalarDof SDCol = ReadScalarDof(pDM, HP, 1);
       if ((SDCol.iOrder != 0) 
	   || (SDCol.pNode->SetDof(0) != DofOrder::DIFFERENTIAL)) {
	  std::cerr << "Sorry, a spring support must be linked to the algebraic value of a differential node" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       }	         
       
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);

       switch (ConstLawType) {
	case DefHingeType::ELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelCrossSpringSupport,	
				  GenelCrossSpringSupport(uLabel, pDO, pCL, 
							  SDRow, SDCol, fOut));
	   break;
	}
	case DefHingeType::VISCOUS:
	case DefHingeType::VISCOELASTIC: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  GenelCrossSpringDamperSupport,
				  GenelCrossSpringDamperSupport(uLabel, pDO, pCL, 
								SDRow, SDCol, fOut));
	   break;
	}
	default: {
	   std::cerr << "You shouldn't be here!" << std::endl;
	   THROW(DataManager::ErrGeneric());
	}	  
       }
       
       break;
    }      
      
    case MASS: {
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel mass" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  std::cerr << "Sorry, only differential dofs are allowed for genel mass" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }

       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      GenelMass,
			      GenelMass(uLabel, pDO, pDC, SD, fOut));
       
       
       break;
    }

#if 0
    case SCALARFILTER: {
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  Node::PARAMETER
	   || SD_u.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel scalar filter" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD_y.pNode->SetDof(0) != DofOrder::DIFFERENTIAL
		  || SD_u.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  std::cerr << "Sorry, only differential dofs are allowed for genel scalar filter" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       unsigned int na = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << na << " order denominator" << std::endl);
       doublereal* pdP = NULL;
       SAFENEWARR(pdP, doublereal, na+1);
       
       pdP[0] = 1.;
       
       if (na > 0) {	
	  for (unsigned int iCnt = 1; iCnt <= na; iCnt++) {
	     pdP[iCnt] = HP.GetReal();
	     DEBUGCOUT("a(" << iCnt << ") = " << pdP[iCnt] << std::endl);
	  }
       }
	     
       unsigned int nb = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << nb << " order numerator" << std::endl);
       doublereal* pdTau = NULL;
       SAFENEWARR(pdTau, doublereal, nb+1);
             
       for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	  pdTau[iCnt] = HP.GetReal();
	  DEBUGCOUT("b(" << iCnt << ") = " << pdTau[iCnt] << std::endl);
       }
            
       if (HP.IsKeyWord("gain")) {
	  doublereal gain = HP.GetReal();
	  DEBUGCOUT("Gain is: " << gain << std::endl);
	  for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	     pdTau[iCnt] *= gain;
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelFilter,
			      GenelFilter(uLabel, pDO, SD_y, SD_u,
					  na, nb, pdP, pdTau, fOut));
              
       break;
    }
#endif /* 0 */
      
    case SCALARFILTER: {
       /* output */
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel scalar filter output" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       /* input */
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);

       /* ordine del denominatore */
       unsigned int na = HP.GetInt();
       DEBUGCOUT("GenelFilter " << uLabel << " has a " << na 
		       << " order denominator" << std::endl);
       doublereal* pdP = NULL;
       SAFENEWARR(pdP, doublereal, na+1);
       pdP[0] = 1.;
       
       /* i coefficienti vengono letti a partire da a1, perche' il polinomio 
	* e' assunto monico (a0 = 1) */
       if (na > 0) {	
	  for (unsigned int iCnt = 1; iCnt <= na; iCnt++) {
	     pdP[iCnt] = HP.GetReal();
	     DEBUGCOUT("a(" << iCnt << ") = " << pdP[iCnt] << std::endl);
	  }
       }
	     
       /* ordine del numeratore */
       unsigned int nb = HP.GetInt();
       if (nb > na) {
	  std::cerr << "illegal (non proper) transfer function for scalar filter " 
	    << uLabel << " at line " << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }

       DEBUGCOUT("GenelFilter " << uLabel << " has a " << nb 
		       << " order numerator" << std::endl);
       doublereal* pdTau = NULL;
       SAFENEWARR(pdTau, doublereal, nb+1);
             
       /* i coefficienti vengono letti a partire da b0; e' possibile inserire 
	* il polinomio in forma monica e modificare direttamente il guadagno */
       for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	  pdTau[iCnt] = HP.GetReal();
	  DEBUGCOUT("b(" << iCnt << ") = " << pdTau[iCnt] << std::endl);
       }
            
       /* guadagno (opzionale) */
       if (HP.IsKeyWord("gain")) {
	  doublereal gain = HP.GetReal();
	  DEBUGCOUT("Gain is: " << gain << std::endl);
	  for (unsigned int iCnt = 0; iCnt <= nb; iCnt++) {
	     pdTau[iCnt] *= gain;
	  }
       }
       
       flag fState(0);
       if (HP.IsKeyWord("state")) {
	  if (HP.IsKeyWord("steady")) {
	     fState = 1;
	  } else {
	     std::cerr << "unknown state option at line " << HP.GetLineData() << std::endl;
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelFilterEq,
			      GenelFilterEq(uLabel, pDO, SD_y, SD_u,
					    na, nb, pdP, pdTau, 
					    fState, fOut));
              
       break;
    }
      
    case STATESPACESISO: {
       ScalarDof SD_y = ReadScalarDof(pDM, HP, 1);
       ScalarDof SD_u = ReadScalarDof(pDM, HP, 1);
       if (SD_y.pNode->GetNodeType() ==  Node::PARAMETER) {
	  std::cerr << "Sorry, parameters are not allowed for genel state space SISO output" << std::endl;
	  THROW(DataManager::ErrGeneric());	      
       } else if (SD_y.pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	  std::cerr << "Sorry, only differential dofs are allowed for genel state space SISO output" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }
       
       unsigned int Order = HP.GetInt();       
       DEBUGCOUT("State Space SISO " << uLabel << " is of order " << Order << std::endl);
       
       if (!HP.IsKeyWord("matrixA")) {
	  std::cerr << "matrix A expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdA = NULL;
       SAFENEWARR(pdA, doublereal, Order*Order);
       doublereal* pd = pdA;
       for (unsigned int i = 0; i < Order*Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixB")) {
	  std::cerr << "matrix B expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdB = NULL;
       SAFENEWARR(pdB, doublereal, Order);
       pd = pdB;
       for (unsigned int i = 0; i < Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixC")) {
	  std::cerr << "matrix C expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdC = NULL;
       SAFENEWARR(pdC, doublereal, Order);
       pd = pdC;
       for (unsigned int i = 0; i < Order; i++) {
	  *pd++ = HP.GetReal();
       }

       doublereal dD = 0.;
       if (HP.IsKeyWord("matrixD")) {
	  dD = HP.GetReal();
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelStateSpaceSISO,
			      GenelStateSpaceSISO(uLabel, pDO, SD_y, SD_u,
						  Order, 
						  pdA, pdB, pdC, dD, fOut));
              
       break;
    }
      
      
    case STATESPACEMIMO: {
       int iNumOutputs = HP.GetInt();
       if (iNumOutputs <= 0) {
	  std::cerr << "illegal number of outputs for state space MIMO " << uLabel
	    << " at line " << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDof* pvSD_y = NULL;
       SAFENEWARR(pvSD_y, ScalarDof, iNumOutputs);
       for (int i = 0; i < iNumOutputs; i++) {	  
	  pvSD_y[i] = ReadScalarDof(pDM, HP, 1);
	  if (pvSD_y[i].pNode->GetNodeType() ==  Node::PARAMETER) {
	     std::cerr << "line " << HP.GetLineData()
	       << ": sorry, parameters are not allowed for genel state space MIMO output" << std::endl;
	     THROW(DataManager::ErrGeneric());	      
	  } else if (pvSD_y[i].pNode->SetDof(0) != DofOrder::DIFFERENTIAL) {
	     std::cerr << "line " << HP.GetLineData() 
	       << ": sorry, only differential dofs are allowed for genel state space MIMO output" << std::endl;
	     THROW(DataManager::ErrGeneric());
	  }
       }
       
       int iNumInputs = HP.GetInt();
       if (iNumInputs <= 0) {
	  std::cerr << "illegal number of inputs for state space MIMO " << uLabel
	    << " at line " << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }       
       
       ScalarDof* pvSD_u = NULL;
       SAFENEWARR(pvSD_u, ScalarDof, iNumInputs);
       for (int i = 0; i < iNumInputs; i++) {	  
	  pvSD_u[i] = ReadScalarDof(pDM, HP, 1);
       }
       
       unsigned int Order = HP.GetInt();       
       DEBUGCOUT("State Space MIMO " << uLabel << " is of order " << Order << std::endl);
       
       if (!HP.IsKeyWord("matrixA")) {
	  std::cerr << "matrix A expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdA = NULL;
       SAFENEWARR(pdA, doublereal, Order*Order);
       doublereal* pd = pdA;
       for (unsigned int i = 0; i < Order*Order; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixB")) {
	  std::cerr << "matrix B expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdB = NULL;
       SAFENEWARR(pdB, doublereal, Order*iNumInputs);
       pd = pdB;
       for (unsigned int i = 0; i < Order*iNumInputs; i++) {
	  *pd++ = HP.GetReal();
       }
              
       if (!HP.IsKeyWord("matrixC")) {
	  std::cerr << "matrix C expected at line " << HP.GetLineNumber() << std::endl;
	  THROW(ErrGeneric());
       }
       doublereal* pdC = NULL;
       SAFENEWARR(pdC, doublereal, iNumOutputs*Order);
       pd = pdC;
       for (unsigned int i = 0; i < iNumOutputs*Order; i++) {
	  *pd++ = HP.GetReal();
       }

       doublereal* pdD = NULL;
       if (HP.IsKeyWord("matrixD")) {
	  SAFENEWARR(pdD, doublereal, iNumOutputs*iNumInputs);
	  pd = pdD;
	  for (int i = 0; i < iNumOutputs*iNumInputs; i++) {
	     *pd++ = HP.GetReal();
	  }
       }       
       
       flag fOut = pDM->fReadOutput(HP, Elem::GENEL);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      GenelStateSpaceMIMO,
			      GenelStateSpaceMIMO(uLabel, pDO, 
						  iNumOutputs, pvSD_y, 
						  iNumInputs, pvSD_u,
						  Order, 
						  pdA, pdB, pdC, pdD, fOut));

       break;
    }
      
      
      
      
      
      
      /* Aggiungere altri genel */
      
    default: {
       std::cerr << "unknown genel type in genel " << uLabel
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
} /* ReadGenel() */

#endif /* USE_ELECTRIC_NODES */

