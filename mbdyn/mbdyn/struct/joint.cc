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

/* joint */

#include <mbconfig.h>

extern "C" {
#include <strings.h>
}

#include <joint.h>
#include <dataman.h>

#include <drvhinge.h>
#include <drvj.h>      /* Vincoli di velocita' imposta */
#include <accj.h>      /* Vincoli di accelerazione imposta */
#include <genj.h>
#include <inplanej.h>  /* Vincoli di giacitura nel piano */
#include <inline.h>
#include <planej.h>
#include <planedj.h>
#include <prismj.h>    /* Vincolo prismatico */
#include <rodj.h>      /* Aste elastiche */
#include <spherj.h>
#include <univj.h>
#include <vehj.h>      /* Giunti deformabili */
#include <vehj2.h>     /* "" */
#include <kinj.h>

/* Provvisorio ?!? */
#include <modal.h>

/* Joint - begin */

Joint::Joint(unsigned int uL, JointType::Type T, const DofOwner* pDO, 
	     flag fOut)
: Elem(uL, ElemType::JOINT, fOut), 
ElemWithDofs(uL, ElemType::JOINT, pDO, fOut), 
InitialAssemblyElem(uL, ElemType::JOINT, fOut), JointT(T)
{ 
   NO_OP;
}
   
Joint::~Joint(void) 
{ 
   NO_OP; 
}
   
void Joint::Output(OutputHandler& OH) const
{
   if(fToBeOutput()) {
#ifdef DEBUG   
      OH.Output() << "Joint " << uLabel 
	<< ", type \"" << psJointNames[GetJointType()] 
	<< "\": sorry, not implemented yet" << endl;
#endif
      
      Output(OH.Joints(), "UnknownJoint", GetLabel(), 
	     Zero3, Zero3, Zero3, Zero3);
   }   
}


/* Output specifico dei vincoli */
ostream& Joint::Output(ostream& out, const char* sJointName, 
		       unsigned int uLabel,
		       const Vec3& FLocal, const Vec3& MLocal,
		       const Vec3& FGlobal, const Vec3& MGlobal) const
{
   /* Modificare le dimensioni del campo per il nome in base 
    * ai nomi dei vincoli futuri */
   ASSERT(strlen(sJointName) <= 16);
   
   /* Nota: non c'e' *endl* perche' i vincoli possono aggiungere outut 
    * ulteriore a quello comune a tutti
   return out << sJointName << setw(16+8-strlen(sJointName)) << uLabel << " "
     << FLocal << " " << MLocal << " " << FGlobal << " " << MGlobal;
    */
   return out << setw(8) << uLabel << " "
     << FLocal << " " << MLocal << " " << FGlobal << " " << MGlobal;
}

/* Joint - end */


/* Legge un vincolo */

Elem* ReadJoint(DataManager* pDM,
		MBDynParser& HP,	       
		const DofOwner* pDO,
		unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadJoint");
   
   const char* sKeyWords[] = {
      "distance",
      "distance" "with" "offset",
      "clamp",
      "coincidence",
      "spherical" "hinge",
      "pin",
      "universal" "hinge",
      "universal" "pin",
      "plane" "hinge",
      "plane" "pin",
      "axial" "rotation",
      "plane" "displacement",
      "plane" "displacement" "pin",
      "in" "plane",
      "in" "line",
      "rod",
      "deformable" "hinge",
      "deformable" "displacement" "hinge",
      "linear" "velocity",
      "angular" "velocity",
      "linear" "acceleration",
      "angular" "acceleration",
      "prismatic",
      "drive" "hinge",
      "kinematic",
      
      "modal",
      
      NULL
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
      
      DISTANCE = 0,
      DISTANCEWITHOFFSET,
      CLAMP,
      COINCIDENCE,
      SPHERICALHINGE,
      PIN,
      UNIVERSALHINGE,
      UNIVERSALPIN,
      PLANEHINGE,
      PLANEPIN,
      AXIALROTATION,
      PLANEDISPLACEMENT,
      PLANEDISPLACEMENTPIN,
      INPLANE,
      INLINE,
      ROD,
      DEFORMABLEHINGE,
      DEFORMABLEDISPHINGE,
      LINEARVELOCITY,
      ANGULARVELOCITY,
      LINEARACCELERATION,
      ANGULARACCELERATION,
      PRISMATIC,
      DRIVEHINGE,
      KINEMATIC,
      
      MODAL,
      
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
      cout << "joint type: " << sKeyWords[CurrKeyWord] << endl;
   }
#endif   

   Joint* pEl = NULL;
   
   switch (CurrKeyWord) {
      
      /* vincolo di distanza */
    case DISTANCE: {
       /* lettura dei dati specifici */
       /* due nodi e tipo di drive, con dati specifici */
       
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData()
	    << ": structural node " << uNode2
	    << " not defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       DriveCaller* pDC = NULL;
       if (HP.IsKeyWord("fromnodes")) {
	  doublereal l = (pNode2->GetXCurr()-pNode1->GetXCurr()).Norm();
	  SAFENEWWITHCONSTRUCTOR(pDC, 
				 ConstDriveCaller,
				 ConstDriveCaller(pDM->pGetDrvHdl(), l),
				 DMmm);
       } else {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       }
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       /* allocazione e costruzione */
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      DistanceJoint,
			      DistanceJoint(uLabel, pDO, pNode1, pNode2, pDC, fOut), 
			      DMmm);
       
       /* scrittura dei dati specifici */	     	     
       break;
    }
      
      /* vincolo di distanza con offset */
    case DISTANCEWITHOFFSET: {
       /* lettura dei dati specifici */
       /* due nodi e tipo di drive, con dati specifici */
       
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));	    				     
       
       DEBUGCOUT("Offset 1: " << f1 << endl);
       
              
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
       
       DEBUGCOUT("Offset 2: " << f2 << endl);
       
              
       /* Legge e costruisce il drivecaller */
       if (!HP.fIsArg()) {
	  cerr << "line " << HP.GetLineData()
	    << ": driver data expected" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }	     
       
       DriveCaller* pDC = NULL;
       if (HP.IsKeyWord("fromnodes")) {
	  doublereal l = (pNode2->GetXCurr()+f2-pNode1->GetXCurr()-f1).Norm();
	  SAFENEWWITHCONSTRUCTOR(pDC, 
				 ConstDriveCaller, 
				 ConstDriveCaller(pDM->pGetDrvHdl(), l),
				 DMmm);
       } else {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       }
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       /* allocazione e costruzione */
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      DistanceJointWithOffset,
			      DistanceJointWithOffset(uLabel, pDO, pNode1, pNode2,
						      f1, f2, pDC, fOut), DMmm);
       
       /* scrittura dei dati specifici */	     	     
       break;
    }
      
      /* vincolo di incastro */
    case CLAMP: {	     
       /* lettura dei dati specifici */
       
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* posizione (vettore di 3 elementi) */
       ReferenceFrame RF(pNode);
       Vec3 X0;
       if (HP.IsKeyWord("node")) { /* stessa posizione del nodo */
	  X0 = pNode->GetXCurr();
       } else {                     /* posizione arbitraria */
	  X0 = HP.GetPosAbs(RF);
       }	   	       	     	       
       
       DEBUGCOUT("X0 =" << endl << X0 << endl);
       
       /* sistema di riferimento (trucco dei due vettori) */
       Mat3x3 R0;
       if (HP.IsKeyWord("node")) { /* stessa giacitura del nodo */
	  R0 = pNode->GetRCurr();
       } else {                     /* giacitura arbitraria */
	  R0 = HP.GetRotAbs(RF);
       }	   
       
       DEBUGCOUT("R0 =" << endl << R0 << endl);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       /* allocazione e costruzione */
       SAFENEWWITHCONSTRUCTOR(pEl,
			      ClampJoint,
			      ClampJoint(uLabel, pDO, pNode, X0, R0, fOut), DMmm);
       break;
    }      	   	            
      
    case PIN: {
       /* lettura dei dati specifici */
       
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 d(HP.GetPosRel(ReferenceFrame(pNode)));
       
       DEBUGCOUT("Node reference frame d:" << endl << d << endl);
       
       
       /* posizione (vettore di 3 elementi) */	 
       Vec3 X0(HP.GetPosAbs(AbsRefFrame));
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pEl,
			      PinJoint,
			      PinJoint(uLabel, pDO, pNode, X0, d, fOut), DMmm);
       
       break;
    }	
      
      /* vincolo di cerniera piana (PLANEHINGE)
	* eventualmente con velocita' di rotazione imposta (AXIALROTATION) */
    case SPHERICALHINGE:
    case PLANEHINGE:
    case UNIVERSALHINGE:
    case AXIALROTATION:
    case PLANEDISPLACEMENT: {
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       ReferenceFrame RF(pNode1);
       Vec3 d1(HP.GetPosRel(RF));	   
       
       DEBUGCOUT("Node 1 reference frame d1:" << endl << d1 << endl);
       
       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);
	  R1h = HP.GetRotRel(RF);
       }
       
       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* Stessa cosa per il nodo 2 */
       
       RF = ReferenceFrame(pNode2);
       Vec3 d2(HP.GetPosRel(RF));	   
       
       DEBUGCOUT("Node 2 reference frame d2:" << endl << d2 << endl);
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);
	  R2h = HP.GetRotRel(RF);
       }
       
       
       DriveCaller* pDC = NULL;
       if (CurrKeyWord == AXIALROTATION) {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	  HP.PutKeyTable(K);
       }	   
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       switch (CurrKeyWord) {
	  
	  
	  /* allocazione e creazione cerniera sferica */
	case SPHERICALHINGE: {		   
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  SphericalHingeJoint,
				  SphericalHingeJoint(uLabel, pDO, 
						      pNode1, pNode2, 
						      d1, R1h,
						      d2, R2h, fOut),
				  DMmm);
	   break;
	}
	  
	  /* allocazione e creazione cerniera piana */
	case PLANEHINGE: {	
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlaneHingeJoint,
				  PlaneHingeJoint(uLabel, pDO, pNode1, pNode2, 
						  d1, d2, R1h, R2h, fOut), 
				  DMmm);
	   break;
	}
	  
	  /* allocazione e creazione cerniera universale */
	case UNIVERSALHINGE: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  UniversalHingeJoint,
				  UniversalHingeJoint(uLabel, pDO, 
						      pNode1, pNode2,
						      d1, d2, R1h, R2h, fOut), 
				  DMmm);
	   break;
	}
	  
	  /* allocazione e creazione cerniera piana 
	   * con velocita' di rotazione imposta */
	case AXIALROTATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  AxialRotationJoint,
				  AxialRotationJoint(uLabel, pDO, 
						     pNode1, pNode2,
						     d1, d2, R1h, R2h, pDC, 
						     fOut), 
				  DMmm);
	   break;
	}
	   
	  /* allocazione e creazione pattino */
	case PLANEDISPLACEMENT: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlaneDispJoint,
				  PlaneDispJoint(uLabel, pDO, pNode1, pNode2, 
						 d1, d2, R1h, R2h, fOut),
				  DMmm);
	   
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}	      
       }	     
       
       break;
    }
      
    case UNIVERSALPIN:
    case PLANEPIN:
    case PLANEDISPLACEMENTPIN: {
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;
	  THROW(DataManager::ErrGeneric());
       }		  
       
       ReferenceFrame RF(pNode);
       Vec3 d(HP.GetPosRel(RF));
       
       DEBUGCOUT("Node reference frame d:" << endl << d << endl);
       
       Mat3x3 Rh(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);
	  Rh = HP.GetRotRel(RF);
	  DEBUGCOUT("Hinge Rotation matrix Rh:" << endl << Rh << endl);
       }

       Vec3 X0(HP.GetPosAbs(AbsRefFrame));
       DEBUGCOUT("Absolute X:" << endl << X0 << endl);

       Mat3x3 R0(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  R0 = HP.GetRotAbs(AbsRefFrame);
       }
       DEBUGCOUT("Absolute R:" << endl << R0 << endl);
       
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       switch (CurrKeyWord) {
	  
	  /* allocazione e creazione cerniera piana */	     
	case PLANEPIN: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlanePinJoint,
				  PlanePinJoint(uLabel, pDO, pNode,
						X0, R0, d, Rh, fOut), 
				  DMmm);
	   break;
	}
	  
	case UNIVERSALPIN: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  UniversalPinJoint,
				  UniversalPinJoint(uLabel, pDO, pNode,
						    X0, R0, d, Rh, fOut), 
				  DMmm);
	   break;
	}
	  
	  /* allocazione e creazione cerniera piana */	     
	case PLANEDISPLACEMENTPIN: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlaneDispPinJoint,
				  PlaneDispPinJoint(uLabel, pDO, pNode,
						    X0, R0, d, Rh, fOut), 
				  DMmm);
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}
       }	     	 
       
       break;
    }
      
    case INPLANE: {
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       ReferenceFrame RF(pNode1);
       Vec3 p(HP.GetPosRel(RF));
       
       DEBUGCOUT("Node 1 reference frame p:" << endl << p << endl);
       
       Vec3 v(HP.GetVecRel(RF));
       doublereal d = v.Dot();
       if (d <= DBL_EPSILON) {
	  cerr << "null direction at line " << HP.GetLineData() << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       v /= sqrt(d);
       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
              
       Vec3 q;
       flag fOffset(0);
       flag fOut(0);
       
       if (HP.fIsArg()) {
	  if (HP.IsKeyWord("offset")) {
	     fOffset = 1;	     
	     q = HP.GetPosRel(ReferenceFrame(pNode2));	     
	     DEBUGCOUT("Node 2 reference frame q:" << endl << p << endl);
	  }
       }
       
       fOut = pDM->fReadOutput(HP, ElemType::JOINT);      	   
       
       if (fOffset) {	      
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 InPlaneWithOffsetJoint,
				 InPlaneWithOffsetJoint(uLabel, pDO,
							pNode1, pNode2, 
							v, p, q, fOut), DMmm);
       } else {	      	   
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 InPlaneJoint,
				 InPlaneJoint(uLabel, pDO, pNode1, pNode2, 
					      v, p, fOut), DMmm);
       }
       
       break;
    }	
      
    case INLINE: {
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       ReferenceFrame RF(pNode1);
       Vec3 p(HP.GetPosRel(RF));
       
       DEBUGCOUT("Node 1 reference frame p:" << endl << p << endl);
       
       Mat3x3 R(HP.GetRotRel(RF));       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
              
       Vec3 q;
       flag fOffset(0);
       flag fOut(0);
       
       if (HP.fIsArg()) {
	  if (HP.IsKeyWord("offset")) {
	     fOffset = 1;	     
	     q = HP.GetPosRel(ReferenceFrame(pNode2));	     
	     DEBUGCOUT("Node 2 reference frame q:" << endl << p << endl);
	  }
       }
       
       fOut = pDM->fReadOutput(HP, ElemType::JOINT);      	   
       
       if (fOffset) {
	  SAFENEWWITHCONSTRUCTOR(pEl,
				 InLineWithOffsetJoint,
				 InLineWithOffsetJoint(uLabel, pDO, 
						       pNode1, pNode2,
						       R, p, q, fOut), 
				 DMmm);
       } else {	      	   
	  SAFENEWWITHCONSTRUCTOR(pEl,
				 InLineJoint,
				 InLineJoint(uLabel, pDO, 
					     pNode1, pNode2,
					     R, p, fOut), 
				 DMmm);
       }
       
       break;
    }	
      
      /* vincolo di distanza */
    case ROD: {
       
       /* lettura dei dati specifici */
       /* due nodi, lunghezza iniziale e tipo di legame elastico */
       
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* Lunghezza iniziale */
       doublereal dL = 0.;
       flag fFromNodes(0);
       if (HP.IsKeyWord("fromnodes")) {
	  fFromNodes = 1;
	  DEBUGCOUT("Initial length will be computed from nodes position" << endl);	  
       } else {
	  dL = HP.GetReal();
	  DEBUGCOUT("Initial length = " << dL << endl);
       }
       
       flag fOffset(0);
       Vec3 f1;
       Vec3 f2;
       
       /* Se si tratta di Rod con Offset, legge gli offset e poi passa 
	* al tipo di legame costitutivo */
       if (HP.IsKeyWord("offset")) {
	  fOffset = 1;	  
	  f1 = HP.GetPosRel(ReferenceFrame(pNode1));	  
	  DEBUGCOUT("Offset 1: " << f1 << endl);

	  f2 = HP.GetPosRel(ReferenceFrame(pNode2));	  
	  DEBUGCOUT("Offset 2: " << f2 << endl);
       }
       
       if (fFromNodes) {
	  Vec3 v = pNode2->GetXCurr()-pNode1->GetXCurr();
	  if (fOffset) {
	     v += pNode2->GetRCurr()*f2-pNode1->GetRCurr()*f1;
	  }
	  dL = v.Norm();
	  DEBUGCOUT("Initial length = " << dL << endl);	  
       }
       
       /* Legame costitutivo */
       DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, ConstLawType);
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       if (fOffset == 1) {	      
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 RodWithOffsetJoint,
				 RodWithOffsetJoint(uLabel, pDO, pCL,
						    pNode1, pNode2, 
						    f1, f2, dL, fOut),
				 DMmm);
       } else {	  
	  if (ConstLawType == DefHingeType::VISCOUS || ConstLawType == DefHingeType::VISCOELASTIC) {
	     SAFENEWWITHCONSTRUCTOR(pEl, 
				    ViscoElasticRodJoint,
				    ViscoElasticRodJoint(uLabel, pDO, pCL,
							 pNode1, pNode2, 
							 dL, fOut), 
				    DMmm);		 
	  } else {		 
	     SAFENEWWITHCONSTRUCTOR(pEl, 
				    RodJoint,
				    RodJoint(uLabel, pDO, pCL, pNode1, pNode2, dL, fOut), 
				    DMmm);
	  }
       }	   
       
       break;
    }
      
    case DEFORMABLEHINGE:
    case DEFORMABLEDISPHINGE: {
       /* lettura dei dati specifici */	  
       
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Offset if displacement hinge */
       Vec3 f1;
       if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	  f1 = HP.GetPosRel(ReferenceFrame(pNode1));
       }	   
       
       Mat3x3 R1(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  R1 = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  

	   
       /* Offset if displacement hinge */
       Vec3 f2;
       if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	  f2 = HP.GetPosRel(ReferenceFrame(pNode2));
       }	   
       
       Mat3x3 R2(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  R2 = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       
       
       /* Legame costitutivo */
       DefHingeType::Type ConstLawType;
       ConstitutiveLaw3D* pCL = pDM->ReadConstLaw3D(HP, ConstLawType);
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       switch (ConstLawType) {
	case DefHingeType::ELASTIC: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ElasticHingeJoint,
				     ElasticHingeJoint(uLabel, pDO, pCL,
						       pNode1, pNode2, 
						       R1, R2, fOut),
				     DMmm);
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ElasticDispHingeJoint,
				     ElasticDispHingeJoint(uLabel, pDO, pCL,
							   pNode1, pNode2, 
							   f1, f2, R1, R2, 
							   fOut),
				     DMmm);
	   }
	   
	   break;
	}
	  
	case DefHingeType::VISCOUS: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscousHingeJoint,
				     ViscousHingeJoint(uLabel, pDO, pCL,
						       pNode1, pNode2, 
						       R1, R2, fOut),
				     DMmm);
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscousDispHingeJoint,
				     ViscousDispHingeJoint(uLabel, pDO, pCL,
							   pNode1, pNode2,
							   f1, f2, R1, R2, 
							   fOut),
				     DMmm);
	   }
	   
	   break;
	}
	  
	case DefHingeType::VISCOELASTIC: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscoElasticHingeJoint,
				     ViscoElasticHingeJoint(uLabel, pDO, pCL,
							    pNode1, pNode2, 
							    R1, R2, fOut),
				     DMmm);
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ViscoElasticDispHingeJoint,
				     ViscoElasticDispHingeJoint(uLabel, 
								pDO, pCL,
								pNode1, pNode2,
								f1, f2, R1, R2,
								fOut),
				     DMmm);
	   }
	   
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}
       }	   	   
       
       break;
    }            
      
    case LINEARVELOCITY:
    case ANGULARVELOCITY: {
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 Dir(HP.GetVecRel(ReferenceFrame(pNode)));
       doublereal d = Dir.Dot();
       ASSERT(d > 0.);	 
       if (d > 0.) {	      
	  Dir /= sqrt(d);
       } else {
	  cerr << "warning, direction vector is null" << endl;
       }	     
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       switch (CurrKeyWord) {
	case LINEARVELOCITY: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  LinearVelocityJoint,
				  LinearVelocityJoint(uLabel, pDO,
						      pNode, Dir, pDC, fOut), 
				  DMmm);
	   break;
	}
	  
	case ANGULARVELOCITY: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  AngularVelocityJoint,
				  AngularVelocityJoint(uLabel, pDO,
						       pNode, Dir, pDC, fOut), 
				  DMmm);
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}
       }
       
       break;
    }
      
    case LINEARACCELERATION:
    case ANGULARACCELERATION: {
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 Dir(HP.GetVecRel(ReferenceFrame(pNode)));
       doublereal d = Dir.Dot();
       ASSERT(d > 0.);	 
       if (d > 0.) {	      
	  Dir /= sqrt(d);
       } else {
	  cerr << "warning, direction vector is null" << endl;
       }	     
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       HP.PutKeyTable(K);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       switch (CurrKeyWord) {
	case LINEARACCELERATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  LinearAccelerationJoint,
				  LinearAccelerationJoint(uLabel, pDO,
							  pNode, Dir, pDC, 
							  fOut), 
				  DMmm);
	   break;
	}
	  
	case ANGULARACCELERATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  AngularAccelerationJoint,
				  AngularAccelerationJoint(uLabel, pDO,
							   pNode, Dir, pDC, 
							   fOut), 
				  DMmm);
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}
       }
       
       break;
    }
      
    case PRISMATIC: {	     
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);	   
	  R1h = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* Stessa cosa per il nodo 2 */
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);
	  R2h = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      PrismaticJoint,
			      PrismaticJoint(uLabel, pDO, pNode1, pNode2, 
					     R1h, R2h, fOut), 
			      DMmm);
       
       break;
    }
      

    case DRIVEHINGE: {
       /* nodo collegato 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode1;
       if ((pNode1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);	   
	  R1h = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode2;
       if ((pNode2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* Stessa cosa per il nodo 2 */
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << endl);
	  R2h = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       

       TplDriveCaller<Vec3>* pDC 
	 = ReadTplDrive(pDM, HP, pDM->pGetDrvHdl(), Vec3(0.));
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      DriveHingeJoint,
			      DriveHingeJoint(uLabel, pDO, pDC, 
					      pNode1, pNode2, 
					      R1h, R2h, fOut), 
			      DMmm);
       
       
       break;
    }
      
      
    case KINEMATIC: {
       /* nodo collegato */
       unsigned int uNode = (unsigned int)HP.GetInt();	    
       DEBUGCOUT("Linked to Node " << uNode << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNode;
       if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       Kinematics* pK = NULL;
       SAFENEWWITHCONSTRUCTOR(pK, 
			      KinematicsTest,
			      KinematicsTest(pDC),
			      DMmm);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      KinJoint,
			      KinJoint(uLabel, pDO, pNode, pK, fOut),
			      DMmm);
       
       break;
    }

      
    case MODAL: {
       pEl = ReadModal(pDM, HP, pDO, uLabel);
       break;
    }
      
      
      /* Aggiungere altri vincoli */
      
    default: {
       cerr << "unknown joint type in joint " << uLabel
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
} /* ReadJoint() */
