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

/* joint */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <strings.h>

#include <joint.h>
#include <dataman.h>

#include <tpldrive.h>
#include <tpldrive_.h>

#include <distance.h>

#include <drvhinge.h>
#include <drvj.h>      /* Vincoli di velocita' imposta */
#include <accj.h>      /* Vincoli di accelerazione imposta */
#include <genj.h>
#include <inplanej.h>  /* Vincoli di giacitura nel piano */
#include <inline.h>
#include <planej.h>
/* No longer supported 
#include <planedj.h>
 */
#include <prismj.h>    /* Vincolo prismatico */
#include <rodj.h>      /* Aste elastiche */
#include <spherj.h>
#include <univj.h>
#include <vehj.h>      /* Giunti deformabili */
#include <vehj2.h>     /* "" */
#include <kinj.h>
#include <beamslider.h>

/* Provvisorio ?!? */
#include <modal.h>

#ifdef __HACK_MODULES_AS_JOINTS__
/* To allow direct loading of wheel2 (and other modules, in the future)
 * as joints ... */
#include <modules.h>
#endif /* __HACK_MODULES_AS_JOINTS__ */

#define MBDYN_X_COMPATIBLE_INPUT

/* Joint - begin */

Joint::Joint(unsigned int uL, Joint::Type T, const DofOwner* pDO, 
	     flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
ElemWithDofs(uL, Elem::JOINT, pDO, fOut), 
InitialAssemblyElem(uL, Elem::JOINT, fOut), JointT(T)
{ 
   NO_OP;
}
   
Joint::~Joint(void) 
{ 
   NO_OP; 
}
   

/* Output specifico dei vincoli */
std::ostream& Joint::Output(std::ostream& out, const char* sJointName, 
		       unsigned int uLabel,
		       const Vec3& FLocal, const Vec3& MLocal,
		       const Vec3& FGlobal, const Vec3& MGlobal) const
{
   /* Modificare le dimensioni del campo per il nome in base 
    * ai nomi dei vincoli futuri */
   ASSERT(strlen(sJointName) <= 16);
   
   /* Nota: non c'e' *std::endl* perche' i vincoli possono aggiungere outut 
    * ulteriore a quello comune a tutti
   return out << sJointName << std::setw(16+8-strlen(sJointName)) << uLabel << " "
     << FLocal << " " << MLocal << " " << FGlobal << " " << MGlobal;
    */
   return out << std::setw(8) << uLabel << " "
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
      "universal" "rotation",
      "universal" "pin",
      "plane" "hinge",
      "revolute" "hinge",
      "revolute" "rotation",
      "plane" "pin",
      "revolute" "pin",
      "axial" "rotation",
      "plane" "displacement",
      "plane" "displacement" "pin",
      "in" "plane",
      "in" "line",
      "rod",
      "rod" "with" "offset",
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

      "beam" "slider",
      "wheel2",
      
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
      UNIVERSALROTATION,
      UNIVERSALPIN,
      PLANEHINGE,
      REVOLUTEHINGE,
      REVOLUTEROTATION,
      PLANEPIN,
      REVOLUTEPIN,
      AXIALROTATION,
      PLANEDISPLACEMENT,
      PLANEDISPLACEMENTPIN,
      INPLANE,
      INLINE,
      ROD,
      RODWITHOFFSET,
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

      BEAMSLIDER,

      WHEEL2,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);
   
   /* lettura del tipo di vincolo */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {
      std::cout << "joint type: " << sKeyWords[CurrKeyWord] << std::endl;
   }
#endif   

   Joint* pEl = NULL;
   
   switch (CurrKeyWord) {
      
      /* vincolo di distanza */
    case DISTANCE: {
       /* lettura dei dati specifici */
       /* due nodi e tipo di drive, con dati specifici */
       
       flag fOffset(0);

       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

       Vec3 f1(0.);
       if (HP.IsKeyWord("position")) {
	       f1 = HP.GetPosRel(ReferenceFrame(pNode1));
	       fOffset = 1;
       }
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 f2(0.);
       if (HP.IsKeyWord("position")) {
	       f2 = HP.GetPosRel(ReferenceFrame(pNode2));
	       fOffset = 1;
       }
       
       DriveCaller* pDC = NULL;
       if (HP.IsKeyWord("from" "nodes")) {
	  doublereal l = (pNode2->GetXCurr()
			  -pNode1->GetXCurr()
			  +pNode2->GetRCurr()*f2
			  -pNode1->GetRCurr()*f1).Norm();

	  pedantic_cout("Distance(" << uLabel << "): "
			  "length from nodes = " << l << std::endl);

	  SAFENEWWITHCONSTRUCTOR(pDC, 
				 ConstDriveCaller,
				 ConstDriveCaller(pDM->pGetDrvHdl(), l));
       } else {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       /* allocazione e costruzione */
       if (fOffset) {
	       SAFENEWWITHCONSTRUCTOR(pEl, 
			       DistanceJointWithOffset,
			       DistanceJointWithOffset(uLabel, pDO, 
				       pNode1, pNode2, f1, f2, pDC, fOut));
       } else {
     	       SAFENEWWITHCONSTRUCTOR(pEl, 
 			       DistanceJoint,
 			       DistanceJoint(uLabel, pDO,
				       pNode1, pNode2, pDC, fOut));
       }
       
       break;
    }
      
      /* vincolo di distanza con offset */
    case DISTANCEWITHOFFSET: {
       /* lettura dei dati specifici */
       /* due nodi e tipo di drive, con dati specifici */

#ifdef MBDYN_X_COMPATIBLE_INPUT
       pedantic_cerr("Joint(" << uLabel << "): \"distance with offset\" "
		       "is deprecated; use \"distance\" instead" << std::endl);
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 f1(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          f1 = HP.GetPosRel(ReferenceFrame(pNode1));
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* !MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("Offset 1: " << f1 << std::endl);
       
              
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 f2(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          f2 = HP.GetPosRel(ReferenceFrame(pNode2));
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* !MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("Offset 2: " << f2 << std::endl);
       
 
       /* Legge e costruisce il drivecaller */
       if (!HP.IsArg()) {
	  std::cerr << "line " << HP.GetLineData()
	    << ": driver data expected" << std::endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       
       DriveCaller* pDC = NULL;
       if (HP.IsKeyWord("from" "nodes")) {
	  doublereal l = (pNode2->GetXCurr()
			  +pNode2->GetRCurr()*f2
			  -pNode1->GetXCurr()
			  -pNode1->GetRCurr()*f1).Norm();

	  pedantic_cout("DistanceWithOffset(" << uLabel << "): "
			  "length from nodes = " << l << std::endl);

	  SAFENEWWITHCONSTRUCTOR(pDC, 
				 ConstDriveCaller, 
				 ConstDriveCaller(pDM->pGetDrvHdl(), l));
       } else {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       /* allocazione e costruzione */
       SAFENEWWITHCONSTRUCTOR(pEl, 
		       DistanceJointWithOffset,
		       DistanceJointWithOffset(uLabel, pDO, pNode1, pNode2,
			       f1, f2, pDC, fOut));
       
       /* scrittura dei dati specifici */
       break;
    }
      
      /* vincolo di incastro */
    case CLAMP: {	     
       /* lettura dei dati specifici */
       
       /* nodo collegato */
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       /* posizione (vettore di 3 elementi) */
       ReferenceFrame RF(pNode);
       /* stessa posizione del nodo */
       Vec3 X0(pNode->GetXCurr());
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          if (!HP.IsKeyWord("node")) {
	     /* posizione arbitraria */
	     X0 = HP.GetPosAbs(RF);
          }	   	       	     	       
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* !MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("X0 =" << std::endl << X0 << std::endl);
       
       /* sistema di riferimento (trucco dei due vettori) */
       /* stessa giacitura del nodo */
       Mat3x3 R0(pNode->GetRCurr());
       if (HP.IsKeyWord("orientation")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"orientation\" at line "
			       << HP.GetLineData());
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          if (!HP.IsKeyWord("node")) {
	     /* giacitura arbitraria */
	     R0 = HP.GetRotAbs(RF);
          }
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* !MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("R0 =" << std::endl << R0 << std::endl);
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       /* allocazione e costruzione */
       SAFENEWWITHCONSTRUCTOR(pEl,
			      ClampJoint,
			      ClampJoint(uLabel, pDO, pNode, X0, R0, fOut));
       break;
    }      	   	            
      
    case PIN: {
       /* lettura dei dati specifici */
       
       /* nodo collegato */
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode);

       Vec3 d(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	  pedantic_cerr("Joint(" << uLabel 
			  << "): missing keyword \"position\" at line "
			  << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          d = HP.GetPosRel(RF);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */

       /* currently unused */
       if (HP.IsKeyWord("orientation")) {
	  (void)HP.GetRotRel(RF);
       }
       
       DEBUGCOUT("Node reference frame d:" << std::endl << d << std::endl);
       
       /* posizione (vettore di 3 elementi) */	 
       Vec3 X0(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	  pedantic_cerr("Joint(" << uLabel 
			  << "): missing keyword \"position\" at line "
			  << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
	  X0 = HP.GetPosAbs(AbsRefFrame);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */

       /* currently unused */
       if (HP.IsKeyWord("orientation")) {
	  (void)HP.GetRotRel(ReferenceFrame(pNode));
       }
       
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pEl,
			      PinJoint,
			      PinJoint(uLabel, pDO, pNode, X0, d, fOut));
       
       break;
    }	
      
      /* vincolo di cerniera piana (PLANEHINGE)
	* eventualmente con velocita' di rotazione imposta (AXIALROTATION) */
    case SPHERICALHINGE:
    case PLANEHINGE:
    case REVOLUTEHINGE:
    case REVOLUTEROTATION:
    case UNIVERSALHINGE:
    case UNIVERSALROTATION:
    case AXIALROTATION:
    case PLANEDISPLACEMENT: {
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode1);
       Vec3 d1(0.);

#ifndef MBDYN_X_COMPATIBLE_INPUT
       if (HP.IsKeyWord("position")) {
	  d1 = HP.GetPosRel(RF);
       }
#else /* MBDYN_X_COMPATIBLE_INPUT */
       switch (CurrKeyWord) {
       case REVOLUTEROTATION:
       case UNIVERSALROTATION:
	  if (HP.IsKeyWord("position")) {
	     /* currently ignored */
	     (void)HP.GetPosRel(RF);
	  } else {
	     pedantic_cerr("Joint(" << uLabel 
	    	     << "): missing keyword \"position\" at line " 
	    	     << HP.GetLineData() << std::endl);
	  }
          break;

       default:
	  if (!HP.IsKeyWord("position")) {
	     pedantic_cerr("Joint(" << uLabel 
	    	     << "): missing keyword \"position\" at line " 
	    	     << HP.GetLineData() << std::endl);
	  }
          d1 = HP.GetPosRel(RF);
   	  DEBUGCOUT("Node 1 reference frame d1:" << std::endl 
			  << d1 << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */

       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("orientation")) {
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R1h = HP.GetRotRel(RF);
#ifdef MBDYN_X_COMPATIBLE_INPUT
       } else if (HP.IsKeyWord("hinge")) {
	  pedantic_cerr("Joint(" << uLabel 
		 << "): keyword \"hinge\" at line " << HP.GetLineData()
		 << " is deprecated; use \"orientation\" instead" << std::endl);
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R1h = HP.GetRotRel(RF);
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       }
       
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       /* Stessa cosa per il nodo 2 */
       
       RF = ReferenceFrame(pNode2);
       Vec3 d2(0.);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       if (HP.IsKeyWord("position")) {
	  d2 = HP.GetPosRel(RF);
       }
#else /* MBDYN_X_COMPATIBLE_INPUT */
       switch (CurrKeyWord) {
       case REVOLUTEROTATION:
       case UNIVERSALROTATION:
	  if (HP.IsKeyWord("position")) {
	     /* currently ignored */
	     (void)HP.GetPosRel(RF);
	  } else {
	     pedantic_cerr("Joint(" << uLabel 
	    	     << "): missing keyword \"position\" at line " 
	    	     << HP.GetLineData() << std::endl);
	  }
          break;

       default:
	  if (!HP.IsKeyWord("position")) {
	     pedantic_cerr("Joint(" << uLabel 
	    	     << "): missing keyword \"position\" at line " 
	    	     << HP.GetLineData() << std::endl);
	  }
          d2 = HP.GetPosRel(RF);
   	  DEBUGCOUT("Node 2 reference frame d2:" << std::endl 
			  << d2 << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("orientation")) {
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R2h = HP.GetRotRel(RF);
#ifdef MBDYN_X_COMPATIBLE_INPUT
       } else if (HP.IsKeyWord("hinge")) {
	  pedantic_cerr("Joint(" << uLabel 
		 << "): keyword \"hinge\" at line " << HP.GetLineData()
		 << " is deprecated; use \"orientation\" instead" << std::endl);
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R2h = HP.GetRotRel(RF);
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       }
       
       
       DriveCaller* pDC = NULL;
       if (CurrKeyWord == AXIALROTATION) {
	  pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       }	   
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       switch (CurrKeyWord) {
	  
	  /* allocazione e creazione cerniera sferica */
	case SPHERICALHINGE: {		   
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  SphericalHingeJoint,
				  SphericalHingeJoint(uLabel, pDO, 
						      pNode1, pNode2, 
						      d1, R1h,
						      d2, R2h, fOut));
	   break;
	}
	  
	  /* allocazione e creazione cerniera piana */
	case PLANEHINGE:
	   std::cerr << "line " << HP.GetLineData() 
		   << ": deprecated 'plane hinge' name;"
		   << " use 'revolute hinge' instead" << std::endl;
	case REVOLUTEHINGE: {
	   doublereal r = 0.;
	   BasicFriction * bf = 0;
	   BasicShapeCoefficient * bsh = 0;
	   if (HP.IsKeyWord("friction")) {
	       r = HP.GetReal();
	       bf = ParseFriction(HP,pDM);
	       bsh = ParseShapeCoefficient(HP);
	   }	
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlaneHingeJoint,
				  PlaneHingeJoint(uLabel, pDO, pNode1, pNode2, 
						  d1, d2, R1h, R2h, fOut,
						  r,bsh,bf));
	   break;
	}

	  /* allocazione e creazione cerniera piana senza vincolo in pos. */
	case REVOLUTEROTATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlaneRotationJoint,
				  PlaneRotationJoint(uLabel, pDO, 
					  pNode1, pNode2, R1h, R2h, fOut));
	   break;
	}
	  
	  /* allocazione e creazione cerniera universale */
	case UNIVERSALHINGE: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  UniversalHingeJoint,
				  UniversalHingeJoint(uLabel, pDO, 
						      pNode1, pNode2,
						      d1, d2, R1h, R2h, fOut));
	   break;
	}
	  
	  /* allocazione e creazione cerniera universale senza vincolo pos. */
	case UNIVERSALROTATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  UniversalRotationJoint,
				  UniversalRotationJoint(uLabel, pDO, 
						      pNode1, pNode2,
						      R1h, R2h, fOut));
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
						     fOut));
	   break;
	}
	   
	  /* allocazione e creazione pattino */
	case PLANEDISPLACEMENT:
	   std::cerr << "PlaneDispJoint(" << uLabel << "): "
		   "no longer supported; "
		   "use an InPlane and a RevoluteRotation"
		   << std::endl;
	   THROW(ErrGeneric());
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   THROW(DataManager::ErrGeneric());
	}	      
       }	     
       
       break;
    }
      
    case UNIVERSALPIN:
    case PLANEPIN:
    case REVOLUTEPIN:
    case PLANEDISPLACEMENTPIN: {
       /* nodo collegato */
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode);
       Vec3 d(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel 
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
	  d = HP.GetPosRel(RF);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("Node reference frame d:" << std::endl << d << std::endl);
       
       Mat3x3 Rh(Eye3);
       if (HP.IsKeyWord("orientation")) {
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  Rh = HP.GetRotRel(RF);
#ifdef MBDYN_X_COMPATIBLE_INPUT
       } else if (HP.IsKeyWord("hinge")) {
	  pedantic_cerr("Joint(" << uLabel 
		 << "): keyword \"hinge\" at line " << HP.GetLineData()
		 << " is deprecated; use \"orientation\" instead" << std::endl);
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  Rh = HP.GetRotRel(RF);
	  DEBUGCOUT("Hinge Rotation matrix Rh:" << std::endl << Rh << std::endl);
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       }
       

       Vec3 X0(0.);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       if (HP.IsKeyWord("position")) {
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          X0 = HP.GetPosAbs(AbsRefFrame);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       DEBUGCOUT("Absolute X:" << std::endl << X0 << std::endl);

       Mat3x3 R0(Eye3);
       if (HP.IsKeyWord("orientation")) {
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R0 = HP.GetRotAbs(AbsRefFrame);
#ifdef MBDYN_X_COMPATIBLE_INPUT
       } else if (HP.IsKeyWord("hinge")) {
	  pedantic_cerr("Joint(" << uLabel 
		 << "): keyword \"hinge\" at line " << HP.GetLineData()
		 << " is deprecated; use \"orientation\" instead" << std::endl);
	  DEBUGCOUT("Hinge orientation matrix is supplied" << std::endl);
	  R0 = HP.GetRotAbs(AbsRefFrame);
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       }
       
       DEBUGCOUT("Absolute R:" << std::endl << R0 << std::endl);
       
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       switch (CurrKeyWord) {
	  
	  /* allocazione e creazione cerniera piana */	     
	case PLANEPIN:
		std::cerr << "deprecated 'plane pin' name;" 
			<< " use 'revolute pin' instead" << std::endl;
	case REVOLUTEPIN: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  PlanePinJoint,
				  PlanePinJoint(uLabel, pDO, pNode,
						X0, R0, d, Rh, fOut));
	   break;
	}
	  
	case UNIVERSALPIN: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  UniversalPinJoint,
				  UniversalPinJoint(uLabel, pDO, pNode,
						    X0, R0, d, Rh, fOut));
	   break;
	}
	  
	  /* allocazione e creazione cerniera piana */	     
	case PLANEDISPLACEMENTPIN: {
	   std::cerr << "PlaneDispJoint(" << uLabel << "): "
		   "no longer supported; "
		   "use an InPlane and a RevoluteRotation"
		   << std::endl;
	   THROW(ErrGeneric());
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
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode1);
       Vec3 p(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          p = HP.GetPosRel(RF);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       DEBUGCOUT("Node 1 reference frame p:" << std::endl << p << std::endl);
       
       Vec3 v(HP.GetVecRel(RF));
       doublereal d = v.Dot();
       if (d <= DBL_EPSILON) {
	  std::cerr << "null direction at line " << HP.GetLineData() << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       v /= sqrt(d);
       
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
              
       Vec3 q(0.);
       flag fOffset(0);
       
       if (HP.IsArg()) {
	  if (HP.IsKeyWord("offset")) {
	     fOffset = 1;	     
	     q = HP.GetPosRel(ReferenceFrame(pNode2));	     
	     DEBUGCOUT("Node 2 reference frame q:" << std::endl
			     << p << std::endl);
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);      	   
       
       if (fOffset) {	      
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 InPlaneWithOffsetJoint,
				 InPlaneWithOffsetJoint(uLabel, pDO,
							pNode1, pNode2, 
							v, p, q, fOut));
       } else {	      	   
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 InPlaneJoint,
				 InPlaneJoint(uLabel, pDO, pNode1, pNode2, 
					      v, p, fOut));
       }
       
       break;
    }	
      
    case INLINE: {
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode1);
       Vec3 p(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          p = HP.GetPosRel(RF);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("Node 1 reference frame p:" << std::endl << p << std::endl);
       
       Mat3x3 R(Eye3);
       if (HP.IsKeyWord("orientation")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"orientation\" at line "
			       << HP.GetLineData());
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          R = HP.GetRotRel(RF);
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
              
       Vec3 q(0.);
       flag fOffset(0);
       
       if (HP.IsArg()) {
	  if (HP.IsKeyWord("offset")) {
	     fOffset = 1;	     
	     q = HP.GetPosRel(ReferenceFrame(pNode2));	     
	     DEBUGCOUT("Node 2 reference frame q:" << std::endl << p << std::endl);
	  }
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);      	   
       
       if (fOffset) {
	  SAFENEWWITHCONSTRUCTOR(pEl,
				 InLineWithOffsetJoint,
				 InLineWithOffsetJoint(uLabel, pDO, 
						       pNode1, pNode2,
						       R, p, q, fOut));
       } else {	      	   
	  SAFENEWWITHCONSTRUCTOR(pEl,
				 InLineJoint,
				 InLineJoint(uLabel, pDO, 
					     pNode1, pNode2,
					     R, p, fOut));
       }
       
       break;
    }	
      
      /* asta con pin alle estremita' */
    case ROD: {
       
       /*
	* lettura dei dati specifici:
        * due nodi, lunghezza iniziale e tipo di legame elastico
	*/
       
       flag fOffset(0);
       
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       DEBUGCOUT("Linked to Node " << pNode1->GetLabel() << std::endl);

       Vec3 f1;
       if (HP.IsKeyWord("position")) {
          f1 = HP.GetPosRel(ReferenceFrame(pNode1));
	  fOffset = 1;
       }
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       DEBUGCOUT("Linked to Node " << pNode2->GetLabel() << std::endl);
       
       Vec3 f2;
       if (HP.IsKeyWord("position")) {
          f2 = HP.GetPosRel(ReferenceFrame(pNode2));
	  fOffset = 1;
       }

       /* Lunghezza iniziale */
       doublereal dL = 0.;
       flag fFromNodes(0);
       if (HP.IsKeyWord("from" "nodes")) {
	  fFromNodes = 1;
	  DEBUGCOUT("Initial length will be computed from nodes position" << std::endl);	  
       } else {
	  dL = HP.GetReal();
	  DEBUGCOUT("Initial length = " << dL << std::endl);
       }
       
       /* Se si tratta di Rod con Offset, legge gli offset e poi passa 
	* al tipo di legame costitutivo */
#ifdef MBDYN_X_COMPATIBLE_INPUT
       if (HP.IsKeyWord("offset")) {
	  pedantic_cerr("Joint(" << uLabel 
		 << "): keyword \"offset\" at line " << HP.GetLineData()
		 << " is deprecated; use \"position\" instead" << std::endl);
	  if (fOffset) {
	     std::cerr << "Joint(" << uLabel << "): offsets already defined "
		     "at line " << HP.GetLineData() << std::endl;
	     THROW(ErrGeneric());
	  }
	  fOffset = 1;
	  f1 = HP.GetPosRel(ReferenceFrame(pNode1));	  
	  DEBUGCOUT("Offset 1: " << f1 << std::endl);

	  f2 = HP.GetPosRel(ReferenceFrame(pNode2));	  
	  DEBUGCOUT("Offset 2: " << f2 << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       if (fFromNodes) {
	  Vec3 v = pNode2->GetXCurr()-pNode1->GetXCurr();
	  if (fOffset) {
	     v += pNode2->GetRCurr()*f2-pNode1->GetRCurr()*f1;
	  }
	  dL = v.Norm();
	  
	  if (fOffset) {
	  	pedantic_cout("RodWithOffset(" << uLabel << "): "
				  "length from nodes = " << dL << std::endl);
	  } else {
	  	pedantic_cout("Rod(" << uLabel << "): "
				  "length from nodes = " << dL << std::endl);
	  }

	  DEBUGCOUT("Initial length = " << dL << std::endl);	  
       }
       
       /* Legame costitutivo */
       ConstLawType::Type CLType = ConstLawType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, CLType);
       
       if (pCL->iGetNumDof() != 0) {
	  std::cerr << "line " << HP.GetLineData()
		  << ": rod joint does not support "
		  "dynamic constitutive laws yet"
		  << std::endl;
	  THROW(ErrGeneric());
       }
	
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       if (fOffset == 1) {	      
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 RodWithOffset,
				 RodWithOffset(uLabel, pDO, pCL,
						    pNode1, pNode2, 
						    f1, f2, dL, fOut));
       } else {	  
	  if (CLType == ConstLawType::VISCOUS
			  || CLType == ConstLawType::VISCOELASTIC) {
	     SAFENEWWITHCONSTRUCTOR(pEl, 
				    ViscoElasticRod,
				    ViscoElasticRod(uLabel, pDO, pCL,
							 pNode1, pNode2, 
							 dL, fOut));
	  } else {		 
	     SAFENEWWITHCONSTRUCTOR(pEl, 
				    Rod,
				    Rod(uLabel, pDO, pCL,
					pNode1, pNode2, dL, fOut));
	  }
       }	   
       
       break;
    }

    case RODWITHOFFSET: {

       /*
	* lettura dei dati specifici:
	* due nodi, lunghezza iniziale e tipo di legame elastico.
	* Nota: prima gli offset erano un optional dei dati del rod normale;
	* questo e' lo stesso rod with offset, ma con un input piu' simile
	* a quello standard.
	*/

#ifdef MBDYN_X_COMPATIBLE_INPUT
       pedantic_cerr("Joint(" << uLabel << "): \"rod with offset\" "
		       "is deprecated; use \"rod\" instead" << std::endl);
#endif /* MBDYN_X_COMPATIBLE_INPUT */

       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 f1(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
          f1 = HP.GetPosRel(ReferenceFrame(pNode1));
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
       
       DEBUGCOUT("Linked to Node " << pNode1->GetLabel()
		       << "with offset " << f1 << std::endl);

       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

       Vec3 f2(0.);
       if (HP.IsKeyWord("position")) {
#ifdef MBDYN_X_COMPATIBLE_INPUT
	       NO_OP;
       } else {
	       pedantic_cerr("Joint(" << uLabel
			       << "): missing keyword \"position\" at line "
			       << HP.GetLineData() << std::endl);
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */
	  f2 = HP.GetPosRel(ReferenceFrame(pNode2));
#ifndef MBDYN_X_COMPATIBLE_INPUT
       }
#endif /* MBDYN_X_COMPATIBLE_INPUT */

       DEBUGCOUT("Linked to Node " << pNode2->GetLabel()
		       << "with offset " << f2 << std::endl);

       /* Lunghezza iniziale */
       doublereal dL = 0.;
       if (HP.IsKeyWord("from" "nodes")) {
          dL = (pNode2->GetXCurr()-pNode1->GetXCurr()
		+pNode2->GetRCurr()*f2-pNode1->GetRCurr()*f1).Norm();

	  pedantic_cout("RodWithOffset(" << uLabel << "): "
			  "length from nodes = " << dL << std::endl);

       } else {
          dL = HP.GetReal();
       }
       DEBUGCOUT("Initial length = " << dL << std::endl);

       /* Legame costitutivo */
       ConstLawType::Type CLType = ConstLawType::UNKNOWN;
       ConstitutiveLaw1D* pCL = pDM->ReadConstLaw1D(HP, CLType);

       if (pCL->iGetNumDof() != 0) {
	  std::cerr << "line " << HP.GetLineData()
		  << ": rod with offset joint does not support "
		  "dynamic constitutive laws yet"
		  << std::endl;
	  THROW(ErrGeneric());
       }
	
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);

       SAFENEWWITHCONSTRUCTOR(pEl,
                              RodWithOffset,
                              RodWithOffset(uLabel, pDO, pCL,
                                            pNode1, pNode2,
                                            f1, f2, dL, fOut));

       break;
    }
      
    case DEFORMABLEHINGE:
    case DEFORMABLEDISPHINGE: {
       /* lettura dei dati specifici */	  
       
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       /* Offset if displacement hinge */
       Vec3 f1(0.);
       if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	  f1 = HP.GetPosRel(ReferenceFrame(pNode1));
       }	   
       
       Mat3x3 R1(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  R1 = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

	   
       /* Offset if displacement hinge */
       Vec3 f2(Zero3);
       if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	  f2 = HP.GetPosRel(ReferenceFrame(pNode2));
       }	   
       
       Mat3x3 R2(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  R2 = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       
       
       /* Legame costitutivo */
       ConstLawType::Type CLType;
       ConstitutiveLaw3D* pCL = pDM->ReadConstLaw3D(HP, CLType);
       
       if (pCL->iGetNumDof() != 0) {
	  std::cerr << "line " << HP.GetLineData()
		  << ": deformable hinge joint does not support "
		  "dynamic constitutive laws yet"
		  << std::endl;
	  THROW(ErrGeneric());
       }
	
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       switch (CLType) {
	case ConstLawType::ELASTIC: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ElasticHingeJoint,
				     ElasticHingeJoint(uLabel, pDO, pCL,
						       pNode1, pNode2, 
						       R1, R2, fOut));
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ElasticDispHingeJoint,
				     ElasticDispHingeJoint(uLabel, pDO, pCL,
							   pNode1, pNode2, 
							   f1, f2, R1, R2, 
							   fOut));
	   }
	   
	   break;
	}
	  
	case ConstLawType::VISCOUS: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscousHingeJoint,
				     ViscousHingeJoint(uLabel, pDO, pCL,
						       pNode1, pNode2, 
						       R1, R2, fOut));
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscousDispHingeJoint,
				     ViscousDispHingeJoint(uLabel, pDO, pCL,
							   pNode1, pNode2,
							   f1, f2, R1, R2, 
							   fOut));
	   }
	   
	   break;
	}
	  
	case ConstLawType::VISCOELASTIC: {
	   if (CurrKeyWord == DEFORMABLEHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     ViscoElasticHingeJoint,
				     ViscoElasticHingeJoint(uLabel, pDO, pCL,
							    pNode1, pNode2, 
							    R1, R2, fOut));
	   } else if (CurrKeyWord == DEFORMABLEDISPHINGE) {
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ViscoElasticDispHingeJoint,
				     ViscoElasticDispHingeJoint(uLabel, 
								pDO, pCL,
								pNode1, pNode2,
								f1, f2, R1, R2,
								fOut));
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
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 Dir(HP.GetVecRel(ReferenceFrame(pNode)));
       doublereal d = Dir.Dot();
       ASSERT(d > 0.);	 
       if (d > 0.) {	      
	  Dir /= sqrt(d);
       } else {
	  std::cerr << "direction is null" << std::endl;
	  THROW(ErrGeneric());
       }	     
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       switch (CurrKeyWord) {
	case LINEARVELOCITY: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  LinearVelocityJoint,
				  LinearVelocityJoint(uLabel, pDO,
						      pNode, Dir, pDC, fOut));
	   break;
	}
	  
	case ANGULARVELOCITY: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  AngularVelocityJoint,
				  AngularVelocityJoint(uLabel, pDO,
						       pNode, Dir, pDC, fOut));
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
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Vec3 Dir(HP.GetVecRel(ReferenceFrame(pNode)));
       doublereal d = Dir.Dot();
       ASSERT(d > 0.);	 
       if (d > 0.) {	      
	  Dir /= sqrt(d);
       } else {
	  std::cerr << "direction is null" << std::endl;
	  THROW(ErrGeneric());
       }	     
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       switch (CurrKeyWord) {
	case LINEARACCELERATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl, 
				  LinearAccelerationJoint,
				  LinearAccelerationJoint(uLabel, pDO,
							  pNode, Dir, pDC, 
							  fOut));
	   break;
	}
	  
	case ANGULARACCELERATION: {
	   SAFENEWWITHCONSTRUCTOR(pEl,
				  AngularAccelerationJoint,
				  AngularAccelerationJoint(uLabel, pDO,
							   pNode, Dir, pDC, 
							   fOut));
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
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << std::endl);	   
	  R1h = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       /* Stessa cosa per il nodo 2 */
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << std::endl);
	  R2h = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      PrismaticJoint,
			      PrismaticJoint(uLabel, pDO, pNode1, pNode2, 
					     R1h, R2h, fOut));
       
       break;
    }
      

    case DRIVEHINGE: {
       /* nodo collegato 1 */
       StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       Mat3x3 R1h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << std::endl);	   
	  R1h = HP.GetRotRel(ReferenceFrame(pNode1));
       }
       
       
       /* nodo collegato 2 */
       StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       /* Stessa cosa per il nodo 2 */
       
       Mat3x3 R2h(Eye3);
       if (HP.IsKeyWord("hinge")) {
	  DEBUGCOUT("Hinge Rotation matrix is supplied" << std::endl);
	  R2h = HP.GetRotRel(ReferenceFrame(pNode2));
       }
       

       TplDriveCaller<Vec3>* pDC 
	 = ReadTplDrive(pDM, HP, pDM->pGetDrvHdl(), Vec3(0.));
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      DriveHingeJoint,
			      DriveHingeJoint(uLabel, pDO, pDC, 
					      pNode1, pNode2, 
					      R1h, R2h, fOut));
       
       
       break;
    }
      
      
    case KINEMATIC: {
       /* nodo collegato */
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       
       Kinematics* pK = NULL;
       SAFENEWWITHCONSTRUCTOR(pK, KinematicsTest, KinematicsTest(pDC));
       
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      KinJoint,
			      KinJoint(uLabel, pDO, pNode, pK, fOut));
       
       break;
    }

      
    case MODAL: {
       pEl = ReadModal(pDM, HP, pDO, uLabel);
       break;
    }

    case BEAMSLIDER: {
       /* Corpo slittante */
       StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
       
       ReferenceFrame RF(pNode);
       Vec3 f(HP.GetPosRel(RF));
       DEBUGCOUT("Linked to Node " << pNode->GetLabel()
		       << "with offset " << f << std::endl);
       
       Mat3x3 R = Eye3;
       if (HP.IsKeyWord("hinge")) {
	       R = HP.GetRotRel(RF);
       }
       DEBUGLCOUT(MYDEBUG_INPUT,
		       "Slider rotation matrix: " << std::endl << R << std::endl);

       /* Slider type */
       BeamSliderJoint::Type sliderType = BeamSliderJoint::SPHERICAL;
       if (HP.IsKeyWord("type")) {
	       if (HP.IsKeyWord("spherical")) {
		       sliderType = BeamSliderJoint::SPHERICAL;
	       } else if (HP.IsKeyWord("classic")) {
		       sliderType = BeamSliderJoint::CLASSIC;
	       } else if (HP.IsKeyWord("spline")) {
		       sliderType = BeamSliderJoint::SPLINE;
	       } else {
		       std::cerr << "unknown slider type at line " 
			       << HP.GetLineData() << std::endl;
		       THROW(ErrGeneric());
	       }
       }

       unsigned int nB = HP.GetInt();
       if (nB < 1) {
	       std::cerr << "At least one beam is required at line " 
		       << HP.GetLineData() << std::endl;
	       THROW(ErrGeneric());
       }
       
       BeamConn **bc = NULL;
       SAFENEWARR(bc, BeamConn *, nB);

       const StructNode* pLastNode = NULL;
       for (unsigned int i = 0; i < nB; i++) {
          /* trave di appoggio */
	  unsigned int uBeam = HP.GetInt();
       
	  Beam* pBeam;
	  Elem* p = (Elem*)(pDM->pFindElem(Elem::BEAM, uBeam));
	  if (p == NULL) {
		  std::cerr << " at line " << HP.GetLineData()
			  << ": Beam(" << uBeam
			  << ") not defined" << std::endl;
		  THROW(DataManager::ErrGeneric());
	  }
	  pBeam = (Beam*)p->pGet();
       
	  /* Nodo 1: */
	  
	  /* Offset del punto rispetto al nodo */
	  const StructNode* pNode1 = pBeam->pGetNode(1);
	  RF = ReferenceFrame(pNode1);
	  Vec3 f1;
	  Mat3x3 R1 = Eye3;
	  if (i != 0) {
		  if (pNode1 != pLastNode) {
			std::cerr 
				<< "line " << HP.GetLineData()
				<< ": Beam(" << pBeam->GetLabel() 
				<< ").Node1(" << pNode1->GetLabel()
 				<< ") and Beam(" 
				<< bc[i-1]->pGetBeam()->GetLabel() 
				<< ").Node3(" << pLastNode->GetLabel()
				<< ") do not match" << std::endl;
			THROW(ErrGeneric());
		  }
		  
  		  if (HP.IsKeyWord("same")) {
			  f1 = bc[i-1]->Getf(3);
		  } else {
			  f1 = HP.GetPosRel(RF);
			  if (f1 != bc[i-1]->Getf(3)) {
				std::cerr 
					<< "line " << HP.GetLineData()
					<< ": Beam(" << pBeam->GetLabel()
					<< ").f1 and Beam("
					<< bc[i-1]->pGetBeam()->GetLabel()
					<< ").f3 do not match" << std::endl;
				THROW(ErrGeneric());
			  }
		  }
	  	  DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

	  	  if (HP.IsKeyWord("hinge")) {
  			  if (HP.IsKeyWord("same")) {
				  R1 = bc[i-1]->GetR(3);
			  } else {
				  R1 = HP.GetRotRel(RF);
				  if (R1 != bc[i-1]->GetR(3)) {
					std::cerr 
						<< "line " << HP.GetLineData()
						<< ": Beam(" << pBeam->GetLabel()
						<< ").R1 and Beam("
						<< bc[i-1]->pGetBeam()->GetLabel()
						<< ").R3 do not match" << std::endl;
					THROW(ErrGeneric());
				  }
			  }
		  	  DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 rotation matrix: " 
					  << std::endl << R1 << std::endl);
		  }

	  } else {
		f1 = HP.GetPosRel(RF);
		if (HP.IsKeyWord("hinge")) {
			R1 = HP.GetRotRel(RF);
		  	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 rotation matrix: " 
					<< std::endl << R1 << std::endl);
		}
	  }
       
	  /* Nodo 2: */
       
	  /* Offset del punto rispetto al nodo */
	  const StructNode* pNode2 = pBeam->pGetNode(2);
       
	  RF = ReferenceFrame(pNode2);
	  Vec3 f2(HP.GetPosRel(RF));
	  DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);
       
	  Mat3x3 R2(Eye3);
	  if (HP.IsKeyWord("hinge")) {
		  R2 = HP.GetRotRel(RF);
		  DEBUGLCOUT(MYDEBUG_INPUT,
				  "Node 2 rotation matrix: " << std::endl 
				  << R2 << std::endl);
	  }
       
	  /* Nodo 3: */
       
	  /* Offset del punto rispetto al nodo */
	  const StructNode* pNode3 = pBeam->pGetNode(3);
       
	  RF = ReferenceFrame(pNode3);
	  Vec3 f3(HP.GetPosRel(RF));
	  DEBUGLCOUT(MYDEBUG_INPUT, "Node 3 offset: " << f3 << std::endl);
	  
	  Mat3x3 R3(Eye3);
	  if (HP.IsKeyWord("hinge")) {
		  R3 = HP.GetRotRel(RF);
		  DEBUGLCOUT(MYDEBUG_INPUT,
				  "Node 3 rotation matrix: " << std::endl 
				  << R3 << std::endl);
	  }
	  
	  pLastNode = pNode3;

	  bc[i] = NULL;
	  SAFENEWWITHCONSTRUCTOR(bc[i], BeamConn, 
			  BeamConn(pBeam, f1, f2, f3, R1, R2, R3));
       }

       unsigned uIB = 1;
       if (HP.IsKeyWord("initial" "beam")) {
	       uIB = HP.GetInt();

	       if (uIB < 1 || uIB > nB) {
		       std::cerr << "illegal initial beam " << uIB
			       << " at line " << HP.GetLineData() << std::endl;
		       THROW(ErrGeneric());
	       }
       }
      
       unsigned uIN = 2;
       if (HP.IsKeyWord("initial" "node")) {
	       uIN = HP.GetInt();

	       if (uIN < 1 || uIN > 3) {
		       std::cerr << "illegal initial node " << uIN
			       << " at line " << HP.GetLineData() << std::endl;
		       THROW(ErrGeneric());
	       }
       }

       doublereal dL = 0.;
       if (HP.IsKeyWord("smearing")) {
	       dL = HP.GetReal();
	       if (dL < 0. || dL > .4) {
		       std::cerr << "illegal smearing factor " << dL 
			       << "; using default" << std::endl;
		       dL = 0.;
	       }
       }
      
       flag fOut = pDM->fReadOutput(HP, Elem::JOINT);
       SAFENEWWITHCONSTRUCTOR(pEl, BeamSliderJoint,
		       BeamSliderJoint(uLabel, pDO, 
			       pNode, 
			       sliderType,
			       nB, bc,
			       uIB, uIN, dL,
			       f, R, fOut));
       break;
    }


    case WHEEL2: {
#ifdef __HACK_MODULES_AS_JOINTS__
       LoadableElem *pLE = NULL;
       SAFENEWWITHCONSTRUCTOR(pLE, LoadableElem,
		       LoadableElem(uLabel, pDO, &module_wheel2_lc, pDM, HP));
       pEl = (Joint *)pLE->pGet();
       pDM->OutputOpen(OutputHandler::LOADABLE);
#endif /* __HACK_MODULES_AS_JOINTS__ */
       break;
    } 
      /* Aggiungere altri vincoli */
      
    default: {
       std::cerr << "unknown joint type in joint " << uLabel
	 << " at line " << HP.GetLineData() << std::endl;       
       THROW(DataManager::ErrGeneric());
    }
   }

   /* Se non c'e' il punto e virgola finale */
   if (HP.IsArg()) {
      std::cerr << "semicolon expected at line " << HP.GetLineData() << std::endl;
      THROW(DataManager::ErrGeneric());
   }

   return pEl;
} /* ReadJoint() */
