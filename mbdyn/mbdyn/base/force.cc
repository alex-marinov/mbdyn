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

/* Forze */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/float.h>

#include <force.h>
#ifdef USE_STRUCT_NODES
#include <strforce.h>
#endif /* USE_STRUCT_NODES */
#include <dataman.h>

/* Force - begin */

std::ostream&
Force::Output(unsigned int NodeLabel, std::ostream& out) const
{
   return out << std::setw(8) << GetLabel() << " "
     << std::setw(8) << NodeLabel << " " 
     << (pGetDriveCaller()->dGet());
}


std::ostream&
Force::Restart(std::ostream& out) const
{
   return out << "  force: " << GetLabel();
}

/* Force - end */


/* AbstractForce - begin */

/* Costruttore non banale */

AbstractForce::AbstractForce(unsigned int uL, const Node* pN, 
		const DriveCaller* pDC, flag fOut)
: Elem(uL, Elem::FORCE, fOut),
Force(uL, Force::ABSTRACTFORCE, pDC, fOut),
pNode(pN)
{
   NO_OP;
}


/* Contributo al file di restart */
std::ostream&
AbstractForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", abstract, " 
     << pNode->GetLabel() << ", " 
     << psReadNodesNodes[pNode->GetNodeType()] << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


/* Assembla il residuo */
SubVectorHandler&
AbstractForce::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering AbstractForce::AssRes()" << std::endl);

   WorkVec.ResizeReset(1);

   /* Dati */   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstIndex = pNode->iGetFirstRowIndex();
   WorkVec.PutRowIndex(1, iFirstIndex+1);

   WorkVec.PutCoef(1, dAmplitude);

   return WorkVec;
}

void
AbstractForce::Output(OutputHandler& OH) const 
{
   if (fToBeOutput()) {
      Force::Output(pNode->GetLabel(), OH.Forces()) << std::endl;
   }
}

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AbstractForce::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering AbstractForce::InitialAssRes()" << std::endl);

   return AssRes(WorkVec, 1., XCurr, XCurr);
}

/* AbstractForce - end */


/* AbstractInternalForce - begin */

/* Costruttore non banale */

AbstractInternalForce::AbstractInternalForce(unsigned int uL,
		const Node* pN1, const Node *pN2, 
		const DriveCaller* pDC, flag fOut)
: Elem(uL, Elem::FORCE, fOut),
Force(uL, Force::ABSTRACTFORCE, pDC, fOut),
pNode1(pN1), pNode2(pN2)
{
   NO_OP;
}


/* Contributo al file di restart */
std::ostream&
AbstractInternalForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", abstract internal, " 
     << pNode1->GetLabel() << ", " 
     << psReadNodesNodes[pNode1->GetNodeType()] << ", "
     << pNode1->GetLabel() << ", " 
     << psReadNodesNodes[pNode1->GetNodeType()] << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


/* Assembla il residuo */
SubVectorHandler&
AbstractInternalForce::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering AbstractInternalForce::AssRes()" << std::endl);

   WorkVec.ResizeReset(2);

   /* Dati */   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstIndex1 = pNode1->iGetFirstRowIndex();
   integer iFirstIndex2 = pNode2->iGetFirstRowIndex();
   WorkVec.PutRowIndex(1, iFirstIndex1+1);
   WorkVec.PutRowIndex(2, iFirstIndex2+1);

   WorkVec.PutCoef(1, dAmplitude);
   WorkVec.PutCoef(2, -dAmplitude);

   return WorkVec;
}

void
AbstractInternalForce::Output(OutputHandler& OH) const 
{
   if (fToBeOutput()) {
      Force::Output(pNode1->GetLabel(), OH.Forces()) << std::endl;
   }
}

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AbstractInternalForce::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering AbstractInternalForce::InitialAssRes()" << std::endl);

   return AssRes(WorkVec, 1., XCurr, XCurr);
}

/* AbstractInternalForce - end */


/* Legge una forza */

Elem* ReadForce(DataManager* pDM, 
		MBDynParser& HP, 
		unsigned int uLabel, 
		flag fCouple)
{
   const char sFuncName[] = "ReadForce()";
   DEBUGCOUT("Entering " << sFuncName << std::endl);
   
   const char* sKeyWords[] = {

#if defined(USE_STRUCT_NODES)
      "conservative",
      "follower",
      
      "conservative" "internal",
      "follower" "internal",
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
      "abstract",
      "abstract" "internal",
#endif /* USE_ELECTRIC_NODES */

      NULL
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,

#if defined(USE_STRUCT_NODES)
	CONSERVATIVE,
	FOLLOWER,
	
	CONSERVATIVEINTERNAL,
	FOLLOWERINTERNAL,
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
	ABSTRACT,
	ABSTRACTINTERNAL,
#endif /* USE_ELECTRIC_NODES */

	LASTKEYWORD 
   };
   
   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);
   
   /* tipo di forza */
   KeyWords CurrType = KeyWords(HP.GetWord());
   if (CurrType == UNKNOWN) {
      silent_cerr(sFuncName << " at line " << HP.GetLineData() 
	<< ": unknown force type" << std::endl);
      throw DataManager::ErrGeneric();
   }   
   
#ifdef DEBUG
   switch(CurrType) {
#if defined(USE_STRUCT_NODES)
    case CONSERVATIVE:
      std::cout << "Force type: \"Conservative\"" << std::endl;
      break;

    case FOLLOWER:
      std::cout << "Force type: \"Follower\"" << std::endl;
      break;

    case CONSERVATIVEINTERNAL:
      std::cout << "Force type: \"Conservative Internal\"" << std::endl;
      break;

    case FOLLOWERINTERNAL:
      std::cout << "Force type: \"Follower Internal\"" << std::endl;
      break;

#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
    case ABSTRACT:
      std::cout << "Force type: \"Abstract\"" << std::endl;
      break;

    case ABSTRACTINTERNAL:
      std::cout << "Force type: \"Abstract Internal\"" << std::endl;
      break;

#endif /* USE_ELECTRIC_NODES */
      
    default:
      throw ErrGeneric();
   }
#endif /* DEBUG */

   Elem* pEl = NULL;
   
#if defined(USE_ELECTRIC_NODES)
   if (CurrType == ABSTRACT) {
      
      /* tabella delle parole chiave */
      KeyTable KDof(HP, psReadNodesNodes);
      
      ScalarDof SD = ReadScalarDof(pDM, HP, 0);                          
      
      DriveCaller* pDC = HP.GetDriveCaller();

      flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
      
      SAFENEWWITHCONSTRUCTOR(pEl,
			     AbstractForce,
			     AbstractForce(uLabel, SD.pNode, pDC, fOut));
            
   } else if (CurrType == ABSTRACTINTERNAL) {
      
      /* tabella delle parole chiave */
      KeyTable KDof(HP, psReadNodesNodes);
      
      ScalarDof SD1 = ReadScalarDof(pDM, HP, 0);                          
      
      ScalarDof SD2 = ReadScalarDof(pDM, HP, 0);                          
      
      DriveCaller* pDC = HP.GetDriveCaller();

      flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
      
      SAFENEWWITHCONSTRUCTOR(pEl,
		      AbstractInternalForce,
		      AbstractInternalForce(uLabel, SD1.pNode, SD2.pNode,
			      pDC, fOut));
 
   } else
     
#endif /* USE_ELECTRIC_NODES */
     
#if defined(USE_STRUCT_NODES)
     { /* CurrType e' <structural> */
         
      /* nodo collegato */
      StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
      
      /* direzione della forza */   	     
      Mat3x3 RNode(pNode->GetRCurr());     
      ReferenceFrame RF(pNode);
      Vec3 Dir(HP.GetVecRel(RF));

      /* Se la forza e' conservativa, viene passata nel sistema globale */
      if (CurrType == CONSERVATIVE || CurrType == CONSERVATIVEINTERNAL) {
	 Dir = RNode*Dir;
      }      	           
      
      /* Normalizza la direzione */
      ASSERT(Dir.Dot() > DBL_EPSILON);
      doublereal d = Dir.Dot();
      if (d > DBL_EPSILON) {      
	 Dir /= sqrt(d);
      } else {      
	 silent_cerr("Force(" << uLabel 
	   << ") has null direction" << std::endl);
	 throw ErrGeneric();
      }        
      
      /* distanza dal nodo (vettore di 3 elementi) ( solo se e' una forza) */
      Vec3 Arm(0.);
      if (fCouple == 0) {	
	 if (!HP.IsKeyWord("null")) {	  
	    Arm = HP.GetPosRel(RF);	    
	    DEBUGCOUT("Distance is supplied" << std::endl);
	 }
      }

      StructNode *pNode2 = NULL;
      Vec3 Arm2(0.);
      if (CurrType == CONSERVATIVEINTERNAL || CurrType == FOLLOWERINTERNAL) {
	      
         /* nodo collegato */
	 pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
      
   	 /* distanza dal nodo (vettore di 3 elementi) ( solo se e' una forza) */
   	 if (fCouple == 0) {	
	    if (!HP.IsKeyWord("null")) {	  
	       Arm2 = HP.GetPosRel(RF);	    
	       DEBUGCOUT("Distance is supplied" << std::endl);
   	    }
   	 }
      }
 
#ifdef DEBUG
      if (CurrType == CONSERVATIVE || CurrType == CONSERVATIVEINTERNAL) {      
	 std::cout << "Global reference frame direction: " << std::endl
		 << Dir << std::endl;
      } else if (CurrType == FOLLOWER || CurrType == FOLLOWERINTERNAL) {      
	 std::cout << "Node reference frame direction: " << std::endl
		 << Dir << std::endl;
      }
      
      if (!fCouple) {      
	 std::cout << "Node reference frame arm: " << std::endl
		 << Arm << std::endl;
      }   
#endif /* DEBUG */
            
      DriveCaller* pDC = HP.GetDriveCaller();
      
      flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
      
      /* Alloca la forza */
      if (fCouple == 0) {
	 switch (CurrType) {
	  case CONSERVATIVE:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    ConservativeForce,
			    ConservativeForce(uLabel, pNode, pDC, 
				    Dir, Arm, fOut));
	    break;

	  case CONSERVATIVEINTERNAL:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    ConservativeInternalForce,
			    ConservativeInternalForce(uLabel, pNode, pNode2,
				    pDC, Dir, Arm, Arm2, fOut));
	    break;

	  case FOLLOWER:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    FollowerForce,
			    FollowerForce(uLabel, pNode, pDC, 
				    Dir, Arm, fOut));
	    break;

	  case FOLLOWERINTERNAL:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    FollowerInternalForce,
			    FollowerInternalForce(uLabel, pNode, pNode2,
				    pDC, Dir, Arm, Arm2, fOut));
	    break;

	  default:
	    ASSERT(0);
	 }

      } else if (fCouple == 1) {
	 switch (CurrType) {
	  case  CONSERVATIVE:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    ConservativeCouple,
			    ConservativeCouple(uLabel, pNode, pDC, 
				    Dir, fOut));
	    break;

	  case CONSERVATIVEINTERNAL:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    ConservativeInternalCouple,
			    ConservativeInternalCouple(uLabel, pNode, pNode2,
				    pDC, Dir, fOut));
	    break;

	  case FOLLOWER:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    FollowerCouple,
			    FollowerCouple(uLabel, pNode, pDC, 
				    Dir, fOut));
	    break;

	  case FOLLOWERINTERNAL:
	    SAFENEWWITHCONSTRUCTOR(pEl, 
			    FollowerInternalCouple,
			    FollowerInternalCouple(uLabel, pNode, pNode2,
				    pDC, Dir, fOut));
	    break;

	  default:
	    ASSERT(0);
	 }
      }
   }
#endif /* USE_STRUCT_NODES */
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.IsArg()) {
      silent_cerr(sFuncName
	<< ": semicolon expected at line " << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric();
   }      
   
   return pEl;
} /* End of ReadForce() */

