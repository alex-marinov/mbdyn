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

/* Forza */

#ifndef FORCE_H
#define FORCE_H

#include "elem.h"
#include "drive.h"
#include "strnode.h"
#include "elecnode.h"

extern const char* psForceNames[];


/* Force - begin */

class Force 
: virtual public Elem, public InitialAssemblyElem, public DriveOwner {
 public:
   /* Tipi di Force */
   enum Type {
      UNKNOWN = -1,
	ABSTRACTFORCE = 0,
	
	CONSERVATIVEFORCE,
	FOLLOWERFORCE,
	CONSERVATIVECOUPLE,
	FOLLOWERCOUPLE,
	
	LASTFORCETYPE
   };

 private:
   Force::Type ForceT;
   
 public:
   /* Costruttore banale */
   Force(unsigned int uL, Force::Type T, 
	 const DriveCaller* pDC, flag fOut)
     : Elem(uL, Elem::FORCE, fOut), 
     InitialAssemblyElem(uL, Elem::FORCE, fOut), 
     DriveOwner(pDC), ForceT(T) { 
	NO_OP; 
     };
      
   virtual ~Force(void) { 
      NO_OP; 
   };
         
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const { 
      return Elem::FORCE; 
   };   
   
   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const = 0;
   
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ ,
	    const VectorHandler& /* XCurr */ , 
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUT("Entering Force::AssJac()" << endl);
	
	WorkMat.SetNullMatrix();
	return WorkMat;
     };

   /* Output comune a tutti i tipi di forza;
    * scrive su file il valore del drive (temporaneo) */
   ostream& Output(unsigned int NodeLabel, ostream& out) const;
   
   virtual ostream& Restart(ostream& out) const;

   virtual unsigned int iGetInitialNumDof(void) const { 
      return 0;
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& /* XCurr */ ) {
	WorkMat.SetNullMatrix();
	return WorkMat;
     };
};

/* Force - end */


/* StructuralForce - begin */

class StructuralForce : virtual public Elem, public Force {
 protected:
   const StructNode* pNode;
   const Vec3 Dir;
   
 public:
   /* Costruttore */
   StructuralForce(unsigned int uL, Force::Type T, 
		   const StructNode* pN,
		   const DriveCaller* pDC, const Vec3& TmpDir,
		   flag fOut);

   virtual ~StructuralForce(void);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();
   };
   /* ************************************************ */
};

/* StructuralForce - end */


/* AbstractForce - begin */

class AbstractForce : virtual public Elem, public Force {
 protected:
   const Node* pNode;
   // const integer iDofNumber;
   
 public:
   /* Costruttore banale */
   AbstractForce(unsigned int uL, const Node* pN, 
		 const DriveCaller* pDC, flag fOut);
      
   virtual ~AbstractForce(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
         
   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const { 
      return Force::ABSTRACTFORCE;
   };

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 1; 
      *piNumCols = 1; 
   };

   /* Contributo al residuo */
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);

   virtual void Output(OutputHandler& OH) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows, 
				    integer* piNumCols) const {
      *piNumRows = 1;
      *piNumCols = 1;
   };

   /* Contributo al residuo nell'assemblaggio iniziale */
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();
   };
   /* ************************************************ */
};

/* AbstractForce - end */


/* ConservativeForce - begin */

class ConservativeForce : virtual public Elem, public StructuralForce {
 protected:
   const Vec3 Arm;
   
 public:
   /* Costruttore non banale */
   ConservativeForce(unsigned int uL, const StructNode* pN, 
		     const DriveCaller* pDC, 
		     const Vec3& TmpDir, const Vec3& TmpArm, 
		     flag fOut);
      
   ~ConservativeForce(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
     
   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const { 
      return Force::CONSERVATIVEFORCE; 
   };   
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 6; 
      *piNumCols = 3; 
   };

   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);     

   virtual void Output(OutputHandler& OH) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 12; 
      *piNumCols = 6; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);

   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);   
};

/* ConservativeForce - end */


/* FollowerForce - begin */

class FollowerForce : virtual public Elem, public StructuralForce {
 protected:
   const Vec3 Arm;
   
 public:
   /* Costruttore banale */
   FollowerForce(unsigned int uL, const StructNode* pN, 
		 const DriveCaller* pDC, 
		 const Vec3& TmpDir, const Vec3& TmpArm,
		 flag fOut);
      
   ~FollowerForce(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { return (void*)this; };
     
   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const {
      return Force::FOLLOWERFORCE; 
   };

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 6;
      *piNumCols = 3;
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);     

   virtual void Output(OutputHandler& OH) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 12; 
      *piNumCols = 6; 
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);
};

/* FollowerForce - end */


/* ConservativeCouple - begin */

class ConservativeCouple : virtual public Elem, public StructuralForce {
   
 public:
   /* Costruttore banale */
   ConservativeCouple(unsigned int uL, const StructNode* pN, 
		      const DriveCaller* pDC, 
		      const Vec3& TmpDir,
		      flag fOut);

   ~ConservativeCouple(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { return (void*)this; };
     
   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const {
      return Force::CONSERVATIVECOUPLE; 
   };

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 3;
      *piNumCols = 1;
   };

   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);     

   virtual void Output(OutputHandler& OH) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows,integer* piNumCols) const { 
      *piNumRows = 3;
      *piNumCols = 1; 
   };

   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);
};

/* ConservativeCouple - end */


/* FollowerCouple - begin */

class FollowerCouple : virtual public Elem, public StructuralForce {
   
 public:
   /* Costruttore banale */
   FollowerCouple(unsigned int uL, const StructNode* pN, 
		  const DriveCaller* pDC, 
		  const Vec3& TmpDir,
		  flag fOut);
      
   ~FollowerCouple(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { return (void*)this; };

   /* Tipo di forza */
   virtual Force::Type GetForceType(void) const {
      return Force::FOLLOWERCOUPLE; 
   };   

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 3; 
      *piNumCols = 3; 
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);

   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);     

   virtual void Output(OutputHandler& OH) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 6;
      *piNumCols = 6;
   };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);
};

/* FollowerCouple - end */

class DataManager;
class MBDynParser;

extern Elem* ReadForce(DataManager* pDM, 
		       MBDynParser& HP, 
		       unsigned int uLabel, 
		       flag fCouple);

#endif
