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

/* Vincoli generali */


#ifndef GENJ_H
#define GENJ_H

#include "joint.h"
#include "drive.h"


/* DistanceJoint - begin */

class DistanceJoint : virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 v;
   doublereal dAlpha;
   
 public:
   /* Costruttore non banale */
   DistanceJoint(unsigned int uL, const DofOwner* pDO,
		 const StructNode* pN1, const StructNode* pN2,
		 const DriveCaller* pDC, flag fOut);
   
   ~DistanceJoint(void);
   
   virtual inline void* pGet(void) const { return (void*)this; };

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::DISTANCE; 
   };
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 4;
   };

#ifdef DEBUG
   virtual DofOrder::Order SetDof(unsigned int i) const
#else
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i >= 0 && i < 4);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 10; *piNumCols = 10; };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void Output(OutputHandler& OH) const;

   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { return 8; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const 
     { *piNumRows = 20; *piNumCols = 20; };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X) const;
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const ;

   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const {
      return 1;
   };   

#ifdef DEBUG
   virtual doublereal dGetPrivData(unsigned int i) const
#else
   virtual doublereal dGetPrivData(unsigned int /* i */ = 0) const
#endif
   {
      ASSERT(i == 1);
      return dGet();
   };

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */
   
   /* Adams output stuff */
   virtual unsigned int iGetNumAdamsDummyParts(void) const {
      return 1;
   };
   virtual void GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& R) const;
   virtual ostream& WriteAdamsDummyPartCmd(ostream& out, unsigned int part, unsigned int firstId) const;
};

/* DistanceJoint - end */


/* DistanceJointWithOffset - begin */

class DistanceJointWithOffset : 
virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 f1;
   Vec3 f2;
   Vec3 v;
   doublereal dAlpha;
   
 public:
   /* Costruttore non banale */
   DistanceJointWithOffset(unsigned int uL, const DofOwner* pDO,
			   const StructNode* pN1, const StructNode* pN2,
			   const Vec3& f1Tmp, const Vec3& f2Tmp,
			   const DriveCaller* pDC, flag fOut);
   
   ~DistanceJointWithOffset(void);
   
   virtual inline void* pGet(void) const { return (void*)this; };

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::DISTANCEWITHOFFSET; 
   };
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 4;
   };
#ifdef DEBUG      
   virtual DofOrder::Order SetDof(unsigned int i) const
#else
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i >= 0 && i < 4);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 16;
      *piNumCols = 16; 
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

   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { return 8; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const 
     { *piNumRows = 32; *piNumCols = 32; };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X) const;
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const ;

   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const {
      return 1;
   };   

#ifdef DEBUG
   virtual doublereal dGetPrivData(unsigned int i = 0) const
#else
   virtual doublereal dGetPrivData(unsigned int /* i */ = 0) const
#endif
   {
      ASSERT(i == 1);
      return dGet();
   };   

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */
   
   /* Adams output stuff */
   virtual unsigned int iGetNumAdamsDummyParts(void) const {
      return 1;
   };
   virtual void GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& R) const;
   virtual ostream& WriteAdamsDummyPartCmd(ostream& out, unsigned int part, unsigned int firstId) const;
};

/* DistanceJointWithOffset - end */


/* ClampJoint - begin */

class ClampJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode;    /* nodo incastrato */
   Vec3 XClamp;                /* posizione imposta */
   Mat3x3 RClamp;              /* assetto imposto */
   Vec3 F;                     /* forza di reazione */
   Vec3 M;                     /* momento di reazione */
   
 public:
   /* Costruttore definitivo (da mettere a punto) */
   ClampJoint(unsigned int uL, const DofOwner*pD, const StructNode* pN, 
	      const Vec3& X0, const Mat3x3& R0, flag fOut);
   
   /* Distruttore + o - banale */
   virtual ~ClampJoint(void);
   
   virtual inline void* pGet(void) const { return (void*)this; };

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::CLAMP; 
   };

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   /* Funzioni obbligatorie, per la gestione dei dof */
   virtual unsigned int iGetNumDof(void) const {
      return 6; 
   };

#ifdef DEBUG
   virtual DofOrder::Order SetDof(unsigned int i) const
#else
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i >= 0 && i < 6);
      return DofOrder::ALGEBRAIC; 
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
     { *piNumRows = 12; *piNumCols = 12; };
      
   /* Assemblaggio matrice jacobiana */
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);

   /* assemblaggio matrici per autovalori */
   void AssEig(VariableSubMatrixHandler& WorkMatA,
	       VariableSubMatrixHandler& WorkMatB,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
   
   /* Assemblaggio residuo */
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void Output(OutputHandler& OH) const;


   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { return 12; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const 
     { *piNumRows = 24; *piNumCols = 24; };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Metodi per l'estrazione di dati "privati".
    * Si suppone che l'estrattore li sappia interpretare.
    * Come default non ci sono dati privati estraibili */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual doublereal dGetPrivData(unsigned int i) const;

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

/* ClampJoint - end */

#endif
