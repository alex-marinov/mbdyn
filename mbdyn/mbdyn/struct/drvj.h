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

/* Giunti di velocita' imposta */


#ifndef DRVJ_H
#define DRVJ_H

#include "joint.h"
#include "drive.h"

/* LinearVelocityJoint - begin */

class LinearVelocityJoint 
: virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode;
   Vec3 Dir;
   doublereal dF;
   
 public:
   /* Costruttore non banale */
   LinearVelocityJoint(unsigned int uL, 
		       const DofOwner* pDO, 
		       const StructNode* pN,
		       const Vec3& TmpDir,
		       const DriveCaller* pDC,
		       flag fOut);
   
   /* Distruttore */
   ~LinearVelocityJoint(void);
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const 
     { return Joint::LINEARVELOCITY; };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
   
#ifdef DEBUG
   DofOrder::Order SetDof(unsigned int i) const
#else
   DofOrder::Order SetDof(unsigned int /* i */ ) const     
#endif     
   {
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC; 
   };
   
   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 4; *piNumCols = 4; };
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void Output(OutputHandler& OH) const;
 

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { return 1; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const { 
      *piNumRows = 4; 
      *piNumCols = 4; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);   

   /* dati privati */
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
     NumNodes = 1;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();
   };
   /* ************************************************ */

};

/* LinearVelocityJoint - end */


/* AngularVelocityJoint - begin */

class AngularVelocityJoint 
: virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode;
   Vec3 Dir;
   doublereal dM;
   
 public:
   /* Costruttore non banale */
   AngularVelocityJoint(unsigned int uL, 
			const DofOwner* pDO,
			const StructNode* pN,
			const Vec3& TmpDir,
			const DriveCaller* pDC,
			flag fOut);
   
   /* Distruttore */
   ~AngularVelocityJoint(void);
   
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const 
     { return Joint::ANGULARVELOCITY; };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
   
#ifdef DEBUG
   DofOrder::Order SetDof(unsigned int i) const
#else
   DofOrder::Order SetDof(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC; 
   };

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 4; *piNumCols = 4; };
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void Output(OutputHandler& OH) const;
 

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 1;
   };
   
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const { 
      *piNumRows = 4; 
      *piNumCols = 7; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);   

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
     NumNodes = 1;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();
   };
   /* ************************************************ */ 
};

/* AngularVelocityJoint - end */

#endif
