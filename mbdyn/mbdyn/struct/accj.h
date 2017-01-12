/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* Giunti di accelerazione imposta */


#ifndef ACCJ_H
#define ACCJ_H

#include "joint.h"
#include "drive.h"

/* LinearAccelerationJoint - begin */

class LinearAccelerationJoint 
: virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode;
   Vec3 Dir;
   doublereal dF;
   
 public:
   /* Costruttore non banale */
   LinearAccelerationJoint(unsigned int uL, 
			   const DofOwner* pDO, 
			   const StructNode* pN,
			   const Vec3& TmpDir,
			   const DriveCaller* pDC,
			   flag fOut);
   
   /* Distruttore */
   ~LinearAccelerationJoint(void);
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const;
   DofOrder::Order GetDofType(unsigned int i) const;
   
   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void Output(OutputHandler& OH) const;
 
	virtual void SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* Xp */ ,
		SimulationEntity::Hints *ph = 0);

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const;
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const;
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual unsigned int iGetPrivDataIdx(const char *s) const;
   virtual doublereal dGetPrivData(unsigned int i) const;

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode;
     
   };
   /* ************************************************ */   
};

/* LinearAccelerationJoint - end */


/* AngularAccelerationJoint - begin */

class AngularAccelerationJoint 
: virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode;
   Vec3 Dir;
   doublereal dM;
   
 public:
   /* Costruttore non banale */
   AngularAccelerationJoint(unsigned int uL, 
			    const DofOwner* pDO,
			    const StructNode* pN,
			    const Vec3& TmpDir,
			    const DriveCaller* pDC,
			    flag fOut);
   
   /* Distruttore */
   ~AngularAccelerationJoint(void);
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const;
   DofOrder::Order GetDofType(unsigned int i) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void Output(OutputHandler& OH) const;
 
	virtual void SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* Xp */ ,
		SimulationEntity::Hints *ph = 0);

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const;
   
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const;
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual unsigned int iGetPrivDataIdx(const char *s) const;
   virtual doublereal dGetPrivData(unsigned int i) const;

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode;
     
   };
   /* ************************************************ */   
};

/* AngularAccelerationJoint - end */

#endif /* ACCJ_H */

