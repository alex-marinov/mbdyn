/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Giunti sferici */

#ifndef SPHERJ_H
#define SPHERJ_H

#include "joint.h"


/* SphericalHingeJoint - begin */

class SphericalHingeJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
#ifdef USE_NETCDF
	NcVar *Var_Phi;
#endif // USE_NETCDF
   Vec3 d1;
   Mat3x3 R1h;
   Vec3 d2;
   Mat3x3 R2h;
   Vec3 F;

 protected:
	OrientationDescription od;

 public:
   /* Costruttore non banale */
   SphericalHingeJoint(unsigned int uL, const DofOwner* pDO,
		       const StructNode* pN1, const StructNode* pN2,
		       const Vec3& dTmp1, const Mat3x3& RTmp1h,
		       const Vec3& dTmp2, const Mat3x3& RTmp2h,
		       const OrientationDescription& od,
		       flag fOut);
   
   ~SphericalHingeJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::SPHERICALHINGE; 
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 3;
   };
   
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 3);
      return DofOrder::ALGEBRAIC;
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 15; 
      *piNumCols = 15; 
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
			    
   DofOrder::Order GetEqType(unsigned int i) const;
   
   void OutputPrepare(OutputHandler &OH);
   virtual void Output(OutputHandler& OH) const;
 
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const {
      return 6;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const { 
      *piNumRows = 30; 
      *piNumCols = 30; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "SphericalHingeJoint";
   };
#endif   

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */
};

/* SphericalHingeJoint - end */


/* PinJoint - begin */

/* Incastro con liberta' di rotazione sui tre assi */

class PinJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode;
   Vec3 X0;
   Vec3 d;
   Vec3 F;
   
 public:
   /* Costruttore non banale */
   PinJoint(unsigned int uL, const DofOwner* pDO,
	    const StructNode* pN, 
	    const Vec3& X0Tmp, const Vec3& dTmp, flag fOut);
   
   ~PinJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const {
      return Joint::PIN; 
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 3;
   };
   
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 3);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 9; 
      *piNumCols = 9; 
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
			    
   DofOrder::Order GetEqType(unsigned int i) const;
   
   virtual void Output(OutputHandler& OH) const;
 
   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const {
      return 6;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const  { 
      *piNumRows = 18; 
      *piNumCols = 18; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "PinJoint";
   };
#endif

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode;
   };
   /* ************************************************ */ 
};

/* PinJoint - end */

#endif
