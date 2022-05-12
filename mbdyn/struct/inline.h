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

#ifndef INLINE_H
#define INLINE_H

#include <joint.h>
#include "friction.h"

/* InLineJoint - begin */

class InLineJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   
   const Mat3x3 Rv;
   const Vec3 p;
   
   Vec3 F;

#ifdef USE_NETCDF
   size_t Var_FF;
   size_t Var_fc; 
   size_t Var_v; 
   size_t Var_displ;
#endif // USE_NETCDF

   /* friction related data */
   BasicShapeCoefficient *const Sh_c;
   BasicFriction *const fc;
   const doublereal preF;
   doublereal F3;
   doublereal v;
   doublereal displ;
   static const unsigned int NumSelfDof;
   static const unsigned int NumDof;
   /* end of friction related data */
   
 public:
   /* Costruttore */
   InLineJoint(unsigned int uL, const DofOwner* pDO,
	       const StructNode* pN1, const StructNode* pN2, 
	       const Mat3x3& RvTmp, const Vec3& pTmp, flag fOut,
	       const doublereal pref = 0.,
	       BasicShapeCoefficient *const sh = 0,
	       BasicFriction *const f = 0);
   
   ~InLineJoint(void);

   /* Tipo di joint */
   virtual Joint::Type GetJointType(void) const {
      return J_INLINE;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const;
      
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = NumDof;
      *piNumCols = NumDof; 
      if (fc) {
          *piNumRows += fc->iGetNumDof();
          *piNumCols += fc->iGetNumDof();
      } 
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

   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 4; 
   };
   
   virtual void 
   InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 24+4; 
      *piNumCols = 24+4; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;

   /* describes the dimension of components of equation */
   virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* InLineJoint - end */


/* InLineWithOffsetJoint - begin */

class InLineWithOffsetJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   
   const Mat3x3 Rv;
   const Vec3 p;
   const Vec3 q;
   
   Vec3 F;
#ifdef USE_NETCDF
   size_t Var_FF;
   size_t Var_fc;
#endif // USE_NETCDF
 public:
   /* Costruttore */
   InLineWithOffsetJoint(unsigned int uL, const DofOwner* pDO,
			 const StructNode* pN1, const StructNode* pN2, 
			 const Mat3x3& RvTmp, 
			 const Vec3& pTmp, const Vec3& qTmp, flag fOut);
   
   ~InLineWithOffsetJoint(void);

   /* Tipo di joint */
   virtual Joint::Type GetJointType(void) const {
      return J_INLINE;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 2;
   };
      
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 2);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 12+2;
      *piNumCols = 12+2; 
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
      return 4; 
   };
   
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 24+4; 
      *piNumCols = 24+4; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;

   /* describes the dimension of components of equation */
   virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* InLineWithOffsetJoint - end */

#endif /* INLINE_H */

