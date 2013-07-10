/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
/*
 * Author: Eduardo Okabe, 2013
 */

// Simple Gear Joint v2 - July 2nd, 2013

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "Rot.hh"

class SimpleGearJoint
: virtual public Elem, public UserDefinedElem {
private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const StructNode* pNodeRef;
   doublereal Gear_r1;
   doublereal Gear_r2;
   Vec3 ThetaOut1;
   Vec3 ThetaOut2;
   Vec3 M1;
   Vec3 M2;
   Vec3 kvector;
   Mat3x3 R1tilde;
   Mat3x3 R2tilde;
   
public:
	SimpleGearJoint(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~SimpleGearJoint(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
};

SimpleGearJoint::SimpleGearJoint(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO), 
pNode1(0), pNode2(0), pNodeRef(0), Gear_r1(0.), Gear_r2(0.)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Simple Gear Joint						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}
      
	// Simple Gear Joint processing .mbd file:

	// Read the node of gear 1 from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
          silent_cerr("Simple Gear Joint (" << GetLabel() << ") - node 1: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the node of gear 2 from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
          silent_cerr("Simple Gear Joint (" << GetLabel() << ") - node 2: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the reference node (support structure) from .mbd file:
   pNodeRef = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef) {
          silent_cerr("Simple Gear Joint (" << GetLabel() << ") - reference node: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the radius of gear 1 from .mbd file:
   Gear_r1 = HP.GetReal();

	// Read the radius of gear 2 from .mbd file:
   Gear_r2 = HP.GetReal();

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

   kvector = Vec3(0., 0., 1.);
  
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Mat3x3 RRf(pNodeRef->GetRCurr());     
  
   R1tilde = RRf.MulTM(R1);  
   R2tilde = RRf.MulTM(R2);  
  
   ThetaOut1 = RotManip::VecRot(R1.MulTM(RRf)*R1tilde);
   ThetaOut2 = RotManip::VecRot(R2.MulTM(RRf)*R2tilde);

}

SimpleGearJoint::~SimpleGearJoint(void)
{
	// destroy private data
	NO_OP;
}

void
SimpleGearJoint::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << M1
      << " " << M2
      << " " << -M1-M2
      << " " << ThetaOut1
      << " " << ThetaOut2
      << std::endl;
   }
}

void
SimpleGearJoint::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 10;
	*piNumCols = 10;
}

VariableSubMatrixHandler& 
SimpleGearJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering SimpleGearJoint::AssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   /* Change the dimension of the submatrix based on the constraint demand */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);


   /* Recover the index of the nodes and reaction moment variables */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iNodeRefFirstPosIndex = pNodeRef->iGetFirstPositionIndex()+3;
   integer iNodeRefFirstMomIndex = pNodeRef->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex()+1;


   /* Set the indexes of equation */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
          WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
          WM.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
          WM.PutColIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
          WM.PutRowIndex(6 + iCnt, iNodeRefFirstMomIndex + iCnt);
          WM.PutColIndex(6 + iCnt, iNodeRefFirstPosIndex + iCnt);
   }
   
   /* Set the index of constraint equation */
   WM.PutRowIndex(10, iFirstReactionIndex);
   WM.PutColIndex(10, iFirstReactionIndex);

   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Mat3x3 RRf(pNodeRef->GetRCurr());
   
   Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf)*R1tilde));
   Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf)*R2tilde));

   //unwrapping function:
   theta1 = Unwrap(ThetaOut1, theta1);
   theta2 = Unwrap(ThetaOut2, theta2);
   
   Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
   Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));
   
   Mat3x3 MTmp1 = (-Mat3x3(MatCross, R1*Gamma1IT*kvector) + R1*Gamma1IT*RotManip::Elle(-theta1, Gamma1IT*kvector))*Gear_r1*dCoef;
   Mat3x3 MTmp2 = (-Mat3x3(MatCross, R2*Gamma2IT*kvector) + R2*Gamma2IT*RotManip::Elle(-theta2, Gamma2IT*kvector))*Gear_r2*dCoef;

   Vec3 C1(R1*Gamma1IT*kvector*Gear_r1);   
   Vec3 C2(R2*Gamma2IT*kvector*Gear_r2);

   Vec3 dTmp1 = kvector*RotManip::DRot_I(theta1).MulMT(R1)*Gear_r1*dCoef;
   Vec3 dTmp2 = kvector*RotManip::DRot_I(theta2).MulMT(R2)*Gear_r2*dCoef;
   
   // Perturbation of couple in node 1:
   WM.Add(1, 1, MTmp1);
   WM.Add(1, 10, C1);

   // Perturbation of couple in node 2:
   WM.Add(3 + 1, 3 + 1, MTmp2); 
   WM.Add(3 + 1, 10, C2); 

   // Perturbation of couple in reference node:
   WM.Sub(6 + 1, 1, MTmp1); 
   WM.Sub(6 + 1, 3 + 1, MTmp2); 
   WM.Sub(6 + 1, 10, C1+C2);

   // Perturbation of Constraint Equation:
   WM.SubT(10, 1, dTmp1);
   WM.SubT(10, 3 + 1, dTmp2);
 
   WM.AddT(10, 6 + 1, dTmp1+dTmp2);
 
	return WorkMat;
}

SubVectorHandler& 
SimpleGearJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering SimpleGearJoint::AssRes()" << std::endl);

   /* Change the dimension of the vector based on the constraint demand */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   /* Get the index of the nodes and reaction moment variables */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iNodeRefFirstMomIndex = pNodeRef->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex()+1;

   for (int iCnt = 1; iCnt <= 3; iCnt++) {
          WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(3 + iCnt, iNode2FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(6 + iCnt, iNodeRefFirstMomIndex + iCnt);
   }

   WorkVec.PutRowIndex(10, iFirstReactionIndex);

   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Mat3x3 RRf(pNodeRef->GetRCurr());
   
   Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf)*R1tilde));
   Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf)*R2tilde));

   DEBUGCOUT("SimpleGearJoint::AssRes(), R1r: " << R1.MulTM(RRf)*R1tilde << std::endl);
   DEBUGCOUT("SimpleGearJoint::AssRes(), R2r: " << R2.MulTM(RRf)*R2tilde << std::endl);      
   DEBUGCOUT("SimpleGearJoint::AssRes(), Norm[theta1]: " << theta1.Norm() << std::endl);
   DEBUGCOUT("SimpleGearJoint::AssRes(), Norm[theta2]: " << theta2.Norm() << std::endl);      
   DEBUGCOUT("SimpleGearJoint::AssRes(), theta1: " << theta1 << std::endl);
   DEBUGCOUT("SimpleGearJoint::AssRes(), theta2: " << theta2 << std::endl);   
   DEBUGCOUT("SimpleGearJoint::AssRes(), ThetaOut1: " << ThetaOut1 << std::endl);
   DEBUGCOUT("SimpleGearJoint::AssRes(), thetaOut2: " << ThetaOut2 << std::endl);   
   
   //unwrapping function:
   theta1 = Unwrap(ThetaOut1, theta1);
   theta2 = Unwrap(ThetaOut2, theta2);

   Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
   Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));
   
   doublereal MTmp = XCurr(iFirstReactionIndex);
   Vec3 C1(R1*Gamma1IT*kvector*Gear_r1*MTmp);   
   Vec3 C2(R2*Gamma2IT*kvector*Gear_r2*MTmp);

   WorkVec.Sub(1, C1);
   WorkVec.Sub(3 + 1, C2);
   WorkVec.Add(6 + 1, C1+C2);
    
   DEBUGCOUT("SimpleGearJoint::AssRes(), theta1 unwrapped: " << theta1 << std::endl);
   DEBUGCOUT("SimpleGearJoint::AssRes(), theta2 unwrapped: " << theta2 << std::endl);
   
   doublereal th1z(theta1*kvector);
   doublereal th2z(theta2*kvector);

   doublereal eps(-th1z*Gear_r1 - th2z*Gear_r2);
   WorkVec.PutCoef(10, eps);
      
	return WorkVec;
}

unsigned int
SimpleGearJoint::iGetNumPrivData(void) const
{
	return 0;
}

int
SimpleGearJoint::iGetNumConnectedNodes(void) const
{
	//Useful, but not essential:
	return 3;
}

void
SimpleGearJoint::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	//Useful, but not essential:
	connectedNodes.resize(3);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
	connectedNodes[2] = pNodeRef;
}

void
SimpleGearJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
SimpleGearJoint::Restart(std::ostream& out) const
{
	return out << "# SimpleGearJoint: not implemented" << std::endl;
}

unsigned int
SimpleGearJoint::iGetInitialNumDof(void) const
{
	return 0;
}

void 
SimpleGearJoint::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
SimpleGearJoint::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
SimpleGearJoint::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


unsigned int
SimpleGearJoint::iGetNumDof(void) const
{
   return 1;
}

 
DofOrder::Order
SimpleGearJoint::GetDofType(unsigned int i) const
{
   return DofOrder::ALGEBRAIC;
}


void
SimpleGearJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Mat3x3 RRf(pNodeRef->GetRCurr());
   
   Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf)*R1tilde));
   Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf)*R2tilde));

   // Calculation of couples M1 and M2 (global reference):
   Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
   Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));
   
   doublereal MTmp = X(iGetFirstIndex()+1);

   M1 = R1*Gamma1IT*kvector*Gear_r1*MTmp;   
   M2 = R2*Gamma2IT*kvector*Gear_r2*MTmp;

   DEBUGCOUT("SGJoint::AfterConv(), theta1, ThetaOut1: " << theta1 << ", " << ThetaOut1 << std::endl);
   DEBUGCOUT("SGJoint::AfterConv(), theta2, ThetaOut2: " << theta2 << ", " << ThetaOut2 << std::endl);

   // Actual rotation angles calculation (local reference):  
   ThetaOut1 = Unwrap(ThetaOut1, theta1);
   ThetaOut2 = Unwrap(ThetaOut2, theta2);

   DEBUGCOUT("SGJoint::AfterConv(), theta1, ThetaOut1: " << theta1 << ", " << ThetaOut1 << std::endl);
   DEBUGCOUT("SGJoint::AfterConv(), theta2, ThetaOut2: " << theta2 << ", " << ThetaOut2 << std::endl);

}


extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<SimpleGearJoint>;

	if (!SetUDE("SimpleGearJoint", rf)) {
		delete rf;

		silent_cerr("module-fabricate: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

