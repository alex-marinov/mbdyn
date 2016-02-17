/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/**
    Library of motion transmission components for "digital fabrication" machines (alpha version) [2013]
    Eduardo Okabe (okabe@unicamp.br)
    Postdoc CNPq at Aero/Polimi
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "Rot.hh"

#include "drive.h"

// Gear Joint v2f - July 25th, 2013

class GearJoint
: virtual public Elem, public UserDefinedElem {
private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const StructNode* pNodeRef;
   doublereal Gear_r1;
   doublereal Gear_r2;
   doublereal MTmp;
   doublereal jr_coef;
   Vec3 ThetaOut1;
   Vec3 ThetaOut2;
   Vec3 M1;
   Vec3 M2;
   Mat3x3 R1rtilde;
   Mat3x3 R2rtilde;
   Mat3x3 R1c;
   Mat3x3 R2c;

public:
	GearJoint(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~GearJoint(void);

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

GearJoint::GearJoint(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
pNode1(0), pNode2(0), pNodeRef(0), Gear_r1(1.), Gear_r2(1.)
{
   DEBUGCOUT("Entering GearJoint constructor" << std::endl);

	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Gear Joint						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// Gear Joint processing .mbd file:

	// Read the node of gear 1 from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
          silent_cerr("Gear Joint (" << GetLabel() << ") - node 1: structural node expected at line " << HP.GetLineData() << std::endl);
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   ReferenceFrame RF1(pNode1);

   R1rtilde = Mat3x3(Eye3);

	// Read the relative reference frame of gear 1, if supplied:
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 1 is supplied" << std::endl);
        R1rtilde = HP.GetRotRel(RF1);
   }

	// Read the node of gear 2 from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
      silent_cerr("Gear Joint (" << GetLabel() << ") - node 2: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   ReferenceFrame RF2(pNode2);

   R2rtilde = Mat3x3(Eye3);

	// Read the relative reference frame of gear 2, if supplied:
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 1 is supplied" << std::endl);
        R2rtilde = HP.GetRotRel(RF2);
   }

	// Read the reference node (support structure) from .mbd file:
   pNodeRef = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef) {
      silent_cerr("Gear Joint (" << GetLabel() << ") - reference node: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }


   if (HP.IsKeyWord("ratio")) {
	   // Set radius of gear 1 equal to 1.:
      Gear_r1 = 1.;
	   // Read the ratio (r2/r1) from .mbd file:
      Gear_r2 = HP.GetReal();
   }
   else {
	   // Read the radius of gear 1 from .mbd file:
      Gear_r1 = HP.GetReal();
	   // Read the radius of gear 2 from .mbd file:
      Gear_r2 = HP.GetReal();
   }

	// Verify if both coefficient are zero (null joint):
   if ((Gear_r1 == 0.) && (Gear_r2 == 0.)) {
      silent_cerr("Gear Joint: both coefficients (r1 and r2) cannot be equal to zero at the same time. Error at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Joint regularization coefficient:

   if (HP.IsKeyWord("regularization")) {
      jr_coef = HP.GetReal();
   } else {
      jr_coef = 0.;
   };

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

   Mat3x3 R1(pNode1->GetRCurr()*R1rtilde);
   Mat3x3 R2(pNode2->GetRCurr()*R2rtilde);
   Mat3x3 RRf(pNodeRef->GetRCurr());


	// Disable initial orientation correction (correction is enabled by default):
   if (HP.IsKeyWord("disable" "initial" "correction")) {
      R1c = Mat3x3(Eye3);
      R2c = Mat3x3(Eye3);
   }
   else {
      R1c = RRf.MulTM(R1);
      R2c = RRf.MulTM(R2);
   }

	// If correction is enabled, both angles should be zero:
   ThetaOut1 = RotManip::VecRot(R1.MulTM(RRf)*R1c);
   ThetaOut2 = RotManip::VecRot(R2.MulTM(RRf)*R2c);

}

GearJoint::~GearJoint(void)
{
	// destroy private data
	NO_OP;
}

void
GearJoint::Output(OutputHandler& OH) const
{

   if (bToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << M1         // Moment on node 1
      << " " << M2         // Moment on node 2
      << " " << -M1-M2     // Moment on reference node
      << " " << ThetaOut1  // Actual relative rotation of node 1
      << " " << ThetaOut2  // Actual relative rotation of node 2
      << " " << MTmp  // Lambda
      << std::endl;
   }
}

void
GearJoint::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 10;
	*piNumCols = 10;
}

VariableSubMatrixHandler&
GearJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering GearJoint::AssJac()" << std::endl);

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

    // Retrieve node rotation matrices:
    Mat3x3 R1r(pNode1->GetRCurr()*R1rtilde);
    Mat3x3 R2r(pNode2->GetRCurr()*R2rtilde);
    Mat3x3 RRf(pNodeRef->GetRCurr());

    // Calculate relative angles:
    Vec3 theta1(RotManip::VecRot(R1r.MulTM(RRf)*R1c));
    Vec3 theta2(RotManip::VecRot(R2r.MulTM(RRf)*R2c));

    Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
    Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));

    // MTmp = XCurr(iFirstReactionIndex); // Its value comes from residue calculation
    Mat3x3 MTmp1 = (-Mat3x3(MatCross, R1r*Gamma1IT.GetCol(3)) + R1r*Gamma1IT*RotManip::Elle(-theta1, Gamma1IT.GetCol(3)))*MTmp*Gear_r1*dCoef;
    Mat3x3 MTmp2 = (-Mat3x3(MatCross, R2r*Gamma2IT.GetCol(3)) + R2r*Gamma2IT*RotManip::Elle(-theta2, Gamma2IT.GetCol(3)))*MTmp*Gear_r2*dCoef;

    Vec3 C1(R1r*Gamma1IT.GetCol(3)*Gear_r1);
    Vec3 C2(R2r*Gamma2IT.GetCol(3)*Gear_r2);

    Vec3 dTmp1 = (RotManip::DRot_I(theta1).MulMT(R1r)).GetRow(3)*Gear_r1*dCoef;
    Vec3 dTmp2 = (RotManip::DRot_I(theta2).MulMT(R2r)).GetRow(3)*Gear_r2*dCoef;

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

    // Perturbation of joint regularization term:
    WM.IncCoef(10, 10, jr_coef);

	return WorkMat;
}

SubVectorHandler&
GearJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering GearJoint::AssRes()" << std::endl);

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

   Mat3x3 R1r(pNode1->GetRCurr()*R1rtilde);
   Mat3x3 R2r(pNode2->GetRCurr()*R2rtilde);
   Mat3x3 RRf(pNodeRef->GetRCurr());

   Vec3 theta1(RotManip::VecRot(R1r.MulTM(RRf)*R1c));
   Vec3 theta2(RotManip::VecRot(R2r.MulTM(RRf)*R2c));

   DEBUGCOUT("GearJoint::AssRes(), theta1: " << theta1 << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), theta2: " << theta2 << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), ThetaOut1: " << ThetaOut1 << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), thetaOut2: " << ThetaOut2 << std::endl);

   Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
   Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));

   MTmp = XCurr(iFirstReactionIndex);
   Vec3 C1(R1r*Gamma1IT.GetCol(3)*Gear_r1*MTmp);
   Vec3 C2(R2r*Gamma2IT.GetCol(3)*Gear_r2*MTmp);

   WorkVec.Sub(1, C1);
   WorkVec.Sub(3 + 1, C2);
   WorkVec.Add(6 + 1, C1+C2);

  // Unwrap angle:
   theta1 = Unwrap(ThetaOut1, theta1);
   theta2 = Unwrap(ThetaOut2, theta2);

   DEBUGCOUT("GearJoint::AssRes(), Lambda: " << MTmp << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), theta1 unwrapped: " << theta1 << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), theta2 unwrapped: " << theta2 << std::endl);
   DEBUGCOUT("GearJoint::AssRes(), C1, C2: " << C1 << ", " << C2 << std::endl);

   doublereal eps(-theta1(3)*Gear_r1 - theta2(3)*Gear_r2 - jr_coef*MTmp);
   WorkVec.PutCoef(10, eps);

	return WorkVec;
}

unsigned int
GearJoint::iGetNumPrivData(void) const
{
	return 0;
}

int
GearJoint::iGetNumConnectedNodes(void) const
{
	//Useful, but not essential:
	return 3;
}

void
GearJoint::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	//Useful, but not essential:
	connectedNodes.resize(3);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
	connectedNodes[2] = pNodeRef;
}

void
GearJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
GearJoint::Restart(std::ostream& out) const
{
	return out << "# GearJoint: not implemented" << std::endl;
}

unsigned int
GearJoint::iGetInitialNumDof(void) const
{
	return 0;
}

void
GearJoint::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
GearJoint::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
GearJoint::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


unsigned int
GearJoint::iGetNumDof(void) const
{
   return 1;
}


DofOrder::Order
GearJoint::GetDofType(unsigned int i) const
{
   return DofOrder::ALGEBRAIC;
}


void
GearJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   Mat3x3 R1r(pNode1->GetRCurr()*R1rtilde);
   Mat3x3 R2r(pNode2->GetRCurr()*R2rtilde);
   Mat3x3 RRf(pNodeRef->GetRCurr());

   Vec3 theta1(RotManip::VecRot(R1r.MulTM(RRf)*R1c));
   Vec3 theta2(RotManip::VecRot(R2r.MulTM(RRf)*R2c));

   // Calculation of couples M1 and M2 (global reference):
   Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));
   Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));

   // Lagrange multiplier:
   doublereal MTmp = X(iGetFirstIndex()+1);

   M1 = R1r*Gamma1IT.GetCol(3)*Gear_r1*MTmp;
   M2 = R2r*Gamma2IT.GetCol(3)*Gear_r2*MTmp;

   // Actual rotation angles calculation (local reference):
   ThetaOut1 = Unwrap(ThetaOut1, theta1);
   ThetaOut2 = Unwrap(ThetaOut2, theta2);

}


// Linear to Linear Transmission Joint v0a - July 16th, 2013

class LinearTransmissionJoint
: virtual public Elem, public UserDefinedElem {
private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const StructNode* pNodeRef1;
   const StructNode* pNodeRef2;
   Vec3 x1_off;
   Vec3 x1ref_off;
   Vec3 x2_off;
   Vec3 x2ref_off;
   doublereal z_corr;
   doublereal Coef_f1;
   doublereal Coef_f2;
   doublereal LambdaTmp;
   Vec3 F1;
   Vec3 M1;
   Vec3 M1ref;
   Vec3 F2;
   Vec3 M2;
   Vec3 M2ref;
   Mat3x3 R1tilde;
   Mat3x3 R2tilde;

public:
	LinearTransmissionJoint(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~LinearTransmissionJoint(void);

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

LinearTransmissionJoint::LinearTransmissionJoint(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
pNode1(0), pNodeRef1(0), pNode2(0), pNodeRef2(0), Coef_f1(1.), Coef_f2(1.)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Linear Transmission Joint						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// Linear Transmission Joint processing .mbd file:

	// Read node 1 from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
      silent_cerr("Linear Transmission Joint (" << GetLabel() << ") - node 1: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of node 1, if supplied:
   ReferenceFrame RF1(pNode1);
   x1_off = Vec3(Zero3);
   if (HP.IsKeyWord("position")) {
      x1_off = HP.GetPosRel(RF1);
   }

	// Read the relative reference of node 1, if supplied:
   R1tilde = Mat3x3(Eye3);
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 1 is supplied" << std::endl);
        R1tilde = HP.GetRotRel(RF1);
   }

	// Read reference node 1 from .mbd file:
   pNodeRef1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef1) {
      silent_cerr("Linear Transmission Joint (" << GetLabel() << ") - reference node 1: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of reference node 1, if supplied:
   ReferenceFrame RF1ref(pNode1);
   x1ref_off = Vec3(Zero3);
   if (HP.IsKeyWord("position")) {
      x1ref_off = HP.GetPosRel(RF1ref);
   }


	// Read node 2 from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
      silent_cerr("Linear Transmission Joint (" << GetLabel() << ") - node 2: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of node 2, if supplied:
   ReferenceFrame RF2(pNode2);
   x2_off = Vec3(Zero3);
   if (HP.IsKeyWord("position")) {
      DEBUGCOUT("Position offset of node 2 is supplied" << std::endl);
      x2_off = HP.GetPosRel(RF2);
   }

	// Read the relative reference of node 2, if supplied:
   R2tilde = Mat3x3(Eye3);
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 2 is supplied" << std::endl);
        R2tilde = HP.GetRotRel(RF2);
   }

	// Read reference node 2 from .mbd file:
   pNodeRef2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef2) {
      silent_cerr("Linear Transmission Joint (" << GetLabel() << ") - reference node 2: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of reference node 2, if supplied:
   ReferenceFrame RF2ref(pNode2);
   x2ref_off = Vec3(Zero3);
   if (HP.IsKeyWord("position")) {
      x2ref_off = HP.GetPosRel(RF2ref);
   }

   if (HP.IsKeyWord("ratio")) {
	   // Set coefficient of node 1 equal to 1.:
      Coef_f1 = 1.;
	   // Read the ratio (f2/f1) from .mbd file:
      Coef_f2 = HP.GetReal();
   }
   else {
	   // Read the coefficient of node 1 from .mbd file:
      Coef_f1 = HP.GetReal();
	   // Read the coefficient of node 2 from .mbd file:
      Coef_f2 = HP.GetReal();
   }

	// Verify if both coefficient are zero (null joint):
   if ((Coef_f1 == 0.) && (Coef_f2 == 0.)) {
      silent_cerr("Linear Transmission Joint: both coefficients (f1 and f2) cannot be equal to zero at the same time. Error at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	// Disable initial orientation correction (correction is enabled by default):
   if (HP.IsKeyWord("disable" "initial" "correction")) {
      z_corr = 0.;
   }
   else {
      Mat3x3 R1h(pNode1->GetRCurr()*R1tilde);
      Mat3x3 R2h(pNode2->GetRCurr()*R2tilde);

      Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - pNodeRef1->GetRCurr()*x1ref_off;
      Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - pNodeRef2->GetRCurr()*x2ref_off;

      z_corr = (R1h.MulTV(b1)+R1tilde.MulTV(x1_off))(3)*Coef_f1 -
               (R2h.MulTV(b2)+R2tilde.MulTV(x2_off))(3)*Coef_f2;
   }

}

LinearTransmissionJoint::~LinearTransmissionJoint(void)
{
	// destroy private data
	NO_OP;
}

void
LinearTransmissionJoint::Output(OutputHandler& OH) const
{

   if (bToBeOutput()) {
   std::ostream& out = OH.Loadable();
   out << std::setw(8) << GetLabel()
      << " " << F1         // Force on node 1
      << " " << M1         // Moment on node 1
      << " " << -F1        // Force on reference node 1
      << " " << M1ref      // Moment on reference node 1
      << " " << F2         // Force on node 2
      << " " << M2         // Moment on node 2
      << " " << -F2        // Force on reference node 2
      << " " << M2ref      // Moment on reference node 2
      << std::endl;
   }
}

void
LinearTransmissionJoint::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	//Space Dimension = F1(3)+M1(3)+F1ref(3)+M1ref(3) + F2(3)+M2(3)+F2ref(3)+M2ref(3) + Constraint(1) = 12 + 12 + 1 = 25
	*piNumRows = 25;
	*piNumCols = 25;
}

VariableSubMatrixHandler&
LinearTransmissionJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering LinearTransmissionJoint::AssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   /* Change the dimension of the submatrix based on the constraint demand */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);


   /* Recover the index of the nodes and reaction moment variables */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNodeRef1FirstPosIndex = pNodeRef1->iGetFirstPositionIndex();
   integer iNodeRef1FirstMomIndex = pNodeRef1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iNodeRef2FirstPosIndex = pNodeRef2->iGetFirstPositionIndex();
   integer iNodeRef2FirstMomIndex = pNodeRef2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex()+1;


   /* Set the indexes of equation */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
          WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
          WM.PutRowIndex(6 + iCnt, iNodeRef1FirstMomIndex + iCnt);
          WM.PutColIndex(6 + iCnt, iNodeRef1FirstPosIndex + iCnt);
          WM.PutRowIndex(12 + iCnt, iNode2FirstMomIndex + iCnt);
          WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
          WM.PutRowIndex(18 + iCnt, iNodeRef2FirstMomIndex + iCnt);
          WM.PutColIndex(18 + iCnt, iNodeRef2FirstPosIndex + iCnt);
   }

   /* Set the index of constraint equation */
   WM.PutRowIndex(25, iFirstReactionIndex);
   WM.PutColIndex(25, iFirstReactionIndex);

   // Retrieve node rotation matrices:
   Mat3x3 R1h(pNode1->GetRCurr()*R1tilde);
   Mat3x3 R2h(pNode2->GetRCurr()*R2tilde);

   // Calculate vectors b1 and b2:
   Vec3 RXref1(pNodeRef1->GetRCurr()*x1ref_off);
   Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RXref1;

   Vec3 RXref2(pNodeRef2->GetRCurr()*x2ref_off);
   Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RXref2;

   // Partial calculation of forces and moments:
   Vec3 f1kR(R1h.GetCol(3)*Coef_f1);
   Vec3 f2kR(R2h.GetCol(3)*Coef_f2);

   // LambdaTmp comes from the AssRes():

   // Calculate constraint force/lambda on node 1:
   F1 = -R1h.GetCol(3)*Coef_f1;

   // Calculate constraint moment/lambda on node 1:
   M1 = (b1.Cross(R1h)).GetCol(3)*Coef_f1;

   // Calculate constraint moment/lambda on reference node 1:
   M1ref = (RXref1.Cross(R1h)).GetCol(3)*Coef_f1;


   // Calculate constraint force/lambda on node 2:
   F2 = R2h.GetCol(3)*Coef_f2;

   // Calculate constraint moment/lambda on node 2:
   M2 = -(b2.Cross(R2h)).GetCol(3)*Coef_f2;

   // Calculate constraint moment/lambda on reference node 2:
   M2ref = -(RXref2.Cross(R2h)).GetCol(3)*Coef_f2;


   // Partial calculations of force/moments perturbation:
   Mat3x3 R1hk_Cross(Mat3x3(MatCross, R1h.GetCol(3)*dCoef*Coef_f1*LambdaTmp));
   Mat3x3 R2hk_Cross(Mat3x3(MatCross, R2h.GetCol(3)*dCoef*Coef_f2*LambdaTmp));
   Mat3x3 RfXo1_Cross(Mat3x3(MatCross, pNodeRef1->GetRCurr()*x1ref_off));
   Mat3x3 RfXo2_Cross(Mat3x3(MatCross, pNodeRef2->GetRCurr()*x2ref_off));


   // Perturbation of force on node 1:
   WM.Add(1, 3 + 1, R1hk_Cross);
   WM.Add(1, 25, F1);

   // Perturbation of moment on node 1:
   WM.Sub(3 + 1, 1, R1hk_Cross);
   WM.Sub(3 + 1, 3 + 1, b1.Cross(R1hk_Cross));
   WM.Add(3 + 1, 6 + 1, R1hk_Cross);
   WM.Sub(3 + 1, 9 + 1, R1hk_Cross*RfXo1_Cross);
   WM.Add(3 + 1, 25, M1);

   // Perturbation of force on reference node 1:
   WM.Sub(6 + 1, 3 + 1, R1hk_Cross);
   WM.Sub(6 + 1, 25, F1);

   // Perturbation of moment on reference node 1:
   WM.Sub(9 + 1, 3 + 1, RfXo1_Cross*R1hk_Cross);
   WM.Add(9 + 1, 9 + 1, R1hk_Cross*RfXo1_Cross);
   WM.Add(9 + 1, 25, M1ref);


   // Perturbation of force on node 2:
   WM.Sub(12 + 1, 15 + 1, R2hk_Cross);
   WM.Add(12 + 1, 25, F2);

   // Perturbation of moment on node 2:
   WM.Add(15 + 1, 12 + 1, R2hk_Cross);
   WM.Add(15 + 1, 15 + 1, b2.Cross(R2hk_Cross));
   WM.Sub(15 + 1, 18 + 1, R2hk_Cross);
   WM.Add(15 + 1, 21 + 1, R2hk_Cross*RfXo2_Cross);
   WM.Add(15 + 1, 25, M2);

   // Perturbation of force on reference node 2:
   WM.Add(18 + 1, 15 + 1, R2hk_Cross);
   WM.Sub(18 + 1, 25, F2);

   // Perturbation of moment on reference node 2:
   WM.Add(21 + 1, 15 + 1, RfXo2_Cross*R2hk_Cross);
   WM.Sub(21 + 1, 21 + 1, R2hk_Cross*RfXo2_Cross);
   WM.Add(21 + 1, 25, M2ref);


   // Perturbation of Constraint Equation:
   WM.AddT(25, 1, f1kR*dCoef);
   WM.AddT(25, 3 + 1, f1kR.Cross(b1)*dCoef);
   WM.SubT(25, 6 + 1, f1kR*dCoef);
   WM.AddT(25, 9 + 1, f1kR.Cross(RXref1)*dCoef);

   WM.SubT(25, 12 + 1, f2kR*dCoef);
   WM.SubT(25, 15 + 1, f2kR.Cross(b2)*dCoef);
   WM.AddT(25, 18 + 1, f2kR*dCoef);
   WM.SubT(25, 21 + 1, f2kR.Cross(RXref2)*dCoef);

	return WorkMat;
}

SubVectorHandler&
LinearTransmissionJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering LinearTransmissionJoint::AssRes()" << std::endl);

   /* Change the dimension of the vector based on the constraint demand */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   /* Get the index of the nodes and reaction moment variables */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNodeRef1FirstMomIndex = pNodeRef1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iNodeRef2FirstMomIndex = pNodeRef2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex()+1;

   /* Set the indexes of vector */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
          WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(6 + iCnt, iNodeRef1FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(12 + iCnt, iNode2FirstMomIndex + iCnt);
          WorkVec.PutRowIndex(18 + iCnt, iNodeRef2FirstMomIndex + iCnt);
   }

   WorkVec.PutRowIndex(25, iFirstReactionIndex);

   // Retrieve Lagrange multiplier:
   LambdaTmp = XCurr(iFirstReactionIndex);


   // Retrieve node rotation matrices:
   Mat3x3 R1h(pNode1->GetRCurr()*R1tilde);
   Mat3x3 R2h(pNode2->GetRCurr()*R2tilde);

   // Calculate vectors b1 and b2:
   Vec3 RXref1(pNodeRef1->GetRCurr()*x1ref_off);
   Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RXref1;

   Vec3 RXref2(pNodeRef2->GetRCurr()*x2ref_off);
   Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RXref2;


   // Calculate constraint force/lambda on node 1:
   F1 = -R1h.GetCol(3)*Coef_f1;

   // Calculate constraint moment/lambda on node 1:
   M1 = (b1.Cross(R1h)).GetCol(3)*Coef_f1;

   // Calculate constraint moment/lambda on reference node 1:
   M1ref = (RXref1.Cross(R1h)).GetCol(3)*Coef_f1;


   // Calculate constraint force/lambda on node 2:
   F2 = R2h.GetCol(3)*Coef_f2;

   // Calculate constraint moment/lambda on node 2:
   M2 = -(b2.Cross(R2h)).GetCol(3)*Coef_f2;

   // Calculate constraint moment/lambda on reference node 2:
   M2ref = -(RXref2.Cross(R2h)).GetCol(3)*Coef_f2;


   // Calculate constraint equation:
   doublereal eps((R1h.MulTV(b1)+R1tilde.MulTV(x1_off))(3)*Coef_f1 -
                  (R2h.MulTV(b2)+R2tilde.MulTV(x2_off))(3)*Coef_f2 - z_corr);

   // DEBUGCOUT("LinearTransmissionJoint::AssRes(), Coef_f1: " << Coef_f1 << std::endl);

   // Put (or add) calculated values in the WorkVec:

   // node 1:
   WorkVec.Sub(1, F1*LambdaTmp);
   WorkVec.Sub(3 + 1, M1*LambdaTmp);

   // reference node 1:
   WorkVec.Add(6 + 1, F1*LambdaTmp);
   WorkVec.Sub(9 + 1, M1ref*LambdaTmp);

   // node 2:
   WorkVec.Sub(12 + 1, F2*LambdaTmp);
   WorkVec.Sub(15 + 1, M2*LambdaTmp);

   // reference node 2:
   WorkVec.Add(18 + 1, F2*LambdaTmp);
   WorkVec.Sub(21 + 1, M2ref*LambdaTmp);

   // constraint:
   WorkVec.PutCoef(25, -eps);

	return WorkVec;
}

unsigned int
LinearTransmissionJoint::iGetNumPrivData(void) const
{
	return 0;
}

int
LinearTransmissionJoint::iGetNumConnectedNodes(void) const
{
	//Useful, but not essential:
	return 4;
}

void
LinearTransmissionJoint::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	//Useful, but not essential:
	connectedNodes.resize(4);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNodeRef1;
	connectedNodes[2] = pNode2;
	connectedNodes[3] = pNodeRef2;
}

void
LinearTransmissionJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
LinearTransmissionJoint::Restart(std::ostream& out) const
{
	return out << "# LinearTransmissionJoint: not implemented" << std::endl;
}

unsigned int
LinearTransmissionJoint::iGetInitialNumDof(void) const
{
	return 0;
}

void
LinearTransmissionJoint::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
LinearTransmissionJoint::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
LinearTransmissionJoint::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


unsigned int
LinearTransmissionJoint::iGetNumDof(void) const
{
   return 1;
}


DofOrder::Order
LinearTransmissionJoint::GetDofType(unsigned int i) const
{
   return DofOrder::ALGEBRAIC;
}


void
LinearTransmissionJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
   /*
   // Lagrange multiplier:
   doublereal LambdaTmp = X(iGetFirstIndex()+1);

   // Retrieve node rotation matrices:
   Mat3x3 R1h(pNode1->GetRCurr()*R1tilde);
   Mat3x3 R2h(pNode2->GetRCurr()*R2tilde);

   // Calculate vectors b1 and b2:
   Vec3 RXref1(pNodeRef1->GetRCurr()*x1ref_off);
   Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RXref1;

   Vec3 RXref2(pNodeRef2->GetRCurr()*x2ref_off);
   Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RXref2;

   // Calculate constraint force on node 1:
   F1 = -R1h.GetCol(3)*Coef_f1*LambdaTmp;

   // Calculate constraint moment on node 1:
   M1 = (b1.Cross(R1h)).GetCol(3)*Coef_f1*LambdaTmp;

   // Calculate constraint moment on reference node 1:
   M1ref = (RXref1.Cross(R1h)).GetCol(3)*Coef_f1*LambdaTmp;


   // Calculate constraint force on node 2:
   F2 = R2h.GetCol(3)*Coef_f2*LambdaTmp;

   // Calculate constraint moment on node 2:
   M2 = -(b2.Cross(R2h)).GetCol(3)*Coef_f2*LambdaTmp;

   // Calculate constraint moment on reference node 2:
   M2ref = -(RXref2.Cross(R2h)).GetCol(3)*Coef_f2*LambdaTmp; */

}



// Motion Transmission Joint v0a - July 17th, 2013

class MotionTransmissionJoint
: virtual public Elem, public UserDefinedElem {
private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   const StructNode* pNodeRef1;
   const StructNode* pNodeRef2;
   unsigned int mtype1;
   unsigned int mtype2;
   integer nWDim;
   Vec3 x1_off;
   Vec3 x1ref_off;
   Vec3 x2_off;
   Vec3 x2ref_off;
   Vec3 ThetaOut1;
   Vec3 ThetaOut2;
   doublereal z_corr;
   doublereal Coef_f1;
   doublereal Coef_f2;
   doublereal LambdaTmp;
   Vec3 F1;
   Vec3 M1;
   Vec3 M1ref;
   Vec3 F2;
   Vec3 M2;
   Vec3 M2ref;
   Mat3x3 R1tilde;
   Mat3x3 R2tilde;
   Mat3x3 R1c;
   Mat3x3 R2c;


public:
	MotionTransmissionJoint(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~MotionTransmissionJoint(void);

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

MotionTransmissionJoint::MotionTransmissionJoint(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
pNode1(0), pNodeRef1(0), pNode2(0), pNodeRef2(0), Coef_f1(1.), Coef_f2(1.)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
      "									\n"
      "Module: 	Motion Transmission Joint						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

   // Dimension of workspace:
   nWDim = 1;

	// Motion Transmission Joint processing .mbd file:

	// Read type of movement (linear/angular) of node 1:
   if (HP.IsKeyWord("linear")) {
      mtype1 = 1;
      nWDim += 12;
   } else if (HP.IsKeyWord("angular")) {
      mtype1 = 2;
      nWDim += 6;
   } else {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << "): "
         "invalid type for node 1 (should be linear or angular)"
         << " at line " << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read node 1 from .mbd file:
   pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode1) {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << ") - node 1: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of node 1, if supplied:
   ReferenceFrame RF1(pNode1);
   if (mtype1==1) {
      x1_off = Vec3(Zero3);
      if (HP.IsKeyWord("position")) x1_off = HP.GetPosRel(RF1);
   }

	// Read the relative reference of node 1, if supplied:
   R1tilde = Mat3x3(Eye3);
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 1 is supplied" << std::endl);
        R1tilde = HP.GetRotRel(RF1);
   }

	// Read reference node 1 from .mbd file:
   pNodeRef1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef1) {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << ") - reference node 1: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of reference node 1, if supplied:
   ReferenceFrame RF1ref(pNode1);
   if (mtype1==1) {
      x1ref_off = Vec3(Zero3);
      if (HP.IsKeyWord("position")) x1ref_off = HP.GetPosRel(RF1ref);
   }


	// Read type of movement (linear/angular) of node 2:
   if (HP.IsKeyWord("linear")) {
      mtype2 = 1;
      nWDim += 12;
   } else if (HP.IsKeyWord("angular")) {
      mtype2 = 2;
      nWDim += 6;
   } else {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << "): "
         "invalid type for node 2 (should be linear or angular)"
         << "at line " << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read node 2 from .mbd file:
   pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNode2) {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << ") - node 2: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of node 2, if supplied:
   ReferenceFrame RF2(pNode2);
   if (mtype2==1) {
      x2_off = Vec3(Zero3);
      if (HP.IsKeyWord("position")) {
         DEBUGCOUT("Position offset of node 2 is supplied" << std::endl);
         x2_off = HP.GetPosRel(RF2);
      }
   }

	// Read the relative reference of node 2, if supplied:
   R2tilde = Mat3x3(Eye3);
   if (HP.IsKeyWord("orientation")) {
        DEBUGCOUT("Rotation orientation matrix of node 2 is supplied" << std::endl);
        R2tilde = HP.GetRotRel(RF2);
   }

	// Read reference node 2 from .mbd file:
   pNodeRef2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

   if (!pNodeRef2) {
      silent_cerr("Motion Transmission Joint (" << GetLabel() << ") - reference node 2: structural node expected at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Read the relative offset of reference node 2, if supplied:
   ReferenceFrame RF2ref(pNode2);
   if (mtype2==1) {
      x2ref_off = Vec3(Zero3);
      if (HP.IsKeyWord("position")) x2ref_off = HP.GetPosRel(RF2ref);
   }

   if (HP.IsKeyWord("ratio")) {
	   // Set coefficient of node 1 equal to 1.:
      Coef_f1 = 1.;
	   // Read the ratio (f2/f1) from .mbd file:
      Coef_f2 = HP.GetReal();
   }
   else {
	   // Read the coefficient of node 1 from .mbd file:
      Coef_f1 = HP.GetReal();
	   // Read the coefficient of node 2 from .mbd file:
      Coef_f2 = HP.GetReal();
   }

	// Verify if both coefficient are zero (null joint):
   if ((Coef_f1 == 0.) && (Coef_f2 == 0.)) {
      silent_cerr("Motion Transmission Joint: coefficients (f1 and f2) cannot be equal to zero at the same time. Error at line " << HP.GetLineData() << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

	// Activate element output:
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

   Mat3x3 R1(pNode1->GetRCurr()*R1tilde);
   Mat3x3 R2(pNode2->GetRCurr()*R2tilde);
   Mat3x3 RRf1(pNodeRef1->GetRCurr());
   Mat3x3 RRf2(pNodeRef2->GetRCurr());

	// Disable initial orientation correction (correction is enabled by default):
   z_corr = 0.;

   if (HP.IsKeyWord("disable" "initial" "correction")) {
      R1c = Mat3x3(Eye3);
      R2c = Mat3x3(Eye3);
   }
   else {
      if (mtype1==1) {
         Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RRf1*x1ref_off;
         z_corr += (R1.MulTV(b1)+R1tilde.MulTV(x1_off))(3)*Coef_f1;
      } else R1c = RRf1.MulTM(R1);

      if (mtype2==1) {
         Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RRf2*x2ref_off;
         z_corr -= (R2.MulTV(b2)+R2tilde.MulTV(x2_off))(3)*Coef_f2;
      } else R2c = RRf2.MulTM(R2);
   }

	// Calculate initial unwrapped angles (zero, if correction is enabled):
   if (mtype1==2) ThetaOut1 = RotManip::VecRot(R1.MulTM(RRf1)*R1c); else ThetaOut1 = Vec3(Zero3);
   if (mtype2==2) ThetaOut2 = RotManip::VecRot(R2.MulTM(RRf2)*R2c); else ThetaOut2 = Vec3(Zero3);

}

MotionTransmissionJoint::~MotionTransmissionJoint(void)
{
	// destroy private data
	NO_OP;
}

void
MotionTransmissionJoint::Output(OutputHandler& OH) const
{

   if (bToBeOutput()) {
      std::ostream& out = OH.Loadable();
      out << std::setw(8) << GetLabel();
      if (mtype1==1) {
         out << " " << F1         // Force on node 1
             << " " << M1         // Moment on node 1
             << " " << -F1        // Force on reference node 1
             << " " << M1ref;     // Moment on reference node 1
      } else {
         out << " " << M1         // Moment on node 1
             << " " << ThetaOut1; // Actual relative rotation of node 1
      }

      if (mtype2==1) {
         out << " " << F2         // Force on node 2
             << " " << M2         // Moment on node 2
             << " " << -F2        // Force on reference node 2
             << " " << M2ref;     // Moment on reference node 2
      } else {
         out << " " << M2         // Moment on node 2
             << " " << ThetaOut2; // Actual relative rotation of node 2
      }
      out << std::endl;
   }
}

void
MotionTransmissionJoint::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	//Space Dimension = F1(3)+M1(3)+F1ref(3)+M1ref(3) + F2(3)+M2(3)+F2ref(3)+M2ref(3) + Constraint(1) = 12 + 12 + 1 = 25
	*piNumRows = *piNumCols = nWDim;
}


VariableSubMatrixHandler&
MotionTransmissionJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    DEBUGCOUT("Entering MotionTransmissionJoint::AssJac()" << std::endl);

    FullSubMatrixHandler& WM = WorkMat.SetFull();
    /* Change the dimension of the submatrix based on the constraint demand */
    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WM.ResizeReset(iNumRows, iNumCols);


    /* Recover the index of the nodes and reaction moment variables */
    integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
    integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
    integer iNodeRef1FirstPosIndex = pNodeRef1->iGetFirstPositionIndex();
    integer iNodeRef1FirstMomIndex = pNodeRef1->iGetFirstMomentumIndex();

    integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
    integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
    integer iNodeRef2FirstPosIndex = pNodeRef2->iGetFirstPositionIndex();
    integer iNodeRef2FirstMomIndex = pNodeRef2->iGetFirstMomentumIndex();

    integer iFirstReactionIndex = iGetFirstIndex()+1;


    /* Set the indexes of equation relative to node 1 */
    int nCnt1 = 3;
    int aCnt = 0;

    if (mtype1==1) nCnt1 += 3; else aCnt += 3;
    for (int iCnt = 1; iCnt <= nCnt1; iCnt++) {
          WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt + aCnt);
          WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt + aCnt);
          WM.PutRowIndex(nCnt1 + iCnt, iNodeRef1FirstMomIndex + iCnt + aCnt);
          WM.PutColIndex(nCnt1 + iCnt, iNodeRef1FirstPosIndex + iCnt + aCnt);
    }


    /* Set the indexes of equation relative to node 2 */
    nCnt1 += nCnt1; // 2 nodes = (2*nCnt1)
    aCnt = 0;
    int nCnt2 = 3;
    if (mtype2==1) nCnt2 += 3; else aCnt += 3;
    for (int iCnt = 1; iCnt <= nCnt2; iCnt++) {
          WM.PutRowIndex(nCnt1 + iCnt, iNode2FirstMomIndex + iCnt + aCnt);
          WM.PutColIndex(nCnt1 + iCnt, iNode2FirstPosIndex + iCnt + aCnt);
          WM.PutRowIndex(nCnt1 + nCnt2 + iCnt, iNodeRef2FirstMomIndex + iCnt + aCnt);
          WM.PutColIndex(nCnt1 + nCnt2 + iCnt, iNodeRef2FirstPosIndex + iCnt + aCnt);
    }

    /* Set the index of constraint equation */
    int nCnt3 = nCnt1 + 2*nCnt2 + 1;
    WM.PutRowIndex(nCnt3, iFirstReactionIndex);
    WM.PutColIndex(nCnt3, iFirstReactionIndex);

    // Retrieve node rotation matrices:
    Mat3x3 R1(pNode1->GetRCurr()*R1tilde);
    Mat3x3 R2(pNode2->GetRCurr()*R2tilde);
    Mat3x3 RRf1(pNodeRef1->GetRCurr());
    Mat3x3 RRf2(pNodeRef2->GetRCurr());

    // LambdaTmp comes from the AssRes():

    if (mtype1==1) {

      // Linear Node 1:

      // Calculate vector b1 and partial calculation of forces and moments:
      Vec3 RXref1(RRf1*x1ref_off);
      Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RXref1;
      Vec3 f1kR(R1.GetCol(3)*Coef_f1);

      // Calculate constraint force/lambda on node 1:
      F1 = -R1.GetCol(3)*Coef_f1;

      // Calculate constraint moment/lambda on node 1:
      M1 = (b1.Cross(R1)).GetCol(3)*Coef_f1;

      // Calculate constraint moment/lambda on reference node 1:
      M1ref = (RXref1.Cross(R1)).GetCol(3)*Coef_f1;

      // Partial calculations of force/moments perturbation of node 1:
      Mat3x3 R1hk_Cross(Mat3x3(MatCross, R1.GetCol(3)*dCoef*Coef_f1*LambdaTmp));
      Mat3x3 RfXo1_Cross(Mat3x3(MatCross, RRf1*x1ref_off));

      // Perturbation of force on node 1:
      WM.Add(1, 3 + 1, R1hk_Cross);
      WM.Add(1, nCnt3, F1);

      // Perturbation of moment on node 1:
      WM.Sub(3 + 1, 1, R1hk_Cross);
      WM.Sub(3 + 1, 3 + 1, b1.Cross(R1hk_Cross));
      WM.Add(3 + 1, 6 + 1, R1hk_Cross);
      WM.Sub(3 + 1, 9 + 1, R1hk_Cross*RfXo1_Cross);
      WM.Add(3 + 1, nCnt3, M1);

      // Perturbation of force on reference node 1:
      WM.Sub(6 + 1, 3 + 1, R1hk_Cross);
      WM.Sub(6 + 1, nCnt3, F1);

      // Perturbation of moment on reference node 1:
      WM.Sub(9 + 1, 3 + 1, RfXo1_Cross*R1hk_Cross);
      WM.Add(9 + 1, 9 + 1, R1hk_Cross*RfXo1_Cross);
      WM.Add(9 + 1, nCnt3, M1ref);

      // Perturbation of Constraint Equation (node 1):
      WM.AddT(nCnt3, 1, f1kR*dCoef);
      WM.AddT(nCnt3, 3 + 1, f1kR.Cross(b1)*dCoef);
      WM.SubT(nCnt3, 6 + 1, f1kR*dCoef);
      WM.AddT(nCnt3, 9 + 1, f1kR.Cross(RXref1)*dCoef);

    } else {

      // Angular node 1:

      // Calculate relative angles:
      Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf1)*R1c));
      Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));

      // LambdaTmp = XCurr(iFirstReactionIndex); // Its value comes from residue calculation
      Mat3x3 MTmp1 = (-Mat3x3(MatCross, R1*Gamma1IT.GetCol(3)) + R1*Gamma1IT*RotManip::Elle(-theta1, Gamma1IT.GetCol(3)))*LambdaTmp*Coef_f1*dCoef;

      // Couple/Lambda in node 1:
      Vec3 C1(R1*Gamma1IT.GetCol(3)*Coef_f1);

      // Perturbation of couple in node 1:
      WM.Add(1, 1, MTmp1);
      WM.Add(1, nCnt3, C1);

      // Perturbation of couple in reference node:
      WM.Sub(3 + 1, 1, MTmp1);
      WM.Sub(3 + 1, nCnt3, C1);

      // Perturbation of Constraint Equation:
      Vec3 dTmp1 = (RotManip::DRot_I(theta1).MulMT(R1)).GetRow(3)*Coef_f1*dCoef;
      WM.SubT(nCnt3, 1, dTmp1);
      WM.AddT(nCnt3, 3 + 1, dTmp1);

    }


    if (mtype2==1) {

      // Linear Node 2:

      // Calculate vector b2 and partial calculation of forces and moments:
      Vec3 RXref2(RRf2*x2ref_off);
      Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RXref2;
      Vec3 f2kR(R2.GetCol(3)*Coef_f2);

      // Calculate constraint force/lambda on node 2:
      F2 = R2.GetCol(3)*Coef_f2;

      // Calculate constraint moment/lambda on node 2:
      M2 = -(b2.Cross(R2)).GetCol(3)*Coef_f2;

      // Calculate constraint moment/lambda on reference node 2:
      M2ref = -(RXref2.Cross(R2)).GetCol(3)*Coef_f2;


      // Partial calculations of force/moments perturbation of node 2:
      Mat3x3 R2hk_Cross(Mat3x3(MatCross, R2.GetCol(3)*dCoef*Coef_f2*LambdaTmp));
      Mat3x3 RfXo2_Cross(Mat3x3(MatCross, RRf2*x2ref_off));

      // Perturbation of force on node 2:
      WM.Sub(nCnt1 + 1, nCnt1 + 4, R2hk_Cross);
      WM.Add(nCnt1 + 1, nCnt3, F2);

      // Perturbation of moment on node 2:
      WM.Add(nCnt1 + 4, nCnt1 + 1, R2hk_Cross);
      WM.Add(nCnt1 + 4, nCnt1 + 4, b2.Cross(R2hk_Cross));
      WM.Sub(nCnt1 + 4, nCnt1 + 7, R2hk_Cross);
      WM.Add(nCnt1 + 4, nCnt1 + 10, R2hk_Cross*RfXo2_Cross);
      WM.Add(nCnt1 + 4, nCnt3, M2);

      // Perturbation of force on reference node 2:
      WM.Add(nCnt1 + 7, nCnt1 + 4, R2hk_Cross);
      WM.Sub(nCnt1 + 7, nCnt3, F2);

      // Perturbation of moment on reference node 2:
      WM.Add(nCnt1 + 10, nCnt1 + 4, RfXo2_Cross*R2hk_Cross);
      WM.Sub(nCnt1 + 10, nCnt1 + 10, R2hk_Cross*RfXo2_Cross);
      WM.Add(nCnt1 + 10, nCnt3, M2ref);


      // Perturbation of Constraint Equation (node 2):
      WM.SubT(nCnt3, nCnt1 + 1, f2kR*dCoef);
      WM.SubT(nCnt3, nCnt1 + 4, f2kR.Cross(b2)*dCoef);
      WM.AddT(nCnt3, nCnt1 + 7, f2kR*dCoef);
      WM.SubT(nCnt3, nCnt1 + 10, f2kR.Cross(RXref2)*dCoef);

    } else {

      // Angular node 2:

      // Calculate relative angles:
      Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf2)*R2c));

      // Couple/Lambda in node 2:
      Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));

      // MTmp = XCurr(iFirstReactionIndex); // Its value comes from residue calculation
      Mat3x3 MTmp2 = (-Mat3x3(MatCross, R2*Gamma2IT.GetCol(3)) + R2*Gamma2IT*RotManip::Elle(-theta2, Gamma2IT.GetCol(3)))*LambdaTmp*Coef_f2*dCoef;

      // Couple/Lambda in node 2:
      Vec3 C2(R2*Gamma2IT.GetCol(3)*Coef_f2);

      // Perturbation of couple in node 2:
      WM.Sub(nCnt1 + 1, nCnt1 + 1, MTmp2);
      WM.Sub(nCnt1 + 1, nCnt3, C2);

      // Perturbation of couple in reference node:
      WM.Add(nCnt1 + 4, nCnt1 + 1, MTmp2);
      WM.Add(nCnt1 + 4, nCnt3, C2);

      // Perturbation of Constraint Equation:
      Vec3 dTmp2 = (RotManip::DRot_I(theta2).MulMT(R2)).GetRow(3)*Coef_f2*dCoef;
      WM.AddT(nCnt3, nCnt1 + 1, dTmp2);
      WM.SubT(nCnt3, nCnt1 + 4, dTmp2);

    }

    return WorkMat;
}

SubVectorHandler&
MotionTransmissionJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering MotionTransmissionJoint::AssRes()" << std::endl);

   /* Change the dimension of the vector based on the constraint demand */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

   /* Get the index of the nodes and reaction moment variables */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNodeRef1FirstMomIndex = pNodeRef1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iNodeRef2FirstMomIndex = pNodeRef2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex()+1;

   /* Set the indexes of equation relative to node 1 */
   int nCnt1 = 3;
   int aCnt = 0;
   if (mtype1==1) nCnt1 += 3; else aCnt += 3;
   for (int iCnt = 1; iCnt <= nCnt1; iCnt++) {
          WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt + aCnt);
          WorkVec.PutRowIndex(nCnt1 + iCnt, iNodeRef1FirstMomIndex + iCnt + aCnt);
   }

   /* Set the indexes of equation relative to node 2 */
   nCnt1 += nCnt1; // 2 nodes = (2*nCnt1)
   aCnt = 0;
   int nCnt2 = 3;
   if (mtype2==1) nCnt2 += 3; else aCnt += 3;
   for (int iCnt = 1; iCnt <= nCnt2; iCnt++) {
          WorkVec.PutRowIndex(nCnt1 + iCnt, iNode2FirstMomIndex + iCnt + aCnt);
          WorkVec.PutRowIndex(nCnt1 + nCnt2 + iCnt, iNodeRef2FirstMomIndex + iCnt + aCnt);
   }

   /* Set the index of constraint equation */
   int nCnt3 = nCnt1 + 2*nCnt2 + 1;

   WorkVec.PutRowIndex(nCnt3, iFirstReactionIndex);

   // Retrieve Lagrange multiplier:
   LambdaTmp = XCurr(iFirstReactionIndex);

   // Retrieve node rotation matrices:
   Mat3x3 R1(pNode1->GetRCurr()*R1tilde);
   Mat3x3 R2(pNode2->GetRCurr()*R2tilde);
   Mat3x3 RRf1(pNodeRef1->GetRCurr());
   Mat3x3 RRf2(pNodeRef2->GetRCurr());

   doublereal eps(0.);

   if (mtype1==1) {

      // Linear node 1:

      // Calculate vectors b1:
      Vec3 RXref1(RRf1*x1ref_off);
      Vec3 b1 = pNode1->GetXCurr() - pNodeRef1->GetXCurr() - RXref1;

      // Calculate constraint force/lambda on node 1:
      F1 = -R1.GetCol(3)*Coef_f1*LambdaTmp;

      // Calculate constraint moment/lambda on node 1:
      M1 = (b1.Cross(R1)).GetCol(3)*Coef_f1*LambdaTmp;

      // Calculate constraint moment/lambda on reference node 1:
      M1ref = (RXref1.Cross(R1)).GetCol(3)*Coef_f1*LambdaTmp;

      // Calculate constraint equation:
      eps += (R1.MulTV(b1)+R1tilde.MulTV(x1_off))(3)*Coef_f1;

      // Put (or add) calculated values in the WorkVec:

      // node 1:
      WorkVec.Sub(1, F1);
      WorkVec.Sub(3 + 1, M1);

      // reference node 1:
      WorkVec.Add(6 + 1, F1);
      WorkVec.Sub(9 + 1, M1ref);

   } else {

      // Angular node 1:
      Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf1)*R1c));

      Mat3x3 Gamma1IT(RotManip::DRot_IT(theta1));

     // Couple :
      Vec3 C1(R1*Gamma1IT.GetCol(3)*Coef_f1*LambdaTmp);

      WorkVec.Sub(1, C1);
      WorkVec.Add(3 + 1, C1);

     // Unwrap angle:
      theta1 = Unwrap(ThetaOut1, theta1);

      eps += theta1(3)*Coef_f1;

   }


   if (mtype2==1) {

      // Linear node 2:

      // Calculate vectors b1 and b2:
      Vec3 RXref2(RRf2*x2ref_off);
      Vec3 b2 = pNode2->GetXCurr() - pNodeRef2->GetXCurr() - RXref2;

      // Calculate constraint force/lambda on node 2:
      F2 = R2.GetCol(3)*Coef_f2*LambdaTmp;

      // Calculate constraint moment/lambda on node 2:
      M2 = -(b2.Cross(R2)).GetCol(3)*Coef_f2*LambdaTmp;

      // Calculate constraint moment/lambda on reference node 2:
      M2ref = -(RXref2.Cross(R2)).GetCol(3)*Coef_f2*LambdaTmp;

      // Calculate constraint equation:
      eps -= (R2.MulTV(b2)+R2tilde.MulTV(x2_off))(3)*Coef_f2;

      // Put (or add) calculated values in the WorkVec:

      // node 2:
      WorkVec.Sub(nCnt1 + 1, F2);
      WorkVec.Sub(nCnt1 + 4, M2);

      // reference node 2:
      WorkVec.Add(nCnt1 + 7, F2);
      WorkVec.Sub(nCnt1 + 10, M2ref);

   } else {

      // Angular node 2:

      Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf2)*R2c));

      Mat3x3 Gamma2IT(RotManip::DRot_IT(theta2));

      Vec3 C2(R2*Gamma2IT.GetCol(3)*Coef_f2*LambdaTmp);

      WorkVec.Add(nCnt1 + 1, C2);
      WorkVec.Sub(nCnt1 + 4, C2);

      // Unwrap angle:
      theta2 = Unwrap(ThetaOut2, theta2);

      eps -= theta2(3)*Coef_f2;

   }

   // constraint:
   WorkVec.PutCoef(nCnt3, -eps + z_corr);

	return WorkVec;
}

unsigned int
MotionTransmissionJoint::iGetNumPrivData(void) const
{
	return 0;
}

int
MotionTransmissionJoint::iGetNumConnectedNodes(void) const
{
	//Useful, but not essential:
	return 4;
}

void
MotionTransmissionJoint::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	//Useful, but not essential:
	connectedNodes.resize(4);
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNodeRef1;
	connectedNodes[2] = pNode2;
	connectedNodes[3] = pNodeRef2;
}

void
MotionTransmissionJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
MotionTransmissionJoint::Restart(std::ostream& out) const
{
	return out << "# MotionTransmissionJoint: not implemented" << std::endl;
}

unsigned int
MotionTransmissionJoint::iGetInitialNumDof(void) const
{
	return 0;
}

void
MotionTransmissionJoint::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
MotionTransmissionJoint::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
MotionTransmissionJoint::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


unsigned int
MotionTransmissionJoint::iGetNumDof(void) const
{
   return 1;
}


DofOrder::Order
MotionTransmissionJoint::GetDofType(unsigned int i) const
{
   return DofOrder::ALGEBRAIC;
}


void
MotionTransmissionJoint::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{

   // Actual rotation angles calculation (local reference):
   if (mtype1==2) {
      Mat3x3 R1(pNode1->GetRCurr()*R1tilde);
      Mat3x3 RRf1(pNodeRef1->GetRCurr());
      Vec3 theta1(RotManip::VecRot(R1.MulTM(RRf1)*R1c));
      ThetaOut1 = Unwrap(ThetaOut1, theta1);
   }

   if (mtype2==2) {
      Mat3x3 R2(pNode2->GetRCurr()*R2tilde);
      Mat3x3 RRf2(pNodeRef2->GetRCurr());
      Vec3 theta2(RotManip::VecRot(R2.MulTM(RRf2)*R2c));
      ThetaOut2 = Unwrap(ThetaOut2, theta2);
   }
}

// Smooth step driver
// based on http://en.wikipedia.org/wiki/Smoothstep - smootherstep (Ken Perlin's version)

class SmoothStepDriveCaller : public DriveCaller {
private:
	doublereal dStartTime;
	doublereal dEndTime;
	doublereal dInitialValue;
	doublereal dFinalValue;

public:
	SmoothStepDriveCaller(const DriveHandler *pDH, const doublereal &d1, const doublereal &d2,
        const doublereal &d3, const doublereal &d4);

	virtual ~SmoothStepDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	inline doublereal dGet(const doublereal& dVar) const;

#if 0
	inline doublereal dGet(void) const;
#endif

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;

#if 0
	virtual inline doublereal dGetP(void) const;
#endif

};

SmoothStepDriveCaller::SmoothStepDriveCaller(const DriveHandler *pDH, const doublereal &d1,
    const doublereal &d2, const doublereal &d3, const doublereal &d4)
    : DriveCaller(pDH), dStartTime(d1), dEndTime(d2), dInitialValue(d3), dFinalValue(d4)
{
	NO_OP;
}

SmoothStepDriveCaller::~SmoothStepDriveCaller(void)
{
	NO_OP;
}

DriveCaller *
SmoothStepDriveCaller::pCopy(void) const
{
    DriveCaller* pDC = 0;
    SAFENEWWITHCONSTRUCTOR(pDC,
        SmoothStepDriveCaller,
        SmoothStepDriveCaller(pDrvHdl, dStartTime, dEndTime, dInitialValue, dFinalValue));
    return pDC;
}

std::ostream&
SmoothStepDriveCaller::Restart(std::ostream& out) const
{
	return out << "SmoothStep, " <<
        dStartTime <<
        dEndTime <<
        dInitialValue <<
        dFinalValue;
}

inline doublereal
SmoothStepDriveCaller::dGet(const doublereal& dVar) const
{
    doublereal dVal = dInitialValue;
	if (dVar>dStartTime) {
        if (dVar<dEndTime) {
            doublereal t1 = (dVar-dStartTime)/(dEndTime-dStartTime);
            dVal += (dFinalValue-dInitialValue)*((t1*t1*t1*(t1*(t1*6. - 15.) + 10.)));
        }
        else dVal = dFinalValue;
    };

	return dVal;
}

/*
inline doublereal
SmoothStepDriveCaller::dGet(void) const
{
	return dConst;
}
*/

inline bool
SmoothStepDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

inline doublereal
SmoothStepDriveCaller::dGetP(const doublereal& dVar) const
{
    doublereal dVal = 0.;

	if (dVar>dStartTime) {
        if (dVar<dEndTime) {
            doublereal t1 = (dVar-dStartTime)/(dEndTime-dStartTime);
            dVal += (dFinalValue-dInitialValue)*((30*t1*t1*(t1*(t1 - 2.) + 1.)));
        }
    };

	return dVal;
}

/*
inline doublereal
SmoothStepDriveCaller::dGetP(void) const
{
	return 0.;
}

/*

/* prototype of the functional object: reads a drive caller */
struct SmoothStepDCR : public DriveCallerRead {
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) {
        const DriveHandler* pDrvHdl = 0;
        if (pDM != 0) {
            pDrvHdl = pDM->pGetDrvHdl();
        }
		// Get initial time:
		doublereal d1 = HP.GetReal(0.);
        DEBUGCOUT("Initial time: " << d1 << std::endl);

		// Get final time:
		doublereal d2 = HP.GetReal(1.);
        DEBUGCOUT("Final time: " << d2 << std::endl);

		// Get initial value:
		doublereal d3 = HP.GetReal(0.);
        DEBUGCOUT("Initial value: " << d3 << std::endl);

		// Get final value:
		doublereal d4 = HP.GetReal(1.);
        DEBUGCOUT("Final value: " << d4 << std::endl);

        DriveCaller* pDC = 0;

        SAFENEWWITHCONSTRUCTOR(pDC,
            SmoothStepDriveCaller,
            SmoothStepDriveCaller(pDrvHdl, d1, d2, d3, d4));

        return pDC;

	};
};



extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<GearJoint>;

	if (!SetUDE("gear" "joint", rf1)) {
		delete rf1;

		silent_cerr("module-fabricate (gear joint): "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf2 = new UDERead<LinearTransmissionJoint>;

	if (!SetUDE("linear" "transmission" "joint", rf2)) {
		delete rf2;

		silent_cerr("module-fabricate (linear joint): "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf3 = new UDERead<MotionTransmissionJoint>;

	if (!SetUDE("motion" "transmission" "joint", rf3)) {
		delete rf3;

		silent_cerr("module-fabricate (motion joint): "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}


	DriveCallerRead	*rf4 = new SmoothStepDCR;

	if (!SetDriveCallerData("smooth" "step", rf4)) {
		delete rf4;

		silent_cerr("SmoothStepDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}


	return 0;
}

