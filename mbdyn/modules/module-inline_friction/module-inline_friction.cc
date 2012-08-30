/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/*
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2012) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-inline_friction.h"

#define FORCE_UPDATE_ASSRES 0
#define FORCE_UPDATE_ASSJAC 0
#define FORCE_UPDATE_AFTERCONVERGENCE 0
#define FORCE_UPDATE_AFTERPREDICT 1

class InlineFriction: virtual public Elem, public UserDefinedElem
{
public:
	InlineFriction(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~InlineFriction(void);
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;
	virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
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
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
	virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

  private:
   	inline doublereal SlidingVelocity(const Vec3& X1, const Mat3x3& R1, const Vec3& XP1, const Vec3& omega1, const Vec3& X2, const Mat3x3& R2, const Vec3& XP2, const Vec3& omega2) const;
   	inline doublereal NormalForceMagnitude() const;
   	inline doublereal FrictionCoefficient(doublereal DeltaXP) const;
   	inline doublereal FrictionForce(doublereal DeltaXP) const;

  private:
   	StructNode* pNode1;
   	Vec3 o1;
   	Mat3x3 e;
   	StructNode* pNode2;
   	Vec3 o2;

	enum LagrangeMultiplierIndex_t
	{
		L1 = 0,
		L2 = 1
	};

	doublereal lambda[2];
	doublereal z;
	doublereal zP;

	doublereal mus;
	doublereal muc;
	doublereal vs;
	doublereal iv;
	doublereal kv;
	doublereal delta;
	doublereal PhiScale;
};

InlineFriction::InlineFriction(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	pNode1(0),
	o1(Zero3),
	e(Eye3),
	pNode2(0),
	o2(Zero3),
	z(0.),
	zP(0.),
	mus(0.),
	muc(0.),
	vs(0.),
	iv(1.),
	kv(0.),
	delta(0.),
	PhiScale(1.)
{
	for (int i = 0; i < 2; ++i) {
		lambda[i] = 0;
	}

	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"\n"
			"Module: 	InlineFriction\n"
			"\n"
			"	This element implements a inline joint with friction\n"
			"\n"
			"	inline friction,\n"
			"		node1, (label) <node1>,\n"
			"			[ offset, (Vec3) <offset>, ]\n"
			"			[ hinge, (Mat3x3) <orientation>, ]\n"
			"		node2, (label) <node2>,\n"
			"			[ offset, (Vec3) <offset>, ]\n"
			"		coulomb friction coefficient, (real) <muc>,\n"
			"		static friction coefficient, (real) <mus>,\n"
			"		sliding velocity coefficient, (real) <vs>,\n"
			"		[ sliding velocity exponent, (real) <i>, ]\n"
			"		micro slip displacement, (real) <delta>,\n"
			"   	[ initial stiction state, (real) <z0>, ]\n"
			"   	[ initial stiction derivative, (real) <zP0>, ]\n"
			"		[ viscous friction coefficient, (real) <kv>, ]\n"
			"		[ stiction state equation scale, (real) <PhiScale> ]\n"
			"\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if ( !HP.IsKeyWord("node1") ) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"node1\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pNode1 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

	if (!pNode1) {
		silent_cerr("inline friction(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const ReferenceFrame refNode1(pNode1);

	if (HP.IsKeyWord("offset")) {
		o1 = HP.GetPosRel(refNode1);
	}

	if (HP.IsKeyWord("hinge") || HP.IsKeyWord("orientation")) {
		e = HP.GetRotRel(refNode1);
	}

	if (!HP.IsKeyWord("node2")) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"node2\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pNode2 = dynamic_cast<StructNode*>(pDM->ReadNode(HP,Node::STRUCTURAL));

	if (!pNode2) {
		silent_cerr("inline friction(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("offset")) {
		const ReferenceFrame refNode2(pNode2);

		o2 = HP.GetPosRel(refNode2);
	}

	if (!HP.IsKeyWord("coulomb" "friction" "coefficient")) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"coulomb friction coefficient\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	muc = HP.GetReal();

	if (muc <= 0) {
		silent_cerr("inline friction(" << GetLabel() << "): coulomb friction coefficient must be greater than zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HP.IsKeyWord("static" "friction" "coefficient")) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"static friction coefficient\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	mus = HP.GetReal();

	if (mus < 0) {
		silent_cerr("inline friction(" << GetLabel() << "): static friction coefficient must be greater than or equal to zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HP.IsKeyWord("sliding" "velocity" "coefficient")) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"sliding velocity coefficient\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	vs = HP.GetReal();

	if (vs <= 0.) {
		silent_cerr("inline friction(" << GetLabel() << "): sliding velocity coefficient must be greater than zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("sliding" "velocity" "exponent")) {
		iv = HP.GetReal();
	}

	if (!HP.IsKeyWord("micro" "slip" "displacement")) {
		silent_cerr("inline friction(" << GetLabel() << "): keyword \"micro slip displacement\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	delta = HP.GetReal();

	if (delta <= 0) {
		silent_cerr("inline friction(" << GetLabel() << "): micro slip displacement must be greater than zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("initial" "stiction" "state")) {
		z = HP.GetReal();
	}

	if (std::abs(z) > 1.) {
		silent_cerr("inline friction(" << GetLabel() << "): initial stiction state must be between -1 and 1 at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("initial" "stiction" "derivative")) {
		zP = HP.GetReal();
	}

	z *= delta;
	zP *= delta;

	if (HP.IsKeyWord("viscous" "friction" "coefficient")) {
		kv = HP.GetReal();
	}

	if (kv < 0) {
		silent_cerr("inline friction(" << GetLabel() << "): viscous friction coefficient must be greater than or equal to zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("stiction" "state" "equation" "scale")) {
		PhiScale = HP.GetReal();
	}

	if (PhiScale == 0.) {
		silent_cerr("inline friction(" << GetLabel() << "): stiction state equation scale must not be equal to zero at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	std::ostream& out = pDM->GetLogFile();

	out << "inline friction: " << GetLabel() << " "
		<< pNode1->GetLabel() << " "
		<< o1 << " "
		<< e << " "
		<< pNode2->GetLabel() << " " << o2 << " "
		<< muc << " "
		<< mus << " "
		<< vs << " "
		<< iv << " "
		<< kv << " "
		<< delta << " "
		<< z << " "
		<< zP
		<< std::endl;
}

InlineFriction::~InlineFriction(void)
{
	// destroy private data
}

unsigned int InlineFriction::iGetNumDof(void) const
{
	return 3u;
}

DofOrder::Order InlineFriction::GetDofType(unsigned int i) const
{
	switch (i) {
		case 0:
		case 1:
			return DofOrder::ALGEBRAIC;
		case 2:
			return DofOrder::DIFFERENTIAL;
		default:
			ASSERT(0);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

DofOrder::Order InlineFriction::GetEqType(unsigned int i) const
{
	switch (i) {
		case 0:
		case 1:
			return DofOrder::ALGEBRAIC;
		case 2:
			return DofOrder::DIFFERENTIAL;
		default:
			ASSERT(0);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

std::ostream& InlineFriction::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	const integer iFirstIndex = iGetFirstIndex();

	out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": reaction forces [lambda1, lambda2]" << std::endl;

	if (bInitial) {
		out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": reaction force derivatives [lambdaP1, lambdaP2]" << std::endl;
	} else {
		out << prefix << iFirstIndex + 3 << ": stiction state [z]" << std::endl;
	}

	return out;
}

std::ostream& InlineFriction::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	const integer iFirstIndex = iGetFirstIndex();

	out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 2 << ": position constraints [c1, c2]" << std::endl;

	if (bInitial) {
		out << prefix << iFirstIndex + 3 << "->" << iFirstIndex + 4 << ": velocity constraints [cP1, cP2]" << std::endl;
	} else {
		out << prefix << iFirstIndex + 3 << ": stick slip transition [Phi]" << std::endl;
	}

	return out;
}

unsigned int InlineFriction::iGetNumPrivData(void) const
{
	return 7;
}

unsigned int InlineFriction::iGetPrivDataIdx(const char *s) const
{
	static const struct {
		int index;
		char name[8];
	} data[] = {
			{ 1, "lambda1" },
			{ 2, "lambda2" },
			{ 3, "tau" },
			{ 4, "mu" },
			{ 5, "z" },
			{ 6, "zP" },
			{ 7, "v" }
	};

	const int N = sizeof(data) / sizeof(data[0]);

	for (int i = 0; i < N; ++i) {
		if (0 == strcmp(data[i].name, s)) {
			return data[i].index;
		}
	}

	return 0;
}

doublereal InlineFriction::dGetPrivData(unsigned int i) const
{
	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Vec3& XP1 = pNode1->GetVCurr();
	const Vec3& omega1 = pNode1->GetWCurr();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Vec3& XP2 = pNode2->GetVCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	switch (i) {
		case 1:
		case 2:
			return lambda[i - 1];
		case 3:
			return FrictionForce(SlidingVelocity(X1, R1, XP1, omega1, X2, R2, XP2, omega2));
		case 4:
			return FrictionCoefficient(SlidingVelocity(X1, R1, XP1, omega1, X2, R2, XP2, omega2));
		case 5:
			return z;
		case 6:
			return zP;
		case 7:
			return SlidingVelocity(X1, R1, XP1, omega1, X2, R2, XP2, omega2);
		default:
			silent_cerr("inline friction(" << GetLabel() << "): invalid private data index " << i << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
InlineFriction::Output(OutputHandler& OH) const
{
	if ( fToBeOutput() )
	{
		if ( OH.UseText(OutputHandler::LOADABLE) )
		{
			std::ostream& os = OH.Loadable();

			os << std::setw(8) << GetLabel() << " " << lambda[L1] << " " << lambda[L2] << " " << z << " " << zP << std::endl;
		}
	}
}

void
InlineFriction::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = *piNumCols = 15;
}

VariableSubMatrixHandler&
InlineFriction::AssJac(VariableSubMatrixHandler& WorkMatVar,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
#if FORCE_UPDATE_ASSJAC == 1
	Update(XCurr, XPrimeCurr);
#endif

	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const integer iFirstMomentumIndexNode1 = pNode1->iGetFirstMomentumIndex();
	const integer iFirstPositionIndexNode1 = pNode1->iGetFirstPositionIndex();
	const integer iFirstMomentumIndexNode2 = pNode2->iGetFirstMomentumIndex();
	const integer iFirstPositionIndexNode2 = pNode2->iGetFirstPositionIndex();
	const integer iFirstIndex = iGetFirstIndex();

	for (integer i = 1; i <= 6; ++i) {
		WorkMat.PutRowIndex(i, iFirstMomentumIndexNode1 + i);
		WorkMat.PutColIndex(i, iFirstPositionIndexNode1 + i);
		WorkMat.PutRowIndex(i + 6, iFirstMomentumIndexNode2 + i);
		WorkMat.PutColIndex(i + 6, iFirstPositionIndexNode2 + i);
	}

	for (integer i = 1; i <= 3; ++i) {
		WorkMat.PutRowIndex(12 + i, iFirstIndex + i);
		WorkMat.PutColIndex(12 + i, iFirstIndex + i);
	}

	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Mat3x3& R1_0 = pNode1->GetRRef();
	const Vec3& XP1 = pNode1->GetVCurr();
	const Vec3& omega1 = pNode1->GetWCurr();
	const Vec3& omega1_0 = pNode1->GetWRef();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Mat3x3& R2_0 = pNode2->GetRRef();
	const Vec3& XP2 = pNode2->GetVCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	const doublereal DeltaXP = SlidingVelocity(X1, R1, XP1, omega1, X2, R2, XP2, omega2);

	// common subexpressions
	const Vec3 R1e1 = R1 * e.GetCol(1);
	const Vec3 R1e2 = R1 * e.GetCol(2);
	const Vec3 R1e3 = R1 * e.GetCol(3);
	const Vec3 R2o2 = R2 * o2;
	const Vec3 R2_0o2 = R2_0 * o2;
	const Vec3 l1 = X2 + R2o2 - X1;

	const Vec3 dDeltaXP_dX1_T = -omega1.Cross(R1e1);
	const Vec3 dDeltaXP_dg1_T = -(XP2 + omega2.Cross(R2o2) - XP1 - omega1.Cross(l1)).Cross(R1_0 * e.GetCol(1))
								- omega1_0.Cross((l1).Cross(R1e1));
	const Vec3 dDeltaXP_dXP1_T = -R1e1;
	const Vec3 dDeltaXP_dgP1_T = -l1.Cross(R1e1);
	const Vec3 dDeltaXP_dX2_T = omega1.Cross(R1e1);
	const Vec3 dDeltaXP_dg2_T = R2_0o2.Cross((omega1 - omega2).Cross(R1e1)) + omega2.Cross(R2o2.Cross(R1e1));
	const Vec3& dDeltaXP_dXP2_T = R1e1;
	const Vec3 dDeltaXP_dgP2_T = R2o2.Cross(R1e1);

	const doublereal tau = FrictionForce(DeltaXP);
	const doublereal mu = FrictionCoefficient(DeltaXP);
	const doublereal lambda_res = NormalForceMagnitude();

	const doublereal dtau_dDeltaXP = kv - muc * lambda_res * z / delta * (mus / muc - 1) * exp(-std::pow(std::abs(DeltaXP) / vs, iv))
									 * iv * std::pow(std::abs(DeltaXP) / vs, iv - 1) * copysign(1., DeltaXP) / vs;
	doublereal dtau_dlambda[2];

	if (lambda_res != 0) {
		for (int i = 0; i < 2; ++i) {
			dtau_dlambda[i] = mu * lambda[i] * z / (lambda_res * delta);
		}
	} else { // avoid division by zero
		for (int i = 0; i < 2; ++i) {
			dtau_dlambda[i] = mu * copysign(1., lambda[i]) * z / delta;
		}
	}

	const doublereal dtau_dz = mu * lambda_res / delta;

	const Vec3 dtau_dX1_T = dDeltaXP_dX1_T * dtau_dDeltaXP;
	const Vec3 dtau_dg1_T = dDeltaXP_dg1_T * dtau_dDeltaXP;
	const Vec3 dtau_dXP1_T = dDeltaXP_dXP1_T * dtau_dDeltaXP;
	const Vec3 dtau_dgP1_T = dDeltaXP_dgP1_T * dtau_dDeltaXP;
	const Vec3 dtau_dX2_T = dDeltaXP_dX2_T * dtau_dDeltaXP;
	const Vec3 dtau_dg2_T = dDeltaXP_dg2_T * dtau_dDeltaXP;
	const Vec3 dtau_dXP2_T = dDeltaXP_dXP2_T * dtau_dDeltaXP;
	const Vec3 dtau_dgP2_T = dDeltaXP_dgP2_T * dtau_dDeltaXP;

	const Vec3 F1 = R1e1 * tau + R1e2 * lambda[L1] + R1e3 * lambda[L2];
	const Vec3 M1 = l1.Cross(F1);
	// const Vec3 F2 = -F1;
	// const Vec3 M2 = R2o2.Cross(F2);

	const Mat3x3 dF1_dX1 = R1e1.Tens(dtau_dX1_T);
	const Mat3x3 dF1_dg1 = R1e1.Tens(dtau_dg1_T) - Mat3x3(MatCross, R1_0 * (e.GetCol(1) * tau + e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]));
	const Mat3x3 dF1_dXP1 = R1e1.Tens(dtau_dXP1_T);
	const Mat3x3 dF1_dgP1 = R1e1.Tens(dtau_dgP1_T);
	const Mat3x3 dF1_dX2 = R1e1.Tens(dtau_dX2_T);
	const Mat3x3 dF1_dg2 = R1e1.Tens(dtau_dg2_T);
	const Mat3x3 dF1_dXP2 = R1e1.Tens(dtau_dXP2_T);
	const Mat3x3 dF1_dgP2 = R1e1.Tens(dtau_dgP2_T);
	const Vec3 dF1_dlambda1 = R1e1 * dtau_dlambda[L1] + R1e2;
	const Vec3 dF1_dlambda2 = R1e1 * dtau_dlambda[L2] + R1e3;
	const Vec3 dF1_dz = R1e1 * dtau_dz;

	const Mat3x3 dM1_dX1 = l1.Cross(dF1_dX1) + Mat3x3(MatCross, F1);
	const Mat3x3 dM1_dg1 = l1.Cross(dF1_dg1);
	const Mat3x3 dM1_dXP1 = l1.Cross(dF1_dXP1);
	const Mat3x3 dM1_dgP1 = l1.Cross(dF1_dgP1);
	const Mat3x3 dM1_dX2 = l1.Cross(dF1_dX2) - Mat3x3(MatCross, F1);
	const Mat3x3 dM1_dg2 = l1.Cross(dF1_dg2) + Mat3x3(MatCrossCross, F1, R2_0o2);
	const Mat3x3 dM1_dXP2 = l1.Cross(dF1_dXP2);
	const Mat3x3 dM1_dgP2 = l1.Cross(dF1_dgP2);
	const Vec3 dM1_dlambda1 = l1.Cross(dF1_dlambda1);
	const Vec3 dM1_dlambda2 = l1.Cross(dF1_dlambda2);
	const Vec3 dM1_dz = l1.Cross(dF1_dz);

	const Mat3x3 dM2_dX1 = (-R2o2).Cross(dF1_dX1);
	const Mat3x3 dM2_dg1 = (-R2o2).Cross(dF1_dg1);
	const Mat3x3 dM2_dXP1 = (-R2o2).Cross(dF1_dXP1);
	const Mat3x3 dM2_dgP1 = (-R2o2).Cross(dF1_dgP1);
	const Mat3x3 dM2_dX2 = (-R2o2).Cross(dF1_dX2);
	const Mat3x3 dM2_dg2 = Mat3x3(MatCrossCross, -F1, R2_0o2) - R2o2.Cross(dF1_dg2);
	const Mat3x3 dM2_dXP2 = (-R2o2).Cross(dF1_dXP2);
	const Mat3x3 dM2_dgP2 = (-R2o2).Cross(dF1_dgP2);
	const Vec3 dM2_dlambda1 = (-R2o2).Cross(dF1_dlambda1);
	const Vec3 dM2_dlambda2 = (-R2o2).Cross(dF1_dlambda2);
	const Vec3 dM2_dz = (-R2o2).Cross(dF1_dz);

	const Vec3 dc1_dX1_T = -R1e2;
	const Vec3 dc1_dg1_T = -l1.Cross(R1_0 * e.GetCol(2));
	const Vec3& dc1_dX2_T = R1e2;
	const Vec3 dc1_dg2_T = R2_0o2.Cross(R1e2);

	const Vec3 dc2_dX1_T = -R1e3;
	const Vec3 dc2_dg1_T = -l1.Cross(R1_0 * e.GetCol(3));
	const Vec3& dc2_dX2_T = R1e3;
	const Vec3 dc2_dg2_T = R2_0o2.Cross(R1e3);

	const doublereal alpha = z / delta * copysign(1., DeltaXP) - 1;

	const Vec3 dPhi_dX1_T = dDeltaXP_dX1_T * alpha;
	const Vec3 dPhi_dg1_T = dDeltaXP_dg1_T * alpha;
	const Vec3 dPhi_dXP1_T = dDeltaXP_dXP1_T * alpha;
	const Vec3 dPhi_dgP1_T = dDeltaXP_dgP1_T * alpha;
	const Vec3 dPhi_dX2_T = dDeltaXP_dX2_T * alpha;
	const Vec3 dPhi_dg2_T = dDeltaXP_dg2_T * alpha;
	const Vec3 dPhi_dXP2_T = dDeltaXP_dXP2_T * alpha;
	const Vec3 dPhi_dgP2_T = dDeltaXP_dgP2_T * alpha;

	const doublereal dPhi_dz = std::abs(DeltaXP) / delta;
	const doublereal dPhi_dzP = 1.;

	WorkMat.Put(1,  1, -dF1_dXP1 - dF1_dX1 * dCoef);
	WorkMat.Put(1,  4, -dF1_dgP1 - dF1_dg1 * dCoef);
	WorkMat.Put(1,  7, -dF1_dXP2 - dF1_dX2 * dCoef);
	WorkMat.Put(1, 10, -dF1_dgP2 - dF1_dg2 * dCoef);
	WorkMat.Put(1, 13, -dF1_dlambda1);
	WorkMat.Put(1, 14, -dF1_dlambda2);
	WorkMat.Put(1, 15, -dF1_dz * dCoef);

	WorkMat.Put(4,  1, -dM1_dXP1 - dM1_dX1 * dCoef);
	WorkMat.Put(4,  4, -dM1_dgP1 - dM1_dg1 * dCoef);
	WorkMat.Put(4,  7, -dM1_dXP2 - dM1_dX2 * dCoef);
	WorkMat.Put(4, 10, -dM1_dgP2 - dM1_dg2 * dCoef);
	WorkMat.Put(4, 13, -dM1_dlambda1);
	WorkMat.Put(4, 14, -dM1_dlambda2);
	WorkMat.Put(4, 15, -dM1_dz * dCoef);

	WorkMat.Put(7,  1, dF1_dXP1 + dF1_dX1 * dCoef);
	WorkMat.Put(7,  4, dF1_dgP1 + dF1_dg1 * dCoef);
	WorkMat.Put(7,  7, dF1_dXP2 + dF1_dX2 * dCoef);
	WorkMat.Put(7, 10, dF1_dgP2 + dF1_dg2 * dCoef);
	WorkMat.Put(7, 13, dF1_dlambda1);
	WorkMat.Put(7, 14, dF1_dlambda2);
	WorkMat.Put(7, 15, dF1_dz * dCoef);

	WorkMat.Put(10,  1, -dM2_dXP1 - dM2_dX1 * dCoef);
	WorkMat.Put(10,  4, -dM2_dgP1 - dM2_dg1 * dCoef);
	WorkMat.Put(10,  7, -dM2_dXP2 - dM2_dX2 * dCoef);
	WorkMat.Put(10, 10, -dM2_dgP2 - dM2_dg2 * dCoef);
	WorkMat.Put(10, 13, -dM2_dlambda1);
	WorkMat.Put(10, 14, -dM2_dlambda2);
	WorkMat.Put(10, 15, -dM2_dz * dCoef);

	WorkMat.PutT(13,  1, -dc1_dX1_T);
	WorkMat.PutT(13,  4, -dc1_dg1_T);
	WorkMat.PutT(13,  7, -dc1_dX2_T);
	WorkMat.PutT(13, 10, -dc1_dg2_T);

	WorkMat.PutT(14,  1, -dc2_dX1_T);
	WorkMat.PutT(14,  4, -dc2_dg1_T);
	WorkMat.PutT(14,  7, -dc2_dX2_T);
	WorkMat.PutT(14, 10, -dc2_dg2_T);

	WorkMat.PutT(15,  1, (dPhi_dXP1_T + dPhi_dX1_T * dCoef) * (-PhiScale));
	WorkMat.PutT(15,  4, (dPhi_dgP1_T + dPhi_dg1_T * dCoef) * (-PhiScale));
	WorkMat.PutT(15,  7, (dPhi_dXP2_T + dPhi_dX2_T * dCoef) * (-PhiScale));
	WorkMat.PutT(15, 10, (dPhi_dgP2_T + dPhi_dg2_T * dCoef) * (-PhiScale));
	WorkMat.PutCoef(15, 15, (dPhi_dzP + dPhi_dz * dCoef) * (-PhiScale));

	return WorkMatVar;
}


SubVectorHandler&
InlineFriction::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
#if	FORCE_UPDATE_ASSRES == 1
	Update(XCurr, XPrimeCurr);
#endif

	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	const integer iFirstMomentumIndexNode1 = pNode1->iGetFirstMomentumIndex();
	const integer iFirstMomentumIndexNode2 = pNode2->iGetFirstMomentumIndex();
	const integer iFirstIndex = iGetFirstIndex();

	for (integer i = 1; i <= 6; ++i) {
		WorkVec.PutRowIndex(i, iFirstMomentumIndexNode1 + i);
		WorkVec.PutRowIndex(i + 6, iFirstMomentumIndexNode2 + i);
	}

	for (integer i = 1; i <= 3; ++i) {
		WorkVec.PutRowIndex(i + 12, iFirstIndex + i);
	}

	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Vec3& XP1 = pNode1->GetVCurr();
	const Vec3& omega1 = pNode1->GetWCurr();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Vec3& XP2 = pNode2->GetVCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	const Vec3 R2o2 = R2 * o2;
	const Vec3 l1 = X2 + R2o2 - X1;

	const doublereal DeltaXP = SlidingVelocity(X1, R1, XP1, omega1, X2, R2, XP2, omega2);
	const doublereal tau = FrictionForce(DeltaXP);

	const Vec3 F1 = R1 * ( e.GetCol(1) * tau + e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2] );
	const Vec3 M1 = l1.Cross(F1);
	const Vec3 F2 = -F1;
	const Vec3 M2 = R2o2.Cross(F2);

	const Vec3 a = R1.MulTV(l1) - o1;
	const doublereal Phi = zP - DeltaXP * (1. - z / delta * copysign(1., DeltaXP));

	WorkVec.Put(1, F1);
	WorkVec.Put(4, M1);
	WorkVec.Put(7, F2);
	WorkVec.Put(10, M2);

	for (int i = 1; i <= 2; ++i) {
		WorkVec.PutCoef(12 + i, e.GetCol(i + 1).Dot(a) / dCoef);
	}

	WorkVec.PutCoef(15, Phi * PhiScale);

	return WorkVec;
}

int
InlineFriction::iGetNumConnectedNodes(void) const
{
	return 2;
}

void
InlineFriction::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(iGetNumConnectedNodes());
	connectedNodes[0] = pNode1;
	connectedNodes[1] = pNode2;
}

void
InlineFriction::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	const integer iFirstIndex = iGetFirstIndex();

	for (int i = 1; i <= 2; ++i)
		X.PutCoef(iFirstIndex + i, lambda[i - 1]);

	X.PutCoef(iFirstIndex + 3, z);
	XP.PutCoef(iFirstIndex + 3, zP);
}

void InlineFriction::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
#if	FORCE_UPDATE_AFTERPREDICT == 1
	Update(X, XP);
#endif
}

void InlineFriction::Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr)
{
	const integer iFirstIndex = iGetFirstIndex();

	for (int i = 1; i <= 2; ++i)
		lambda[i - 1] = XCurr(iFirstIndex + i);

	z = XCurr(iFirstIndex + 3);
	zP = XPrimeCurr(iFirstIndex + 3);
}

void InlineFriction::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
#if FORCE_UPDATE_AFTERCONVERGENCE == 1
	Update(X, XP);
#endif
}

std::ostream&
InlineFriction::Restart(std::ostream& out) const
{
	return out;
}

unsigned int
InlineFriction::iGetInitialNumDof(void) const
{
	return 4;
}

void
InlineFriction::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = *piNumCols = 28;
}

VariableSubMatrixHandler&
InlineFriction::InitialAssJac(
	VariableSubMatrixHandler& WorkMatVar,
	const VectorHandler& XCurr)
{
	integer iNumRows, iNumCols;

	InitialWorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Mat3x3& R1_0 = pNode1->GetRRef();
	const Vec3& XP1 = pNode1->GetVCurr();
	const Vec3& omega1 = pNode1->GetWCurr();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Mat3x3& R2_0 = pNode2->GetRRef();
	const Vec3& XP2 = pNode2->GetVCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	const integer iFirstIndexNode1 = pNode1->iGetFirstIndex();
	const integer iFirstIndexNode2 = pNode2->iGetFirstIndex();
	const integer iFirstIndex = iGetFirstIndex();

	for (integer i = 1; i <= 12; ++i) {
		WorkMat.PutRowIndex(i, iFirstIndexNode1 + i);
		WorkMat.PutColIndex(i, iFirstIndexNode1 + i);
		WorkMat.PutRowIndex(i + 12, iFirstIndexNode2 + i);
		WorkMat.PutColIndex(i + 12, iFirstIndexNode2 + i);
	}

	for (integer i = 1; i <= 4; ++i) {
		WorkMat.PutRowIndex(i + 24, iFirstIndex + i);
		WorkMat.PutColIndex(i + 24, iFirstIndex + i);
	}

	doublereal lambdaP[2];

	for (int i = 1; i <= 2; ++i) {
		lambda[i - 1] = XCurr.dGetCoef(iFirstIndex + i);
		lambdaP[i - 1] = XCurr.dGetCoef(iFirstIndex + i + 2);
	}

	const Vec3 R1e2 = R1 * e.GetCol(2);
	const Vec3 R1e3 = R1 * e.GetCol(3);
	const Vec3 R2o2 = R2 * o2;
	const Vec3 R2_0o2 = R2_0 * o2;
	const Vec3 l1 = X2 + R2o2 - X1;
	const Vec3 lP1 = XP2 + omega2.Cross(R2o2) - XP1;

	const Vec3 F1 = R1e2 * lambda[L1] + R1e3 * lambda[L2];
	const Vec3 FP1 = omega1.Cross(R1 * (e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]))
					+ R1 * (e.GetCol(2) * lambdaP[L1] + e.GetCol(3) * lambdaP[L2]);
	// const Vec3 F2 = -F1;

	const Mat3x3 dF1_dg1 = Mat3x3(MatCross, R1_0 * (e.GetCol(2) * (-lambda[L1]) + e.GetCol(3) * (-lambda[L2])));
	const Vec3& dF1_dlambda1 = R1e2;
	const Vec3& dF1_dlambda2 = R1e3;

	const Mat3x3 dM1_dX1 = Mat3x3(MatCross, F1);
	const Mat3x3 dM1_dg1 = l1.Cross(dF1_dg1);
	const Mat3x3 dM1_dX2 = -dM1_dX1;
	const Mat3x3 dM1_dg2 = Mat3x3(MatCrossCross, F1, R2_0o2);
	const Vec3 dM1_dlambda1 = l1.Cross(dF1_dlambda1);
	const Vec3 dM1_dlambda2 = l1.Cross(dF1_dlambda2);

	const Mat3x3 dFP1_dg1 = Mat3x3(MatCrossCross, -omega1, R1_0 * ( e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]))
							- Mat3x3(MatCross, R1_0 * ( e.GetCol(2) * lambdaP[L1] + e.GetCol(3) * lambdaP[L2]));
	const Mat3x3 dFP1_domega1 = Mat3x3(MatCross, R1 * (e.GetCol(2) * (-lambda[L1]) + e.GetCol(3) * (-lambda[L2])));
	const Vec3 dFP1_dlambda1 = omega1.Cross(R1e2);
	const Vec3 dFP1_dlambda2 = omega1.Cross(R1e3);
	const Vec3& dFP1_dlambdaP1 = R1e2;
	const Vec3& dFP1_dlambdaP2 = R1e3;

	const Mat3x3 dMP1_dX1 = Mat3x3(MatCross, FP1);
	const Mat3x3 dMP1_dg1 = lP1.Cross(dF1_dg1) + l1.Cross(dFP1_dg1);
	const Mat3x3 dMP1_dXP1 = Mat3x3(MatCross, F1);
	const Mat3x3 dMP1_domega1 = l1.Cross(dFP1_domega1);
	const Mat3x3 dMP1_dX2 = Mat3x3(MatCross, -FP1);
	const Mat3x3 dMP1_dg2 = (Mat3x3(MatCrossCross, F1, omega2) + Mat3x3(MatCross, FP1)).MulVCross(R2_0o2);
	const Mat3x3 dMP1_dXP2 = Mat3x3(MatCross, -F1);
	const Mat3x3 dMP1_domega2 = Mat3x3(MatCrossCross, F1, R2o2);
	const Vec3 dMP1_dlambda1 = lP1.Cross(dF1_dlambda1) + l1.Cross(dFP1_dlambda1);
	const Vec3 dMP1_dlambda2 = lP1.Cross(dF1_dlambda2) + l1.Cross(dFP1_dlambda2);
	const Vec3 dMP1_dlambdaP1 = l1.Cross(dFP1_dlambdaP1);
	const Vec3 dMP1_dlambdaP2 = l1.Cross(dFP1_dlambdaP2);

	const Mat3x3 dM2_dg1 = (-R2o2).Cross(dF1_dg1);
	const Mat3x3 dM2_dg2 = Mat3x3(MatCrossCross, -F1, R2_0o2);
	const Vec3 dM2_dlambda1 = R2o2.Cross(-dF1_dlambda1);
	const Vec3 dM2_dlambda2 = R2o2.Cross(-dF1_dlambda2);

	const Mat3x3 dMP2_dg1 = (-omega2.Cross(R2o2)).Cross(dF1_dg1) - R2o2.Cross(dFP1_dg1);
	const Mat3x3 dMP2_domega1 = (-R2o2).Cross(dFP1_domega1);
	const Mat3x3 dMP2_dg2 = (Mat3x3(MatCrossCross, -F1, omega2) - Mat3x3(MatCross, FP1)).MulVCross(R2_0o2);
	const Mat3x3 dMP2_domega2 = Mat3x3(MatCrossCross, -F1, R2o2);
	const Vec3 dMP2_dlambda1 = (-omega2.Cross(R2o2)).Cross(dF1_dlambda1) - R2o2.Cross(dFP1_dlambda1);
	const Vec3 dMP2_dlambda2 = (-omega2.Cross(R2o2)).Cross(dF1_dlambda2) - R2o2.Cross(dFP1_dlambda2);
	const Vec3 dMP2_dlambdaP1 = -R2o2.Cross(dFP1_dlambdaP1);
	const Vec3 dMP2_dlambdaP2 = -R2o2.Cross(dFP1_dlambdaP2);

	const Vec3 dc1_dX1_T = -R1e2;
	const Vec3 dc1_dg1_T = -l1.Cross(R1_0 * e.GetCol(2));
	const Vec3& dc1_dX2_T = R1e2;
	const Vec3 dc1_dg2_T = R2_0o2.Cross(R1e2);

	const Vec3 dc2_dX1_T = -R1e3;
	const Vec3 dc2_dg1_T = -l1.Cross(R1_0 * e.GetCol(3));
	const Vec3& dc2_dX2_T = R1e3;
	const Vec3 dc2_dg2_T = R2_0o2.Cross(R1e3);

	const Vec3 dcP1_dX1_T = -omega1.Cross(R1e2);
	const Vec3 dcP1_dg1_T = (omega1.Cross(l1) - lP1).Cross(R1_0 * e.GetCol(2));
	const Vec3 dcP1_dXP1_T = -R1e2;
	const Vec3 dcP1_domega1_T = -l1.Cross(R1e2);
	const Vec3 dcP1_dX2_T = omega1.Cross(R1e2);
	const Vec3 dcP1_dg2_T = R2_0o2.Cross((omega1 - omega2).Cross(R1e2));
	const Vec3& dcP1_dXP2_T = R1e2;
	const Vec3 dcP1_domega2_T = R2o2.Cross(R1e2);

	const Vec3 dcP2_dX1_T = -omega1.Cross(R1e3);
	const Vec3 dcP2_dg1_T = (omega1.Cross(l1) - lP1).Cross(R1_0 * e.GetCol(3));
	const Vec3 dcP2_dXP1_T = -R1e3;
	const Vec3 dcP2_domega1_T = -l1.Cross(R1e3);
	const Vec3 dcP2_dX2_T = omega1.Cross(R1e3);
	const Vec3 dcP2_dg2_T = R2_0o2.Cross((omega1 - omega2).Cross(R1e3));
	const Vec3& dcP2_dXP2_T = R1e3;
	const Vec3 dcP2_domega2_T = R2o2.Cross(R1e3);

	WorkMat.Put( 1,  4, -dF1_dg1);
	WorkMat.Put( 1, 25, -dF1_dlambda1);
	WorkMat.Put( 1, 26, -dF1_dlambda2);

	WorkMat.Put( 4,  1, -dM1_dX1);
	WorkMat.Put( 4,  4, -dM1_dg1);
	WorkMat.Put( 4, 13, -dM1_dX2);
	WorkMat.Put( 4, 16, -dM1_dg2);
	WorkMat.Put( 4, 25, -dM1_dlambda1);
	WorkMat.Put( 4, 26, -dM1_dlambda2);

	WorkMat.Put( 7,  4, -dFP1_dg1);
	WorkMat.Put( 7, 10, -dFP1_domega1);
	WorkMat.Put( 7, 25, -dFP1_dlambda1);
	WorkMat.Put( 7, 26, -dFP1_dlambda2);
	WorkMat.Put( 7, 27, -dFP1_dlambdaP1);
	WorkMat.Put( 7, 28, -dFP1_dlambdaP2);

	WorkMat.Put(10,  1, -dMP1_dX1);
	WorkMat.Put(10,  4, -dMP1_dg1);
	WorkMat.Put(10,  7, -dMP1_dXP1);
	WorkMat.Put(10, 10, -dMP1_domega1);
	WorkMat.Put(10, 13, -dMP1_dX2);
	WorkMat.Put(10, 16, -dMP1_dg2);
	WorkMat.Put(10, 19, -dMP1_dXP2);
	WorkMat.Put(10, 22, -dMP1_domega2);
	WorkMat.Put(10, 25, -dMP1_dlambda1);
	WorkMat.Put(10, 26, -dMP1_dlambda2);
	WorkMat.Put(10, 27, -dMP1_dlambdaP1);
	WorkMat.Put(10, 28, -dMP1_dlambdaP2);

	WorkMat.Put(13,  4, dF1_dg1);   	// dF2_dg1 = -dF1_dg1
	WorkMat.Put(13, 25, dF1_dlambda1);
	WorkMat.Put(13, 26, dF1_dlambda2);

	WorkMat.Put(16,  4, -dM2_dg1);
	WorkMat.Put(16, 16, -dM2_dg2);
	WorkMat.Put(16, 25, -dM2_dlambda1);
	WorkMat.Put(16, 26, -dM2_dlambda2);

	WorkMat.Put(19,  4, dFP1_dg1);
	WorkMat.Put(19, 10, dFP1_domega1);
	WorkMat.Put(19, 25, dFP1_dlambda1);
	WorkMat.Put(19, 26, dFP1_dlambda2);
	WorkMat.Put(19, 27, dFP1_dlambdaP1);
	WorkMat.Put(19, 28, dFP1_dlambdaP2);

	WorkMat.Put(22,  4, -dMP2_dg1);
	WorkMat.Put(22, 10, -dMP2_domega1);
	WorkMat.Put(22, 16, -dMP2_dg2);
	WorkMat.Put(22, 22, -dMP2_domega2);
	WorkMat.Put(22, 25, -dMP2_dlambda1);
	WorkMat.Put(22, 26, -dMP2_dlambda2);
	WorkMat.Put(22, 27, -dMP2_dlambdaP1);
	WorkMat.Put(22, 28, -dMP2_dlambdaP2);

	WorkMat.PutT(25,  1, -dc1_dX1_T);
	WorkMat.PutT(25,  4, -dc1_dg1_T);
	WorkMat.PutT(25, 13, -dc1_dX2_T);
	WorkMat.PutT(25, 16, -dc1_dg2_T);

	WorkMat.PutT(26,  1, -dc2_dX1_T);
	WorkMat.PutT(26,  4, -dc2_dg1_T);
	WorkMat.PutT(26, 13, -dc2_dX2_T);
	WorkMat.PutT(26, 16, -dc2_dg2_T);

	WorkMat.PutT(27,  1, -dcP1_dX1_T);
	WorkMat.PutT(27,  4, -dcP1_dg1_T);
	WorkMat.PutT(27,  7, -dcP1_dXP1_T);
	WorkMat.PutT(27, 10, -dcP1_domega1_T);
	WorkMat.PutT(27, 13, -dcP1_dX2_T);
	WorkMat.PutT(27, 16, -dcP1_dg2_T);
	WorkMat.PutT(27, 19, -dcP1_dXP2_T);
	WorkMat.PutT(27, 22, -dcP1_domega2_T);

	WorkMat.PutT(28,  1, -dcP2_dX1_T);
	WorkMat.PutT(28,  4, -dcP2_dg1_T);
	WorkMat.PutT(28,  7, -dcP2_dXP1_T);
	WorkMat.PutT(28, 10, -dcP2_domega1_T);
	WorkMat.PutT(28, 13, -dcP2_dX2_T);
	WorkMat.PutT(28, 16, -dcP2_dg2_T);
	WorkMat.PutT(28, 19, -dcP2_dXP2_T);
	WorkMat.PutT(28, 22, -dcP2_domega2_T);

	return WorkMatVar;
}

SubVectorHandler&
InlineFriction::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	integer iNumRows, iNumCols;

	InitialWorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	const Vec3& X1 = pNode1->GetXCurr();
	const Mat3x3& R1 = pNode1->GetRCurr();
	const Vec3& XP1 = pNode1->GetVCurr();
	const Vec3& omega1 = pNode1->GetWCurr();
	const Vec3& X2 = pNode2->GetXCurr();
	const Mat3x3& R2 = pNode2->GetRCurr();
	const Vec3& XP2 = pNode2->GetVCurr();
	const Vec3& omega2 = pNode2->GetWCurr();

	const integer iFirstIndexNode1 = pNode1->iGetFirstIndex();
	const integer iFirstIndexNode2 = pNode2->iGetFirstIndex();
	const integer iFirstIndex = iGetFirstIndex();

	for (integer i = 1; i <= 12; ++i) {
		WorkVec.PutRowIndex(i, iFirstIndexNode1 + i);
		WorkVec.PutRowIndex(i + 12, iFirstIndexNode2 + i);
	}

	for (integer i = 1; i <= 4; ++i) {
		WorkVec.PutRowIndex(i + 24, iFirstIndex + i);
	}

	doublereal lambdaP[2];

	for (int i = 1; i <= 2; ++i) {
		lambda[i - 1] = XCurr.dGetCoef(iFirstIndex + i);
		lambdaP[i - 1] = XCurr.dGetCoef(iFirstIndex + i + 2);
	}

	const Vec3 R2o2 = R2 * o2;
	const Vec3 l1 = X2 + R2o2 - X1;

	const Vec3 F1 = R1 * (e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]);
	const Vec3 M1 = l1.Cross(F1);
	const Vec3 FP1 = omega1.Cross(R1 * (e.GetCol(2) * lambda[L1] + e.GetCol(3) * lambda[L2]))
					+ R1 * (e.GetCol(2) * lambdaP[L1] + e.GetCol(3) * lambdaP[L2]);
	const Vec3 MP1 = -F1.Cross(XP2 + omega2.Cross(R2o2) - XP1) + l1.Cross(FP1);
	const Vec3 F2 = -F1;
	const Vec3 M2 = R2o2.Cross(F2);
	const Vec3 FP2 = -FP1;
	const Vec3 MP2 = (omega2.Cross(R2o2)).Cross(F2) + R2o2.Cross(FP2);

	const Vec3 a = R1.MulTV(l1) - o1;
	const Vec3 aP = R1.MulTV(l1.Cross(omega1) + XP2 + omega2.Cross(R2o2) - XP1);

	WorkVec.Put( 1, F1);
	WorkVec.Put( 4, M1);
	WorkVec.Put( 7, FP1);
	WorkVec.Put(10, MP1);
	WorkVec.Put(13, F2);
	WorkVec.Put(16, M2);
	WorkVec.Put(19, FP2);
	WorkVec.Put(22, MP2);

	for (int i = 1; i <= 2; ++i) {
		WorkVec.PutCoef(24 + i, e.GetCol(i + 1).Dot(a));
		WorkVec.PutCoef(26 + i, e.GetCol(i + 1).Dot(aP));
	}

	return WorkVec;
}

doublereal InlineFriction::SlidingVelocity(const Vec3& X1, const Mat3x3& R1, const Vec3& XP1, const Vec3& omega1, const Vec3& X2, const Mat3x3& R2, const Vec3& XP2, const Vec3& omega2) const
{
	return e.GetCol(1).Dot(R1.MulTV(XP2 - (R2 * o2).Cross(omega2) - XP1 + (X2 + R2 * o2 - X1).Cross(omega1)));
}

doublereal InlineFriction::NormalForceMagnitude() const
{
	doublereal lambda_res = 0;

	for (int i = 0; i < 2; ++i)
		lambda_res += std::pow(lambda[i], 2);

	lambda_res = sqrt(lambda_res);

	return lambda_res;
}

doublereal InlineFriction::FrictionCoefficient(doublereal DeltaXP) const
{
	const doublereal mu = (1. + (mus / muc - 1.) * exp(-std::pow(std::abs(DeltaXP) / vs, iv))) * muc;
	return mu;
}

doublereal InlineFriction::FrictionForce(doublereal DeltaXP) const
{
	const doublereal mu = FrictionCoefficient(DeltaXP);
	const doublereal lambda_res = NormalForceMagnitude();
	const doublereal tau = mu * lambda_res * z / delta + kv * DeltaXP;

	return tau;
}

bool inline_friction_set(void)
{
	UserDefinedElemRead *rf = new UDERead<InlineFriction>;

	if (!SetUDE("inline" "friction", rf))
	{
		delete rf;
		return false;
	}

	return true;
}

//#ifndef STATIC_MODULES

extern "C"
{

int module_init(const char *module_name, void *pdm, void *php)
{
	if (!inline_friction_set())
	{
		silent_cerr("inline friction: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

}

//#endif // ! STATIC_MODULE
