/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

// ModuleIMU: begin

class ModuleIMU
: virtual public Elem, public UserDefinedElem {
private:
	// add private data
	const StructNode *m_pNode;
	Vec3 m_tilde_f;
	Mat3x3 m_tilde_Rh;

	Vec3 m_overline_f;

	Vec3 m_overline_Omega;
	Vec3 m_overline_Acceleration;

public:
	ModuleIMU(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleIMU(void);

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
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
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
};

ModuleIMU::ModuleIMU(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pNode(0),
m_tilde_f(::Zero3),
m_tilde_Rh(::Eye3),
m_overline_f(::Zero3),
m_overline_Omega(::Zero3),
m_overline_Acceleration(::Zero3)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	imu							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"All rights reserved							\n"
"									\n"
"Syntax:								\n"
"	user defined : <label> , imu ,					\n"
"		<node_label>						\n"
"		[ , position , (Vec3)<offset> ]				\n"
"		[ , orientation , (OrientationMatrix)<orientation> ]	\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	if (!m_pNode->bComputeAccelerations()) {
		const_cast<StructNode *>(m_pNode)->ComputeAccelerations(true);
	}

	ReferenceFrame RF(m_pNode);
	if (HP.IsKeyWord("position")) {
		m_tilde_f = HP.GetPosRel(RF);
	}

	if (HP.IsKeyWord("orientation")) {
		m_tilde_Rh = HP.GetRotRel(RF);
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	m_overline_f = m_tilde_Rh.MulTV(m_tilde_f);
}

ModuleIMU::~ModuleIMU(void)
{
	NO_OP;
}

void
ModuleIMU::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel()
			<< " " << m_overline_Omega
			<< " " << m_overline_Acceleration
			<< std::endl;
	}
}

void
ModuleIMU::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleIMU::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleIMU::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);

	Mat3x3 R(m_pNode->GetRCurr()*m_tilde_Rh);
	m_overline_Omega = R.MulTV(m_pNode->GetWCurr());

	m_overline_Acceleration = R.MulTV(m_pNode->GetXPPCurr());
	if (!m_overline_f.IsNull()) {
		Vec3 OmegaP(R.MulTV(m_pNode->GetWPCurr()));
		m_overline_Acceleration += OmegaP.Cross(m_overline_f);
		m_overline_Acceleration += m_overline_Omega.Cross(m_overline_Omega.Cross(m_overline_f));
	}

	return WorkVec;
}

unsigned int
ModuleIMU::iGetNumPrivData(void) const
{
	return 3 + 3;
}

unsigned int
ModuleIMU::iGetPrivDataIdx(const char *s) const
{
	unsigned idx = 0;

	switch (s[0]) {
	case 'w':
		break;

	case 'a':
		idx += 3;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		idx += 1;
		break;

	case 'y':
		idx += 2;
		break;

	case 'z':
		idx += 3;
		break;

	default:
		return 0;
	}

	return idx;
}

doublereal
ModuleIMU::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 1 && i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
		return m_overline_Omega(i);

	case 4:
	case 5:
	case 6:
		return m_overline_Acceleration(i - 3);
	}

	return 0.;
}

int
ModuleIMU::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
ModuleIMU::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = m_pNode;
}

void
ModuleIMU::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleIMU::Restart(std::ostream& out) const
{
	return out << "# ModuleIMU: not implemented" << std::endl;
}

unsigned int
ModuleIMU::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleIMU::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleIMU::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleIMU::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

// ModuleIMU: end

// ModuleIMUConstraint: begin

class ModuleIMUConstraint
: virtual public Elem, public UserDefinedElem {
private:
	// add private data
	const StructNode *m_pNode;
	Vec3 m_tilde_f;
	Mat3x3 m_tilde_Rh;

	Vec3 m_overline_f;

	TplDriveOwner<Vec3> m_OmegaDrv;
	TplDriveOwner<Vec3> m_AccelerationDrv;

	Vec3 m_LambdaOmega;
	Vec3 m_LambdaAcceleration;
	Vec3 m_v;
	Vec3 m_vP;

public:
	ModuleIMUConstraint(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleIMUConstraint(void);

	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;

#if 0
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
#endif

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
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
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
};

ModuleIMUConstraint::ModuleIMUConstraint(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pNode(0),
m_tilde_f(::Zero3),
m_tilde_Rh(::Eye3),
m_overline_f(::Zero3),
m_LambdaOmega(::Zero3),
m_LambdaAcceleration(::Zero3),
m_v(::Zero3),
m_vP(::Zero3)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	imu							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"All rights reserved							\n"
"									\n"
"Syntax:								\n"
"	user defined : <label> , imu constraint,			\n"
"		<node_label> ,						\n"
"		[ position , (Vec3)<offset> , ]				\n"
"		[ orientation , (OrientationMatrix)<orientation> , ]	\n"
"		(TplDriveCaller<Vec3>) <omega> ,			\n"
"		(TplDriveCaller<Vec3>) <acceleration> 			\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	flag fOut = m_pNode->fToBeOutput();
	if (!(fOut & StructNode::OUTPUT_ACCELERATIONS)) {
		const_cast<StructNode *>(m_pNode)->SetOutputFlag(fOut | StructNode::OUTPUT_ACCELERATIONS);
	}

	ReferenceFrame RF(m_pNode);
	if (HP.IsKeyWord("position")) {
		m_tilde_f = HP.GetPosRel(RF);
	}

	if (HP.IsKeyWord("orientation")) {
		m_tilde_Rh = HP.GetRotRel(RF);
	}

	ReferenceFrame RFh(0, ::Zero3,
		m_pNode->GetRCurr()*m_tilde_Rh,
		::Zero3, ::Zero3,
		pDM->GetOrientationDescription());

	TplDriveCaller<Vec3> *pOmegaDC = ReadDCVecRel(pDM, HP, RFh);
	m_OmegaDrv.Set(pOmegaDC);

	TplDriveCaller<Vec3> *pAccelerationDC = ReadDCVecRel(pDM, HP, RFh);
	m_AccelerationDrv.Set(pAccelerationDC);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	m_overline_f = m_tilde_Rh.MulTV(m_tilde_f);
}

ModuleIMUConstraint::~ModuleIMUConstraint(void)
{
	NO_OP;
}

unsigned int
ModuleIMUConstraint::iGetNumDof(void) const
{
	return 9;
}

DofOrder::Order
ModuleIMUConstraint::GetDofType(unsigned int i) const
{
	if (i >= 0 && i < 6) {
		return DofOrder::ALGEBRAIC;
	}

	if (i >= 6 && i < 9) {
		return DofOrder::DIFFERENTIAL;
	}

	ASSERT(0);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
ModuleIMUConstraint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel()
			<< " " << m_LambdaOmega
			<< " " << m_LambdaAcceleration
			<< " " << m_v
			<< " " << m_vP
			<< std::endl;
	}
}

void
ModuleIMUConstraint::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6 + 9;
	*piNumCols = 6 + 9;
}

VariableSubMatrixHandler& 
ModuleIMUConstraint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.ResizeReset(15, 15);

	integer iFirstPositionIndex = m_pNode->iGetFirstPositionIndex();
	integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex();
	integer iFirstIndex = iGetFirstIndex();

	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
	}

	for (integer iCnt = 1; iCnt <= 9; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstIndex + iCnt);
	}

	Mat3x3 R(m_pNode->GetRCurr()*m_tilde_Rh);
	Vec3 f(m_pNode->GetRCurr()*m_tilde_f);

	Vec3 F(R*(m_LambdaAcceleration*dCoef));
	Vec3 M(R*(m_LambdaOmega*dCoef) + f.Cross(F));

	WM.Sub(1, 4, Mat3x3(MatCross, F));
	WM.Sub(4, 4, Mat3x3(MatCross, M));

	WM.Add(1, 7, R);
	WM.Add(4, 10, R);

	WM.AddT(7, 1, R);
	WM.AddT(10, 4, R);

	Mat3x3 Tmp(f.Cross(R));
	WM.Add(4, 7, Tmp);
	Tmp -= (m_pNode->GetVCurr()*dCoef).Cross(R);
	WM.AddT(7, 4, Tmp);

	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutCoef(6 + iCnt, 12 + iCnt, -dCoef);
		WM.PutCoef(12 + iCnt, 12 + iCnt, 1.);
	}

	return WorkMat;
}

SubVectorHandler& 
ModuleIMUConstraint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(15);

	integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex();
	integer iFirstIndex = iGetFirstIndex();

	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
	}

	for (integer iCnt = 1; iCnt <= 9; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstIndex + iCnt);
	}

	m_LambdaAcceleration = Vec3(XCurr, iFirstIndex + 1);
	m_LambdaOmega = Vec3(XCurr, iFirstIndex + 4);
	m_v = Vec3(XCurr, iFirstIndex + 7);
	m_vP = Vec3(XPrimeCurr, iFirstIndex + 7);

	Mat3x3 R(m_pNode->GetRCurr()*m_tilde_Rh);
	Vec3 f(m_pNode->GetRCurr()*m_tilde_f);

	Vec3 F(R*m_LambdaAcceleration);
	Vec3 M(R*m_LambdaOmega + f.Cross(F));

	Vec3 Omega(R.MulTV(m_pNode->GetWCurr()));

	WorkVec.Sub(1, F);
	WorkVec.Sub(4, M);
	WorkVec.Add(7, m_v - R.MulTV(m_pNode->GetVCurr()) - Omega.Cross(m_overline_f));
	WorkVec.Add(10, m_OmegaDrv.Get() - Omega);
	WorkVec.Add(13, m_AccelerationDrv.Get() - m_vP);

	return WorkVec;
}

unsigned int
ModuleIMUConstraint::iGetNumPrivData(void) const
{
	return 0;
}

unsigned int
ModuleIMUConstraint::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

doublereal
ModuleIMUConstraint::dGetPrivData(unsigned int i) const
{
	return 0.;
}

int
ModuleIMUConstraint::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
ModuleIMUConstraint::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = m_pNode;
}

void
ModuleIMUConstraint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	Mat3x3 R(m_pNode->GetRCurr()*m_tilde_Rh);
	Vec3 f(m_pNode->GetRCurr()*m_tilde_f);

	m_v = R.MulTV(m_pNode->GetVCurr() + m_pNode->GetWCurr().Cross(f));
	X.Put(iGetFirstIndex() + 7, m_v);
}

std::ostream&
ModuleIMUConstraint::Restart(std::ostream& out) const
{
	return out << "# ModuleIMUConstraint: not implemented" << std::endl;
}

unsigned int
ModuleIMUConstraint::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleIMUConstraint::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleIMUConstraint::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleIMUConstraint::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

// ModuleIMUConstraint: end

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<ModuleIMU>;

	if (!SetUDE("imu", rf1)) {
		delete rf1;

		silent_cerr("ModuleIMU: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rf2 = new UDERead<ModuleIMUConstraint>;

	if (!SetUDE("imu" "constraint", rf2)) {
		delete rf2;

		silent_cerr("ModuleIMUConstraint: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

