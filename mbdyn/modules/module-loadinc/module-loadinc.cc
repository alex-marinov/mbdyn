/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
 * Authors:	Yang Ding <dingyang@gatech.edu>
 *		Pierangelo Masarati <masarati@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "solver.h"
#include "dataman.h"
#include "userelem.h"
#include "driven.h"
#include "Rot.hh"

class LoadIncNorm
: virtual public Elem, public UserDefinedElem {
private:
	int m_FirstSteps;
	doublereal m_dP;
	doublereal m_dPPrev;
	DriveOwner m_DeltaS;
	doublereal m_dDeltaS;
	doublereal m_dS;


	// stop at or past max load
	doublereal m_dPMax;

	// TODO: implement other strategies?

	struct NodeData {
		const StructNode *pNode;
		Vec3 DX;
		Vec3 DTheta;
	};
	std::vector<NodeData> m_Nodes;
	doublereal m_dCompliance;
	doublereal m_dRefLen;
	integer m_iDofOffset;
	doublereal m_dDeltaP;

	doublereal m_DeltaS2;

public:
	LoadIncNorm(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~LoadIncNorm(void);

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
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);

	// get private variable
	doublereal dGetP(void) const;
};

LoadIncNorm::LoadIncNorm(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_FirstSteps(2),
m_dS(0.),
m_dPMax(1.),
m_dCompliance(1.),
m_dRefLen(1.),
m_iDofOffset(6)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        loadinc - load increment normalization\n"
"Author:        Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:  Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               http://www.aero.polimi.it/\n"
"\n"
"               All rights reserved\n"
"\n"
"Usage:\n"
"\n"
"user defined : <label> , load increment normalization ,\n"
"        [ max load , <max_load> , ]          # bails out when p >= max_load\n"
"        [ compliance , <compliance> , ]      # compliance*p = length\n"
"        [ reference length, <ref_length> , ] # multiplies DeltaTheta\n"
"        (DriveCaller) <DeltaS> ;             # arc length increment\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("max" "load")) {
		m_dPMax = HP.GetReal();
		if (m_dPMax <= 0.) {
			silent_cerr("LoadIncNorm(" << uLabel << "): invalid \"max load\" at line " << HP.GetLineData() << std::endl);
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("compliance")) {
		m_dCompliance = HP.GetReal();
		if (m_dCompliance <= std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("LoadIncNorm(" << uLabel << "): invalid \"compliance\" at line " << HP.GetLineData() << std::endl);
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("reference" "length")) {
		m_dRefLen = HP.GetReal();
		if (m_dRefLen < 0.) {
			silent_cerr("LoadIncNorm(" << uLabel << "): invalid \"reference length\" at line " << HP.GetLineData() << std::endl);
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}

		if (m_dRefLen == 0.) {
			m_iDofOffset = 3;
		}
	}

	DriveCaller *pDC = HP.GetDriveCaller();
	m_DeltaS.Set(pDC);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	// initialize data structures
	unsigned uCnt;
	DataManager::NodeContainerType::const_iterator i;
	for (i = pDM->begin(Node::STRUCTURAL), uCnt = 0; i != pDM->end(Node::STRUCTURAL); ++i) {
		StructNode *pNode(dynamic_cast<StructNode *>(i->second));
		ASSERT(pNode != 0);
		if (pNode->GetStructNodeType() == StructNode::DUMMY) {
			continue;
		}
		uCnt++;
	}

	m_Nodes.resize(uCnt);

	for (i = pDM->begin(Node::STRUCTURAL), uCnt = 0; i != pDM->end(Node::STRUCTURAL); ++i) {
		StructNode *pNode(dynamic_cast<StructNode *>(i->second));
		ASSERT(pNode != 0);
		if (pNode->GetStructNodeType() == StructNode::DUMMY) {
			continue;
		}
		m_Nodes[uCnt].pNode = pNode;
		uCnt++;
	}

	m_dDeltaS = m_DeltaS.dGet();
	if (m_dDeltaS < 0.) {
		silent_cerr("LoadIncNorm(" << uLabel << "): DeltaS must be positive" << std::endl);
		throw NoErr(MBDYN_EXCEPT_ARGS);
	}
	m_dDeltaP = m_dDeltaS/m_dCompliance;
	m_dPPrev = -m_dDeltaP;

	m_DeltaS2 = m_dDeltaS;

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_dP=" << m_dP << " m_dDeltaP=" << m_dDeltaP << " m_dPPrev=" << m_dPPrev << std::endl;
}

LoadIncNorm::~LoadIncNorm(void)
{
	NO_OP;
}

void
LoadIncNorm::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << GetLabel()
			<< " " << m_dP
			<< " " << m_dS
			<< std::endl;
	}
}

void
LoadIncNorm::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	NO_OP;
}

void
LoadIncNorm::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	if (m_dP >= m_dPMax) {
		mbdyn_set_stop_at_end_of_time_step();
	}

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_FirstSteps=" << m_FirstSteps << std::endl;

	if (m_FirstSteps) {
		m_FirstSteps--;
	}

	m_dS += m_DeltaS.dGet();
	m_dPPrev = m_dP;

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_dP=" << m_dP << " m_dDeltaP=" << m_dDeltaP << " m_dPPrev=" << m_dPPrev << std::endl;
}

void
LoadIncNorm::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 1;
	*piNumCols = m_iDofOffset*m_Nodes.size() + 1;
}

VariableSubMatrixHandler& 
LoadIncNorm::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iIndex = iGetFirstIndex() + 1;

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_FirstSteps=" << m_FirstSteps << std::endl;

	if (m_FirstSteps || m_dDeltaS == 0.) {
		WM.ResizeReset(1, 1);
		WM.PutRowIndex(1, iIndex);
		WM.PutColIndex(1, iIndex);
		WM.PutCoef(1, 1, m_dCompliance);

	} else {
		integer iNumRows, iNumCols;
		WorkSpaceDim(&iNumRows, &iNumCols);
		WM.ResizeReset(iNumRows, iNumCols);

		WM.PutRowIndex(1, iIndex);
		WM.PutColIndex(iNumCols, iIndex);
		doublereal d = (m_dDeltaP*m_dCompliance*m_dCompliance)/m_DeltaS2;
		if (std::abs(d) <= 10.*std::numeric_limits<doublereal>::epsilon()) {
			d = m_DeltaS.dGet()*m_dCompliance;
		}
		WM.PutCoef(1, iNumCols, d);


		// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_dP=" << m_dP << " m_dDeltaP=" << m_dDeltaP << " m_dPPrev=" << m_dPPrev << std::endl;

		doublereal dCoefX = dCoef/m_DeltaS2;
		doublereal dCoefTheta = dCoefX*m_dRefLen*m_dRefLen;

		integer iNodeOffset = 0;
		for (std::vector<NodeData>::const_iterator i = m_Nodes.begin(); i != m_Nodes.end(); ++i) {
			integer iFirstPositionIndex = i->pNode->iGetFirstIndex();

			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.PutColIndex(iNodeOffset + iCnt, iFirstPositionIndex + iCnt);
				WM.PutCoef(1, iNodeOffset + iCnt, dCoefX*i->DX(iCnt));
			}

			if (m_dRefLen > 0.) {
				// FIXME: check linearization
				Mat3x3 Gm1 = RotManip::DRot_I(i->DTheta);
				Vec3 VTmp = Gm1.MulTV(i->DTheta);

				for (integer iCnt = 1; iCnt <= 3; iCnt++) {
					WM.PutColIndex(iNodeOffset + 3 + iCnt, iFirstPositionIndex + 3 + iCnt);
					WM.PutCoef(1, iNodeOffset + 3 + iCnt, dCoefTheta*VTmp(iCnt));
				}
			}

			iNodeOffset += m_iDofOffset;
		}
	}

	return WorkMat;
}

SubVectorHandler& 
LoadIncNorm::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iIndex = iGetFirstIndex() + 1;

	m_dP = XCurr(iIndex);
	m_dDeltaP = m_dP - m_dPPrev;

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_dP=" << m_dP << " m_dDeltaP=" << m_dDeltaP << " m_dPPrev=" << m_dPPrev << std::endl;

	m_dDeltaS = m_DeltaS.dGet();
	if (m_dDeltaS < 0.) {
		silent_cerr("LoadIncNorm(" << uLabel << ")::AssRes(): DeltaS must be positive" << std::endl);
		throw NoErr(MBDYN_EXCEPT_ARGS);
	}
	doublereal d;

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " m_FirstSteps=" << m_FirstSteps << std::endl;

	if (m_FirstSteps == 2) {
		d = -m_dP*m_dCompliance;

	} else if (m_FirstSteps == 1 || m_dDeltaS == 0.) {
		d = m_dDeltaS - m_dP*m_dCompliance;

	} else {
		d = 0.;

		for (std::vector<NodeData>::iterator i = m_Nodes.begin(); i != m_Nodes.end(); ++i) {
			i->DX = i->pNode->GetXCurr() - i->pNode->GetXPrev();
			d += i->DX.Dot();

			if (m_dRefLen > 0.) {
				i->DTheta = RotManip::VecRot(i->pNode->GetRCurr().MulMT(i->pNode->GetRPrev()));
				d += (m_dRefLen*m_dRefLen)*i->DTheta.Dot();
			}
		}

		d += (m_dDeltaP*m_dCompliance)*(m_dDeltaP*m_dCompliance);
		m_DeltaS2 = std::sqrt(d);

		d = m_dDeltaS - m_DeltaS2;
	}

	WorkVec.PutItem(1, iIndex, d);

	return WorkVec;
}

unsigned int
LoadIncNorm::iGetNumPrivData(void) const
{
	return 3;
}

unsigned int
LoadIncNorm::iGetPrivDataIdx(const char *s) const
{
	if (strcmp(s, "p") == 0) {
		return 1;
	}

	if (strcmp(s, "S") == 0) {
		return 2;
	}

	if (strcmp(s, "DeltaS") == 0) {
		return 3;
	}

	return 0;
}

doublereal
LoadIncNorm::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
		return m_dP;

	case 2:
		return m_dS;

	case 3:
		return m_dDeltaS;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

int
LoadIncNorm::iGetNumConnectedNodes(void) const
{
	return m_Nodes.size();
}

void
LoadIncNorm::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(m_Nodes.size());

	for (unsigned uCnt = 0; uCnt < m_Nodes.size(); uCnt++) {
		connectedNodes[uCnt] = m_Nodes[uCnt].pNode;
	}
}

void
LoadIncNorm::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iIndex = iGetFirstIndex() + 1;
	X(iIndex) = 0.;
	XP(iIndex) = m_dPPrev;
}

unsigned int
LoadIncNorm::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
LoadIncNorm::GetDofType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::ALGEBRAIC;
}

DofOrder::Order
LoadIncNorm::GetEqType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::ALGEBRAIC;
}

std::ostream&
LoadIncNorm::Restart(std::ostream& out) const
{
	return out << "# LoadIncNorm: not implemented" << std::endl;
}

unsigned int
LoadIncNorm::iGetInitialNumDof(void) const
{
	return 0;
}

void 
LoadIncNorm::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
LoadIncNorm::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
LoadIncNorm::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

doublereal
LoadIncNorm::dGetP(void) const
{
	return m_dP;
}

class LoadIncForce
: virtual public Elem, public UserDefinedElem {
private:
	bool m_bCouple;
	bool m_bFollower;

	LoadIncNorm *m_pLoadIncNorm;
	DrivenElem *m_pDrivenLoadIncNorm;
	StructNode *m_pNode;
	Vec3 m_b;
	Vec3 m_Dir;
	Vec3 m_F;
	Vec3 m_M;

public:
	LoadIncForce(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~LoadIncForce(void);

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
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
};

LoadIncForce::LoadIncForce(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pDrivenLoadIncNorm(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        loadinc - load increment force\n"
"Author:        Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:  Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               http://www.aero.polimi.it/\n"
"\n"
"               All rights reserved\n"
"\n"
"Usage:\n"
"\n"
"user defined : <label> , load increment force ,\n"
"        { force | couple } ,\n"
"        { absolute | follower } ,\n"
"        <node_label> ,\n"
"        [ position , (Vec3) <position> , ] # meaningless for couple\n"
"        (Vec3) <direction> , \n"
"        load increment normalization , <lin_label> ;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("force")) {
		m_bCouple = false;

	} else if (HP.IsKeyWord("couple")) {
		m_bCouple = true;

	} else {
		silent_cerr("LoadIncForce(" << uLabel << "): missing { force | couple } or unexpected type at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("absolute")) {
		m_bFollower = false;

	} else if (HP.IsKeyWord("follower")) {
		m_bFollower = true;

	} else {
		silent_cerr("LoadIncForce(" << uLabel << "): missing { absolute | follower } or unexpected type at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pNode = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));
	ASSERT(m_pNode != 0);

	ReferenceFrame rf(m_pNode);

	m_b = Zero3;
	if (HP.IsKeyWord("position")) {
		m_b = HP.GetPosRel(rf);
	}

	if (m_bFollower) {
		m_Dir = HP.GetVecRel(rf);

	} else {
		m_Dir = HP.GetVecAbs(rf);
	}

	if (!HP.IsKeyWord("load" "increment" "normalization")) {
		silent_cerr("LoadIncForce(" << uLabel << "): missing \"load increment normalization\" keyword at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// if m_pLoadIncNorm is driven, make sure the LoadIncForce
	// and the LoadIncNorm are simultaneously active
	Elem *pEl = pDM->ReadElem(HP, Elem::LOADABLE);
	m_pLoadIncNorm = dynamic_cast<LoadIncNorm *>(pEl);
	if (m_pLoadIncNorm == 0) {
		m_pDrivenLoadIncNorm = dynamic_cast<DrivenElem *>(pEl);
		if (m_pDrivenLoadIncNorm == 0) {
			silent_cerr("LoadIncForce(" << uLabel << "): invalid \"load increment normalization\" UseDefined(" << pEl->GetLabel() << ") element at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		m_pLoadIncNorm = dynamic_cast<LoadIncNorm *>(m_pDrivenLoadIncNorm->pGetElem());
		if (m_pLoadIncNorm == 0) {
			silent_cerr("LoadIncForce(" << uLabel << "): invalid \"load increment normalization\" UseDefined(" << pEl->GetLabel() << ") element at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

LoadIncForce::~LoadIncForce(void)
{
	NO_OP;
}

void
LoadIncForce::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << GetLabel() << "@" << m_pNode->GetLabel();
		if (!m_bCouple) {
			out << " " << m_F;
		}
		out << " " << m_M << std::endl;
	}
}

void
LoadIncForce::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	NO_OP;
}

void
LoadIncForce::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	if (m_bCouple) {
		*piNumRows = 3;
	} else {
		*piNumRows = 6;
	}

	if (m_bFollower || !m_bCouple) {
		*piNumCols = 3 + 1;

	} else {
		*piNumCols = 1;
	}
}

VariableSubMatrixHandler& 
LoadIncForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// if m_pLoadIncNorm is driven, make sure the LoadIncForce
	// and the LoadIncNorm are simultaneously active
	if (m_pDrivenLoadIncNorm && !m_pDrivenLoadIncNorm->bIsActive()) {
		silent_cerr("LoadIncForce(" << GetLabel() << "): LoadIncNorm(" << m_pLoadIncNorm->GetLabel() << ") inactive" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	// std::cerr << "### " << __PRETTY_FUNCTION__ << " iNumRows=" << iNumRows << " iNumCols=" << iNumCols << std::endl;

	if (m_bCouple) {
		integer iFirstPositionIndex = m_pNode->iGetFirstPositionIndex() + 3;
		integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex() + 3;
		integer iNormIndex = m_pLoadIncNorm->iGetFirstIndex() + 1;

		for (integer iCnt = 1; iCnt <= 3; iCnt++) {
			WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		}

		if (m_bFollower) {
			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
			}
			WM.PutColIndex(4, iNormIndex);

			WM.Add(1, 1, m_M*dCoef);

			const Mat3x3& R(m_pNode->GetRRef());
			Vec3 TmpDir(R*m_Dir);
			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.DecCoef(iCnt, 4, TmpDir(iCnt));
			}

		} else {
			WM.PutColIndex(1, iNormIndex);

			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.DecCoef(iCnt, 1, m_Dir(iCnt));
			}
		}

	} else {
		integer iFirstPositionIndex = m_pNode->iGetFirstPositionIndex() + 3;
		integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex();
		integer iNormIndex = m_pLoadIncNorm->iGetFirstIndex() + 1;

		for (integer iCnt = 1; iCnt <= 3; iCnt++) {
			WM.PutColIndex(iCnt, iFirstPositionIndex + iCnt);
		}
		WM.PutColIndex(4, iNormIndex);

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		}

		const Mat3x3& R(m_pNode->GetRRef());

		if (m_bFollower) {
			WM.Add(1, 1, m_F*dCoef);
			WM.Add(3 + 1, 1, m_M*dCoef);

			Vec3 TmpDir(R*m_Dir);
			Vec3 TmpM(R*m_b.Cross(m_Dir));
			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.DecCoef(iCnt, 4, TmpDir(iCnt));
				WM.DecCoef(3 + iCnt, 4, TmpM(iCnt));
			}

		} else {
			Vec3 TmpArm(R*m_b);

			WM.Sub(3 + 1, 1, Mat3x3(MatCrossCross, m_F, TmpArm*dCoef));

			Vec3 TmpM(TmpArm.Cross(m_Dir));

			// std::cerr << "### " << __PRETTY_FUNCTION__ << " TmpM=" << TmpM << std::endl;

			for (integer iCnt = 1; iCnt <= 3; iCnt++) {
				WM.DecCoef(iCnt, 4, m_Dir(iCnt));
				WM.DecCoef(3 + iCnt, 4, TmpM(iCnt));
			}
		}
	}

	return WorkMat;
}

SubVectorHandler& 
LoadIncForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// if m_pLoadIncNorm is driven, make sure the LoadIncForce
	// and the LoadIncNorm are simultaneously active
	if (m_pDrivenLoadIncNorm && !m_pDrivenLoadIncNorm->bIsActive()) {
		silent_cerr("LoadIncForce(" << GetLabel() << "): LoadIncNorm(" << m_pLoadIncNorm->GetLabel() << ") inactive" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	doublereal dp = m_pLoadIncNorm->dGetP();

	if (m_bCouple) {
		integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex() + 3;
		for (integer iCnt = 1; iCnt <= 3; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		}

		if (m_bFollower) {
			const Mat3x3& R(m_pNode->GetRCurr());
			m_M = R*(m_Dir*dp);

		} else {
			m_M = m_Dir*dp;
		}

		WorkVec.Add(1, m_M);

	} else {
		integer iFirstMomentumIndex = m_pNode->iGetFirstMomentumIndex();
		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
		}

		const Mat3x3& R(m_pNode->GetRCurr());

		if (m_bFollower) {
			Vec3 FTmp(m_Dir*dp);

			m_F = R*FTmp;
			m_M = R*m_b.Cross(FTmp);

		} else {
			m_F = m_Dir*dp;
			m_M = (R*m_b).Cross(m_F);
		}

		WorkVec.Add(1, m_F);
		WorkVec.Add(4, m_M);
	}

	return WorkVec;
}

unsigned int
LoadIncForce::iGetNumPrivData(void) const
{
	return 0;
}

int
LoadIncForce::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
LoadIncForce::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);

	connectedNodes[0] = m_pNode;
}

void
LoadIncForce::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
LoadIncForce::Restart(std::ostream& out) const
{
	return out << "# LoadIncForce: not implemented" << std::endl;
}

unsigned int
LoadIncForce::iGetInitialNumDof(void) const
{
	return 0;
}

void 
LoadIncForce::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
LoadIncForce::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
LoadIncForce::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf;

	rf = new UDERead<LoadIncForce>;
	if (!SetUDE("load" "increment" "force", rf)) {
		delete rf;

		silent_cerr("module-loadinc: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	rf = new UDERead<LoadIncNorm>;
	if (!SetUDE("load" "increment" "normalization", rf)) {
		delete rf;

		silent_cerr("module-loadinc: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

