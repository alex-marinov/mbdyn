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

/*
 * Copyright (C) 2011-2015
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * All rights reserved.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstdlib>
#include <cerrno>
#include <cfloat>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "dataman.h"
#include "userelem.h"
#include "rotor.h"
#include "drive.h"

class RotorLooseCoupling
: virtual public Elem,
	public UserDefinedElem,
	public UniformRotor
{
private:
	const DataManager* m_pDM;
	integer m_nCols;

	struct NodalForce {
		const StructNode *pNode;
		Vec3 F;
		Vec3 M;
	};

	unsigned m_stepsPerRevolution;
	std::vector<NodalForce>::iterator m_i_curr;
	std::deque<std::vector<NodalForce> *> m_qF;
	bool m_bFirstAssembly;

	DriveOwner m_Record;
	DriveOwner m_Apply;

	bool bRecord(void) const {
		return (m_Record.dGet() != 0.);
	}

	bool bApply(void) const {
		return (m_Apply.dGet() != 0.);
	}

public:
	RotorLooseCoupling(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~RotorLooseCoupling(void);

	Elem::Type GetElemType(void) const;
	InducedVelocity::Type GetInducedVelocityType(void) const;

	virtual void
	SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints* pHints);

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	virtual void
	AddForce(const Elem *pEl,
		const StructNode *pNode,
		const Vec3& F,
		const Vec3& M,
		const Vec3& X);

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

RotorLooseCoupling::RotorLooseCoupling(unsigned uLabel,
	const DofOwner *pDO,
	DataManager* pDM,
	MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
UniformRotor(uLabel, pDO, 0, Eye3, 0, 0, 0, -1., -1., 0, 0, 0., 0., 1., 1., flag(0)),
m_pDM(pDM),
m_nCols(0),
m_bFirstAssembly(true)
{
	if (HP.IsKeyWord("nodes")) {
		integer nNodes = HP.GetInt();
		if (nNodes < 0) {
			// error
		}
		m_nCols = 6*nNodes;
	}

	if (!HP.IsKeyWord("steps" "per" "revolution")) {
		// error
	}
	integer iSPR = HP.GetInt();
	if (iSPR <= 1) {
		// error
	}
	m_stepsPerRevolution = unsigned(iSPR);

	if (HP.IsKeyWord("record")) {
		m_Record.Set(HP.GetDriveCaller());

	} else {
		m_Record.Set(new OneDriveCaller);
	}

	if (HP.IsKeyWord("apply")) {
		m_Apply.Set(HP.GetDriveCaller());

	} else {
		m_Apply.Set(new OneDriveCaller);
	}

	m_qF.push_back(new std::vector<NodalForce>);

	// read rotor data
}

RotorLooseCoupling::~RotorLooseCoupling(void)
{
	while (!m_qF.empty()) {
		std::vector<NodalForce> *pV = m_qF.back();
		delete pV;
		m_qF.pop_back();
	}
}

Elem::Type
RotorLooseCoupling::GetElemType(void) const
{
	return Elem::LOADABLE;
}

InducedVelocity::Type
RotorLooseCoupling::GetInducedVelocityType(void) const
{
	return InducedVelocity::USER_DEFINED;
}

void
RotorLooseCoupling::SetValue(DataManager *pDM,
	VectorHandler& X,
	VectorHandler& XP,
	SimulationEntity::Hints* pHints)
{
	NO_OP;
}

void
RotorLooseCoupling::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	// pessimistic
	*piNumRows = m_nCols;
	*piNumCols = 1;
}

SubVectorHandler&
RotorLooseCoupling::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	UniformRotor::AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
	if (m_bFirstAssembly) {
		m_bFirstAssembly = false;
	}

	if (bApply() && m_qF.size() == m_stepsPerRevolution) {
		std::vector<NodalForce>* pF = m_qF.back();
		unsigned nNodes = pF->size();
		WorkVec.ResizeReset(6*nNodes);
		for (unsigned iNode = 0; iNode < nNodes; iNode++) {
			const NodalForce& nf = (*pF)[iNode];
			integer iFirstIndex = nf.pNode->iGetFirstMomentumIndex();
			for (integer iCnt = 1; iCnt <= 6; iCnt++) {
				WorkVec.PutRowIndex(6*iNode + iCnt, iFirstIndex + iCnt);
			}
			WorkVec.Put(6*iNode + 1, nf.F);
			WorkVec.Put(6*iNode + 4, nf.M);
		}

	} else {
		WorkVec.Resize(0);
	}

	m_i_curr = m_qF[0]->begin();

	return WorkVec;
}

void
RotorLooseCoupling::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	if (m_qF.size() < m_stepsPerRevolution) {
		m_qF.push_front(new std::vector<NodalForce>[m_qF[0]->size()]);

	} else {
		m_qF.push_front(m_qF.back()); 
		m_qF.pop_back();
	}

	m_i_curr = m_qF[0]->begin();
}

void
RotorLooseCoupling::AddForce(const Elem *pEl,
	const StructNode *pNode,
	const Vec3& F,
	const Vec3& M,
	const Vec3& X)
{
	if (m_bFirstAssembly) {
		NodalForce n;
		n.pNode = pNode;
		m_qF[0]->push_back(n);
		m_i_curr++;
	}

	if (bRecord()) {
		ASSERT(m_i_curr->pNode == pNode);
		m_i_curr->F = F;
		m_i_curr->M = M;
	}
}

unsigned int
RotorLooseCoupling::iGetInitialNumDof(void) const
{
	return 0;
}

void 
RotorLooseCoupling::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
RotorLooseCoupling::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
RotorLooseCoupling::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool
mb_rotor_loose_coupling_set(void)
{
	UserDefinedElemRead *rf = new UDERead<RotorLooseCoupling>;

	bool b = SetUDE("rotor" "loose" "coupling", rf);
	if (!b) {
		delete rf;
	}

	return b;
}

// #ifdef STATIC_MODULES, the function is registered by InitUDE()
#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!mb_rotor_loose_coupling_set()) {
		silent_cerr("module-rotor-loose-coupling: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}
#endif // ! STATIC_MODULES
