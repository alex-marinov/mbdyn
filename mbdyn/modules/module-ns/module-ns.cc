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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

class ModuleNS1
: virtual public Elem, public UserDefinedElem {
private:
	const StructNode *m_pNode;
	const Vec3 m_X0;
	const Mat3x3 m_R0;
	const Vec3 m_N;
	const doublereal m_dE;

public:
	ModuleNS1(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleNS1(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
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
};

ModuleNS1::ModuleNS1(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_dE(0.)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	ns1							\n"
"Author: 	Matteo Fancello <>		\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_pNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	// assume that normal is direction 3 of m_R0
	const_cast<Vec3&>(m_X0) = HP.GetPosAbs(AbsRefFrame);
	const_cast<Mat3x3&>(m_R0) = HP.GetRotAbs(AbsRefFrame);
	const_cast<Vec3&>(m_N) = m_R0.GetVec(3);

	const_cast<doublereal&>(m_dE) = HP.GetReal();
}

ModuleNS1::~ModuleNS1(void)
{
	// destroy private data
	NO_OP;
}

void
ModuleNS1::Output(OutputHandler& OH) const
{
	// should do something useful
	NO_OP;
}

void
ModuleNS1::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 1;
}

void
ModuleNS1::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	// called right after prediction, before any iteration within the time step
}

void
ModuleNS1::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	// called after a step converged
}

VariableSubMatrixHandler& 
ModuleNS1::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleNS1::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkVec.ResizeReset(3);

	integer iGetFirstIndex = m_pNode->iGetFirstMomentumIndex();
	const Vec3& X(m_pNode->GetXCurr());
	Vec3 DX = X - m_X0;

	doublereal dZ = m_N*DX;
	doublereal dZP = m_N*m_pNode->GetVCurr();
	// dZ and dZP are positive along m_N

	doublereal dF = 0.;
	// compute force and assign to dF
	// force is positive along m_N

	WorkVec.PutItem(1, iGetFirstIndex + 1, m_N(1)*dF);
	WorkVec.PutItem(2, iGetFirstIndex + 2, m_N(2)*dF);
	WorkVec.PutItem(3, iGetFirstIndex + 3, m_N(3)*dF);

	return WorkVec;
}

unsigned int
ModuleNS1::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleNS1::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ModuleNS1::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ModuleNS1::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleNS1::Restart(std::ostream& out) const
{
	return out << "# ModuleNS1: not implemented" << std::endl;
}

unsigned int
ModuleNS1::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleNS1::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleNS1::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleNS1::InitialAssRes(
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
	UserDefinedElemRead *rf = new UDERead<ModuleNS1>;

	if (!SetUDE("ns1", rf)) {
		delete rf;

		silent_cerr("module-ns1: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

