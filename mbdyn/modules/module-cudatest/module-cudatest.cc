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
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Tingnan Zhang <tingnan1986@gatech.edu>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "dataman.h"
#include "userelem.h"

#include "cudatest.h"

class MBDynCUDATest
: virtual public Elem, public UserDefinedElem {
private:
	const StructDispNode *m_pNode;
	CUDATest *m_pCUDATest;
	Vec3 m_F;

public:
	MBDynCUDATest(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~MBDynCUDATest(void);

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
	virtual void Output(OutputHandler& OH) const;
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

MBDynCUDATest::MBDynCUDATest(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pNode(0),
m_pCUDATest(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        cudatest\n"
"Author:        Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               <http://www.aero.polimi.it/>\n"
"\n"
"               All rights reserved\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	m_pNode = pDM->ReadNode<const StructDispNode, Node::STRUCTURAL>(HP);

	int nPoints = HP.GetInt();
	if (nPoints <= 0) {
		silent_cerr("MBDynCUDATest(" << uLabel << "): "
			"invalid points number " << nPoints
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<double> k(nPoints), r(nPoints);
	for (int i = 0; i < nPoints; i++) {
		k[i] = HP.GetReal();
		r[i] = HP.GetReal();
	}

	m_pCUDATest = (CUDATest *)mbdyn_CUDATest_init(k.size(), &k[0], &r[0]);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

}

MBDynCUDATest::~MBDynCUDATest(void)
{
	delete m_pCUDATest;
}

void
MBDynCUDATest::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << GetLabel() << "@" << m_pNode->GetLabel()
			<< " " << m_F
			<< std::endl;
	}
}

void
MBDynCUDATest::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
MBDynCUDATest::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// no contribution to Jacobian matrix
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBDynCUDATest::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(3);

	integer iFirstIndex = m_pNode->iGetFirstMomentumIndex();
	for (integer iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	m_pCUDATest->GetForce(m_F, m_pNode->GetXCurr(), m_pNode->GetVCurr());

	WorkVec.Add(1, m_F);

	return WorkVec;
}

unsigned int
MBDynCUDATest::iGetNumPrivData(void) const
{
	return 0;
}

int
MBDynCUDATest::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
MBDynCUDATest::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = m_pNode;
}

void
MBDynCUDATest::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	// initialize X and XP according to the initial state of MD as needed
}

std::ostream&
MBDynCUDATest::Restart(std::ostream& out) const
{
	// don't worry about "soft" restart by now
	return out << "# MBDynCUDATest: not implemented" << std::endl;
}

unsigned int
MBDynCUDATest::iGetInitialNumDof(void) const
{
	return 0;
}

void 
MBDynCUDATest::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
MBDynCUDATest::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBDynCUDATest::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<MBDynCUDATest>;

	if (!SetUDE("cuda" "test", rf)) {
		delete rf;

		silent_cerr("module-cudatest: "
			"mbdyn_CUDAtest_set() "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

