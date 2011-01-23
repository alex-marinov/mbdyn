/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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

#include "dataman.h"
#include "userelem.h"

class MD
: virtual public Elem, public UserDefinedElem {
private:
	// add private data
	struct NodeData {
		StructNode *pNode;
		Vec3 F;
		Vec3 M;
	};

	int iCoupling;
	int iCouplingCounter;

	std::vector<NodeData> m_nodeData;

public:
	MD(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~MD(void);

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

MD::MD(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
iCoupling(0),
iCouplingCounter(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module:        md\n"
"Author:        Yang Ding <dingyang@gatech.edu>\n"
"               Pierangelo Masarati <masarati@aero.polimi.it>\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               http://www.aero.polimi.it/\n"
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

	// do something useful
	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("tight")) {
			iCoupling = 1;

		} else if (HP.IsKeyWord("loose")) {
			iCoupling = 0;

		} else {
			iCoupling = HP.GetInt();
			if (iCoupling < 0) {
				silent_cerr("MD(" << uLabel << "): invalid coupling "
					" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("MD(" << uLabel << "): invalid node number " << n
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_nodeData.resize(n);
	for (unsigned i = 0; i < unsigned(n); i++) {
		m_nodeData[i].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		ASSERT(m_nodeData[i].pNode != 0);

		for (unsigned c = 0; c < i; c++) {
			if (m_nodeData[c].pNode == m_nodeData[i].pNode) {
				silent_cerr("MD(" << uLabel << "): "
					"StructNode(" << m_nodeData[i].pNode->GetLabel() << "), "
					"#" << i << "/" << n << ", "
					"already provided as #" << c << "/" << n
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

MD::~MD(void)
{
	// destroy private data
	NO_OP;
}

void
MD::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (std::vector<NodeData>::const_iterator n = m_nodeData.begin();
			n != m_nodeData.end(); n++)
		{
			out << GetLabel() << "@" << n->pNode->GetLabel()
				<< " " << n->F
				<< " " << n->M
				<< std::endl;
		}
	}
}

void
MD::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	iCouplingCounter = 0;
}

void
MD::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*m_nodeData.size();
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
MD::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MD::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkVec.ResizeReset(6*m_nodeData.size());

	/*
	 * two possible approaches:
	 *   1) loose coupling: pass predicted kinematics,
	 *	receive forces once, reuse them for all iterations
	 *   2) tight coupling: pass actual kinematics,
	 *      receive forces at each iteration
	 */

	if ((iCoupling == 0 && iCouplingCounter == 0) || (iCouplingCounter%iCoupling) == 0) {
		// get forces from MD
	}
	iCouplingCounter++;

	integer r = 0;
	for (std::vector<NodeData>::const_iterator n = m_nodeData.begin();
		n != m_nodeData.end(); n++)
	{
		integer iFirstIndex = n->pNode->iGetFirstMomentumIndex();

		for (integer iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(r + iCnt, iFirstIndex + iCnt);
		}

		WorkVec.Add(r + 1, n->F);
		WorkVec.Add(r + 4, n->M);

		r += 6;
	}

	return WorkVec;
}

unsigned int
MD::iGetNumPrivData(void) const
{
	return 0;
}

int
MD::iGetNumConnectedNodes(void) const
{
	return m_nodeData.size();
}

void
MD::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(m_nodeData.size());

	for (unsigned n = 0; n < m_nodeData.size(); n++)
	{
		connectedNodes[n] = m_nodeData[n].pNode;
	}
}

void
MD::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
MD::Restart(std::ostream& out) const
{
	return out << "# MD: not implemented" << std::endl;
}

unsigned int
MD::iGetInitialNumDof(void) const
{
	return 0;
}

void 
MD::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
MD::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MD::InitialAssRes(
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
	UserDefinedElemRead *rf = new UDERead<MD>;

	if (!SetUDE("md", rf)) {
		delete rf;

		silent_cerr("module-md: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

