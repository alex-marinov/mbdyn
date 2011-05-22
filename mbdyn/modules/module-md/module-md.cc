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
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Tingnan Zhang <tingnan1986@gatech.edu>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "dataman.h"
#include "userelem.h"

// data structure and container declarations
struct NodeData {
	StructNode *pNode;
	Vec3 F;
	Vec3 M;
};

typedef std::vector<NodeData> NodeContainer;

// examples of MD-specific handlers (they could also be member functions
// of the MBMD module class
static int
MD_init(void)
{
	// do something useful
	return 0;
}

static int
MD_get_loads(NodeContainer& nodeData)
{
	// do something useful
	for (NodeContainer::iterator n = nodeData.begin();
		n != nodeData.end(); ++n)
	{
		std::cout << "Node(" << n->pNode->GetLabel() << ") X={" << n->pNode->GetXCurr() << "} V={" << n->pNode->GetVCurr() << "}" << std::endl;

		n->F = Vec3(0., 0., 0.);
		n->M = Vec3(0., 0., 0.);
	}

	return 0;
}

static void
MD_destroy(void)
{
	return;
}

class MBMD
: virtual public Elem, public UserDefinedElem {
private:
	int iCoupling;
	int iCouplingCounter;

	NodeContainer m_nodeData;

public:
	MBMD(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~MBMD(void);

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
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
};

MBMD::MBMD(
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
"Author:        Pierangelo Masarati <masarati@aero.polimi.it>\n"
"		Tingnan Zhang <tingnan1986@gatech.edu>\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale\n"
"               Politecnico di Milano\n"
"               <http://www.aero.polimi.it/>\n"
"		\"Crab Lab\"\n"
"		Georgia Institute of Technology\n"
"		<http://www.physics.gatech.edu/research/goldman/>\n"
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

	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("tight")) {
			iCoupling = 1;

		} else if (HP.IsKeyWord("loose")) {
			iCoupling = 0;

		} else {
			iCoupling = HP.GetInt();
			if (iCoupling < 0) {
				silent_cerr("MBMD(" << uLabel << "): invalid coupling "
					" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("MBMD(" << uLabel << "): invalid node number " << n
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_nodeData.resize(n);
	for (unsigned i = 0; i < unsigned(n); i++) {
		m_nodeData[i].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

		for (unsigned c = 0; c < i; c++) {
			if (m_nodeData[c].pNode == m_nodeData[i].pNode) {
				silent_cerr("MBMD(" << uLabel << "): "
					"StructNode(" << m_nodeData[i].pNode->GetLabel() << "), "
					"#" << i << "/" << n << ", "
					"already provided as #" << c << "/" << n
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	if (MD_init() != 0) {
		silent_cerr("MBMD(" << uLabel << "): "
			"MD_init() failed"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

MBMD::~MBMD(void)
{
	// destroy private data
	MD_destroy();
}

void
MBMD::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (NodeContainer::const_iterator n = m_nodeData.begin();
			n != m_nodeData.end(); ++n)
		{
			// format:
			// - for each node
			//   - element label "@" node label
			//   - three components of force
			//   - three components of moment
			out << GetLabel() << "@" << n->pNode->GetLabel()
				<< " " << n->F
				<< " " << n->M
				<< std::endl;
		}
	}
}

void
MBMD::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	iCouplingCounter = 0;
}

void
MBMD::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	// 6 rows for each node; 1 column since only residual is assembled
	*piNumRows = 6*m_nodeData.size();
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
MBMD::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// no contribution to Jacobian matrix
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBMD::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(6*m_nodeData.size());

	// two possible approaches:
	//   1) loose coupling: pass predicted kinematics,
	//	receive forces once, reuse them for all iterations
	//   2) tight coupling: pass actual kinematics,
	//      receive forces at each iteration

	if ((iCoupling == 0 && iCouplingCounter == 0)
		|| (iCoupling > 0 && (iCouplingCounter%iCoupling) == 0))
	{
		// get loads from MD
		// this helper fills F and M of each node
		if (MD_get_loads(m_nodeData) != 0) {
			silent_cerr("MBMD(" << uLabel << "): "
				"MD_get_loads() failed"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	iCouplingCounter++;

	integer r = 0;
	for (NodeContainer::const_iterator n = m_nodeData.begin();
		n != m_nodeData.end(); ++n)
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
MBMD::iGetNumPrivData(void) const
{
	return 0;
}

int
MBMD::iGetNumConnectedNodes(void) const
{
	return m_nodeData.size();
}

void
MBMD::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(m_nodeData.size());

	for (unsigned n = 0; n < m_nodeData.size(); n++) {
		connectedNodes[n] = m_nodeData[n].pNode;
	}
}

void
MBMD::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	// initialize X and XP according to the initial state of MD as needed
}

std::ostream&
MBMD::Restart(std::ostream& out) const
{
	// don't worry about "soft" restart by now
	return out << "# MBMD: not implemented" << std::endl;
}

unsigned int
MBMD::iGetInitialNumDof(void) const
{
	return 0;
}

void 
MBMD::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
MBMD::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
MBMD::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<MBMD>;

	if (!SetUDE("md", rf)) {
		delete rf;

		silent_cerr("module-md: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

