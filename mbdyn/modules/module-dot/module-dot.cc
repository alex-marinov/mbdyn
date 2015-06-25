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

class ModuleDOT
: virtual public Elem, public UserDefinedElem {
private:
	const DataManager *m_pDM;
	mutable bool m_bDone;
	std::string m_fname;
	bool m_bStructOnly;

public:
	ModuleDOT(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleDOT(void);

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
};

ModuleDOT::ModuleDOT(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pDM(pDM),
m_bDone(false),
m_bStructOnly(false)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	dot      						\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
"\n"
"user defined : <label> , dot\n"
"        [ , help ]\n"
"        [ , file name , \"filename\" ]\n"
"        [ , structure only , { yes | no | <bool> } ]\n"
";\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// do something useful
	if (HP.IsKeyWord("file" "name")) {
		const char *s = HP.GetFileName();
		m_fname = s;
	}

	if (HP.IsKeyWord("structure" "only")) {
		m_bStructOnly = HP.GetYesNoOrBool();
	}

	if (m_fname.empty()) {
		m_fname = "mbdyn.dot";
	}
}

ModuleDOT::~ModuleDOT(void)
{
	// destroy private data
	NO_OP;
}

void
ModuleDOT::Output(OutputHandler& OH) const
{
	// should do something useful
	if (m_bDone) {
		return;
	}

	std::ofstream out(m_fname.c_str());

	std::string props;

	out << "graph mbdyn {\n"
		"layout=neato;\n";

	// Iteration on elements
	out << "node [shape=box];\n";
	for (int Type = 0; Type < Elem::LASTELEMTYPE; Type++) {
		if (m_bStructOnly) {
			switch (Elem::Type(Type)) {
			case Elem::JOINT:
				props = " [color=red]";
				break;

			case Elem::BEAM:
				props = " [color=green]";
				break;

			case Elem::PLATE:
				props = " [color=blue]";
				break;

			case Elem::BODY:
				props = " [color=cyan]";
				break;

			default:
				continue;
			}
		}

		DataManager::ElemContainerType::const_iterator e = m_pDM->begin(Elem::Type(Type));
		DataManager::ElemContainerType::const_iterator e_last = m_pDM->end(Elem::Type(Type));
		for (; e != e_last; ++e) {
			out << "\"" << psElemNames[e->second->GetElemType()] << "(" << e->second->GetLabel() << ")\"" << props << ";\n";
		}
	}

	// Iteration on nodes
	out << "node [shape=ellipse];\n";
	for (int Type = 0; Type < Node::LASTNODETYPE; Type++) {
		if (m_bStructOnly) {
			switch (Node::Type(Type)) {
			case Node::STRUCTURAL:
				break;

			default:
				continue;
			}
		}

		DataManager::NodeContainerType::const_iterator n = m_pDM->begin(Node::Type(Type));
		DataManager::NodeContainerType::const_iterator n_last = m_pDM->end(Node::Type(Type));
		for (; n != n_last; ++n) {
			out << "\"" << psNodeNames[n->second->GetNodeType()] << "(" << n->second->GetLabel() << ")\";\n";
		}
	}

	// Iteration on elements
	out << "edge [style=solid];\n";
	std::vector<const Node *> connectedNodes;
	for (int Type = 0; Type < Elem::LASTELEMTYPE; Type++) {
		if (m_bStructOnly) {
			switch (Elem::Type(Type)) {
			case Elem::JOINT:
			case Elem::BEAM:
			case Elem::PLATE:
			case Elem::BODY:
				break;

			default:
				continue;
			}
		}

		DataManager::ElemContainerType::const_iterator e = m_pDM->begin(Elem::Type(Type));
		DataManager::ElemContainerType::const_iterator e_last = m_pDM->end(Elem::Type(Type));
		for (; e != e_last; ++e) {
			out << "\"" << psElemNames[e->second->GetElemType()] << "(" << e->second->GetLabel() << ")\" --";

			e->second->GetConnectedNodes(connectedNodes);

			if (connectedNodes.size() == 1) {
				out << " \"" << psNodeNames[connectedNodes[0]->GetNodeType()]
					<< "(" << connectedNodes[0]->GetLabel() << ")\";\n";

			} else if (!connectedNodes.empty()) {
				std::vector<const Node *>::const_iterator n = connectedNodes.begin();
				out << " {\"" << psNodeNames[(*n)->GetNodeType()]
					<< "(" << (*n)->GetLabel() << ")\"";
				for ( ++n; n != connectedNodes.end(); ++n) {
					out << " \"" << psNodeNames[(*n)->GetNodeType()] << "(" << (*n)->GetLabel() << ")\"";
				}
				out << "};\n";
			}
		}
	}

	out << "}";

	m_bDone = true;
}

void
ModuleDOT::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleDOT::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleDOT::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkVec.ResizeReset(0);

	return WorkVec;
}

unsigned int
ModuleDOT::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleDOT::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ModuleDOT::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ModuleDOT::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleDOT::Restart(std::ostream& out) const
{
	return out << "# ModuleDOT: not implemented" << std::endl;
}

unsigned int
ModuleDOT::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleDOT::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleDOT::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleDOT::InitialAssRes(
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
	UserDefinedElemRead *rf = new UDERead<ModuleDOT>;

	if (!SetUDE("dot", rf)) {
		delete rf;

		silent_cerr("module-dot: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

