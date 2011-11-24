/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

/* node manager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "search.h"

/* DataManager - begin */

const char sNMClassName[] = "DataManager";

void
DataManager::NodeManager(void)
{
	for (int i = 0; i < Node::LASTNODETYPE; i++) {
#if 0
		NodeData[i].ppFirstNode = NULL;
		NodeData[i].iNum = 0;
#endif
		NodeData[i].iExpectedNum = 0;
		NodeData[i].uFlags = 0U;
		NodeData[i].DefaultOut(::fDefaultOut == 1); /* Da "output.h" */
		NodeData[i].OutFile = OutputHandler::UNKNOWN; /* Da "output.h" */
	}

	/* Se un tipo scrive su un file di output, aggiungere qui il tipo di file */
	NodeData[Node::ABSTRACT].OutFile = OutputHandler::ABSTRACT;
	NodeData[Node::ABSTRACT].Desc = "Abstract";
	NodeData[Node::ABSTRACT].ShortDesc = "abstr";

	NodeData[Node::STRUCTURAL].OutFile = OutputHandler::STRNODES;
	NodeData[Node::STRUCTURAL].Desc = "Structural";
	NodeData[Node::STRUCTURAL].ShortDesc = "struct";

	NodeData[Node::ELECTRIC].OutFile = OutputHandler::ELECTRIC;
	NodeData[Node::ELECTRIC].Desc = "Electric";
	NodeData[Node::ELECTRIC].ShortDesc = "elec";

	NodeData[Node::THERMAL].OutFile = OutputHandler::THERMALNODES;   
	NodeData[Node::THERMAL].Desc = "Thermal";
	NodeData[Node::THERMAL].ShortDesc = "thermal";

	NodeData[Node::HYDRAULIC].OutFile = OutputHandler::PRESNODES;
	NodeData[Node::HYDRAULIC].Desc = "Pressure";
	NodeData[Node::HYDRAULIC].ShortDesc = "pres";

	NodeData[Node::PARAMETER].OutFile = OutputHandler::PARAMETERS;
	NodeData[Node::PARAMETER].Desc = "Parameter";
	NodeData[Node::PARAMETER].ShortDesc = "param";
}

void
DataManager::NodeManagerDestructor(void)
{
	DEBUGCOUT("Entering DataManager::NodeManagerDestructor()" << std::endl);

	for (NodeVecType::iterator p = Nodes.begin(); p != Nodes.end(); ++p) {
		DEBUGCOUT("deleting node "
			<< psNodeNames[(*p)->GetNodeType()]
			<< "(" << (*p)->GetLabel() << ")"
			<< std::endl);
		SAFEDELETE(*p);
	}
}

void
DataManager::NodeDataInit(void)
{
	for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		iTotNodes += NodeData[iCnt].iExpectedNum;
	}

	DEBUGCOUT("iTotNodes = " << iTotNodes << std::endl);

	if (iTotNodes > 0) {
		Nodes.resize(iTotNodes);

		NodeIter.Init(&Nodes[0], iTotNodes);

		for (NodeVecType::iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
			*i = 0;
		}

	} else {
		silent_cerr("warning, no nodes are defined" << std::endl);
	}
}

void
DataManager::NodeOutputPrepare(OutputHandler& OH)
{
#ifdef USE_NETCDF
	for (unsigned nt = 0; nt < Node::LASTNODETYPE; nt++) {
		if (!NodeData[nt].NodeContainer.empty() && OH.UseNetCDF(NodeData[nt].OutFile)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			integer iNumNodes = NodeData[nt].NodeContainer.size();

			OutputHandler::AttrValVec attrs(1);
			attrs[0] = OutputHandler::AttrVal("description", std::string(NodeData[nt].Desc) + " nodes labels");

			OutputHandler::NcDimVec dim(1);
			dim[0] = OH.CreateDim(std::string(NodeData[nt].ShortDesc) + "_node_labels_dim", iNumNodes);

			NcVar *VarLabels = OH.CreateVar(std::string("node.") + NodeData[nt].ShortDesc, ncInt, attrs, dim);

			NodeContainerType::const_iterator p = NodeData[nt].NodeContainer.begin();
			for (unsigned i = 0; i < unsigned(iNumNodes); i++, p++) {
				VarLabels->set_cur(i);
				const long l = p->second->GetLabel();
				VarLabels->put(&l, 1);
			}
		}
	}
#endif // USE_NETCDF

	for (NodeVecType::iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->OutputPrepare(OH);
	}
}

void
DataManager::NodeOutput(OutputHandler& OH) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->Output(OH);
	}
}

void
DataManager::NodeOutput(
	OutputHandler& OH,
	const VectorHandler& X,
	const VectorHandler& XP) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->Output(OH, X, XP);
	}
}

flag
DataManager::fGetDefaultOutputFlag(const Node::Type& t) const
{
	return NodeData[t].bDefaultOut();
}


/* cerca un nodo qualsiasi */
Node*
DataManager::pFindNode(Node::Type Typ, unsigned int uL) const
{
	NodeMapToListType::const_iterator p = NodeData[Typ].NodeMapToList.find(uL);
	if (p == NodeData[Typ].NodeMapToList.end()) {
		return 0;
	}

	return p->second->second;
}

/* DataManager - end */
