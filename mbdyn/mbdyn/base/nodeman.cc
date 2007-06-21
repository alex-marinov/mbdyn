/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <search.h>

/* DataManager - begin */

const char sNMClassName[] = "DataManager";

void
DataManager::NodeManager(void)
{
	for (int i = 0; i < Node::LASTNODETYPE; i++) {
		NodeData[i].ppFirstNode = NULL;
		NodeData[i].iNum = 0;
		NodeData[i].uFlags = 0U;
		NodeData[i].DefaultOut(::fDefaultOut == 1); /* Da "output.h" */
		NodeData[i].OutFile = OutputHandler::UNKNOWN; /* Da "output.h" */
	}

	/* Se un tipo scrive su un file di output, aggiungere qui il tipo di file */
	NodeData[Node::ABSTRACT].OutFile = OutputHandler::ABSTRACT;
	NodeData[Node::STRUCTURAL].OutFile = OutputHandler::STRNODES;
	NodeData[Node::ELECTRIC].OutFile = OutputHandler::ELECTRIC;
   NodeData[Node::THERMAL].OutFile = OutputHandler::THERMALNODES;   
	NodeData[Node::HYDRAULIC].OutFile = OutputHandler::PRESNODES;
	NodeData[Node::PARAMETER].OutFile = OutputHandler::PARAMETERS;
}

void
DataManager::NodeManagerDestructor(void)
{
	DEBUGCOUT("Entering DataManager::NodeManagerDestructor()" << std::endl);

	ASSERT(ppNodes != NULL);

	if (ppNodes != NULL) {
		Node** pp = ppNodes;
		while (pp < ppNodes+iTotNodes) {
			ASSERT(*pp != NULL);
			if (*pp != NULL) {
				DEBUGCOUT("deleting "
					<< psNodeNames[(*pp)->GetNodeType()] << "(" << (*pp)->GetLabel() << ")"
					<< std::endl);
				SAFEDELETE(*pp);
			}
			pp++;
		}

		DEBUGCOUT("deleting node structure" << std::endl);
		SAFEDELETEARR(ppNodes);
	}
}

void
DataManager::NodeDataInit(void)
{
	for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		iTotNodes += NodeData[iCnt].iNum;
	}

	DEBUGCOUT("iTotNodes = " << iTotNodes << std::endl);

	if (iTotNodes > 0) {
		SAFENEWARR(ppNodes, Node*, iTotNodes);

		NodeIter.Init(ppNodes, iTotNodes);

		Node** ppTmp = ppNodes;
		while (ppTmp < ppNodes + iTotNodes) {
			*ppTmp++ = NULL;
		}

		NodeData[0].ppFirstNode = ppNodes;
		for (int iCnt = 0; iCnt < Node::LASTNODETYPE-1; iCnt++) {
			NodeData[iCnt+1].ppFirstNode =
				NodeData[iCnt].ppFirstNode + NodeData[iCnt].iNum;
		}
	} else {
		silent_cerr("warning, no nodes are defined" << std::endl);
	}
}

void
DataManager::NodeOutputPrepare(OutputHandler& OH)
{
	Node** ppTmpNode = ppNodes;
	for (; ppTmpNode < ppNodes + iTotNodes; ppTmpNode++) {
		(*ppTmpNode)->OutputPrepare(OH);
	}
}

void
DataManager::NodeOutput(OutputHandler& OH) const
{
	Node** ppTmpNode = ppNodes;
	for (; ppTmpNode < ppNodes + iTotNodes; ppTmpNode++) {
		(*ppTmpNode)->Output(OH);
	}
}

void
DataManager::NodeOutput(
	OutputHandler& OH,
	const VectorHandler& X,
	const VectorHandler& XP) const
{
	Node** ppTmpNode = ppNodes;
	for (; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {
		(*ppTmpNode)->Output(OH, X, XP);
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
	ASSERT(NodeData[Typ].ppFirstNode != NULL);
	ASSERT(NodeData[Typ].iNum > 0);
	ASSERT(uL > 0);

	return pLabelSearch(NodeData[Typ].ppFirstNode, NodeData[Typ].iNum, uL);
}

/* cerca un nodo strutturale*/
StructNode*
DataManager::pFindStructNode(unsigned int uL) const
{
	ASSERT(NodeData[Node::STRUCTURAL].ppFirstNode != NULL);
	ASSERT(NodeData[Node::STRUCTURAL].iNum > 0);
	ASSERT(uL > 0);

	return (StructNode*)pLabelSearch(NodeData[Node::STRUCTURAL].ppFirstNode, NodeData[Node::STRUCTURAL].iNum, uL);
}

/* cerca un nodo elettrico */
ElectricNode*
DataManager::pFindElectricNode(unsigned int uL) const
{
	ASSERT(NodeData[Node::ELECTRIC].ppFirstNode != NULL);
	ASSERT(NodeData[Node::ELECTRIC].iNum > 0);
	ASSERT(uL > 0);

	return (ElectricNode*)pLabelSearch(NodeData[Node::ELECTRIC].ppFirstNode, NodeData[Node::ELECTRIC].iNum, uL);
}

/* DataManager - end */
