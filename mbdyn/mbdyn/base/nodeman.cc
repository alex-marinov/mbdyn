/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

void DataManager::NodeManager(void)
{
   for(int i = 0; i < Node::LASTNODETYPE; i++) {
      NodeData[i].ppFirstNode = NULL;
      NodeData[i].iNum = 0;
      NodeData[i].fDefaultOut = fDefaultOut; /* Da "output.h" */
      NodeData[i].OutFile = OutputHandler::UNKNOWN; /* Da "output.h" */
   }
   
   /* Se un tipo scrive su un file di output, aggiungere qui il tipo di file */
   NodeData[Node::ABSTRACT].OutFile = OutputHandler::ABSTRACT;
   NodeData[Node::STRUCTURAL].OutFile = OutputHandler::STRNODES;   
   NodeData[Node::ELECTRIC].OutFile = OutputHandler::ELECTRIC;   
   NodeData[Node::HYDRAULIC].OutFile = OutputHandler::PRESNODES;   
   NodeData[Node::PARAMETER].OutFile = OutputHandler::PARAMETERS;   
}


void DataManager::NodeManagerDestructor(void)
{
   DEBUGCOUT("Entering DataManager::NodeManagerDestructor()" << std::endl);

   ASSERT(ppNodes != NULL);
   
   if(ppNodes != NULL) {
      Node** pp = ppNodes;
      while(pp < ppNodes+iTotNodes) {
	 ASSERT(*pp != NULL);
	 if(*pp != NULL) {		  
	    DEBUGCOUT("deleting node " << (*pp)->GetLabel() 
		      << ", type " << psNodeNames[(*pp)->GetNodeType()]
		      << std::endl);
	    SAFEDELETE(*pp);
	 }
	 
	 pp++;
      }	  
      
      DEBUGCOUT("deleting node structure" << std::endl);
      SAFEDELETEARR(ppNodes);
   }	
}


void DataManager::NodeDataInit(void)
{
   for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {      
      iTotNodes += NodeData[iCnt].iNum;
   }

   DEBUGCOUT("iTotNodes = " << iTotNodes << std::endl);
   
   if (iTotNodes > 0) {	
      SAFENEWARR(ppNodes, Node*, iTotNodes);
      
      NodeIter.Init(ppNodes, iTotNodes);
      
      /* Con iteratore:
       Node* pTmp = NULL;
       if(NodeIter.bGetFirst(pTmp)) {	  
	  do {
	     pTmp = NULL;
	  } while(NodeIter.bGetNext(pTmp));
       }
       */
	
      Node** ppTmp = ppNodes;
      while(ppTmp < ppNodes+iTotNodes) {	 
	 *ppTmp++ = NULL;
      }      
      
	
      NodeData[0].ppFirstNode = ppNodes;
      for(int iCnt = 0; iCnt < Node::LASTNODETYPE-1; iCnt++) {	 
	 NodeData[iCnt+1].ppFirstNode =
	   NodeData[iCnt].ppFirstNode+NodeData[iCnt].iNum;
      }      
   } else {
      silent_cerr("warning, no nodes are defined" << std::endl);
   }
}


void DataManager::NodeOutput(OutputHandler& OH) const
{
   Node** ppTmpNode = ppNodes;
   for(; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {      
     (*ppTmpNode)->Output(OH);
   }   

#if 0
   /* Con iteratore: */
    Node* pTmpNode = NULL;
    if(NodeIter.bGetFirst(pTmpNode)) {      
       do {
	  pTmpNode->Output(OH);
       } while(NodeIter.bGetNext(pTmpNode));
    }
#endif
}


void
DataManager::NodeOutput(
		OutputHandler& OH,
		const VectorHandler& X,
		const VectorHandler& XP
		) const
{
   Node** ppTmpNode = ppNodes;
   for(; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {      
     (*ppTmpNode)->Output(OH, X, XP);
   }   

#if 0
   /* Con iteratore: */
    Node* pTmpNode = NULL;
    if(NodeIter.bGetFirst(pTmpNode)) {      
       do {
	  pTmpNode->Output(OH, X, XP);
       } while(NodeIter.bGetNext(pTmpNode));
    }
#endif
}


void
DataManager::NodeOutput_pch(
		std::ostream& pch
		) const
{
   Node** ppTmpNode = ppNodes;
   for(; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {      
     (*ppTmpNode)->Output_pch(pch);
   }   

#if 0
   /* Con iteratore: */
    Node* pTmpNode = NULL;
    if(NodeIter.bGetFirst(pTmpNode)) {      
       do {
	  pTmpNode->Output_pch(pch);
       } while(NodeIter.bGetNext(pTmpNode));
    }
#endif
}


void
DataManager::NodeOutput_f06(
		std::ostream& f06,
		const VectorHandler& X
		) const
{
   Node** ppTmpNode = ppNodes;
   for(; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {      
     (*ppTmpNode)->Output_f06(f06, X);
   }   

#if 0
   /* Con iteratore: */
    Node* pTmpNode = NULL;
    if(NodeIter.bGetFirst(pTmpNode)) {      
       do {
	  pTmpNode->Output_f06(f06, X);
       } while(NodeIter.bGetNext(pTmpNode));
    }
#endif
}


void
DataManager::NodeOutput_f06(
		std::ostream& f06,
		const VectorHandler& Xr,
		const VectorHandler& Xi
		) const
{
   Node** ppTmpNode = ppNodes;
   for(; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {      
     (*ppTmpNode)->Output_f06(f06, Xr, Xi);
   }   

#if 0
   /* Con iteratore: */
    Node* pTmpNode = NULL;
    if(NodeIter.bGetFirst(pTmpNode)) {      
       do {
	  pTmpNode->Output_f06(f06, Xr, Xi);
       } while(NodeIter.bGetNext(pTmpNode));
    }
#endif
}


flag DataManager::fGetDefaultOutputFlag(const Node::Type& t) const
{
   return NodeData[t].fDefaultOut;
}


/* cerca un nodo qualsiasi */
Node* DataManager::pFindNode(Node::Type Typ, unsigned int uL) const
{
   ASSERT(NodeData[Typ].ppFirstNode != NULL);
   ASSERT(NodeData[Typ].iNum > 0);
   ASSERT(uL > 0);

   return pLabelSearch(NodeData[Typ].ppFirstNode, NodeData[Typ].iNum, uL);
}


/* cerca un nodo strutturale*/
#if defined(USE_STRUCT_NODES)   
StructNode* DataManager::pFindStructNode(unsigned int uL) const
{
   ASSERT(NodeData[Node::STRUCTURAL].ppFirstNode != NULL);
   ASSERT(NodeData[Node::STRUCTURAL].iNum > 0);
   ASSERT(uL > 0);
  
   return (StructNode*)pLabelSearch(NodeData[Node::STRUCTURAL].ppFirstNode, NodeData[Node::STRUCTURAL].iNum, uL);
}
#endif // USE_STRUCT_NODES

/* cerca un nodo elettrico */
#if defined(USE_ELECTRIC_NODES)   
ElectricNode* DataManager::pFindElectricNode(unsigned int uL) const
{
   ASSERT(NodeData[Node::ELECTRIC].ppFirstNode != NULL);
   ASSERT(NodeData[Node::ELECTRIC].iNum > 0);
   ASSERT(uL > 0);
  
   return (ElectricNode*)pLabelSearch(NodeData[Node::ELECTRIC].ppFirstNode, NodeData[Node::ELECTRIC].iNum, uL);
}
#endif // USE_ELECTRIC_NODES

/* DataManager - end */
