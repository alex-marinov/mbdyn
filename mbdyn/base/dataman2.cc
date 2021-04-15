/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* DataManager -
 * continua qui perche' il file dataman.cc sta diventando lungo */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <algorithm>
#include <set>
#include <cmath>
#include <sstream>
// #include <typeinfo>

#include "dataman.h"
#include "dataman_.h"

#include "gravity.h"
#include "aerodyn.h"
#include "solver.h"
#include "ls.h"

#include "Rot.hh"
#include "naivemh.h"
#include "spmapmh.h"

#include "bufferstream_out_elem.h"
#include "bufferstreamdrive.h"

const LoadableCalls *
DataManager::GetLoadableElemModule(std::string name) const
{
	for (int j = 0; name[j]; j++) {
		name[j] = tolower(name[j]);
	}

	typedef std::map<std::string,const LoadableCalls *> mleh;
	mleh::const_iterator i = MapOfLoadableElemHandlers.find(name);
	if (i == MapOfLoadableElemHandlers.end()) {
		return 0;
	}
	return i->second;
}

void
DataManager::SetLoadableElemModule(std::string name,
		const LoadableCalls *calls, ModuleInsertMode mode)
{
	for (int j = 0; name[j]; j++) {
		name[j] = tolower(name[j]);
	}

	const LoadableCalls *tmp = GetLoadableElemModule(name);

	if (tmp != 0) {
		switch (mode) {
		case MIM_FAIL:
		default:
			silent_cerr("DataManager::SetLoadableElemModule(): "
				"loadable element handler \"" << name
				<< "\" already defined" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

		case MIM_IGNORE:
			silent_cout("DataManager::SetLoadableElemModule(): "
				"loadable element handler \"" << name
				<< "\" already defined; "
				"new definition ignored" << std::endl);
			return;

		case MIM_REPLACE:
			silent_cout("DataManager::SetLoadableElemModule(): "
				"loadable element handler \"" << name
				<< "\" already defined; "
				"replaced by new definition" << std::endl);
			break;
		}
	}

	MapOfLoadableElemHandlers[name] = calls;
}

const doublereal&
DataManager::dGetInitialPositionStiffness(void) const
{
	return dInitialPositionStiffness;
}

const doublereal&
DataManager::dGetInitialVelocityStiffness(void) const
{
	return dInitialVelocityStiffness;
}

bool
DataManager::bDoesOmegaRotate(void) const
{
	return bOmegaRotates;
}

void
DataManager::IncElemCount(Elem::Type type)
{
	/* FIXME: assert the data structure has not been allocated yet */
	ElemData[type].iExpectedNum++;
}

/* Setta il valore della variabile Time nel DataManager
 * usato dal metodo numerico all'inizio di ogni step temporale */

void
DataManager::SetTime(const doublereal& dTime, const doublereal& dTimeStep,
	const integer& iStep, bool bServePending)
{
	/* Setta il tempo nel DriveHandler */
	DrvHdl.SetTime(dTime, dTimeStep, iStep);

	DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
			"Global symbol table:" << std::endl
	  		<< MathPar.GetSymbolTable() << std::endl);

	/* serve i drive pending */
	if (bServePending) {
		for (int iType = 0; iType < Drive::LASTDRIVETYPE; iType++) {
			for (unsigned int iCnt = 0; iCnt < DriveData[iType].iNum; iCnt++) {
				DriveData[iType].ppFirstDrive[iCnt]->ServePending(dTime);
			}
		}
	}

	// updates rigid body kinematics, if any
	if (pRBK) {
		pRBK->Update();
	}
} /* End of DataManager::SetTime() */

doublereal
DataManager::dGetTime(void) const
{
	return DrvHdl.dGetTime();
} /* End of DataManager::dGetTime() */

/* Collega il DataManager ed il DriveHandler alla soluzione */
void
DataManager::LinkToSolution(VectorHandler& XCurr,
	VectorHandler& XPrimeCurr)
{
	pXCurr = &XCurr;
	pXPrimeCurr = &XPrimeCurr;
	DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
}

/* Inizializzatore dei dof di ogni elemento */
void
DataManager::DofOwnerInit(void)
{
	DEBUGCOUTFNAME("DataManager::DofOwnerInit");
	ASSERT(!Dofs.empty());
	ASSERT(!Nodes.empty());

	if ( uPrintFlags & PRINT_TO_FILE ) {
		OutHdl.Open(OutputHandler::DOFSTATS);
	}

	std::ostream& out_ds = (uPrintFlags & PRINT_TO_FILE)
		? OutHdl.DofStats()
		: std::cout;

	bool pds =
#ifdef DEBUG
		DEBUG_LEVEL_MATCH(MYDEBUG_INIT|MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
		(!silent_output && (uPrintFlags & PRINT_DOF_STATS));

	/* NOTE: further direct use of std::cout instead
	 * of silent_cout() macro because silent_cout is
	 * tested in "pds".
	 */
	if (pds) {
		out_ds << "Regular steps dof stats" << std::endl;
	}

	/* per ogni nodo strutturale */
	if (!NodeData[Node::STRUCTURAL].NodeContainer.empty()) {

		/*
		 * used by POD stuff: if any, output
		 * the list of the first dof (minus 1)
		 * of each structural node, so it's easy
		 * to get the struct node values
		 * in MATLAB: given a vector "X" with all
		 * the states, and a vector
		 * "v" with the first dof of each
		 * structural node, then the x coordinate
		 * is X(v+1) and so forth
		 */

		OutHdl.Log() << "struct node dofs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode) {
				if (pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pNode->iGetFirstPositionIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct node eqs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode) {
				if (pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pNode->iGetFirstMomentumIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct node momentum dofs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode) {
				switch (pNode->GetStructNodeType()) {
				case StructNode::STATIC:
				case StructNode::DUMMY:
					continue;

				default:
					break;
				}
				OutHdl.Log() << " " << pNode->iGetFirstMomentumIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct node momentum eqs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode) {
				switch (pNode->GetStructNodeType()) {
				case StructNode::STATIC:
				case StructNode::DUMMY:
					continue;

				default:
					break;
				}
				OutHdl.Log() << " " << pNode->iGetFirstIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct displacement node dofs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			if (pDispNode) {
				const StructNode* pNode = dynamic_cast<const StructNode*>(pDispNode);
				if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pDispNode->iGetFirstPositionIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct displacement node eqs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			if (pDispNode) {
				const StructNode* pNode = dynamic_cast<const StructNode*>(pDispNode);
				if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pDispNode->iGetFirstMomentumIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct displacement node momentum dofs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			if (pDispNode) {
				const StructNode* pNode = dynamic_cast<const StructNode*>(pDispNode);
				if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				switch (pDispNode->GetStructDispNodeType()) {
				case StructDispNode::STATIC:
					continue;

				default:
					break;
				}
				OutHdl.Log() << " " << pDispNode->iGetFirstMomentumIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct displacement node momentum eqs:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			if (pDispNode) {
				const StructNode* pNode = dynamic_cast<const StructNode*>(pDispNode);
				if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				switch (pDispNode->GetStructDispNodeType()) {
				case StructNode::STATIC:
					continue;

				default:
					break;
				}
				OutHdl.Log() << " " << pDispNode->iGetFirstIndex();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct node labels:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode) {
				if (pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pNode->GetLabel();
			}
		}
		OutHdl.Log() << std::endl;

		OutHdl.Log() << "struct displacement node labels:";
		for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			if (pDispNode) {
				const StructNode* pNode = dynamic_cast<const StructNode*>(pDispNode);
				if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
					continue;
				}
				OutHdl.Log() << " " << pDispNode->GetLabel();
			}
		}
		OutHdl.Log() << std::endl;
	}

	/* per ogni nodo */
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
				psNodeNames[(*i)->GetNodeType()]
				<< "(" << (*i)->GetLabel() << ")"
				<< std::endl);

		/* chiede al nodo quanti dof possiede */
		unsigned int iNumDof = (*i)->iGetNumDof();
		if (iNumDof > 0) {
			ASSERT((*i)->iGetFirstIndex() >= 0);

			/* si fa passare il primo Dof */
			Dof* pDf = &Dofs[(*i)->iGetFirstIndex()];

#ifdef DEBUG
			DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
					psNodeNames[(*i)->GetNodeType()]
					<< "(" << (*i)->GetLabel() << "): "
					"first dof = " << pDf->iIndex + 1
					<< std::endl);
#endif /* DEBUG */

			if (pds) {
				unsigned int nd = (*i)->iGetNumDof();
				integer fd = pDf->iIndex;

				out_ds << psNodeNames[(*i)->GetNodeType()]
					<< "(" << (*i)->GetLabel() << "): "
					<< nd << " " << fd + 1;
				if (nd > 1) {
					out_ds << "->" << fd + nd;
				}
				out_ds << std::endl;
				if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
					(*i)->DescribeDof(out_ds,
							     "        ");
				}

				if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
					(*i)->DescribeEq(out_ds,
							     "        ");
				}
			}

			/* per ogni Dof, chiede al nodo di che tipo e' e lo
			 * setta nel DofOwner */
			std::vector<std::string> DofDesc;
			(*i)->DescribeDof(DofDesc);
			if (DofDesc.size() == iNumDof) {
				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					pDf[iCnt].Description = DofDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*i)->GetNodeType()]
					<< "(" << (*i)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": dof(" << iCnt + 1 << ")";
					pDf[iCnt].Description = os.str();
				}
			}

			std::vector<std::string> EqDesc;
			(*i)->DescribeEq(EqDesc);
			if (EqDesc.size() == iNumDof) {
				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					pDf[iCnt].EqDescription = EqDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*i)->GetNodeType()]
					<< "(" << (*i)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": equation(" << iCnt + 1 << ")";
					pDf[iCnt].EqDescription = os.str();
				}
			}

			for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
				pDf[iCnt].Order = (*i)->GetDofType(iCnt);
				pDf[iCnt].EqOrder = (*i)->GetEqType(iCnt);
			}
		}
	}

	/* per ogni elemento */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			ASSERT(pEl != NULL);
			DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
					"Elem type " << pEl->GetElemType()
					<< " (" << psElemNames[pEl->GetElemType()]
					<< "(" << pEl->GetLabel() << "))" << std::endl);

			/* chiede all'elemento quanti dof possiede */
			unsigned int iNumDof = pEl->iGetNumDof();
			if (iNumDof > 0) {
				ElemWithDofs* pEWD = Cast<ElemWithDofs>(pEl);

				ASSERT(pEWD->iGetFirstIndex() >= 0);

				/* si fa passare il DofOwner */
				Dof* pDf = &Dofs[pEWD->iGetFirstIndex()];

#ifdef DEBUG
				DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
						psElemNames[pEl->GetElemType()]
						<< "(" << pEWD->GetLabel() << "): "
						"first dof = " << pDf->iIndex + 1
						<< std::endl);
#endif /* DEBUG */

				if (pds) {
					unsigned int nd = pEWD->iGetNumDof();
					integer fd = pDf->iIndex;

					out_ds << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << "): "
						<< nd << " " << fd + 1;
					if (nd > 1) {
						out_ds << "->" << fd + nd;
					}
					out_ds << std::endl;
					if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
						pEWD->DescribeDof(out_ds,
								"        ");
					}

					if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
						pEWD->DescribeEq(out_ds,
								"        ");
					}
				}


				/* per ogni Dof, chiede all'elemento
				 * di che tipo e' e lo setta
				 * nel DofOwner */
				std::vector<std::string> DofDesc;
				pEWD->DescribeDof(DofDesc);
				if (DofDesc.size() == iNumDof) {
					for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
						pDf[iCnt].Description = DofDesc[iCnt];
					}

				} else {
					std::ostringstream os;
					os << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << ")";
					std::string name(os.str());

					for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
						os.str(name);
						os.seekp(0, std::ios_base::end);
						os << ": dof(" << iCnt + 1 << ")";
						pDf[iCnt].Description = os.str();
					}
				}

				std::vector<std::string> EqDesc;
				pEWD->DescribeEq(EqDesc);
				if (EqDesc.size() == iNumDof) {
					for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
						pDf[iCnt].EqDescription = EqDesc[iCnt];
					}

				} else {
					std::ostringstream os;
					os << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << ")";
					std::string name(os.str());

					for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
						os.str(name);
						os.seekp(0, std::ios_base::end);
						os << ": equation(" << iCnt + 1 << ")";
						pDf[iCnt].EqDescription = os.str();
					}
				}

				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					pDf[iCnt].Order = pEWD->GetDofType(iCnt);
					pDf[iCnt].EqOrder = pEWD->GetEqType(iCnt);
				}
			}
		} while (ElemIter.bGetNext(pEl));
	}

	/* FIXME: this should rather go before initial assembly */
	/* NOTE: we run code anyway, but only print if requested/allowed
	 * for consistency checking purposes */

	/* per ogni elemento */
	if (ElemIter.bGetFirst(pEl)) {
		/* create node connection structure */
		typedef std::set<const Elem *> elmap;
		typedef std::map<const Node *, elmap *> nodemap;
		std::vector<nodemap> connectedElems(Node::LASTNODETYPE);

		/* element connections get populated directly by elements */
		std::vector<const Node *> connectedNodes;

		if (uPrintFlags & PRINT_EL_CONNECTION) {
			out_ds << "Element connections" << std::endl;
		}

		do {
			pEl->GetConnectedNodes(connectedNodes);

			if (connectedNodes.size() > 0) {
				if (uPrintFlags & PRINT_EL_CONNECTION) {
					out_ds << psElemNames[pEl->GetElemType()]
						<< "(" << pEl->GetLabel() << ") connecting" << std::endl;
				}
				for (std::vector<const Node *>::const_iterator i = connectedNodes.begin();
					i != connectedNodes.end(); ++i)
				{
					const Node *real_i = (*i)->GetNode();
					if (uPrintFlags & PRINT_EL_CONNECTION) {
						out_ds << "        "
							<< psNodeNames[real_i->GetNodeType()]
							<< "(" << real_i->GetLabel() << ")" << std::endl;
					}

					nodemap::iterator n = connectedElems[real_i->GetNodeType()].find(real_i);
					if (n == connectedElems[real_i->GetNodeType()].end()) {
						connectedElems[real_i->GetNodeType()][real_i] = new elmap;
					}
					connectedElems[real_i->GetNodeType()][real_i]->insert(pEl);
				}

			} else {
				if (uPrintFlags & PRINT_EL_CONNECTION) {
					out_ds << psElemNames[pEl->GetElemType()]
						<< "(" << pEl->GetLabel() << ") not connected" << std::endl;
				}
			}

		} while (ElemIter.bGetNext(pEl));


		if (uPrintFlags & PRINT_NODE_CONNECTION) {
			out_ds << "Node connections" << std::endl;
		}
		for (unsigned t = 0; t < Node::LASTNODETYPE; t++) {
			for (nodemap::iterator n = connectedElems[t].begin();
				n != connectedElems[t].end(); ++n)
			{
				if (uPrintFlags & PRINT_NODE_CONNECTION) {
					out_ds << psNodeNames[n->first->GetNodeType()]
						<< "(" << n->first->GetLabel() << ") connected to" << std::endl;
				}
				for (elmap::const_iterator e = n->second->begin();
					e != n->second->end(); ++e)
				{
					if (uPrintFlags & PRINT_NODE_CONNECTION) {
						out_ds << "        "
							<< psElemNames[(*e)->GetElemType()]
							<< "(" << (*e)->GetLabel() << ")" << std::endl;
					}
				}

				delete n->second;
				n->second = 0;
			}
		}
	}
} /* End of DataManager::DofOwnerInit() */

/* Inizializzazione della struttura dei dof
 * per l'assemblaggio iniziale dei vincoli */
void
DataManager::InitialJointAssembly(void)
{
	if ( uPrintFlags & PRINT_TO_FILE ) {
		OutHdl.Open(OutputHandler::DOFSTATS);
	}

	std::ostream& out_ds = (uPrintFlags & PRINT_TO_FILE)
		? OutHdl.DofStats()
		: std::cout;

	/* Costruisce la struttura temporanea dei Dof */

	ASSERTMSG(DofData[DofOwner::JOINT].iNum > 0,
		"Warning, no joints are defined; "
		"You shouldn't have reached this point");
	ASSERT(DofData[DofOwner::STRUCTURALNODE].iNum > 0);

	/* Nodi strutturali: mette gli indici ai DofOwner */
	bool pds =
#ifdef DEBUG
			DEBUG_LEVEL_MATCH(MYDEBUG_INIT|MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
			(!silent_output && (uPrintFlags & PRINT_DOF_STATS));

	if (pds) {
		out_ds << "Initial assembly dof stats" << std::endl;
	}

	/* Numero totale di Dof durante l'assemblaggio iniziale */
	integer iInitialNumDofs = 0;
	for (NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
		i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
	{
		iInitialNumDofs += dynamic_cast<const StructDispNode*>(i->second)->iGetInitialNumDof();
	}

	/* Elementi: mette gli indici agli eventuali DofOwner */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Per ogni tipo di elemento */
		if (ElemData[iCnt1].bToBeUsedInAssembly() && !ElemData[iCnt1].ElemContainer.empty()) {
			/* Se deve essere usato nell'assemblaggio e ne sono definiti */

			/* Tipo di dof dell'elemento corrente */
			DofOwner::Type CurrDofType =
				ElemData[iCnt1].DofOwnerType;

			if (CurrDofType != DofOwner::UNKNOWN) {
				ASSERT((unsigned)DofData[CurrDofType].iNum == ElemData[iCnt1].ElemContainer.size());

				/* Iterazione sugli Elem */
				for (ElemContainerType::const_iterator e = ElemData[iCnt1].ElemContainer.begin();
					e != ElemData[iCnt1].ElemContainer.end(); ++e)
				{
					InitialAssemblyElem *pEl = dynamic_cast<InitialAssemblyElem *>(e->second);
					if (pEl == 0 || (bNotDeformableInitial && pEl->bIsDeformable())) {
						/* Ignore elements
						 * not subjected
						 * to initial assembly */
						continue;
					}

					ElemWithDofs *pDOEl = dynamic_cast<ElemWithDofs *>(e->second);
					if (pDOEl == 0) {
						/* Ignore elements subjected
						 * to initial assembly
						 * but without dofs */
						continue;
					}

					iInitialNumDofs += pEl->iGetInitialNumDof();
				}
			}
		}
	}

	/*
	 * Alla fine, i DofOwner di nodi e joint contengono gli indici giusti per
	 * l'assemblaggio iniziale. Corrispondono a:
	 * - per ogni nodo:
	 *   - posizione x
	 *   - parametri di rotazione g
	 *   - velocita' xP
	 *   - velocita' angolare omega
	 * - per ogni joint:
	 *   - se vincolo in posizione, reazione e sua derivata
	 *   - se vincolo in velocita', reazione.
	 * - per vincoli misti si hanno reazioni ed eventualmente loro derivate
	 *   in base al tipo
	 */

	/* Creazione e costruzione array Dof */
	Dofs.resize(iInitialNumDofs);

	// just to make sure nothing strange occurs...
	iTotDofs = iInitialNumDofs;

	integer iIndex;    /* Indice dei gradi di liberta' */
	for (iIndex = 0; iIndex < iInitialNumDofs; iIndex++) {
		Dofs[iIndex].iIndex = iIndex;
	}

	/* mette a posto i dof */
	iIndex = 0;
	for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
		n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
	{
		// numero di dof di un owner
		const StructDispNode *pNode = dynamic_cast<const StructDispNode *>(n->second);
		unsigned int iNumDofs = pNode->iGetInitialNumDof();
		if (iNumDofs > 0) {
			DofOwner* pFDO = const_cast<DofOwner *>(pNode->pGetDofOwner());
			pFDO->iNumDofs = iNumDofs;
			pFDO->iFirstIndex = iIndex;
			if (pds) {
				unsigned int nd = iNumDofs;
				integer fd = iIndex;

				out_ds << psNodeNames[pNode->GetNodeType()]
					<< "(" << pNode->GetLabel()
					<< "): " << nd << " " << fd + 1;
				if (nd > 1) {
					out_ds << "->" << fd + nd;
				}
				out_ds << std::endl;
				if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
					pNode->DescribeDof(out_ds,
							     "        ", true);
				}

				if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
					pNode->DescribeEq(out_ds,
							     "        ", true);
				}
			}

			std::vector<std::string> DofDesc;
			pNode->DescribeDof(DofDesc, true);
			if (DofDesc.size() == iNumDofs) {
				for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
					Dofs[iIndex + iCnt].Description = DofDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[pNode->GetNodeType()]
					<< "(" << pNode->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": dof(" << iCnt + 1 << ")";
					Dofs[iIndex + iCnt].Description = os.str();
				}
			}

			std::vector<std::string> EqDesc;
			pNode->DescribeEq(EqDesc, true);
			if (EqDesc.size() == iNumDofs) {
				for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
					Dofs[iIndex + iCnt].EqDescription = EqDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[pNode->GetNodeType()]
					<< "(" << pNode->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": equation(" << iCnt + 1 << ")";
					Dofs[iIndex + iCnt].EqDescription = os.str();
				}
			}

			iIndex += iNumDofs;

		} else {
			pedantic_cerr(psNodeNames[pNode->GetNodeType()]
				<< "(" << n->second->GetLabel() << ") has 0 dofs"
				<< std::endl);
		}
	}

	/* Elementi: mette gli indici agli eventuali DofOwner */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Pre ogni tipo di elemento */
		if (ElemData[iCnt1].bToBeUsedInAssembly() && !ElemData[iCnt1].ElemContainer.empty()) {
			/* Se deve essere usato nell'assemblaggio e ne sono definiti */

			/* Tipo di dof dell'elemento corrente */
			DofOwner::Type CurrDofType =
				ElemData[iCnt1].DofOwnerType;

			if (CurrDofType != DofOwner::UNKNOWN) {
				ASSERT((unsigned)DofData[CurrDofType].iNum == ElemData[iCnt1].ElemContainer.size());

				/* Iterazione sugli Elem */
				for (ElemContainerType::const_iterator p = ElemData[iCnt1].ElemContainer.begin();
					p != ElemData[iCnt1].ElemContainer.end();
					++p)
				{
					InitialAssemblyElem *pEl = dynamic_cast<InitialAssemblyElem *>(p->second);
					if (pEl == 0 || (bNotDeformableInitial && pEl->bIsDeformable())) {
						/* Ignore elements
						 * not subjected
						 * to initial assembly */
						continue;
					}

					ElemWithDofs *pDOEl = dynamic_cast<ElemWithDofs *>(p->second);
					if (pDOEl == 0) {
						/* Ignore elements subjected
						 * to initial assembly
						 * but without dofs */
						continue;
					}

					DofOwner *pDO = const_cast<DofOwner *>(pDOEl->pGetDofOwner());
					ASSERT(pDO != 0);
					// numero di dof di un owner
					unsigned int iNumDofs = pEl->iGetInitialNumDof();
					pDO->iNumDofs = iNumDofs;
					if (iNumDofs > 0) {
						pDO->iFirstIndex = iIndex;
						if (pds) {
							unsigned int nd = iNumDofs;
							integer fd = iIndex;
							ElemWithDofs* pEWD = Cast<ElemWithDofs>(p->second);

							out_ds << psElemNames[pEl->GetElemType()]
								<< "(" << pEl->GetLabel()
								<< "): " << nd << " " << fd + 1;
							if (nd > 1) {
								out_ds << "->" << fd + nd;
							}
							out_ds << std::endl;
							if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
								pEWD->DescribeDof(out_ds,
										"        ", true);
							}

							if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
								pEWD->DescribeEq(out_ds,
										"        ", true);
							}
						}

						std::vector<std::string> DofDesc;
						pEl->DescribeDof(DofDesc, true);
						if (DofDesc.size() == iNumDofs) {
							for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
								Dofs[iIndex + iCnt].Description = DofDesc[iCnt];
							}

						} else {
							std::ostringstream os;
							os << psElemNames[pEl->GetElemType()]
								<< "(" << pEl->GetLabel() << ")";
							std::string name(os.str());

							for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
								os.str(name);
								os.seekp(0, std::ios_base::end);
								os << ": dof(" << iCnt + 1 << ")";
								Dofs[iIndex + iCnt].Description = os.str();
							}
						}

						std::vector<std::string> EqDesc;
						pEl->DescribeEq(EqDesc, true);
						if (EqDesc.size() == iNumDofs) {
							for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
								Dofs[iIndex + iCnt].EqDescription = EqDesc[iCnt];
							}

						} else {
							std::ostringstream os;
							os << psElemNames[pEl->GetElemType()]
								<< "(" << pEl->GetLabel() << ")";
							std::string name(os.str());

							for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
								os.str(name);
								os.seekp(0, std::ios_base::end);
								os << ": equation(" << iCnt + 1 << ")";
								Dofs[iIndex + iCnt].EqDescription = os.str();
							}
						}

						iIndex += iNumDofs;

					} else {
						pedantic_cerr(psElemNames[iCnt1]
							<< "(" << pEl->GetLabel() << ") "
							"has 0 dofs" << std::endl);
					}
				}
			}
		}
	}

	ASSERT(iIndex == iInitialNumDofs);

	/* Trova le massime dimensioni del workspace
	 * per l'assemblaggio iniziale */
	integer iMaxRowsRes = 0;
	integer iMaxRowsJac = 0;
	integer iMaxColsJac = 0;
	integer iMaxItemsJac = 0;

	InitialAssemblyIterator IAIter(&ElemData);
	InitialAssemblyElem* pEl = IAIter.GetFirst();
	while (pEl != NULL) {
		integer iCurrRows = 0;
		integer iCurrCols = 0;
		pEl->InitialWorkSpaceDim(&iCurrRows, &iCurrCols);

		if (iCurrRows >= 0) {
			// Assume a full Jacobian matrix
			iMaxRowsJac = std::max(iMaxRowsJac, iCurrRows);
			iMaxColsJac = std::max(iMaxColsJac, iCurrCols);
		} else {
			// Assume a sparse Jacobian matrix
			iCurrRows = std::abs(iCurrRows);
		}

		iMaxRowsRes = std::max(iMaxRowsRes, iCurrRows);
		iMaxItemsJac = std::max(iMaxItemsJac, iCurrRows * iCurrCols);

		pEl = IAIter.GetNext();
	}

	/* Ciclo di iterazioni fino a convergenza */

	/* Crea il solutore lineare, tenendo conto dei tipi
	 * supportati, di quanto scelto nel file di configurazione
	 * e di eventuali paraametri extra */
	SolutionManager* pSM = CurrSolver.GetSolutionManager(iInitialNumDofs);

#ifdef DEBUG_MEMMANAGER
	DEBUGLCOUT(MYDEBUG_MEM|MYDEBUG_ASSEMBLY,
			"After initialisation in InitialJointAssembly" << std::endl
			<< defaultMemoryManager << std::endl);
#endif /* DEBUG_MEMMANAGER */

	MyVectorHandler X(iInitialNumDofs);
	X.Reset();

	/* Linka il DriveHandler al vettore soluzione */
	LinkToSolution(X, X);

	/* Setta i valori iniziali dei gradi di liberta' dei nodi strutturali
	 * durante l'assemblaggio iniziale */
	for (NodeContainerType::iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
		i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
	{
		dynamic_cast<StructDispNode *>(i->second)->SetInitialValue(X);
	}

	/* Setta i valori iniziali dei gradi di liberta' dei vincoli
	 * durante l'assemblaggio iniziale */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Pre ogni tipo di elemento */
		if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
			ElemData[iCnt1].bToBeUsedInAssembly() &&
			!ElemData[iCnt1].ElemContainer.empty())
		{
			for (ElemContainerType::const_iterator p = ElemData[iCnt1].ElemContainer.begin();
				p != ElemData[iCnt1].ElemContainer.end(); ++p)
			{
				ElemWithDofs *pEWD = Cast<ElemWithDofs>(p->second);
				pEWD->SetInitialValue(X);
			}
		}
	}

	/* Vettore di lavoro */
	VectorHandler* pResHdl = pSM->pResHdl();
	MySubVectorHandler WorkVec(iMaxRowsRes);

	/* Matrice di lavoro */
	MatrixHandler* pMatHdl = pSM->pMatHdl();
	VariableSubMatrixHandler WorkMat(iMaxRowsJac, iMaxColsJac, iMaxItemsJac);

	/* Soluzione */
	VectorHandler* pSolHdl = pSM->pSolHdl();

	if (
#ifdef DEBUG
			DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
			outputIters())
	{
		silent_cout("Assembly Tol=" << dInitialAssemblyTol << std::endl);
	}

	/* Ciclo di assemblaggio */
	for (integer iNumIter = 0; ; iNumIter++) {
		/* Assemblo il residuo */
		pResHdl->Reset();

		/* Contributo dei nodi */
		for (NodeContainerType::iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			int iOffset = 3;
			if (pNode) {
				iOffset = 6;
			}

			integer iFirstIndex = pDispNode->iGetFirstPositionIndex();

			/* Nuova feature: ogni nodo ha la sua stiffness */
			doublereal dPosStiff = pDispNode->dGetPositionStiffness();
			doublereal dVelStiff = pDispNode->dGetVelocityStiffness();

			Vec3 TmpVec;

			/* Posizione: k*Delta_x = k(x_0-x) + F */
			TmpVec = pDispNode->GetXPrev() - pDispNode->GetXCurr();
			pResHdl->Add(iFirstIndex + 1, TmpVec*dPosStiff);

			if (pNode) {
				/* Rotazione: k*Delta_g = -k*g(R_Delta) + M */
				Mat3x3 RDelta = pNode->GetRPrev().MulMT(pNode->GetRCurr());
				TmpVec = Vec3(CGR_Rot::Param, RDelta);
				pResHdl->Add(iFirstIndex + 4, TmpVec*dPosStiff);

				/* Velocita' angolare: k*(Delta_w+(R_Delta*w0)/\Delta_g) =
				 *                                    k*(R_Delta*w0-w) + M */
				const Vec3& wPrev(pNode->GetWPrev());
				const Vec3& wCurr(pNode->GetWCurr());

				if (pNode->bOmegaRotates()) {
					/* con questa la velocita' angolare e' solidale
					 * con il nodo */
					TmpVec = RDelta*wPrev - wCurr;
				} else {
					/* con questa la velocita' angolare e' solidale
					 * col riferimento assoluto */
					TmpVec = wPrev - wCurr;
				}

				pResHdl->Add(iFirstIndex + iOffset + 4, TmpVec*dVelStiff);
			}

			/* Velocita': k*Delta_v = k*(v0-Delta_v) + F */
			TmpVec = pDispNode->GetVPrev() - pDispNode->GetVCurr();
			pResHdl->Add(iFirstIndex + iOffset + 1, TmpVec*dVelStiff);

		}

		/* Elementi (con iteratore): */
		pEl = IAIter.GetFirst();
		while (pEl != NULL && !(bNotDeformableInitial && pEl->bIsDeformable())) {
			try {
				*pResHdl += pEl->InitialAssRes(WorkVec, X);
			}
			catch (Elem::ChangedEquationStructure& e) {
				// do nothing: Jacobian matrix
				// is always recomputed anyway...
			}
			pEl = IAIter.GetNext();
		}

		if (
#ifdef DEBUG
				DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL) ||
#endif /* DEBUG */
				outputRes())
		{
			/* Output del residuo */
			PrintResidual(*pResHdl, iNumIter);
		}

		/* Eseguo il test di convergenza; se e' positivo, esco */
		/* FIXME: why /(1.+X.Dot()) ??? */
		doublereal dTest = pResHdl->Dot()/(1. + X.Dot());
		if (!std::isfinite(dTest)) {
			silent_cerr("Assembly diverged; aborting..." << std::endl);
			throw DataManager::ErrAssemblyDiverged(MBDYN_EXCEPT_ARGS);
		}
		dTest = sqrt(dTest);

		if ((dTest <= dInitialAssemblyTol ||
			iNumIter >= iMaxInitialIterations) && (
#ifdef DEBUG
				DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
			outputIters()))
		{
			silent_cout("\tIteration(" << iNumIter << ") "
				"" << dTest << std::endl);
		}

		/* Se la tolleranza e' raggiunta, esce dal ciclo */
		if (dTest <= dInitialAssemblyTol) {
			DEBUGLCOUT(MYDEBUG_ASSEMBLY, "Initial assembly "
					"performed successfully in "
					<< iNumIter << " iterations"
					<< std::endl);
			goto endofcycle;
		}

		/* Se ho raggiunto il numero massimo di iterazioni */
		if (iNumIter >= iMaxInitialIterations) {
			silent_cerr("Initial assembly iterations "
				"reached maximum number "
				<< iMaxInitialIterations << "; aborting..."
				<< std::endl);
			throw DataManager::ErrAssemblyMaxIters(MBDYN_EXCEPT_ARGS);
		}

		/* Assemblo lo jacobiano e risolvo */
		pSM->MatrInitialize();

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
		NodesUpdateJac(1., NodeIter);
#endif	
		
		/* Contributo dei nodi */
		for (NodeContainerType::iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			const StructDispNode *pDispNode = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pNode = dynamic_cast<const StructNode *>(i->second);
			if (pNode && pNode->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			int iOffset = 3;
			if (pNode) {
				iOffset = 6;
			}

			integer iFirstIndex = pDispNode->iGetFirstPositionIndex();

			/* Nuova feature: ogni nodo ha la sua stiffness */
			doublereal dPosStiff = pDispNode->dGetPositionStiffness();
			doublereal dVelStiff = pDispNode->dGetVelocityStiffness();

			/* NOTE: iOffset equal to number of position/velocity equations */
			for (int iCnt = 1; iCnt <= iOffset; iCnt++) {
				/* Posizione, rotazione */
				integer iTmp = iFirstIndex + iCnt;
#ifdef USE_SPARSE_AUTODIFF
				sp_grad::SpGradient g;
				g.Reset(0., iTmp, dPosStiff);
				pMatHdl->AddItem(iTmp, g);
#else
				pMatHdl->PutCoef(iTmp, iTmp, dPosStiff);
#endif

				/* Velocita', velocita' angolare */
				iTmp += iOffset;
#ifdef USE_SPARSE_AUTODIFF
				g.Reset(0., iTmp, dVelStiff);
				pMatHdl->AddItem(iTmp, g);
#else
				pMatHdl->PutCoef(iTmp, iTmp, dVelStiff);
#endif
			}

			if (pNode && pNode->bOmegaRotates()) {
				/* con questi la velocita' angolare e' solidale con il nodo */

				/* Velocita' angolare - termine di rotazione: R_Delta*w0/\ */
				const Mat3x3& R0 = pNode->GetRPrev();
				const Mat3x3& R = pNode->GetRCurr();
				const Vec3& W0 = pNode->GetWPrev();
				Vec3 TmpVec = R*(R0.MulTV(W0*dVelStiff));

				/* W1 in m(3, 2), -W1 in m(2, 3) */
				doublereal d = TmpVec(1);
#ifdef USE_SPARSE_AUTODIFF
				sp_grad::SpGradient g;
				g.Reset(0., iFirstIndex + 5, d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 6, g);
				g.Reset(0., iFirstIndex + 6, -d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 5, g);
#else
				pMatHdl->PutCoef(iFirstIndex + iOffset + 6, iFirstIndex + 5, d);
				pMatHdl->PutCoef(iFirstIndex + iOffset + 5, iFirstIndex + 6, -d);
#endif

				/* W2 in m(1, 3), -W2 in m(3, 1) */
				d = TmpVec(2);
#ifdef USE_SPARSE_AUTODIFF
				g.Reset(0., iFirstIndex + 6, d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 4, g);
				g.Reset(0., iFirstIndex + 4, -d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 6, g);
#else
				pMatHdl->PutCoef(iFirstIndex + iOffset + 4, iFirstIndex + 6, d);
				pMatHdl->PutCoef(iFirstIndex + iOffset + 6, iFirstIndex + 4, -d);
#endif

				/* W3 in m(2, 1), -W3 in m(1, 2) */
				d = TmpVec(3);

#ifdef USE_SPARSE_AUTODIFF
				g.Reset(0., iFirstIndex + 4, d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 5, g);
				g.Reset(0.,  iFirstIndex + 5, -d);
				pMatHdl->AddItem(iFirstIndex + iOffset + 4, g);
#else
				pMatHdl->PutCoef(iFirstIndex + iOffset + 5, iFirstIndex + 4, d);
				pMatHdl->PutCoef(iFirstIndex + iOffset + 4, iFirstIndex + 5, -d);
#endif
			} /* altrimenti la velocita' angolare e' solidale con il nodo */
		}

		/* Contributo degli elementi */
		pEl = IAIter.GetFirst();
		while (pEl != NULL && !(bNotDeformableInitial && pEl->bIsDeformable())) {
			*pMatHdl += pEl->InitialAssJac(WorkMat, X);
			pEl = IAIter.GetNext();
		}

		if (
#ifdef DEBUG
				DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_JAC) ||
#endif /* DEBUG */
				outputJac())
		{
		     if (silent_out) {
			silent_cout("Jacobian:" << std::endl);
			// Use the same format like the nonlinear solver
			pMatHdl->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
		     }
		}

		/* Fattorizza e risolve con jacobiano e residuo appena calcolati */
		try {
			pSM->Solve();
		}
		catch (LinearSolver::ErrFactor& err) {
			silent_cerr("Initial assembly failed because no pivot element "
				"could be found for column " << err.iCol
				<< " (" << GetDofDescription(err.iCol) << "); "
				"aborting..." << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (
#ifdef DEBUG
				DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY|MYDEBUG_RESIDUAL) ||
#endif /* DEBUG */
				outputSol())
		{
			/* Output della soluzione */
			PrintSolution(*pSolHdl, iNumIter);
		}

		if (outputIters() || outputSolverConditionNumber()) {
			if (outputIters()) {
				silent_cout("\tIteration(" << iNumIter << ") " << std::setw(12) << dTest << " J");
		}

		if (outputSolverConditionNumber()) {
			silent_cout(" cond=");
				doublereal dCond;
				if (pSM->bGetConditionNumber(dCond)) {
					silent_cout(dCond);
				} else {
					silent_cout("NA");
				}
			}

			silent_cout(std::endl);
                }

		/* Aggiorno la soluzione */
		if (dEpsilon != 1.) {
			*pSolHdl *= dEpsilon;
		}
		X += *pSolHdl;

		/* Correggo i nodi */
		for (NodeContainerType::iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			i != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++i)
		{
			dynamic_cast<StructDispNode *>(i->second)->InitialUpdate(X);
		}
	}

endofcycle:
	/* Resetta e distrugge la struttura temporanea dei Dof */

	/* Elementi: rimette a posto il numero di Dof propri dei vincoli */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Per ogni tipo di elemento */
		if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
			ElemData[iCnt1].bToBeUsedInAssembly() &&
			!ElemData[iCnt1].ElemContainer.empty())
		{
			/* Se possiede dofs, se deve essere usato nell'assemblaggio
			 * e se ne sono presenti */
			for (ElemContainerType::const_iterator p = ElemData[iCnt1].ElemContainer.begin();
				p != ElemData[iCnt1].ElemContainer.end();
				++p)
			{
				ElemWithDofs *pEWD = Cast<ElemWithDofs>(p->second);
				DofOwner *pDO = const_cast<DofOwner *>(pEWD->pGetDofOwner());
				pDO->iNumDofs = p->second->iGetNumDof();
			}
		}
	}

	/* Dealloca il vettore dei Dof */
	ASSERT(!Dofs.empty());

	// restore
	iTotDofs = 0;

	SAFEDELETE(pSM);
} /* End of InitialJointAssembly */

/* Aggiorna i DofOwner con il numero di dofs dell'elemento */

void
DataManager::DofOwnerSet(void)
{
	DEBUGCOUTFNAME("DataManager::DofOwnerSet");

	/* Setta i DofOwner dei nodi */
	for (NodeVecType::iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		DofOwner* pDO = const_cast<DofOwner *>((*i)->pGetDofOwner());
		pDO->iNumDofs = (*i)->iGetNumDof();
	}

	/* Setta i DofOwner degli elementi (chi li possiede) */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		DofOwner::Type DT = ElemData[iCnt].DofOwnerType;
		if (DT != DofOwner::UNKNOWN) {
			DEBUGLCOUT(MYDEBUG_INIT, "Elem type " << iCnt
					<< " (" << psElemNames[iCnt] << ")"
					<< std::endl);

			for (ElemContainerType::const_iterator p = ElemData[iCnt].ElemContainer.begin();
				p != ElemData[iCnt].ElemContainer.end(); ++p)
			{
				ElemWithDofs* pEWD = Cast<ElemWithDofs>(p->second);

				DEBUGLCOUT(MYDEBUG_INIT, "    " << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << ")" << std::endl);

				DofOwner* pDO = const_cast<DofOwner *>(pEWD->pGetDofOwner());
				pDO->iNumDofs = pEWD->iGetNumDof();
				DEBUGLCOUT(MYDEBUG_INIT, "    num dofs: " << pDO->iNumDofs << std::endl);
			}
		}
	}
} /* end of DofOwnerSet() */


void
DataManager::SetValue(VectorHandler& X, VectorHandler& XP)
{
	/* Nodi */
	for (NodeVecType::iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->SetValue(this, X, XP);
	}

	/* Elementi */
	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->SetValue(this, X, XP);
		} while (ElemIter.bGetNext(pEl));
	}
	if (solArrFileName != NULL) {
		std::ifstream fp(solArrFileName);
#ifdef HAVE_ISOPEN
   		if (!fp.is_open()) {
			silent_cerr("DataManager::SetValue(): "
				"Cannot open file \"" << solArrFileName << "\""
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#endif /* HAVE_ISOPEN */
		fp.read((char*)X.pdGetVec() , X.iGetSize()*sizeof(double));
		if (fp.gcount() != std::streamsize(X.iGetSize()*sizeof(double))) {
			silent_cerr("DataManager::SetValue(): "
				"File(" << solArrFileName << ") too short!"
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		fp.read((char*)XP.pdGetVec() , XP.iGetSize()*sizeof(double));
		if (fp.gcount() != std::streamsize(XP.iGetSize()*sizeof(double))) {
			silent_cerr("DataManager::SetValue(): "
				"File(" << solArrFileName << ") too short!"
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		SAFEDELETEARR(solArrFileName);
		fp.close();
	}
} /* End of SetValue */


/* Output dati */
void
DataManager::OutputPrepare(void)
{
#ifdef USE_NETCDF
	/* Set up NetCDF stuff if required */
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		OutHdl.Open(OutputHandler::NETCDF);
		ASSERT(OutHdl.IsOpen(OutputHandler::NETCDF));

		Var_Step = OutHdl.CreateVar<integer>("run.step", 
			OutputHandler::Dimensions::Dimensionless, "time step index");
		Var_Time = OutHdl.CreateVar<doublereal>("time", 
			OutputHandler::Dimensions::Time, "simulation time");
		Var_TimeStep = OutHdl.CreateVar<doublereal>("run.timestep", 
			OutputHandler::Dimensions::Time, "integration time step");
	}
#endif /* USE_NETCDF */

	/* Dati dei nodi */
	NodeOutputPrepare(OutHdl);

	/* Dati degli elementi */
	ElemOutputPrepare(OutHdl);
}

/* Output setup for Eigenanalysis parameters */
void
DataManager::OutputEigPrepare(const integer iNumAnalyses, const integer iSize)
{
#ifdef USE_NETCDF
	/* Set up additional NetCDF stuff for eigenanalysis output */
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {

		OutputHandler::NcDimVec dim(1);
		dim[0] = OutHdl.CreateDim("eigensolutions", iNumAnalyses);

		m_Dim_Eig_iSize = OutHdl.CreateDim("eig_iSize", iSize);
		m_Dim_Eig_iComplex = OutHdl.CreateDim("complex_var_dim", 2);

		OutputHandler::AttrValVec attrs2(2);
		attrs2[0] = OutputHandler::AttrVal("type", "integer");
		attrs2[1] = OutputHandler::AttrVal("description",
				"timestep index of eigensolution");

		Var_Eig_lStep = OutHdl.CreateVar("eig.step", MbNcInt, attrs2, dim);

		OutputHandler::AttrValVec attrs3(3);
		attrs3[0] = OutputHandler::AttrVal("units", "s");
		attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
		attrs3[2] = OutputHandler::AttrVal("description",
				"simulation time at which the eigensolution was computed");

		Var_Eig_dTime = OutHdl.CreateVar("eig.time", MbNcDouble, attrs3, dim);

		attrs3[0] = OutputHandler::AttrVal("units", "-");
		attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
		attrs3[2] = OutputHandler::AttrVal("description",
				"coefficient used to build Aplus and Aminus matrices");

		Var_Eig_dCoef = OutHdl.CreateVar("eig.dCoef", MbNcDouble, attrs3, dim);

		OutputHandler::NcDimVec dim2(2);
		integer iNumNodes = NodeData[Node::STRUCTURAL].NodeContainer.size();

		dim2[0] = dim[0];
		dim2[1] = OutHdl.CreateDim("eig_iIdxSize", iNumNodes);

		attrs2[0] = OutputHandler::AttrVal("type", "integer");
		attrs2[1] = OutputHandler::AttrVal("description",
				"structural nodes base index");

		Var_Eig_Idx = OutHdl.CreateVar("eig.idx", MbNcInt, attrs2, dim2);
	}
#endif /* USE_NETCDF */
}

/* Output of Eigenanalysis parameters */
void
DataManager::OutputEigParams(const doublereal& dTime,
		const doublereal& dCoef,
		const unsigned uCurrEigSol,
		const int iResultsPrecision)
{
	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iResultsPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iResultsPrecision + 7;
			out.width(iNewWidth);
			out.precision(iResultsPrecision);
		}

		// header
		out
			<< "% time: " << dTime << std::endl;
		out
			<< "dTime = " << dTime << ';' << std::endl;

		// coefficient
		out
			<< "% coefficient" << std::endl
			<< "dCoef = " << dCoef << ";" << std::endl;
	}
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		long lStep = OutHdl.GetCurrentStep();
		OutHdl.WriteNcVar(Var_Eig_dTime, dTime, uCurrEigSol);
		OutHdl.WriteNcVar(Var_Eig_lStep, lStep, uCurrEigSol);
		OutHdl.WriteNcVar(Var_Eig_dCoef, dCoef, uCurrEigSol);
	}
#endif /* USE_NETCDF */
}

void
DataManager::OutputEigFullMatrices(const MatrixHandler* pMatA,
			const MatrixHandler* pMatB,
			const unsigned uCurrEigSol,
			const int iMatrixPrecision)
{
	const FullMatrixHandler& MatB = dynamic_cast<const FullMatrixHandler &>(*pMatB);
	const FullMatrixHandler& MatA = dynamic_cast<const FullMatrixHandler &>(*pMatA);
	integer nrows = MatB.iGetNumRows();
	integer ncols = MatB.iGetNumCols();

#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		OutputHandler::NcDimVec dim2(2);
		dim2[0] = m_Dim_Eig_iSize;
		dim2[1] = m_Dim_Eig_iSize;

		OutputHandler::AttrValVec attrs3(3);
		attrs3[0] = OutputHandler::AttrVal("units", "-");
		attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
		attrs3[2] = OutputHandler::AttrVal("description", "F/xPrime + dCoef * F/x");

		std::stringstream varname_ss;
		varname_ss << "eig." << uCurrEigSol << ".Aplus";
		Var_Eig_dAplus = OutHdl.CreateVar(varname_ss.str(), MbNcDouble, attrs3, dim2);
		attrs3[2] = OutputHandler::AttrVal("description", "F/xPrime - dCoef * F/x");

		varname_ss.str("");
		varname_ss.clear();
		varname_ss << "eig." << uCurrEigSol << ".Aminus";
		Var_Eig_dAminus = OutHdl.CreateVar(varname_ss.str(), MbNcDouble, attrs3, dim2);


#if defined(USE_NETCDFC)
		Var_Eig_dAplus->put(MatB.pdGetMat(), nrows, ncols);
		Var_Eig_dAminus->put(MatA.pdGetMat(), nrows, ncols);
#elif defined(USE_NETCDF4)  /*! USE_NETCDFC */
		std::vector<size_t> ncStartPos, ncCount;
		ncStartPos.push_back(0); // implicit cast here ok?
		ncStartPos.push_back(0); // implicit cast here ok?
		ncCount.push_back(nrows);
		ncCount.push_back(ncols);
		Var_Eig_dAplus.putVar(MatB.pdGetMat()); // seems that there is no purpose in giving matrix size as for old c++ (legacy) interface...
		Var_Eig_dAminus.putVar(MatA.pdGetMat());
#endif  /* USE_NETCDF4 */

	}
#endif /* USE_NETCDF */

	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iMatrixPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iMatrixPrecision + 7;
			out.width(iNewWidth);
			out.precision(iMatrixPrecision);
		}

		out
			<< "% F/xPrime + dCoef * F/x" << std::endl
			<< "Aplus" << " = [";


		for (integer r = 1; r <= nrows; r++) {
			for (integer c = 1; c <= ncols; c++) {
				out << MatB(r, c) << ' ';
			}

			if (r == nrows) {
				out << "];" << std::endl;

			} else {
				out << ";" << std::endl;
			}
		}

		out
			<< "% F/xPrime - dCoef * F/x" << std::endl
			<< "Aminus" << " = [";

		for (integer r = 1; r <= nrows; r++) {
			for (integer c = 1; c <= ncols; c++) {
				out << MatA(r, c) << ' ';
			}

			if (r == nrows) {
				out << "];" << std::endl;

			} else {
				out << ";" << std::endl;
			}
		}

	}
}

void
DataManager::OutputEigSparseMatrices(const MatrixHandler* pMatA,
	const MatrixHandler* pMatB,
	const unsigned uCurrEigSol,
	const int iMatrixPrecision)
{
	const SpMapMatrixHandler& MatB = dynamic_cast<const SpMapMatrixHandler &>(*pMatB);
	const SpMapMatrixHandler& MatA = dynamic_cast<const SpMapMatrixHandler &>(*pMatA);

	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iMatrixPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iMatrixPrecision + 7;
			out.width(iNewWidth);
			out.precision(iMatrixPrecision);
		}

		out
			<< "% F/xPrime + dCoef *F/x" << std::endl
			<< "Aplus" << " = [";

		for (SpMapMatrixHandler::const_iterator i = MatB.begin();
				i != MatB.end(); ++i)
		{
			if (i->dCoef != 0.) {
				out << i->iRow + 1 << " " << i->iCol + 1 << " " << i->dCoef << ";" << std::endl;
			}
		}

		out << "];" << std::endl
			<< "Aplus = spconvert(Aplus);" << std::endl;

		out
			<< "% F/xPrime - dCoef *F/x" << std::endl
			<< "Aminus" << " = [";

		for (SpMapMatrixHandler::const_iterator i = MatA.begin();
				i != MatA.end(); ++i)
		{
			if (i->dCoef != 0.) {
				out << i->iRow + 1 << " " << i->iCol + 1 << " " << i->dCoef << ";" << std::endl;
			}
		}

		out << "];" << std::endl
			<< "Aminus = spconvert(Aminus);" << std::endl;
	}
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		
		std::stringstream varname_ss;
		varname_ss << "eig." << uCurrEigSol << ".Aplus";
		Var_Eig_dAplus = OutHdl.CreateVar<Vec3>(varname_ss.str(), 
				OutputHandler::Dimensions::Dimensionless, "F/xPrime - dCoef * F/x");

		varname_ss.str("");
		varname_ss.clear();
		varname_ss << "eig." << uCurrEigSol << ".Aminus";
		Var_Eig_dAminus = OutHdl.CreateVar<Vec3>(varname_ss.str(), 
			OutputHandler::Dimensions::Dimensionless, "F/xPrime + dCoef * F/x");

		size_t iCnt = 0;
		Vec3 v;
		for (SpMapMatrixHandler::const_iterator i = MatB.begin();
				i != MatB.end(); ++i)
		{
			if (i->dCoef != 0.) {
				v = Vec3(i->iRow + 1, i->iCol + 1, i->dCoef);
				OutHdl.WriteNcVar(Var_Eig_dAplus, v, iCnt);
				iCnt++;
			}
		}

		iCnt = 0;
		for (SpMapMatrixHandler::const_iterator j = MatA.begin();
				j != MatA.end(); ++j)
		{
			if (j->dCoef != 0.) {
				v = Vec3(j->iRow + 1, j->iCol + 1, j->dCoef);
				OutHdl.WriteNcVar(Var_Eig_dAminus, v, iCnt);
				iCnt++;
			}
		}

	}
#endif
}

void
DataManager::OutputEigNaiveMatrices(const MatrixHandler* pMatA,
	const MatrixHandler* pMatB,
	const unsigned uCurrEigSol,
	const int iMatrixPrecision)
{
	const NaiveMatrixHandler& MatB = dynamic_cast<const NaiveMatrixHandler &>(*pMatB);
	const NaiveMatrixHandler& MatA = dynamic_cast<const NaiveMatrixHandler &>(*pMatA);

	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iMatrixPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iMatrixPrecision + 7;
			out.width(iNewWidth);
			out.precision(iMatrixPrecision);
		}

		out
			<< "% F/xPrime + dCoef *F/x" << std::endl
			<< "Aplus" << " = [";

		for (NaiveMatrixHandler::const_iterator i = MatB.begin();
				i != MatB.end(); ++i)
		{
			if (i->dCoef != 0.) {
				out << i->iRow + 1 << " " << i->iCol + 1 << " " << i->dCoef << ";" << std::endl;
			}
		}

		out << "];" << std::endl
			<< "Aplus = spconvert(Aplus);" << std::endl;

		out
			<< "% F/xPrime - dCoef *F/x" << std::endl
			<< "Aminus" << " = [";

		for (NaiveMatrixHandler::const_iterator j = MatA.begin();
				j != MatA.end(); ++j)
		{
			if (j->dCoef != 0.) {
				out << j->iRow + 1 << " " << j->iCol + 1 << " " << j->dCoef << ";" << std::endl;
			}
		}

		out << "];" << std::endl
			<< "Aminus = spconvert(Aminus);" << std::endl;
	}
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {

		std::stringstream varname_ss;
		varname_ss << "eig." << uCurrEigSol << ".Aplus";
		Var_Eig_dAplus = OutHdl.CreateVar<Vec3>(varname_ss.str(), 
			OutputHandler::Dimensions::Dimensionless, "F/xPrime + dCoef * F/x");

		varname_ss.str("");
		varname_ss.clear();
		varname_ss << "eig." << uCurrEigSol << ".Aminus";
		Var_Eig_dAminus = OutHdl.CreateVar<Vec3>(varname_ss.str(), 
			OutputHandler::Dimensions::Dimensionless, "F/xPrime - dCoef * F/x");

		size_t iCnt = 0;
		Vec3 v;
		for (NaiveMatrixHandler::const_iterator i = MatB.begin();
				i != MatB.end(); ++i)
		{
			if (i->dCoef != 0.) {
				v = Vec3(i->iRow + 1, i->iCol + 1, i->dCoef);
				OutHdl.WriteNcVar(Var_Eig_dAplus, v, iCnt);
				iCnt++;
			}
		}

		iCnt = 0;
		for (NaiveMatrixHandler::const_iterator i = MatA.begin();
				i != MatA.end(); ++i)
		{
			if (i->dCoef != 0.) {
				v = Vec3(i->iRow + 1, i->iCol + 1, i->dCoef);
				OutHdl.WriteNcVar(Var_Eig_dAminus, v, iCnt);
				iCnt++;
			}
		}

	}
#endif
}

void
DataManager::OutputEigGeometry(const unsigned uCurrEigSol, const int iResultsPrecision)
{
	NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
	NodeContainerType::const_iterator e = NodeData[Node::STRUCTURAL].NodeContainer.end();

	// no structural nodes!
	if (i == e) {
		return;
	}

	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iResultsPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iResultsPrecision + 7;
			out.width(iNewWidth);
			out.precision(iResultsPrecision);
		}

		out
			<< "% structural nodes labels" << std::endl
			<< "labels = [" << std::endl;

		for (; i != e; ++i) {
			const StructDispNode *pN = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pSN = dynamic_cast<const StructNode *>(pN);
			ASSERT(pN != 0);

			if (pSN && pSN->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			out << pN->GetLabel() << ";" << std::endl;
		}

		out << "];" << std::endl;

		out
			<< "% structural nodes base index" << std::endl
			<< "idx = [" << std::endl;

		for (i = NodeData[Node::STRUCTURAL].NodeContainer.begin(); i != e; ++i) {
			const StructDispNode *pN = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pSN = dynamic_cast<const StructNode *>(pN);
			ASSERT(pN != 0);

			if (pSN && pSN->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			out << pN->iGetFirstIndex() << ";" << std::endl;
		}

		out << "];" << std::endl;

		out
			<< "% structural nodes reference configuration (X, Phi)" << std::endl
			<< "X0 = [" << std::endl;

		for (i = NodeData[Node::STRUCTURAL].NodeContainer.begin(); i != e; ++i) {
			const StructDispNode *pN = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pSN = dynamic_cast<const StructNode *>(pN);
			ASSERT(pN != 0);

			if (pSN && pSN->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			const Vec3& X(pN->GetX());
			Vec3 Phi(mb_zero<Vec3>());
			if (pSN) {
				Phi = RotManip::VecRot(pSN->GetR());
			}

			out
				<< X(1) << ";" << std::endl
				<< X(2) << ";" << std::endl
				<< X(3) << ";" << std::endl
				<< Phi(1) << ";" << std::endl
				<< Phi(2) << ";" << std::endl
				<< Phi(3) << ";" << std::endl;
		}

		out << "];" << std::endl;
	}
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {

		/* start corner and count vector for NetCDF matrix output.
		 * Since we are writing a matrix element-by-element, count will
		 * always be (1, 1) and start will move to the desired place in the
		 * matrix */
		std::vector<size_t> start (2, 0);	
		const std::vector<size_t> count (2, 1);

		start[0] = uCurrEigSol;
		
		for (i = NodeData[Node::STRUCTURAL].NodeContainer.begin(); i != e; ++i) {
			const StructDispNode *pN = dynamic_cast<const StructDispNode *>(i->second);
			const StructNode *pSN = dynamic_cast<const StructNode *>(pN);
			integer iNodeIndex;
			ASSERT(pN != 0);

			if (pSN && pSN->GetStructNodeType() == StructNode::DUMMY) {
				start[1]++; 	// skip the dummy node
				continue;
			}

			iNodeIndex = pSN->iGetFirstIndex();
			OutHdl.WriteNcVar(Var_Eig_Idx, iNodeIndex, start, count); 
			start[1]++;
		}

	}
#endif // USE_NETCDF
}

void
DataManager::OutputEigenvectors(const VectorHandler *pBeta,
		const VectorHandler& R, const VectorHandler& I,
		const doublereal& dShiftR,
		const MatrixHandler *pVL, const MatrixHandler& VR,
		const std::vector<bool>& vOut,
		const unsigned uCurrEigSol,
		const int iResultsPrecision)
{
	const char signs[] = {'-', '+'};
	int iSign;

	integer iSize = VR.iGetNumRows();
	integer iNVec = VR.iGetNumCols();

	integer iEigenValues = 0;

	for (integer r = 1; r <= iNVec; r++) {
		if (!vOut[r - 1]) {
			continue;
		}
		++iEigenValues;
	}

	if (iEigenValues == 0)
		return; // this allows to load the .m file into Matlab/Octave even
			// if no eigenvalues have converged

	// alphar, alphai, beta
	if (OutHdl.UseText(OutputHandler::EIGENANALYSIS)) {
		std::ostream& out = OutHdl.Eigenanalysis();

		if (iResultsPrecision) {
			// 7 = number of characters requested by scientific notation
			int iNewWidth = iResultsPrecision + 7;
			out.width(iNewWidth);
			out.precision(iResultsPrecision);
		}

		out
			<< "% alphar, alphai, beta" << std::endl
			<< "alpha = [";

		for (integer r = 1; r <= iNVec; r++) {
			if (!vOut[r - 1]) {
				continue;
			}

			out
				<< R(r) + dShiftR << ' '
				<< I(r) << ' '
				<< (pBeta ? (*pBeta)(r) : 1.)
				<< ";" << std::endl;
		}

		out << "];" << std::endl;

		if (pVL) {
			// VL
			out
				<< "% left eigenvectors" << std::endl
				<< "VL = [" << std::endl;
			for (integer r = 1; r <= iSize; r++) {
				for (integer c = 1; c <= iNVec; c++) {
					if (!vOut[c - 1]) {
						continue;
					}

					if (I(c) != 0.) {
						ASSERTMSG(c < iNVec, "partial eigenanalysis output: complex eigenvalue with real part of left eigenvector only");
						ASSERT(I(c) > 0.);

						doublereal re = (*pVL)(r, c);
						// NOTE: we cannot assume that if c == iNVec
						// it corresponds to a real-valued eigenvalue;
						// "im" will be wrong, but at least we no not sigsegv
						doublereal im = (c < iNVec) ? (*pVL)(r, c + 1) : 0.;
						if (im < 0) {
							iSign = 0;
							im = -im;
						} else {
							iSign = 1;
						}

						out
							<< re << signs[iSign] << "i*" << im << ' ';
						if (vOut[c]) {
							out
								<< re << signs[1 - iSign] << "i*" << im << ' ';
						}
						c++;
					} else {
						 out
							<< (*pVL)(r, c) << ' ';
					}
				}

				if (r < iSize) {
					out << ";" << std::endl;
				} else {
					out << "];" << std::endl;
				}
			}
		}

		// VR
		out
			<< "% right eigenvectors" << std::endl
			<< "VR = [" << std::endl;
		for (integer r = 1; r <= iSize; r++) {
			for (integer c = 1; c <= iNVec; c++) {
				if (!vOut[c - 1]) {
					continue;
				}

				if(I(c) != 0.) {
					ASSERTMSG(c < iNVec, "partial eigenanalysis output: complex eigenvalue with real part of right eigenvector only");
					ASSERT(I(c) > 0.);

					doublereal re = VR(r, c);
					// NOTE: we cannote assume that if c == iNVec
					// it corresponds to a real-valued eigenvalue;
					// "im" will be wrong, but at least we do not sigsev
					doublereal im = (c < iNVec) ? VR(r, c + 1) : 0.;
					if (im < 0.) {
						iSign = 0;
						im = -im;
					} else {
						iSign = 1;
					}
					out
						<< re << signs[iSign] << "i*" << im << ' ';
					if (vOut[c]) {
						out
							<< re << signs[1 - iSign] << "i*" << im << ' ';
					}
					c++;
				} else {
					out
						<< VR(r, c) << ' ';
				}
			}

			if (r < iSize) {
				out << ";" << std::endl;
			} else {
				out << "];" << std::endl;
			}
		}
	}

#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		OutputHandler::NcDimVec dim_alpha(2);

		std::stringstream dimname_ss;
		dimname_ss << "eig_" << uCurrEigSol << "_iNVec_out";

		integer iNVecOut = 0;
		for (integer r = 1; r <= iNVec; r++)
		{
			if(vOut[r -1]){
				iNVecOut++;
			}
		}

		dim_alpha[0] = OutHdl.CreateDim(dimname_ss.str(), iNVecOut);
		dim_alpha[1] = OutHdl.DimV3();

		OutputHandler::AttrValVec attrs3(3);
		attrs3[0] = OutputHandler::AttrVal("units", "-");
		attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
		attrs3[2] = OutputHandler::AttrVal("description", "alpha matrix");

		std::stringstream varname_ss;
		varname_ss << "eig." << uCurrEigSol << ".alpha";
		Var_Eig_dAlpha = OutHdl.CreateVar(varname_ss.str(), MbNcDouble, attrs3, dim_alpha);

		Vec3 v;
		size_t uNRec = 0;
		for (integer r = 1; r <= iNVec; r++) {
			if (!vOut[r - 1]) {
				continue;
			}

			v(1) = R(r) + dShiftR;
			v(2) = I(r);
			v(3) = (pBeta ? (*pBeta)(r) : 1.);
			OutHdl.WriteNcVar(Var_Eig_dAlpha, v, uNRec);
			uNRec++;
		}

		OutputHandler::NcDimVec dim_v(3);
		dim_v[0] = m_Dim_Eig_iComplex;
		dim_v[1] = dim_alpha[0];
		dim_v[2] = m_Dim_Eig_iSize;

		/* start corner and count vector for NetCDF matrix output.
		 * Since we are writing a matrix element-by-element, count will
		 * always be (1, 1, 1) and start will move to the desired place in the
		 * matrix */
		std::vector<size_t> start (3, 0);
		const std::vector<size_t> count (3, 1);

		if (pVL) {
			// VL
			OutputHandler::AttrValVec attrs3(3);
			attrs3[0] = OutputHandler::AttrVal("units", "-");
			attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
			attrs3[2] = OutputHandler::AttrVal("description", "VL - Left eigenvectors matrix");

			std::stringstream varname_ss;
			varname_ss << "eig." << uCurrEigSol << ".VL";
			Var_Eig_dVL = OutHdl.CreateVar(varname_ss.str(), MbNcDouble, attrs3, dim_v);

			doublereal re;
			doublereal im;


			for (integer r = 1; r <= iSize; r++) {
				start[1] = 0;
				for (integer c = 1; c <= iNVec; c++) {
					if (!vOut[c - 1]) {
						continue;
					}

					start[2] = r - 1;
					if (I(c) != 0.) {
						ASSERTMSG(c < iNVec, "partial eigenanalysis output: complex eigenvalue with real part of left eigenvector only");
						ASSERT(I(c) > 0.);

						re = (*pVL)(r, c); // see above comment
						im = (c < iNVec) ? (*pVL)(r, c + 1) : 0.;

						// NetCDF indexing is zero-based!
						start[0] = 0;	// real part in first "page" of VL
						OutHdl.WriteNcVar(Var_Eig_dVL, re, start, count);

						start[0] = 1;	// imaginary part in second "page"
						OutHdl.WriteNcVar(Var_Eig_dVL, re, start, count);

						start[1]++;

						if (vOut[c]) {
							im = -im;

							start[0] = 0;
							OutHdl.WriteNcVar(Var_Eig_dVL, re, start, count);

							start[0] = 1;
							OutHdl.WriteNcVar(Var_Eig_dVL, im, start, count);

							start[1]++;
						}
						c++;
					} else {
						re = (*pVL)(r, c);
						im = 0.;

						start[0] = 0;
						OutHdl.WriteNcVar(Var_Eig_dVL, re, start, count);

						start[0] = 1;
						OutHdl.WriteNcVar(Var_Eig_dVL, im, start, count);

						start[1]++;
					}
				}
			}
		}

		// VR
		attrs3[0] = OutputHandler::AttrVal("units", "-");
		attrs3[1] = OutputHandler::AttrVal("type", "doublereal");
		attrs3[2] = OutputHandler::AttrVal("description", "VR - Right eigenvectors matrix");

		varname_ss.str("");
		varname_ss.clear();
		varname_ss << "eig." << uCurrEigSol << ".VR";
		Var_Eig_dVR = OutHdl.CreateVar(varname_ss.str(), MbNcDouble, attrs3, dim_v);

		doublereal re;
		doublereal im;
		for (integer r = 1; r <= iSize; r++) {
			start[1] = 0;
			for (integer c = 1; c <= iNVec; c++) {
				if (!vOut[c - 1]) {
					continue;
				}

				start[2] = r - 1;
				if (I(c) != 0.) {
					ASSERTMSG(c < iNVec, "partial eigenanalysis output: complex eigenvalue with real part of right eigenvector only");
					ASSERT(I(c) > 0.);

					re = VR(r, c); // see above comments
					im = (c < iNVec) ? VR(r, c + 1) : 0.;

					// NetCDF indexing is zero-based!
					start[0] = 0;	// real part in first 'page' of VR
					OutHdl.WriteNcVar(Var_Eig_dVR, re, start, count);

					start[0] = 1;	// imaginary part in second 'page'
					OutHdl.WriteNcVar(Var_Eig_dVR, im, start, count);

					start[1]++;

					if (vOut[c]) {
						im = -im;

						start[0] = 0;
						OutHdl.WriteNcVar(Var_Eig_dVR, re, start, count);

						start[0] = 1;
						OutHdl.WriteNcVar(Var_Eig_dVR, im, start, count);
						start[1]++;
					}
					c++;
				} else {
					re = VR(r, c);
					im = 0.;

					start[0] = 0;
					OutHdl.WriteNcVar(Var_Eig_dVR, re, start, count);

					start[0] = 1;
					OutHdl.WriteNcVar(Var_Eig_dVR, im, start, count);

					start[1]++;
				}
			}
		}

	}
#endif /* USE_NETCDF */
}

bool
DataManager::OutputEigClose(void)
{
	return OutHdl.Close(OutputHandler::EIGENANALYSIS);
}

/* Output dati */
bool
DataManager::Output(long lStep,
	const doublereal& dTime,
	const doublereal& dTimeStep,
	bool force) const
{
	/* Nota: il casting di OutHdl e' necessario in quanto la funzione propria
	 * <void DataManager::Output(void) const> e' dichiarata, appunto, <const>.
	 * Questo fa si' che un oggetto proprio della classe DataManager sia
	 * implicitamente definito come <const> agli occhi della funzione.
	 * Dal momento che le funzioni
	 * <void NodeManager::Output(OutputHandler&) const> e
	 * <void ElemManager::Output(OutputHandler&) const> ricevono come argomento
	 * un oggetto di tipo <OutputHandler&> che non e' <const> in quanto su di
	 * esso si scrive, il casting e' necessario per spiegare alla funzione
	 * <void DataManager::Output(void) const> che le funzioni invocate
	 * modificano si' l'<OutputHandler> passato loro, ma solo nel modo
	 * consentito e quindi la sua dichiarazione come funzione <const> e'
	 * dovuta al fatto che i dati propri non vengono modificati in modo
	 * incontrollabile */

	DriveTrace(OutHdl); // trace output will be written for every time step

	/* output only when allowed by the output meter */
	if (!force && !pOutputMeter->dGet()) {
		return false;
	}

	/*
	 * Write general simulation data to binary NetCDF file
	 *   the current time step index
	 *   the current simulatin time
	 *   the current integration time step
	 */
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		OutHdl.WriteNcVar(Var_Step, lStep);
		OutHdl.WriteNcVar(Var_Time, dTime);
		OutHdl.WriteNcVar(Var_TimeStep, dTimeStep);
	}
#endif /* USE_NETCDF */

	/* Dati dei nodi */
	NodeOutput(OutHdl);

	/* Dati degli elementi */
	ElemOutput(OutHdl);

	DriveOutput(OutHdl);

	OutHdl.IncCurrentStep();
#ifdef USE_NETCDF
	if (bNetCDFsync) {
		OutHdl.pGetBinFile()->sync(); // only works with netcdf-cxx4 >= 4.3.0, check implemented in configure.ac (see also https://github.com/Unidata/netcdf-cxx4/commit/e013ab35f0219fff92ed8237d2f385e89fd1cf77#diff-59778321f93df82ff613be3e32d9a3ec)

	}
#endif /* USE_NETCDF */

	return true;
}

void
DataManager::DriveTrace(OutputHandler& OH) const
{
	if (!OH.IsOpen(OutputHandler::TRACES)) {
		return;
	}

	const MBDynParser::DCType& DC = MBPar.GetDriveCallerContainer();

	for (MBDynParser::DCType::const_iterator i = DC.begin(); i != DC.end(); ++i) {
		i->second->Trace(OH);
	}
}

void
DataManager::DriveOutput(OutputHandler& OH) const
{
	if (!OH.IsOpen(OutputHandler::DRIVECALLERS)) {
		return;
	}

	const MBDynParser::DCType& DC = MBPar.GetDriveCallerContainer();

	for (MBDynParser::DCType::const_iterator i = DC.begin(); i != DC.end(); ++i) {
		i->second->Output(OH);
	}
}

/* Output dati */
void
DataManager::Output(const VectorHandler& X, const VectorHandler& XP) const
{
	/* Dati dei nodi */
	NodeOutput(OutHdl, X, XP);

	/* Dati degli elementi */
	ElemOutput(OutHdl, X, XP);
}

void
DataManager::BeforePredict(VectorHandler& X, VectorHandler& XP,
	VectorHandler& XPrev,
	VectorHandler& XPPrev) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->BeforePredict(X, XP, XPrev, XPPrev);
	}

	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->BeforePredict(X, XP, XPrev, XPPrev);
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::AfterPredict(void) const
{
	/* reset any external convergence requirement before starting
	 * a new step */
	for (Converged_t::iterator i = m_IsConverged.begin();
		i != m_IsConverged.end(); ++i)
	{
		*i = Converged::NOT_CONVERGED;
	}

	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		try {
			(*i)->AfterPredict(*pXCurr, *pXPrimeCurr);
		}
		catch (Elem::ChangedEquationStructure& e) {
			// ignore by now
			silent_cerr("DataManager::AfterPredict: "
				"warning, caught Elem::ChangedEquationStructure while processing "
				<< psNodeNames[(*i)->GetNodeType()] << "(" << (*i)->GetLabel() << ")" << std::endl);
		}
	}

	Elem* pEl = 0;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			try {
				pEl->AfterPredict(*pXCurr, *pXPrimeCurr);
			}
			catch (Elem::ChangedEquationStructure& e) {
				// ignore by now
				silent_cerr("DataManager::AfterPredict: "
					"warning, caught Elem::ChangedEquationStructure while processing "
					<< psElemNames[pEl->GetElemType()] << "(" << pEl->GetLabel() << ")" << std::endl);
			}
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::Update(void) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->Update(*pXCurr, *pXPrimeCurr);
	}

	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->Update(*pXCurr, *pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::AfterConvergence(void) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->AfterConvergence(*pXCurr, *pXPrimeCurr);
	}

	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->AfterConvergence(*pXCurr,
				*pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}

	/* Restart condizionato */
	switch (RestartEvery) {
	case NEVER:
		break;

	case ITERATIONS:
		if (++iCurrRestartIter == iRestartIterations) {
			iCurrRestartIter = 0;
			const_cast<DataManager *>(this)->MakeRestart();
		}
		break;

	case TIME: {
		doublereal dT = DrvHdl.dGetTime();
		if (dT - dLastRestartTime >= dRestartTime) {
			dLastRestartTime = dT;
			const_cast<DataManager *>(this)->MakeRestart();
		}
		break;
	}

	case TIMES: {
		doublereal dT = DrvHdl.dGetTime()
			+ pSolver->dGetInitialTimeStep()/100.;
		if (iCurrRestartTime == iNumRestartTimes) {
			break;
		}

		ASSERT(iCurrRestartTime < iNumRestartTimes);

		if (dT >= pdRestartTimes[iCurrRestartTime]) {
			iCurrRestartTime++;
			const_cast<DataManager *>(this)->MakeRestart();
		}
		break;
	}

	default:
		ASSERT(0);
		break;
	}
}


void
DataManager::DerivativesUpdate(void) const
{
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
	}

	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::PrintResidual(const VectorHandler& Res, integer iIterCnt) const
{
 	silent_cout("Residual(" << DrvHdl.iGetStep() << ":" << iIterCnt << ") "
		"t=" << DrvHdl.dGetTime()
		<< " dt=" << DrvHdl.dGetTimeStep() << std::endl);
	integer iSize = Res.iGetSize();
 	for (int iTmpCnt = 1; iTmpCnt <= iSize; iTmpCnt++) {
    		silent_cout("Eq  " << std::setw(8)
			<< iTmpCnt << ": "
			<< std::setw(20) << Res(iTmpCnt)
			<< " " << Dofs[iTmpCnt - 1].EqDescription
			<< std::endl);
	}
}

void
DataManager::PrintSolution(const VectorHandler& Sol, integer iIterCnt) const
{
 	silent_cout("Solution(" << DrvHdl.iGetStep() << ":" << iIterCnt << ") "
		"t=" << DrvHdl.dGetTime()
		<< " dt=" << DrvHdl.dGetTimeStep() << std::endl);
	integer iSize = Sol.iGetSize();
 	for (integer iTmpCnt = 1; iTmpCnt <= iSize; iTmpCnt++) {
    		silent_cout("Dof " << std::setw(8)
			<< iTmpCnt << ": "
			<< std::setw(20) << Sol(iTmpCnt)
			<< " " << Dofs[iTmpCnt - 1].Description
			<< std::endl);
	}
}

const std::string&
DataManager::GetDofDescription(int i) const
{
        if (i == -1) {
                static const std::string strUnknownDof("unknown");
                return strUnknownDof;
        }
    
	ASSERT(i > 0 && i <= iTotDofs);
	return Dofs[i - 1].Description;
}

const std::string&
DataManager::GetEqDescription(int i) const
{
	ASSERT(i > 0 && i <= iTotDofs);
	return Dofs[i - 1].EqDescription;
}

DofOrder::Order
DataManager::GetDofType(int i) const
{
	ASSERT(i > 0 && i <= iTotDofs);
	return Dofs[i - 1].Order;
}

DofOrder::Order
DataManager::GetEqType(int i) const
{
	ASSERT(i > 0 && i <= iTotDofs);
	return Dofs[i - 1].EqOrder;
}

bool
DataManager::bFDJac(void) const
{
	if (pFDJacMeter) {
		return (pFDJacMeter->dGet() != 0.);
	}

	return false;
}

unsigned
DataManager::ConvergedRegister(void)
{
	unsigned idx = m_IsConverged.size();
	m_IsConverged.resize(idx + 1);
	m_IsConverged[idx] = Converged::CONVERGED;
	return idx;
}

void
DataManager::ConvergedSet(unsigned idx, Converged::State s)
{
	ASSERT(idx < m_IsConverged.size());

	m_IsConverged[idx] = s;
}

bool
DataManager::IsConverged(void) const
{
	for (Converged_t::const_iterator i = m_IsConverged.begin();
		i != m_IsConverged.end(); ++i)
	{
		if (*i == Converged::NOT_CONVERGED) {
			return false;
		}
	}

	return true;
}

bool
DataManager::EndOfSimulation(void) const
{
	for (Converged_t::const_iterator i = m_IsConverged.begin();
		i != m_IsConverged.end(); ++i)
	{
		if (*i == Converged::END_OF_SIMULATION) {
			return true;
		}
	}

	return false;
}

std::vector<doublereal>&
DataManager::GetBufIn(unsigned uL)
{
	Drive* pD = pFindDrive(Drive::FILEDRIVE, uL);
	if (pD == 0) {
		silent_cerr("unable to find FileDrive(" << uL << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	BufferStreamDrive *pBSD = dynamic_cast<BufferStreamDrive *>(pD);
	if (pBSD == 0) {
		silent_cerr("FileDrive(" << uL << ") is not a BufferStreamDrive" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pBSD->GetBuf();
}

const std::vector<doublereal>&
DataManager::GetBufOut(unsigned uL) const
{
	BufferStreamElem *pBSE = pFindElem<BufferStreamElem, StreamOutElem, Elem::SOCKETSTREAM_OUTPUT>(uL);
	if (pBSE == 0) {
		silent_cerr("unable to find StreamOutElem(" << uL << "), or not a BufferStreamElem" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pBSE->GetBuf();
}

doublereal *
DataManager::GetBufInRaw(unsigned uL)
{
	Drive* pD = pFindDrive(Drive::FILEDRIVE, uL);
	if (pD == 0) {
		silent_cerr("unable to find FileDrive(" << uL << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	BufferStreamDrive_base *pBSD = dynamic_cast<BufferStreamDrive_base *>(pD);
	if (pBSD == 0) {
		silent_cerr("FileDrive(" << uL << ") is not a BufferStreamDrive_base" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// NOTE: this cast is needed because the input buf must be writeable.
	return (doublereal *)pBSD->GetBufRaw();
}

void
DataManager::SetBufInRaw(unsigned uL, integer n, const doublereal *p)
{
	Drive* pD = pFindDrive(Drive::FILEDRIVE, uL);
	if (pD == 0) {
		silent_cerr("unable to find FileDrive(" << uL << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	BufferStreamDriveRaw *pBSD = dynamic_cast<BufferStreamDriveRaw *>(pD);
	if (pBSD == 0) {
		silent_cerr("FileDrive(" << uL << ") is not a BufferStreamDriveRaw" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pBSD->bOwnsMemory()) {
		silent_cerr("FileDrive(" << uL << ") owns its memory, unable to set buffer" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pBSD->SetBufRaw(n, p);
}

const doublereal *
DataManager::GetBufOutRaw(unsigned uL) const
{
	BufferStreamElem_base *pBSE = pFindElem<BufferStreamElem_base, StreamOutElem, Elem::SOCKETSTREAM_OUTPUT>(uL);
	if (pBSE == 0) {
		silent_cerr("unable to find StreamOutElem(" << uL << "), or not a BufferStreamElem_base" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pBSE->GetBufRaw();
}

void
DataManager::SetBufOutRaw(unsigned uL, integer n, const doublereal *p)
{
	BufferStreamElemRaw *pBSE = pFindElem<BufferStreamElemRaw, StreamOutElem, Elem::SOCKETSTREAM_OUTPUT>(uL);
	if (pBSE == 0) {
		silent_cerr("unable to find StreamOutElem(" << uL << "), or not a BufferStreamElemRaw" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pBSE->bOwnsMemory()) {
		silent_cerr("StreamOutElem(" << uL << ") owns its memory, unable to set buffer" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pBSE->SetBufRaw(n, p);
}

/* DataManager - end */
