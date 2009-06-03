/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <set>

#include <cmath>
#include <sstream>

#include "dataman.h"
#include "dataman_.h"

#include "gravity.h"
#include "aerodyn.h"
#include "solver.h"
#include "driven.h"
#include "ls.h"

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

ElemWithDofs *
DataManager::CastElemWithDofs(Elem *pEl)
{
	ASSERT(pEl != NULL);

	ElemWithDofs *pEWD = dynamic_cast<ElemWithDofs *>(pEl);

	if (pEWD == 0) {
		DrivenElem *pDE = dynamic_cast<DrivenElem *>(pEl);
		if (pDE == 0) {
			silent_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") "
				"to ElemWithDofs" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pEWD = dynamic_cast<ElemWithDofs *>(pDE->pGetElem());
	}

	return pEWD;
}

ElemGravityOwner *
DataManager::CastElemGravityOwner(Elem *pEl)
{
	ASSERT(pEl != NULL);

	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pEl);

	if (pEGO == 0) {
		DrivenElem *pDE = dynamic_cast<DrivenElem *>(pEl);
		if (pDE == 0) {
			silent_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") "
				"to ElemGravityOwner" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pEGO = dynamic_cast<ElemGravityOwner *>(pDE->pGetElem());
	}

	return pEGO;
}

AerodynamicElem *
DataManager::CastAerodynamicElem(Elem *pEl)
{
	ASSERT(pEl != NULL);

	AerodynamicElem *pAEL = dynamic_cast<AerodynamicElem *>(pEl);

	if (pAEL == 0) {
		DrivenElem *pDE = dynamic_cast<DrivenElem *>(pEl);
		if (pDE == 0) {
			silent_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") "
				"to AerodynamicElem" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pAEL = dynamic_cast<AerodynamicElem *>(pDE->pGetElem());
	}

	return pAEL;
}

/* Collega il DataManager ed il DriveHandler alla soluzione */
void
DataManager::LinkToSolution(const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	pXCurr = const_cast<VectorHandler*>(&XCurr);
	pXPrimeCurr = const_cast<VectorHandler*>(&XPrimeCurr);
	DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
}

/* Inizializzatore dei dof di ogni elemento */
void
DataManager::DofOwnerInit(void)
{
	DEBUGCOUTFNAME("DataManager::DofOwnerInit");
	ASSERT(pDofs != NULL);
	ASSERT(ppNodes != NULL);

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
		std::cout << "Regular steps dof stats" << std::endl;
	}

	/* per ogni nodo strutturale */
	if (NodeData[Node::STRUCTURAL].iNum > 0) {

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

		StructNode** ppNd = (StructNode **)NodeData[Node::STRUCTURAL].ppFirstNode;
		for (unsigned long i = 0; i < NodeData[Node::STRUCTURAL].iNum; i++) {
			ASSERT(ppNd[i] != NULL);
			if (ppNd[i]->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}
			OutHdl.Log() << " " << ppNd[i]->iGetFirstIndex();
		}

		OutHdl.Log() << std::endl;
	}

	/* per ogni nodo */
	Node** ppNd = ppNodes;
	while (ppNd < ppNodes+iTotNodes) {
		ASSERT(*ppNd != NULL);
		DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
				psNodeNames[(*ppNd)->GetNodeType()]
				<< "(" << (*ppNd)->GetLabel() << ")"
				<< std::endl);

		unsigned int iNumDof;

		/* chiede al nodo quanti dof possiede */
		if ((iNumDof = (*ppNd)->iGetNumDof()) > 0) {
			/* si fa passare il primo Dof */
			Dof* pDf = pDofs + (*ppNd)->iGetFirstIndex();

#ifdef DEBUG
			DEBUGLCOUT(MYDEBUG_INIT|MYDEBUG_ASSEMBLY,
					psNodeNames[(*ppNd)->GetNodeType()]
					<< "(" << (*ppNd)->GetLabel() << "): "
					"first dof = " << pDf->iIndex+1
					<< std::endl);
#endif /* DEBUG */

			if (pds) {
				unsigned int nd = (*ppNd)->iGetNumDof();
				integer fd = pDf->iIndex;

				std::cout << psNodeNames[(*ppNd)->GetNodeType()]
					<< "(" << (*ppNd)->GetLabel() << "): "
					<< nd << " " << fd + 1;
				if (nd > 1) {
					std::cout << "->" << fd + nd;
				}
				std::cout << std::endl;
				if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
					(*ppNd)->DescribeDof(std::cout,
							     "        ");
				}

				if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
					(*ppNd)->DescribeEq(std::cout,
							     "        ");
				}
			}

			/* per ogni Dof, chiede al nodo di che tipo e' e lo
			 * setta nel DofOwner */
			std::vector<std::string> DofDesc;
			(*ppNd)->DescribeDof(DofDesc);
			if (DofDesc.size() == iNumDof) {
				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					pDf[iCnt].Description = DofDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*ppNd)->GetNodeType()]
					<< "(" << (*ppNd)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": dof(" << iCnt + 1 << ")";
					pDf[iCnt].Description = os.str();
				}
			}

			std::vector<std::string> EqDesc;
			(*ppNd)->DescribeEq(EqDesc);
			if (EqDesc.size() == iNumDof) {
				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					pDf[iCnt].EqDescription = EqDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*ppNd)->GetNodeType()]
					<< "(" << (*ppNd)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": equation(" << iCnt + 1 << ")";
					pDf[iCnt].EqDescription = os.str();
				}
			}

			for (unsigned int iCnt = 0; iCnt < iNumDof; iCnt++) {
				pDf[iCnt].Order = (*ppNd)->GetDofType(iCnt);
				pDf[iCnt].EqOrder = (*ppNd)->GetEqType(iCnt);
			}
		}
		ppNd++;
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

			unsigned int iNumDof;

			/* chiede all'elemento quanti dof possiede */
			if ((iNumDof = pEl->iGetNumDof()) > 0) {
				ElemWithDofs* pEWD = CastElemWithDofs(pEl);

				/* si fa passare il DofOwner */
				Dof* pDf = pDofs + pEWD->iGetFirstIndex();

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

					std::cout << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << "): "
						<< nd << " " << fd + 1;
					if (nd > 1) {
						std::cout << "->" << fd + nd;
					}
					std::cout << std::endl;
					if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
						pEWD->DescribeDof(std::cout,
								"        ");
					}

					if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
						pEWD->DescribeEq(std::cout,
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
			silent_cout("Element connections" << std::endl);
		}

		do {
			pEl->GetConnectedNodes(connectedNodes);

			if (connectedNodes.size() > 0) {
				if (uPrintFlags & PRINT_EL_CONNECTION) {
					silent_cout(psElemNames[pEl->GetElemType()]
						<< "(" << pEl->GetLabel() << ") connecting" << std::endl);
				}
				for (std::vector<const Node *>::const_iterator i = connectedNodes.begin();
					i != connectedNodes.end();
					i++)
				{
					const Node *real_i = (*i)->GetNode();
					if (uPrintFlags & PRINT_EL_CONNECTION) {
						silent_cout("        "
							<< psNodeNames[real_i->GetNodeType()]
							<< "(" << real_i->GetLabel() << ")" << std::endl);
					}
	
					nodemap::iterator n = connectedElems[real_i->GetNodeType()].find(real_i);
					if (n == connectedElems[real_i->GetNodeType()].end()) {
						connectedElems[real_i->GetNodeType()][real_i] = new elmap;
					}
					connectedElems[real_i->GetNodeType()][real_i]->insert(pEl);
				}

			} else {
				if (uPrintFlags & PRINT_EL_CONNECTION) {
					silent_cout(psElemNames[pEl->GetElemType()]
						<< "(" << pEl->GetLabel() << ") not connected" << std::endl);
				}
			}

		} while (ElemIter.bGetNext(pEl));


		if (uPrintFlags & PRINT_NODE_CONNECTION) {
			silent_cout("Node connections" << std::endl);
		}
		for (unsigned t = 0; t < Node::LASTNODETYPE; t++) {
			for (nodemap::iterator n = connectedElems[t].begin();
				n != connectedElems[t].end();
				n++)
			{
				if (uPrintFlags & PRINT_NODE_CONNECTION) {
					silent_cout(psNodeNames[n->first->GetNodeType()]
						<< "(" << n->first->GetLabel() << ") connected to" << std::endl);
				}
				for (elmap::const_iterator e = n->second->begin();
					e != n->second->end();
					e++)
				{
					if (uPrintFlags & PRINT_NODE_CONNECTION) {
						silent_cout("        "
							<< psElemNames[(*e)->GetElemType()]
							<< "(" << (*e)->GetLabel() << ")" << std::endl);
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
	/* Costruisce la struttura temporanea dei Dof */

	ASSERTMSG(DofData[DofOwner::JOINT].iNum > 0,
			"Warning, no joints are defined; "
			"You shouldn't have reached this point");
	ASSERT(DofData[DofOwner::STRUCTURALNODE].iNum > 0);

	/* Nodi strutturali: mette gli indici ai DofOwner */
	StructNode** ppFirstNode =
		(StructNode**)NodeData[Node::STRUCTURAL].ppFirstNode;
	integer iNumNodes = NodeData[Node::STRUCTURAL].iNum;

	StructNode** ppNode = ppFirstNode;
	DofOwner* pTmp = DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner;

	bool pds =
#ifdef DEBUG
			DEBUG_LEVEL_MATCH(MYDEBUG_INIT|MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
			(!silent_output && (uPrintFlags & PRINT_DOF_STATS));

	if (pds) {
		silent_cout("Initial assembly dof stats" << std::endl);
	}

	/* Numero totale di Dof durante l'assemblaggio iniziale */
	integer iInitialNumDofs = 0;
	for (int iCnt = 0; iCnt < iNumNodes; iCnt++) {
		iInitialNumDofs += ppFirstNode[iCnt]->iGetInitialNumDof();
	}

	/* Elementi: mette gli indici agli eventuali DofOwner */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Per ogni tipo di elemento */
		if (ElemData[iCnt1].bToBeUsedInAssembly() && !ElemData[iCnt1].ElemMap.empty()) {
			/* Se deve essere usato nell'assemblaggio e ne sono definiti */

			/* Tipo di dof dell'elemento corrente */
			DofOwner::Type CurrDofType =
				ElemData[iCnt1].DofOwnerType;

			if (CurrDofType != DofOwner::UNKNOWN) {
				ASSERT((unsigned)DofData[CurrDofType].iNum == ElemData[iCnt1].ElemMap.size());

				/* Iterazione sugli Elem */
				for (ElemMapType::const_iterator p = ElemData[iCnt1].ElemMap.begin();
					p != ElemData[iCnt1].ElemMap.end();
					p++)
				{
					InitialAssemblyElem *pEl = dynamic_cast<InitialAssemblyElem *>(p->second);
					if (pEl == 0) {
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
	SAFENEWARRNOFILL(pDofs, Dof, iInitialNumDofs);

	// just to make sure nothing strange occurs...
	iTotDofs = iInitialNumDofs;

	integer iIndex;    /* Indice dei gradi di liberta' */
	for (iIndex = 0; iIndex < iInitialNumDofs; iIndex++) {
		pDofs[iIndex].iIndex = iIndex;
	}

	/* mette a posto i dof */
	
	iIndex = 0;
	for (int iCnt = 1;
		pTmp < DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner
			+ DofData[DofOwner::STRUCTURALNODE].iNum;
		iCnt++, pTmp++, ppNode++)
	{
		// numero di dof di un owner
		unsigned int iNumDofs = pTmp->iNumDofs = (*ppNode)->iGetInitialNumDof();
		if (iNumDofs > 0) {
			pTmp->iFirstIndex = iIndex;
			if (pds) {
				unsigned int nd = iNumDofs;
				integer fd = iIndex;

				std::cout << psNodeNames[(*ppNode)->GetNodeType()]
					<< "(" << (*ppNode)->GetLabel()
					<< "): " << nd << " " << fd + 1;
				if (nd > 1) {
					std::cout << "->" << fd + nd;
				}
				std::cout << std::endl;
				if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
					(*ppNode)->DescribeDof(std::cout,
							     "        ", true);
				}

				if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
					(*ppNode)->DescribeEq(std::cout,
							     "        ", true);
				}
			}

			std::vector<std::string> DofDesc;
			(*ppNode)->DescribeDof(DofDesc, true);
			if (DofDesc.size() == iNumDofs) {
				for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
					pDofs[iIndex + iCnt].Description = DofDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*ppNode)->GetNodeType()]
					<< "(" << (*ppNode)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": dof(" << iCnt + 1 << ")";
					pDofs[iIndex + iCnt].Description = os.str();
				}
			}

			std::vector<std::string> EqDesc;
			(*ppNode)->DescribeEq(EqDesc, true);
			if (EqDesc.size() == iNumDofs) {
				for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
					pDofs[iIndex + iCnt].EqDescription = EqDesc[iCnt];
				}

			} else {
				std::ostringstream os;
				os << psNodeNames[(*ppNode)->GetNodeType()]
					<< "(" << (*ppNode)->GetLabel() << ")";
				std::string name(os.str());

				for (unsigned int iCnt = 0; iCnt < iNumDofs; iCnt++) {
					os.str(name);
					os.seekp(0, std::ios_base::end);
					os << ": equation(" << iCnt + 1 << ")";
					pDofs[iIndex + iCnt].EqDescription = os.str();
				}
			}

			iIndex += iNumDofs;

		} else {
			pedantic_cerr(psNodeNames[(*ppNode)->GetNodeType()]
				<< "(" << iCnt << ") has 0 dofs"
				<< std::endl);
		}
	}

	/* Elementi: mette gli indici agli eventuali DofOwner */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Pre ogni tipo di elemento */
		if (ElemData[iCnt1].bToBeUsedInAssembly() && !ElemData[iCnt1].ElemMap.empty()) {
			/* Se deve essere usato nell'assemblaggio e ne sono definiti */

			/* Tipo di dof dell'elemento corrente */
			DofOwner::Type CurrDofType =
				ElemData[iCnt1].DofOwnerType;

			if (CurrDofType != DofOwner::UNKNOWN) {
				ASSERT((unsigned)DofData[CurrDofType].iNum == ElemData[iCnt1].ElemMap.size());

				/* Iterazione sugli Elem */
				for (ElemMapType::const_iterator p = ElemData[iCnt1].ElemMap.begin();
					p != ElemData[iCnt1].ElemMap.end();
					p++)
				{
					InitialAssemblyElem *pEl = dynamic_cast<InitialAssemblyElem *>(p->second);
					if (pEl == 0) {
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

					pTmp = const_cast<DofOwner *>(pDOEl->pGetDofOwner());
					ASSERT(pTmp != 0);
					// numero di dof di un owner
					unsigned int iNumDofs = pEl->iGetInitialNumDof();
					pTmp->iNumDofs = iNumDofs;
					if (iNumDofs > 0) {
						pTmp->iFirstIndex = iIndex;
						if (pds) {
							unsigned int nd = iNumDofs;
							integer fd = iIndex;
							ElemWithDofs* pEWD = CastElemWithDofs(p->second);

							silent_cout(psElemNames[pEl->GetElemType()]
								<< "(" << pEl->GetLabel()
								<< "): " << nd << " " << fd + 1);
							if (nd > 1) {
								silent_cout("->" << fd + nd);
							}
							silent_cout(std::endl);
							if (uPrintFlags & PRINT_DOF_DESCRIPTION) {
								pEWD->DescribeDof(std::cout,
										"        ", true);
							}

							if (uPrintFlags & PRINT_EQ_DESCRIPTION) {
								pEWD->DescribeEq(std::cout,
										"        ", true);
							}
						}

						std::vector<std::string> DofDesc;
						pEl->DescribeDof(DofDesc, true);
						if (DofDesc.size() == iNumDofs) {
							for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
								pDofs[iIndex + iCnt].Description = DofDesc[iCnt];
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
								pDofs[iIndex + iCnt].Description = os.str();
							}
						}

						std::vector<std::string> EqDesc;
						pEl->DescribeEq(EqDesc, true);
						if (EqDesc.size() == iNumDofs) {
							for (unsigned iCnt = 0; iCnt < iNumDofs; iCnt++) {
								pDofs[iIndex + iCnt].EqDescription = EqDesc[iCnt];
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
								pDofs[iIndex + iCnt].EqDescription = os.str();
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
	integer iMaxRows = 0;
	integer iMaxCols = 0;

	InitialAssemblyIterator IAIter(&ElemData);
	InitialAssemblyElem* pEl = IAIter.GetFirst();
	while (pEl != NULL) {
		integer iCurrRows = 0;
		integer iCurrCols = 0;
		pEl->InitialWorkSpaceDim(&iCurrRows, &iCurrCols);
		if (iCurrRows > iMaxRows) {
			iMaxRows = iCurrRows;
		}
		if (iCurrCols > iMaxCols) {
			iMaxCols = iCurrCols;
		}
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
	for (StructNode** ppTmpNode = ppFirstNode;
		ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++)
	{
		(*ppTmpNode)->SetInitialValue(X);
	}

	/* Setta i valori iniziali dei gradi di liberta' dei vincoli
	 * durante l'assemblaggio iniziale */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Pre ogni tipo di elemento */
		if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
			ElemData[iCnt1].bToBeUsedInAssembly() &&
			!ElemData[iCnt1].ElemMap.empty())
		{
			for (ElemMapType::const_iterator p = ElemData[iCnt1].ElemMap.begin();
				p != ElemData[iCnt1].ElemMap.end();
				p++)
			{
				ElemWithDofs *pEWD = CastElemWithDofs(p->second);
				pEWD->SetInitialValue(X);
			}
		}
	}

	/* Vettore di lavoro */
	VectorHandler* pResHdl = pSM->pResHdl();
	MySubVectorHandler WorkVec(iMaxRows);

	/* Matrice di lavoro */
	MatrixHandler* pMatHdl = pSM->pMatHdl();
	VariableSubMatrixHandler WorkMat(iMaxRows, iMaxCols);

	/* Soluzione */
	VectorHandler* pSolHdl = pSM->pSolHdl();

	/* Ciclo di assemblaggio */
	for (integer iNumIter = 0; ; iNumIter++) {
		/* Assemblo il residuo */
		pResHdl->Reset();

		/* Contributo dei nodi */
		for (StructNode** ppTmpNode = ppFirstNode;
			ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++)
		{
			if ((*ppTmpNode)->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			integer iFirstIndex = ((*ppTmpNode)->pGetDofOwner())->iFirstIndex;

			/* Nuova feature: ogni nodo ha la sua stiffness */
			doublereal dPosStiff = (*ppTmpNode)->dGetPositionStiffness();
			doublereal dVelStiff = (*ppTmpNode)->dGetVelocityStiffness();

			/* Posizione: k*Delta_x = k(x_0-x) + F */
			Vec3 TmpVec = (*ppTmpNode)->GetXPrev()-(*ppTmpNode)->GetXCurr();
			pResHdl->Add(iFirstIndex+1, TmpVec*dPosStiff);

			/* Rotazione: k*Delta_g = -k*g(R_Delta) + M */
			Mat3x3 R0 = (*ppTmpNode)->GetRPrev();
			Mat3x3 RDelta = (*ppTmpNode)->GetRCurr()*R0.Transpose();
			TmpVec = -Vec3(CGR_Rot::Param, RDelta);
			pResHdl->Add(iFirstIndex+4, TmpVec*dPosStiff);

			/* Velocita': k*Delta_v = k*(v0-Delta_v) + F */
			TmpVec = (*ppTmpNode)->GetVPrev()-(*ppTmpNode)->GetVCurr();
			pResHdl->Add(iFirstIndex+7, TmpVec*dVelStiff);

			/* Velocita' angolare: k*(Delta_w+(R_Delta*w0)/\Delta_g) =
			 *                                    k*(R_Delta*w0-w) + M */
			Vec3 wPrev((*ppTmpNode)->GetWPrev());
			Vec3 wCurr((*ppTmpNode)->GetWCurr());

			if ((*ppTmpNode)->bOmegaRotates()) {
				/* con questa la velocita' angolare e' solidale
				 * con il nodo */
				TmpVec = RDelta*wPrev-wCurr;
			} else {
				/* con questa la velocita' angolare e' solidale
				 * col riferimento assoluto */
				TmpVec = wPrev-wCurr;
			}

			pResHdl->Add(iFirstIndex+10, TmpVec*dVelStiff);
		}

		/* Elementi (con iteratore): */
		pEl = IAIter.GetFirst();
		while (pEl != NULL) {
			try {
				*pResHdl += pEl->InitialAssRes(WorkVec, X);
			}
			catch (Elem::ChangedEquationStructure) {
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

		if (
#ifdef DEBUG
				DEBUG_LEVEL_MATCH(MYDEBUG_ASSEMBLY) ||
#endif /* DEBUG */
				outputIters())
		{
			silent_cout("Iteration(" << iNumIter << ") "
				"" << dTest
				<< " (Tol=" << dInitialAssemblyTol << ")"
				<< std::endl);
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
		pSM->MatrReset();

		/* Contributo dei nodi */
		for (StructNode** ppTmpNode = ppFirstNode;
				ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) 
		{
			if ((*ppTmpNode)->GetStructNodeType() == StructNode::DUMMY) {
				continue;
			}

			integer iFirstIndex = ((*ppTmpNode)->pGetDofOwner())->iFirstIndex;

			/* Nuova feature: ogni nodo ha la sua stiffness */
			doublereal dPosStiff = (*ppTmpNode)->dGetPositionStiffness();
			doublereal dVelStiff = (*ppTmpNode)->dGetVelocityStiffness();

			for (int iCnt = 1; iCnt <= 6; iCnt++) {
				/* Posizione, rotazione */
				integer iTmp = iFirstIndex+iCnt;
				pMatHdl->PutCoef(iTmp, iTmp, dPosStiff);

				/* Velocita', velocita' angolare */
				iTmp += 6;
				pMatHdl->PutCoef(iTmp, iTmp, dVelStiff);
			}

			if ((*ppTmpNode)->bOmegaRotates()) {
				/* con questi la velocita' angolare e' solidale con il nodo */

				/* Velocita' angolare - termine di rotazione: R_Delta*w0/\ */
				Mat3x3 R0 = (*ppTmpNode)->GetRPrev();
				Mat3x3 R = (*ppTmpNode)->GetRCurr();
				Vec3 W0 = (*ppTmpNode)->GetWPrev();
				Vec3 TmpVec = R*(R0.Transpose()*(W0*dVelStiff));

				/* W1 in m(3, 2), -W1 in m(2, 3) */
				doublereal d = TmpVec.dGet(1);
				pMatHdl->PutCoef(iFirstIndex+12, iFirstIndex+5, d);
				pMatHdl->PutCoef(iFirstIndex+11, iFirstIndex+6, -d);

				/* W2 in m(1, 3), -W2 in m(3, 1) */
				d = TmpVec.dGet(2);
				pMatHdl->PutCoef(iFirstIndex+10, iFirstIndex+6, d);
				pMatHdl->PutCoef(iFirstIndex+12, iFirstIndex+4, -d);

				/* W3 in m(2, 1), -W3 in m(1, 2) */
				d = TmpVec.dGet(3);
				pMatHdl->PutCoef(iFirstIndex+11, iFirstIndex+4, d);
				pMatHdl->PutCoef(iFirstIndex+10, iFirstIndex+5, -d);
			} /* altrimenti la velocita' angolare e' solidale con il nodo */
		}

		/* Contributo degli elementi */
		pEl = IAIter.GetFirst();
		while (pEl != NULL) {
			*pMatHdl += pEl->InitialAssJac(WorkMat, X);
			pEl = IAIter.GetNext();
		}

		/* Fattorizza e risolve con jacobiano e residuo appena calcolati */
		try {
			pSM->Solve();
		}
		catch (LinearSolver::ErrFactor err) {
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

		/* Aggiorno la soluzione */
		if (dEpsilon != 1.) {
			*pSolHdl *= dEpsilon;
		}
		X += *pSolHdl;

		/* Correggo i nodi */
		for (StructNode** ppTmpNode = ppFirstNode;
				ppTmpNode < ppFirstNode+iNumNodes; ppTmpNode++) {
			(*ppTmpNode)->InitialUpdate(X);
		}
	}

endofcycle:
	/* Resetta e distrugge la struttura temporanea dei Dof */

	/* Elementi: rimette a posto il numero di Dof propri dei vincoli */
	for (int iCnt1 = 0; iCnt1 < Elem::LASTELEMTYPE; iCnt1++) {
		/* Per ogni tipo di elemento */
		if (ElemData[iCnt1].DofOwnerType != DofOwner::UNKNOWN &&
			ElemData[iCnt1].bToBeUsedInAssembly() &&
			!ElemData[iCnt1].ElemMap.empty())
		{
			/* Se possiede dofs, se deve essere usato nell'assemblaggio
			 * e se ne sono presenti */

			/* Tipo di dof dell'elemento corrente */
			DofOwner::Type CurrDofType = ElemData[iCnt1].DofOwnerType;

			/* Puntatore al primo DofOwner */
			pTmp = DofData[CurrDofType].pFirstDofOwner;

			for (ElemMapType::const_iterator p = ElemData[iCnt1].ElemMap.begin();
				p != ElemData[iCnt1].ElemMap.end();
				p++)
			{
				pTmp->iNumDofs = p->second->iGetNumDof();
			}
		}
	}

	/* Dealloca il vettore dei Dof */
	ASSERT(pDofs != NULL);
	if (pDofs != NULL) {
		SAFEDELETEARR(pDofs);
	}

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
	Node** ppTmpNode = ppNodes;
	for (; ppTmpNode < ppNodes+iTotNodes; ppTmpNode++) {
		DofOwner* pDO = (DofOwner*)(*ppTmpNode)->pGetDofOwner();
		pDO->iNumDofs = (*ppTmpNode)->iGetNumDof();
	}

	/* Setta i DofOwner degli elementi (chi li possiede) */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		DofOwner::Type DT = ElemData[iCnt].DofOwnerType;
		if (DT != DofOwner::UNKNOWN) {
			DEBUGLCOUT(MYDEBUG_INIT, "Elem type " << iCnt
					<< " (" << psElemNames[iCnt] << ")"
					<< std::endl);

			for (ElemMapType::const_iterator p = ElemData[iCnt].ElemMap.begin();
				p != ElemData[iCnt].ElemMap.end();
				p++)
			{
				ElemWithDofs* pEWD = CastElemWithDofs(p->second);

				DEBUGLCOUT(MYDEBUG_INIT, "    " << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << ")" << std::endl);

				DofOwner* pDO = (DofOwner*)pEWD->pGetDofOwner();
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
	for (Node** ppNode = ppNodes; ppNode < ppNodes+iTotNodes; ppNode++) {
		ASSERT(*ppNode != NULL);
		(*ppNode)->SetValue(this, X, XP);
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
#if 0
		int count = 0;
		while (!fp.eof()) {
			char tmp;
			count++;
			fp.get(tmp);
		}
		count--;
		if (count != (X.iGetSize()+XP.iGetSize())*sizeof(double)) {
			silent_cerr("DataManager::SetValue(): "
				"File(" << solArrFileName << ") too short!" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#endif
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

		/* get a pointer to binary NetCDF file  -->  pDM->OutHdl.BinFile */
		NcFile *pBinFile = OutHdl.pGetBinFile();

		/* Add general NetCDF (output) variables to the BinFile object
		 * and save the NcVar* pointer returned from add_var as handle
		 * for later write accesses. Define also variable attributes */
		Var_Step = pBinFile->add_var("run.step", ncLong,
			OutHdl.DimTime(), OutHdl.DimV1());
		if (Var_Step == 0) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_Step->add_att("units", "-")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_Step->add_att("type", "integer")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_Step->add_att("description", "time step index")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		Var_Time = pBinFile->add_var("run.time", ncDouble,
			OutHdl.DimTime(), OutHdl.DimV1());
		if (Var_Time == 0) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* FIXME: needs to be configurable? */
		if (!Var_Time->add_att("units", "s")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_Time->add_att("type", "doublereal")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_Time->add_att("description", "simulation time")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		Var_TimeStep = pBinFile->add_var("run.timestep", ncDouble,
			OutHdl.DimTime(), OutHdl.DimV1());
		if (Var_TimeStep == 0) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* FIXME: needs to be configurable? */
		if (!Var_TimeStep->add_att("units", "s")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_TimeStep->add_att("type", "doublereal")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!Var_TimeStep->add_att("description", "integration time step")) {
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
#endif /* USE_NETCDF */

	/* Dati dei nodi */
	NodeOutputPrepare(OutHdl);

	/* Dati degli elementi */
	ElemOutputPrepare(OutHdl);
}

/* Output dati */
void
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

	/* output only when allowed by the output meter */
	if ((!force) && !pOutputMeter->dGet()) {
		return;
	}

	/*
	 * Write general simulation data to binary NetCDF file
	 *   the current time step index
	 *   the current simulatin time
	 *   the current integration time step
	 */
#ifdef USE_NETCDF
	if (OutHdl.UseNetCDF(OutputHandler::NETCDF)) {
		Var_Step->put_rec(&lStep, OutHdl.GetCurrentStep());
		Var_Time->put_rec(&dTime, OutHdl.GetCurrentStep());
		Var_TimeStep->put_rec(&dTimeStep, OutHdl.GetCurrentStep());
	}
#endif /* USE_NETCDF */

	/* Dati dei nodi */
	NodeOutput(OutHdl);

	/* Dati degli elementi */
	ElemOutput(OutHdl);

#if defined(USE_ADAMS) || defined(USE_MOTIONVIEW)
	iOutputBlock++;
#endif /* defined(USE_ADAMS) || defined(USE_MOTIONVIEW) */

#ifdef USE_ADAMS
	/* Se richiesto, esegue l'output delle condizioni iniziali*/
	if (bAdamsOutput()) {
		AdamsResOutput(iOutputBlock, "DYNAMIC", "MBDyn");
	}
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW
	/* Se richiesto, esegue l'output delle condizioni iniziali*/
	if (bMotionViewOutput()) {
		MotionViewResOutput(iOutputBlock, "DYNAMIC", "MBDyn");
	}
#endif /* USE_MOTIONVIEW */

	OutHdl.IncCurrentStep();
#ifdef USE_NETCDF
	if (bNetCDFsync) {
		OutHdl.pGetBinFile()->sync();
	}
#endif /* USE_NETCDF */
}

/* Output dati */
void
DataManager::Output(const VectorHandler& X, const VectorHandler& XP) const
{
	/* Dati dei nodi */
	NodeOutput((OutputHandler&)OutHdl, X, XP);

	/* Dati degli elementi */
	ElemOutput((OutputHandler&)OutHdl, X, XP);
}

#if 0
/* Output dati pch */
void
DataManager::Output_pch(std::ostream& pch) const
{
	/* Dati dei nodi */
	NodeOutput_pch(pch);

	/* Dati degli elementi */
	ElemOutput_pch(pch);
}

/* Output dati f06 */
void
DataManager::Output_f06(std::ostream& f06, const VectorHandler& X) const
{
	/* Dati dei nodi */
	NodeOutput_f06(f06, X);

	/* Dati degli elementi */
	ElemOutput_f06(f06, X);
}

/* Output dati f06 */
void
DataManager::Output_f06(std::ostream& f06, const VectorHandler& Xr,
		const VectorHandler& Xi) const
{
	/* Dati dei nodi */
	NodeOutput_f06(f06, Xr, Xi);

	/* Dati degli elementi */
	ElemOutput_f06(f06, Xr, Xi);
}
#endif

#if 0
/* OpenDX output (tbi) */
void
DataManager::Output_OpenDX(std::ostream& dx, const VectorHandler& Xr, const VectorHandler& Xi) const
{
	NO_OP;
}
#endif

void
DataManager::BeforePredict(VectorHandler& X, VectorHandler& XP,
	VectorHandler& XPrev,
	VectorHandler& XPPrev) const
{
	Node** ppLastNode = ppNodes+iTotNodes;
	for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
		ASSERT(*ppTmp != NULL);
		(*ppTmp)->BeforePredict(X, XP, XPrev, XPPrev);
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
	for (std::vector<bool>::iterator i = IsConverged.begin();
		i != IsConverged.end();
		i++)
	{
		*i = false;
	}

	Node** ppLastNode = ppNodes+iTotNodes;
	for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
		ASSERT(*ppTmp != NULL);
		(*ppTmp)->AfterPredict(*pXCurr, *pXPrimeCurr);
	}

	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->AfterPredict(*pXCurr, *pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::Update(void) const
{
	Node** ppLastNode = ppNodes+iTotNodes;
	for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
		ASSERT(*ppTmp != NULL);
		(*ppTmp)->Update(*pXCurr, *pXPrimeCurr);
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
	Node** ppLastNode = ppNodes+iTotNodes;
	for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
		ASSERT(*ppTmp != NULL);
		(*ppTmp)->AfterConvergence(*pXCurr,
			*pXPrimeCurr);
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
			((DataManager*)this)->MakeRestart();
		}
		break;

	case TIME: {
		doublereal dT = DrvHdl.dGetTime();
		if (dT - dLastRestartTime >= dRestartTime) {
			dLastRestartTime = dT;
			((DataManager*)this)->MakeRestart();
		}
		break;
	}

	case TIMES: {
		doublereal dT = DrvHdl.dGetTime()
			+ pSolver->GetDInitialTimeStep()/100.;
		if (iCurrRestartTime == iNumRestartTimes) {
			break;
		}

		ASSERT(iCurrRestartTime < iNumRestartTimes);

		if (dT >= pdRestartTimes[iCurrRestartTime]) {
			iCurrRestartTime++;
			((DataManager*)this)->MakeRestart();
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
	Node** ppLastNode = ppNodes+iTotNodes;
	for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
		ASSERT(*ppTmp != NULL);
		(*ppTmp)->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
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
 	silent_cout("Residual(" << iIterCnt << "):"
		<< std::endl);
	integer iSize = Res.iGetSize();
 	for (int iTmpCnt = 1; iTmpCnt <= iSize; iTmpCnt++) {
    		silent_cout("Eq  " << std::setw(8)
			<< iTmpCnt << ": " 
			<< std::setw(20) << Res.dGetCoef(iTmpCnt)
			<< " " << pDofs[iTmpCnt-1].EqDescription
			<< std::endl);
	}
}

void
DataManager::PrintSolution(const VectorHandler& Sol, integer iIterCnt) const
{
 	silent_cout("Solution(" << iIterCnt << "):"
		<< std::endl);
	integer iSize = Sol.iGetSize();
 	for (integer iTmpCnt = 1; iTmpCnt <= iSize; iTmpCnt++) {
    		silent_cout("Dof " << std::setw(8)
			<< iTmpCnt << ": " 
			<< std::setw(20) << Sol.dGetCoef(iTmpCnt)
			<< " " << pDofs[iTmpCnt-1].Description
			<< std::endl);
	}
}

const std::string&
DataManager::GetDofDescription(int i) const
{
	ASSERT(i > 0 && i <= iTotDofs);
	return pDofs[i - 1].Description;
}

const std::string&
DataManager::GetEqDescription(int i) const
{
	ASSERT(i > 0 && i <= iTotDofs);
	return pDofs[i - 1].EqDescription;
}

unsigned
DataManager::ConvergedRegister(void)
{
	unsigned idx = IsConverged.size();
	IsConverged.resize(idx + 1);
	IsConverged[idx] = true;
	return idx;
}

void
DataManager::ConvergedSet(unsigned idx, bool b)
{
	ASSERT(idx < IsConverged.size());

	IsConverged[idx] = b;
}

bool
DataManager::Converged(void) const
{
	for (std::vector<bool>::const_iterator i = IsConverged.begin();
		i != IsConverged.end();
		i++)
	{
		if (!*i) {
			return false;
		}
	}

	return true;
}

/* DataManager - end */
