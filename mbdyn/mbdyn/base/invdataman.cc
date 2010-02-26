/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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

/* Inverse Dynamics DataManager */

/*
 * Copyright 2008 Alessandro Fumagalli <alessandro.fumagalli@polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */


#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <strings.h>
#include <time.h>

#include "dataman.h"
#include "friction.h"

#include "solver.h"
#include "invsolver.h"
#include "constltp.h"
#include "dataman_.h"
/* add-ons for math parser */
#include "dofpgin.h"
#include "privpgin.h"
#include "dummypgin.h"
#ifdef USE_TCL
#include "tclpgin.h"
#endif /* USE_TCL */
#include "modelns.h"

/* temporary? */
#include "beam.h"

/* To allow direct loading of modules */
#include "modules.h"

/* To handle  of Elem2Param */
#include "j2p.h"

/* deformable joint */
#include "joint.h"
#include "nestedelem.h"

void 
DataManager::LinkToSolution(const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr,
                    const VectorHandler& XPrimePrimeCurr,
                    const VectorHandler& LambdaCurr)
{
	pXCurr = const_cast<VectorHandler*>(&XCurr);
	pXPrimeCurr = const_cast<VectorHandler*>(&XPrimeCurr);
	pXPrimePrimeCurr = const_cast<VectorHandler*>(&XPrimePrimeCurr);
	pLambdaCurr = &LambdaCurr;
	DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
}
void
DataManager::AssConstrJac(MatrixHandler& JacHdl)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	ASSERT(pWorkMat != NULL);
	ASSERT(Elems.begin() != Elems.end());

	AssConstrJac(JacHdl, ElemIter, *pWorkMat);
}

void
DataManager::AssConstrJac(MatrixHandler& JacHdl,
		VecIter<Elem *> &Iter,
		VariableSubMatrixHandler& WorkMat)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	JacHdl.Reset();

	for (ElemMapType::iterator j = ElemData[Elem::JOINT].ElemMap.begin();
		j != ElemData[Elem::JOINT].ElemMap.end();
		j++)
	{
		JacHdl += j->second->AssJac(WorkMat, *pXCurr);
	}
}

/* Constraint residual assembly: */
void
DataManager::AssConstrRes(VectorHandler& ResHdl, int iOrder) 
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	AssConstrRes(ResHdl, ElemIter, *pWorkVec, iOrder);
}

void
DataManager::AssConstrRes(VectorHandler& ResHdl,
		VecIter<Elem *> &Iter,
		SubVectorHandler& WorkVec, int iOrder)
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	bool ChangedEqStructure(false);
	for (ElemMapType::iterator j = ElemData[Elem::JOINT].ElemMap.begin();
		j != ElemData[Elem::JOINT].ElemMap.end();
		j++)
	{
		try {
			ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, 
						iOrder);
		}
		catch (Elem::ChangedEquationStructure) {
			ResHdl += WorkVec;
			ChangedEqStructure = true;
		}
	}

	if (ChangedEqStructure) {
		throw ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

/* Equilibrium residual assembly, no constraints */
void
DataManager::AssRes(VectorHandler& ResHdl) 
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	AssRes(ResHdl, ElemIter, *pWorkVec);
}

void
DataManager::AssRes(VectorHandler& ResHdl,
		VecIter<Elem *> &Iter,
		SubVectorHandler& WorkVec)
	throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	const Elem::Type ElemType[] = {
		Elem::BODY,
		Elem::BEAM,
		Elem::FORCE,
		
		Elem::LASTELEMTYPE
	};


	bool ChangedEqStructure(false);
	
	for (int et = 0; ElemType[et] != Elem::LASTELEMTYPE; et++) {
		for (ElemMapType::iterator j = ElemData[ElemType[et]].ElemMap.begin();
			j != ElemData[ElemType[et]].ElemMap.end();
			j++)
		{
			try {
				ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
							*pXPrimeCurr, *pXPrimePrimeCurr, 
							-1);
			}
			catch (Elem::ChangedEquationStructure) {
				ResHdl += WorkVec;
				ChangedEqStructure = true;
			}
		}
	}

	for (ElemMapType::iterator j = ElemData[Elem::JOINT].ElemMap.begin();
		j != ElemData[Elem::JOINT].ElemMap.end();
		j++)
	{
		// cast to Joint *
		Joint *pj = dynamic_cast<Joint *>(j->second);
		if (pj == 0) {
			// In case of failure, must be driven...
			NestedElem *pn = dynamic_cast<NestedElem *>(j->second);
			if (pn == 0) {
				// error
			}
			// ... and the driven element must be a joint!
			pj = dynamic_cast<Joint *>(pn->pGetElem());
			if (pj == 0) {
				// error
			}
		}

		switch (pj->GetJointType()) {
		case Joint::DEFORMABLEHINGE:
		case Joint::DEFORMABLEDISPJOINT:
		case Joint::DEFORMABLEJOINT:
			try {
				ResHdl += pj->AssRes(WorkVec, *pXCurr, 
							*pXPrimeCurr, *pXPrimePrimeCurr, 
							-1);
			}
			catch (Elem::ChangedEquationStructure) {
				ResHdl += WorkVec;
				ChangedEqStructure = true;
			}
			break;

		default:
			continue;
		}

		try {
			ResHdl += pj->AssRes(WorkVec, *pXCurr, 
				*pXPrimeCurr, *pXPrimePrimeCurr, 
				-1);
		}
		catch (Elem::ChangedEquationStructure) {
			ResHdl += WorkVec;
			ChangedEqStructure = true;
		}
	}

	if (ChangedEqStructure) {
		throw ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

void
DataManager::Update(int iOrder) const
{
	/* Nodes: */
	switch(iOrder){
	case 0:
		// Update nodes positions
		for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); i++) {
			(*i)->Update(*pXCurr, iOrder);
		}
		break;

	case 1:
		// Update nodes velocities
		for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); i++) {
			(*i)->Update(*pXPrimeCurr, iOrder);
		}
		break;

	case 2:
		// Update nodes accelerations
		for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); i++) {
			(*i)->Update(*pXPrimePrimeCurr, iOrder);
		}
		break;

	case -1:
		// Update constraints reactions (for output only...)*/
		for (ElemMapType::const_iterator j = ElemData[Elem::JOINT].ElemMap.begin();
			j != ElemData[Elem::JOINT].ElemMap.end(); j++)
		{
			j->second->Update(*pLambdaCurr, iOrder);
		}
		break;

	default:
		break;
	}

#if 0
	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->Update(*pXCurr, *pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
#endif
}

void
DataManager::IDAfterConvergence(void) const
{	
	/* Nodes: */
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); i++) {
		(*i)->AfterConvergence(*pXCurr, *pXPrimeCurr, *pXPrimePrimeCurr);
	}

	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->AfterConvergence(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr, *(VectorHandler*)pXPrimePrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}

	/* Restart condizionato */
#if 0
	/* Versione con iteratore: */
	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->Update(*pXCurr, *pXPrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
#endif
}

bool
DataManager::InverseDofOwnerSet(void)
{
	DEBUGCOUTFNAME("DataManager::DofOwnerSet");
	int iNodeTotNumDofs = 0;
	int iJointTotNumDofs = 0;
	
	/* Setta i DofOwner dei nodi */
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); i++) {
		DofOwner* pDO = const_cast<DofOwner *>((*i)->pGetDofOwner());
		pDO->iNumDofs = (*i)->iGetNumDof();
		iNodeTotNumDofs += pDO->iNumDofs;
	}

	for (ElemMapType::iterator j = ElemData[Elem::JOINT].ElemMap.begin();
                j != ElemData[Elem::JOINT].ElemMap.end(); j++)
        {
		ElemWithDofs* pEWD = CastElemWithDofs(j->second);
		iJointTotNumDofs += pEWD->iGetNumDof();

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
	
	if(iNodeTotNumDofs == iJointTotNumDofs)	{
		return true;
	} else	{
		return false;
	}

}


void DataManager::InverseDofInit(bool bIsSquare)
{  
   	if( iTotDofOwners > 0) {	
      
      		/* Di ogni DofOwner setta il primo indice
      		 * e calcola il numero totale di Dof:
		 * poichè per il problema inverso non si 
		 * possono aggiungere incognite diverse
		 * da posizione (velocità e accelerazione)
		 * dei nodi e reazioni vincolari, viene
		 * controllato che gli elementi che aggiungono 
		 * dof siano solo nodi e vincoli.
		 * 
		 * Per problemi mal posti (DoF nodi != Dof vincoli): 
		 * iTotDofs = (DoF nodi) + (Dof vincoli), altrimenti
		 * 
		 * Per problemi ben posti:
		 * iTotDofs = DoF nodi (= DoF vincoli)
		 * */

		
      		/* Mette gli indici ai DofOwner dei nodi strutturali: */
		/* contatore dei Dof dei nodi */
      		integer iNodeIndex = 0;

		NodeMapType::const_iterator i = NodeData[Node::STRUCTURAL].NodeMap.begin();
		for (int iDOm1 = 0; iDOm1 < DofData[DofOwner::STRUCTURALNODE].iNum;
			iDOm1++, i++)
		{
			DofOwner *pDO = const_cast<DofOwner *>(i->second->pGetDofOwner());
      			unsigned iNumDofs = pDO->iNumDofs = i->second->iGetNumDof();
			if (iNumDofs > 0) {
				pDO->iFirstIndex = iNodeIndex;
				iNodeIndex += iNumDofs;

		 	} else {
		    		pDO->iFirstIndex = -1;
		    		DEBUGCERR("warning, item " << (iDOm1 + 1) << " has 0 dofs" << std::endl);
		 	}
      		}
		
		/* Gli indici dei nodi sono ok */
		
		/* Se il problema è ben posto, gli indici delle equazioni di vincolo
		 * hanno numerazione indipendente dai nodi. Altrimenti la numerazione 
		 * è a partire dagli indici dei nodi (per fare spazio alla matrice 
		 * peso nello jacobiano) */

		/* contatore dei Dof dei joint */
      		integer iJointIndex = 0;
		if (!bIsSquare) {
			iJointIndex = iNodeIndex;
		}

		for (ElemMapType::iterator j = ElemData[Elem::JOINT].ElemMap.begin();
			j != ElemData[Elem::JOINT].ElemMap.end();
			j++)
		{
			DofOwner *pDO = const_cast<DofOwner *>((dynamic_cast<ElemWithDofs *>(j->second))->pGetDofOwner());
			if (pDO) {
				unsigned iNumDofs = pDO->iNumDofs;
				if (iNumDofs > 0) {
					pDO->iFirstIndex = iJointIndex;
					iJointIndex += iNumDofs;

		 		} else {
		    			pDO->iFirstIndex = -1;
		    			DEBUGCERR("warning, Joint(" << j->second->GetLabel() << ") "
						"has 0 dofs" << std::endl);
		 		}
			}
				
		}
		
		if (bIsSquare)	{
			iTotDofs = iNodeIndex;
		} else {
			iTotDofs = iJointIndex;
		}

	      	DEBUGLCOUT(MYDEBUG_INIT, "iTotDofs = " << iTotDofs << std::endl);
   	} else {
     	 	DEBUGCERR("");
     	 	silent_cerr("no dof owners are defined" << std::endl);
    	  
      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}	   	
  	 
  	 	
  	/* Crea la struttura dinamica dei Dof */
  	if (iTotDofs > 0) {
		if (bIsSquare)	{
  	    		SAFENEWARRNOFILL(pDofs, Dof, 2*iTotDofs);
  	      	    	/* Inizializza l'iteratore sui Dof */
  	    		DofIter.Init(pDofs, 2*iTotDofs);
  	    
  	   	 	/* Inizializza la struttura dinamica dei Dof */
		 
		 	/*FIXME:*/
	   	 	Dof* pTmp = pDofs;
	   	 	integer iIndex = pDofOwners[0].iFirstIndex;
     				 while(pTmp < pDofs+2*iTotDofs) {
					 pTmp->iIndex = iIndex++;
					 pTmp->Order = DofOrder::DIFFERENTIAL;
					 pTmp++;
      				}		
		} else {
			SAFENEWARRNOFILL(pDofs, Dof, iTotDofs);/* Inizializza l'iteratore sui Dof */
  	    		DofIter.Init(pDofs, iTotDofs);
  	    
  	   	 	/* Inizializza la struttura dinamica dei Dof */
		 
		 	/*FIXME:*/
	   	 	Dof* pTmp = pDofs;
	   	 	integer iIndex = pDofOwners[0].iFirstIndex;
     				 while(pTmp < pDofs+iTotDofs) {
					 pTmp->iIndex = iIndex++;
					 pTmp->Order = DofOrder::DIFFERENTIAL;
					 pTmp++;
      				}
		}
	
  	} else {
		DEBUGCERR("");
      		silent_cerr("no dofs are defined" << std::endl);
     	 
      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   	}	   	
}  

