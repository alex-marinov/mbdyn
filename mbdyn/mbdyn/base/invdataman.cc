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

/* Inverse Dynamics DataManager */

/*
 * Copyright 2007 Alessandro Fumagalli <alessandro.fumagalli@polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */


#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "friction.h"

extern "C" {
#include <strings.h>
#include <time.h>
}

#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */

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
#include <modules.h>

/* To handle  of Elem2Param */
#include "j2p.h"
void 
DataManager::LinkToSolution(const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr,
                    const VectorHandler& XPrimePrimeCurr)
{
	pXCurr = &XCurr;
	pXPrimeCurr = &XPrimeCurr;
	pXPrimePrimeCurr = &XPrimePrimeCurr;
//TODO: DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
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

/* Assemblaggio del residuo */
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
		catch(Elem::ChangedEquationStructure) {
			ResHdl += WorkVec;
			ChangedEqStructure = true;
		}
	}

	if (ChangedEqStructure) {
		throw ChangedEquationStructure();
	}
}

void
DataManager::Update(int iOrder) const
{	
	/* Nodes: */
	switch(iOrder)	{
		case 0:	{
			Node** ppLastNode = ppNodes+iTotNodes;
			for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
				ASSERT(*ppTmp != NULL);
				(*ppTmp)->Update(*pXCurr, iOrder);
			}
			break;
		}
		case 1:	{
			Node** ppLastNode = ppNodes+iTotNodes;
			for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
				ASSERT(*ppTmp != NULL);
				(*ppTmp)->Update(*pXPrimeCurr, iOrder);
			}
			break;
		}
		case 2:	{
			Node** ppLastNode = ppNodes+iTotNodes;
			for (Node** ppTmp = ppNodes; ppTmp < ppLastNode; ppTmp++) {
				ASSERT(*ppTmp != NULL);
				(*ppTmp)->Update(*pXPrimePrimeCurr, iOrder);
			}
			break;
		}
		default:
			NO_OP;
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

