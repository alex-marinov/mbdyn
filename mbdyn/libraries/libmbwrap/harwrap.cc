/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/*****************************************************************************
 *                                                                           *
 *                            HARWLIB C++ WRAPPER                            *
 *                                                                           *
 *****************************************************************************/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_HARWELL

#include "harwrap.h"

/* HarwellSparseSolutionManager - begin: code */

/* Costruttore */
HarwellSparseSolutionManager::HarwellSparseSolutionManager(integer iSize, 
					 		       integer iWorkSpaceSize,
							       const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
iColStart(iSize + 1),
MH(iSize),
pVH(NULL), pLU(NULL),
fHasBeenReset(1)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
      		iWorkSpaceSize = iSize*iSize;
   	}
      
   	/* Alloca arrays */
   	dVec.resize(iMatSize,0.);
   
   	/* Alloca handlers ecc. */
   	SAFENEWWITHCONSTRUCTOR(pVH,
			       MyVectorHandler,
			       MyVectorHandler(iMatSize, &(dVec[0])));
	iRow.reserve(iWorkSpaceSize);
	iCol.reserve(iWorkSpaceSize);
	dMat.reserve(iWorkSpaceSize);
   	SAFENEWWITHCONSTRUCTOR(pLU, 
			       HarwellSolver,
			       HarwellSolver(iMatSize, iWorkSpaceSize,
			       		       &iRow, &iCol, 
					       &dMat, &(dVec[0]), dPivotFactor));
   
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
HarwellSparseSolutionManager::~HarwellSparseSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	/* Dealloca oggetti strani */
   	if (pLU != NULL) {	
      		SAFEDELETE(pLU);
   	}
   	if (pVH != NULL) {      
      		SAFEDELETE(pVH);
   	}
}

#ifdef DEBUG
/* Test di validita' del manager */
void 
HarwellSparseSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatMaxSize > 0);
   	ASSERT(iMatSize > 0);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pVH));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLU));
#endif /* DEBUG_MEMMANAGER */
   
   	ASSERT((pVH->IsValid(), 1));
   	ASSERT((pLU->IsValid(), 1));
}
#endif /* DEBUG */

/* Prepara i vettori e la matrice per il solutore */
void
HarwellSparseSolutionManager::PacVec(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	ASSERT(fHasBeenReset == 1);
   
   	pLU->iNonZeroes = MH.MakeIndexForm(dMat, iRow, iCol, iColStart, 1);
}

/* Inizializza il gestore delle matrici */
void
HarwellSparseSolutionManager::MatrReset()
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
	/* FIXME: TOTALLY UNTESTED */
	pLS->Reset();
   	fHasBeenReset = flag(1);
}

/* Risolve il problema */
void
HarwellSparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

   	if (fHasBeenReset == 1) {
      		PacVec();
      		fHasBeenReset = flag(0);
      		if (!pLU->bLUFactor()) {
	 		throw HarwellSparseSolutionManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
      		}
   	}
   	pLU->Solve();
}

/* HarwellSparseSolutionManager - end */

#endif /* USE_HARWELL */

