/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_HARWELL

#include <harwrap.h>

/* HarwellLUSolver - begin */

#ifdef DEBUG_MEMMANAGER
clMemMan LUmm(sLUClassName);
#endif

/* HarwellLUSolver - end */


/* HarwellSparseLUSolutionManager - begin: code */

#ifdef DEBUG_MEMMANAGER
clMemMan SMmm("HarwellSparseLUSolutionManager");
#endif

/* Costruttore */
HarwellSparseLUSolutionManager::HarwellSparseLUSolutionManager(integer iSize, 
					 		       integer iWorkSpaceSize,
							       const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
piRow(NULL), piCol(NULL), 
pdMat(NULL), pdVec(NULL),
pMH(NULL), pVH(NULL), pLU(NULL),
fHasBeenReset(1)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
      		iWorkSpaceSize = iSize*iSize;
   	}
      
   	/* Alloca arrays */
   	SAFENEWARR(piRow, integer, iWorkSpaceSize, SMmm);
   	SAFENEWARR(piCol, integer, iWorkSpaceSize, SMmm);
   	SAFENEWARR(pdMat, doublereal, iWorkSpaceSize, SMmm);
   	SAFENEWARR(pdVec, doublereal, iMatSize, SMmm);
   
   	/* Alloca handlers ecc. */
   	SAFENEWWITHCONSTRUCTOR(pMH, 
			       SparseMatrixHandler,
			       SparseMatrixHandler(iMatSize, &piRow, 
			       			   &piCol, &pdMat,
			       			   iWorkSpaceSize), 
			       SMmm);
   	SAFENEWWITHCONSTRUCTOR(pVH,
			       MyVectorHandler,
			       MyVectorHandler(iMatSize, pdVec),
			       SMmm);
   	SAFENEWWITHCONSTRUCTOR(pLU, 
			       HarwellLUSolver,
			       HarwellLUSolver(iMatSize, iWorkSpaceSize,
			       		       &piRow, &piCol, 
					       &pdMat, pdVec, dPivotFactor),
			       SMmm);
   
#ifdef DEBUG
   	IsValid();
#endif      
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
HarwellSparseLUSolutionManager::~HarwellSparseLUSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif      
   
   	/* Dealloca oggetti strani */
   	if (pLU != NULL) {	
      		SAFEDELETE(pLU, SMmm);
   	}
   	if (pVH != NULL) {      
      		SAFEDELETE(pVH, SMmm);
   	}
   	if (pMH != NULL) {
      		SAFEDELETE(pMH, SMmm);
   	}
   
   	/* Dealloca arrays */
   	if (pdVec != NULL) {	
      		SAFEDELETEARR(pdVec, SMmm);
   	}
   	if (pdMat != NULL) {	
      		SAFEDELETEARR(pdMat, SMmm);
   	}
   	if (piCol != NULL) {	
      		SAFEDELETEARR(piCol, SMmm);
   	}
   	if (piRow != NULL) {	
      		SAFEDELETEARR(piRow, SMmm);
   	}
}

/* Test di validita' del manager */
void 
HarwellSparseLUSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatMaxSize > 0);
   	ASSERT(iMatSize > 0);
   	ASSERT(pMH != NULL);
   	ASSERT(pdMat != NULL);
   	ASSERT(piRow != NULL);
   	ASSERT(piCol != NULL);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(SMmm.fIsPointerToBlock((void*)piRow));
   	ASSERT(SMmm.fIsPointerToBlock((void*)piCol));
   	ASSERT(SMmm.fIsPointerToBlock((void*)pdMat));
   	ASSERT(SMmm.fIsPointerToBlock((void*)pMH));
   	ASSERT(SMmm.fIsPointerToBlock((void*)pVH));
   	ASSERT(SMmm.fIsPointerToBlock((void*)pLU));
#endif
   
   	ASSERT((pMH->IsValid(), 1));
   	ASSERT((pVH->IsValid(), 1));
   	ASSERT((pLU->IsValid(), 1));
}

/* Prepara i vettori e la matrice per il solutore */
void
HarwellSparseLUSolutionManager::PacVec(void)
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	ASSERT(fHasBeenReset == 1);
   
   	pLU->iNonZeroes = pMH->iPacVec();
}

/* Inizializza il gestore delle matrici */
void
HarwellSparseLUSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	pMH->Init(dResetVal);
   	fHasBeenReset = flag(1);
}

/* Risolve il problema */
void
HarwellSparseLUSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif

   	if (fHasBeenReset == 1) {
      		this->PacVec();
      		fHasBeenReset = flag(0);
      		flag fReturnFlag = pLU->fLUFactor();
      		if (fReturnFlag < 0) {	 
	 		THROW(HarwellSparseLUSolutionManager::ErrGeneric());
      		}
   	}
   	pLU->Solve();
}

/* HarwellSparseLUSolutionManager - end */

#endif /* USE_HARWELL */

