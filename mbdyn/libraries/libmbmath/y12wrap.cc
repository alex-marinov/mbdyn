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

#include <y12wrap.h>
#include <y12lib.h>

/* Y12LUSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
Y12LUSolver::Y12LUSolver(integer iMatOrd, integer iSize,
			 integer** ppiTmpRow, integer** ppiTmpCol,
			 doublereal** ppdTmpMat,
			 doublereal* pdTmpRhs, integer iPivotParam)
: iMatSize(iSize),
ppiRow(ppiTmpRow),
ppiCol(ppiTmpCol),
ppdMat(ppdTmpMat),
iN(iMatOrd),
iNonZeroes(0),
pdRhs(pdTmpRhs),
piHA(NULL),
pdPIVOT(NULL),
iFirstSol(-1)
{
	ASSERT(iMatSize > 0);
	ASSERT(ppiRow != NULL);
	ASSERT(ppiCol != NULL);
	ASSERT(ppdMat != NULL);
	ASSERT(*ppiRow != NULL);
	ASSERT(*ppiCol != NULL);
	ASSERT(*ppdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(SMmm.fIsValid((void*)*ppiRow, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppiCol, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppdMat, iMatSize*sizeof(doublereal)));
	ASSERT(SMmm.fIsValid((void*)pdRhs, iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

	SAFENEWARR(piHA, integer, 11*iN, LUmm);
	SAFENEWARR(pdPIVOT, doublereal, iN, LUmm);
	
#ifdef DEBUG
	for (int iCnt = 0; iCnt < 11*iN; iCnt++) {
		piHA[iCnt] = 0;
	}
	for (int iCnt = 0; iCnt < iN; iCnt++) {
		pdPIVOT[iCnt] = 0;
	}
#endif /* DEBUG */

	iIFLAG[2] = iPivotParam;
}

/* Distruttore */
Y12LUSolver::~Y12LUSolver(void)
{
	if (pdPIVOT != NULL) {
		SAFEDELETEARR(pdPIVOT, LUmm);
	}
	if (piHA != NULL) {
		SAFEDELETEARR(piHA, LUmm);
	}
}

void 
Y12LUSolver::IsValid(void) const
{
	ASSERT(iMatSize > 0);
	ASSERT(ppiRow != NULL);
	ASSERT(ppiCol != NULL);
	ASSERT(ppdMat != NULL);
	ASSERT(*ppiRow != NULL);
	ASSERT(*ppiCol != NULL);
	ASSERT(*ppdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0); 
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(SMmm.fIsValid((void*)*ppiRow, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppiCol, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppdMat, iMatSize*sizeof(doublereal)));
	
	ASSERT(SMmm.fIsValid((void*)pdRhs, iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

	ASSERT(piHA != NULL);
	ASSERT(piPIVOT != NULL);
	
#ifdef DEBUG_MEMMANAGER
	ASSERT(LUmm.fIsBlock((void*)piHA, 11*iN*sizeof(integer)));
	ASSERT(LUmm.fIsBlock((void*)pdPIVOT, 1*iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}

/* Fattorizza la matrice */
flag
Y12LUSolver::fLUFactor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	ASSERT(iNonZeroes > 0);
	
	/* Sets parameters */
	integer iIFAIL = 0;
	
	iIFLAG[0] = 0;
	iIFLAG[1] = 3;
	/* iIFLAG[2] = iPivotParam; */
	iIFLAG[3] = 0;
	iIFLAG[4] = 2;
	
	iFirstSol = 1;
	
	/* Prepares for fatorization */
	__FC_DECL__(y12mbf)(&iN, &iNonZeroes, *ppdMat,
			    &iMatSize, *ppiCol,
			    &iMatSize, *ppiRow,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);
			    
	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during factorisation, code "
			<< iIFAIL << endl;
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	/* actual factorization */
	__FC_DECL__(y12mcf)(&iN, &iNonZeroes, *ppdMat,
			    &iMatSize, *ppiCol,
			    &iMatSize, *ppiRow,
			    pdPIVOT, pdRhs,
			    piHA, &iN,
			    dAFLAG, iIFLAG, &iIFAIL);

	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during factorisation, code "
			<< iIFAIL << endl;
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	return iIFAIL;
}

/* Risolve */
void
Y12LUSolver::Solve(void)
{
#ifdef DEBUG
	IsValid();
#endif

	integer iIFAIL = 0;
	
	__FC_DECL__(y12mdf)(&iN, *ppdMat, &iMatSize, pdRhs,
			    pdPIVOT, *ppiCol,
			    piHA, &iN,
			    iIFLAG, &iIFAIL);
	
	if (iIFAIL != 0) {
		cerr << "Y12LUSolver: error during factorisation, code "
			<< iIFAIL << endl;
		THROW(Y12LUSolver::ErrFactorisation(iIFAIL));
	}
	
	if (iFirstSol == 1) {
		iIFLAG[4] = 3;
		iFirstSol = 0;
	}
}

/* Y12LUSolver - end */


/* Y12SparseLUSolutionManager - begin: code */

/* Costruttore */
Y12SparseLUSolutionManager::Y12SparseLUSolutionManager(integer iSize, 
						       integer iWorkSpaceSize,
						       const doublereal& dPivotFactor) :
iMatMaxSize(iSize),
iMatSize(iSize), 
piRow(NULL),
piCol(NULL), 
pdMat(NULL),
pdVec(NULL),
pMH(NULL),
pVH(NULL),
pLU(NULL),
fHasBeenReset(1)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
      		iWorkSpaceSize = iSize*iSize;
   	}

	integer iPivot;
	if (dPivotFactor == 0.) {
		iPivot = 0;
	} else {
		iPivot = 1;
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
			       Y12LUSolver,
			       Y12LUSolver(iMatSize, iWorkSpaceSize,
			       		   &piRow, &piCol, 
					   &pdMat, pdVec, iPivot),
			       SMmm);
   
#ifdef DEBUG
   	IsValid();
#endif      
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
Y12SparseLUSolutionManager::~Y12SparseLUSolutionManager(void)
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
Y12SparseLUSolutionManager::IsValid(void) const
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
Y12SparseLUSolutionManager::PacVec(void)
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	ASSERT(pMH->iCurSize > 0);   
   	ASSERT(fHasBeenReset == 1);
   
   	pLU->iNonZeroes = pMH->iPacVec();
}

/* Inizializza il gestore delle matrici */
void
Y12SparseLUSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	pMH->Init(dResetVal);
   	fHasBeenReset = flag(1);
}

/* Risolve il problema */
void
Y12SparseLUSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif

   	if (fHasBeenReset == 1) {
      		this->PacVec();
      		fHasBeenReset = flag(0);
      		flag fReturnFlag = pLU->fLUFactor();
      		if (fReturnFlag < 0) {	 
	 		THROW(Y12SparseLUSolutionManager::ErrGeneric());
      		}
   	}
   	pLU->Solve();
}

/* Y12SparseLUSolutionManager - end */


