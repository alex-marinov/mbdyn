/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_MESCHACH

#include <iostream>
#include <iomanip>

#include "mschwrap.h"

/* MeschachVectorHandler - begin */

MeschachVectorHandler::MeschachVectorHandler(integer iSize)
: pv(VNULL),
pdVecm1(NULL)
{
	/* Note: MeschachVectorHandler owns its workspace memory */
   	if (iSize > 0) {
      		pv = v_get(iSize);
      		if (pv == VNULL) {
	 		silent_cerr("out of memory?" << std::endl);
	 		throw ErrMemory(MBDYN_EXCEPT_ARGS);
      		}
		pdVecm1 = pv->ve - 1;
   	}
}

MeschachVectorHandler::~MeschachVectorHandler(void) 
{
	/* Note: MeschachVectorHandler owns its workspace memory */
   	if (pv != VNULL) {
      		if (v_free(pv) != 0) {
			/* FIXME: hanlde error */
		}
   	}
}

#ifdef DEBUG
/* Usata per il debug */
void
MeschachVectorHandler::IsValid(void) const 
{
   	ASSERT(pv != VNULL);
   	ASSERT(pv->max_dim > 0);
   	ASSERT(pv->dim > 0);
   	ASSERT(pv->max_dim >= pv->dim);
	ASSERT(pv->ve != NULL);
	ASSERT(pdVecm1 != NULL);
	ASSERT(pdVecm1 == pv->ve - 1);
}
#endif /* DEBUG */

void
MeschachVectorHandler::Resize(integer iNewSize) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* Note: MeschachVectorHandler owns its workspace memory */
   	VEC* p = v_resize(pv, iNewSize);
   	if (p == VNULL) {
      		silent_cerr("out of memory?" << std::endl);
      		throw ErrMemory(MBDYN_EXCEPT_ARGS);
   	}
   	pv = p;
	pdVecm1 = pv->ve - 1;
}

void
MeschachVectorHandler::Reset(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
	v_zero(pv);
}

/* MeschachVectorHandler - end */

/* MeschachSparseMatrixHandler - begin */

MeschachSparseMatrixHandler::MeschachSparseMatrixHandler(integer m, 
							 integer n,
							 integer maxlen) 
: mat(SMNULL)
{
   	if (m > 0 && n > 0) {
      		Create(m, n, maxlen);
   	}
}

MeschachSparseMatrixHandler::~MeschachSparseMatrixHandler(void) 
{
   	if (mat != SMNULL) {
	
		/*
		 * Note: MeschachSparseMatrixHandler owns its workspace memory
		 */
      		if (sp_free(mat)) {
	 		/* FIXME: handle error */
      		}
   	}
}

/* costruisce la matrice */
void
MeschachSparseMatrixHandler::Create(integer m,
				    integer n,
				    integer maxlen) 
{
	/*
	 * Note: MeschachSparseMatrixHandler owns its workspace memory 
	 */
   	if (mat != SMNULL) {
      		mat = sp_resize(mat, m, n);
   	} else {
		/* FIXME: in case maxlen == 0, use n */
      		mat = sp_get(m, n, maxlen ? maxlen : n);
   	}
}

#ifdef DEBUG
void
MeschachSparseMatrixHandler::IsValid(void) const 
{
   	ASSERT(mat != SMNULL);
}
#endif /* DEBUG */

void
MeschachSparseMatrixHandler::Reset(void)
{
   	sp_zero(mat);
}

/* MeschachSparseMatrixHandler - end */

/* MeschachSparseSolutionManager - begin */

void
MeschachSparseSolutionManager::Create(unsigned integer iSize,
					unsigned integer iMaxSize) 
{
   	if (prhs == NULL) {
      		SAFENEWWITHCONSTRUCTOR(prhs, 
			     	       MeschachVectorHandler,
				       MeschachVectorHandler(iSize));
   	} else {
      		prhs->Resize(iSize);
   	}
   
   	if (pivot == PNULL || pivot->size < iSize) {
      		PERM* p = px_resize(pivot, iSize);
      		if (p == PNULL) {
	 		silent_cerr("out of memory?" << std::endl);
	 		throw ErrMemory(MBDYN_EXCEPT_ARGS);
      		}
      		pivot = p;
   	}
   
   	if (pmh != NULL
            && (pmh->iGetNumRows() < (integer)iSize
	        || pmh->iGetNumCols() < (integer)iSize)) {
      		SAFEDELETE(pmh);
   	}
	
   	if (pmh == NULL) {
      		SAFENEWWITHCONSTRUCTOR(pmh,
				       MeschachSparseMatrixHandler,
				       MeschachSparseMatrixHandler(iSize,
				       				   iSize, 
								   iMaxSize));
   	}
}

void
MeschachSparseSolutionManager::Factor(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	spLUfactor(pmh->pGetMAT(), pivot, alpha);
   	fStatus = FACTORED;
}

MeschachSparseSolutionManager::MeschachSparseSolutionManager(integer iSize,
		integer iMaxSize, 
		const doublereal& a)
: prhs(NULL), pivot(PNULL), pmh(NULL), fStatus(RESET), alpha (a)
{
   	Create(iSize, iMaxSize);
   	MatrReset();
}

MeschachSparseSolutionManager::~MeschachSparseSolutionManager(void) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	if (prhs != NULL) {
      		SAFEDELETE(prhs);
   	}
	
   	if (pivot != PNULL) {
      		px_free(pivot);
   	}
	
   	if (pmh != NULL) {
      		SAFEDELETE(pmh);
   	}
}

#ifdef DEBUG
void
MeschachSparseSolutionManager::IsValid(void) const 
{
   	ASSERT(prhs != NULL);
   	prhs->IsValid();
   	ASSERT(pivot != PNULL);
   	ASSERT(pmh != NULL);
   	pmh->IsValid();
   	ASSERT(prhs->iGetSize() >= pmh->iGetNumCols());
   	ASSERT(pivot->size >= pmh->iGetNumCols());
   	ASSERT(pivot->size >= pmh->iGetNumRows());
}
#endif /* DEBUG */

void
MeschachSparseSolutionManager::MatrReset(void) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	fStatus = RESET;
	/* FIXME: TOTALLY UNTESTED */
	pLS->Reset();
}

void
MeschachSparseSolutionManager::Solve(void) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
      
   	if (fStatus == RESET) {
		/* Factor() is in charge of switching fStatus to FACTORED */
      		Factor();
   	}
   
   	spLUsolve(pmh->pGetMAT(),
		  pivot,
		  prhs->pGetMeschachVEC(),
		  prhs->pGetMeschachVEC());
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler *
MeschachSparseSolutionManager::pMatHdl(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	return pmh;
}

/* Rende disponibile l'handler per il termine noto */
VectorHandler *
MeschachSparseSolutionManager::pResHdl(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	return prhs;
}

/*
 * Rende disponibile l'handler per la soluzione (e' lo stesso
 * del termine noto, ma concettualmente sono separati)
 */
VectorHandler *
MeschachSparseSolutionManager::pSolHdl(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	return prhs;
}

std::ostream&
operator << (std::ostream& out, const MeschachSparseMatrixHandler& MH)
{
   	SPMAT* p = MH.pGetMAT();
   	for (integer i = 0; i < p->m; i++) {
      		for (integer j = 0; j < p->n; j++) {
	 		silent_cout(std::setw(16) << sp_get_val(p, i, j));
      		}
      		silent_cout(std::endl);
   	}
   
   	return out;
}

#endif /* USE_MESCHACH */

