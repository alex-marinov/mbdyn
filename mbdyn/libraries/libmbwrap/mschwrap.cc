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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_MESCHACH

#include <mschwrap.h>

/* MeschachVectorHandler - begin */

MeschachVectorHandler::MeschachVectorHandler(int iSize)
: pv(VNULL),
pdVecm1(NULL)
{
	/* Note: MeschachVectorHandler owns its workspace memory */
   	if (iSize > 0) {
      		pv = v_get(iSize);
      		if (pv == VNULL) {
	 		cerr << "out of memory?" << endl;
	 		THROW(ErrMemory());
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

void
MeschachVectorHandler::Resize(integer iNewSize) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* Note: MeschachVectorHandler owns its workspace memory */
   	VEC* p = v_resize(pv, iNewSize);
   	if (p == VNULL) {
      		cerr << "out of memory?" << endl;
      		THROW(ErrMemory());
   	}
   	pv = p;
	pdVecm1 = pv->ve - 1;
}

void
MeschachVectorHandler::Reset(doublereal dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   	if (dResetVal == 0.) {
      		v_zero(pv);
   	} else {
      		v_init(pv, dResetVal);
   	}
}

/* MeschachVectorHandler - end */

/* MeschachSparseMatrixHandler - begin */

MeschachSparseMatrixHandler::MeschachSparseMatrixHandler(int m, 
							 int n,
							 int maxlen) 
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
MeschachSparseMatrixHandler::Create(unsigned int m,
				    unsigned int n,
				    unsigned int maxlen) 
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

void
MeschachSparseMatrixHandler::IsValid(void) const 
{
   	ASSERT(mat != SMNULL);
}

void
MeschachSparseMatrixHandler::Init(const doublereal& d)
{
   	/* FIXME: d is always assumed 0. */
   	sp_zero(mat);
}

/* MeschachSparseMatrixHandler - end */

/* MeschachSparseLUSolutionManager - begin */

void
MeschachSparseLUSolutionManager::Create(unsigned int iSize,
					unsigned int iMaxSize) 
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
	 		cerr << "out of memory?" << endl;
	 		THROW(ErrMemory());
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
MeschachSparseLUSolutionManager::Factor(void)
{
#ifdef DEBUG
   	IsValid();
#endif
   
   	spLUfactor(pmh->pGetMAT(), pivot, alpha);
   	fStatus = FACTORED;
}

MeschachSparseLUSolutionManager::MeschachSparseLUSolutionManager(int iSize,
								 int iMaxSize, 
								 double a)
: prhs(NULL), pivot(PNULL), pmh(NULL), fStatus(RESET), alpha (a)
{
   	Create(iSize, iMaxSize);
   	MatrInit(0.);
}

MeschachSparseLUSolutionManager::~MeschachSparseLUSolutionManager(void) 
{
#ifdef DEBUG
   	IsValid();
#endif
   
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

void
MeschachSparseLUSolutionManager::IsValid(void) const 
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

void
MeschachSparseLUSolutionManager::MatrInit(const doublereal& d) 
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	fStatus = RESET;
   	pmh->Init(d);
}

void
MeschachSparseLUSolutionManager::Solve(void) 
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
MeschachSparseLUSolutionManager::pMatHdl(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	return pmh;
}

/* Rende disponibile l'handler per il termine noto */
VectorHandler *
MeschachSparseLUSolutionManager::pResHdl(void) const
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
MeschachSparseLUSolutionManager::pSolHdl(void) const
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	return prhs;
}

ostream&
operator << (ostream& out, const MeschachSparseMatrixHandler& MH)
{
   	SPMAT* p = MH.pGetMAT();
   	for (int i = 0; i < p->m; i++) {
      		for (int j = 0; j < p->n; j++) {
	 		cout << setw(16) << sp_get_val(p, i, j);
      		}
      		cout << endl;
   	}
   
   	return out;
}

#endif /* USE_MESCHACH */

