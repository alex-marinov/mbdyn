/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
 *                          SuperLU C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_SUPERLU_MT

#include <signal.h>

#include "spmh.h"
#include "spmapmh.h"
#include "dirccmh.h"
#include "ccmh.h"

extern "C" {
#include <pdsp_defs.h>
#include <util.h>
extern void *(*pdgstrf_thread)(void *);
extern void pdgstrf_finalize(pdgstrf_options_t *, SuperMatrix*);
}

#include "superluwrap.h"


struct SuperLUSolverData {
	SuperMatrix		A,
				AC,
				L,
				U,
				B;
	Gstat_t			Gstat;
	std::vector<int>	perm_c, /* column permutation vector */
				perm_r; /* row permutations from partial pivoting */
	pdgstrf_options_t	pdgstrf_options;
	pxgstrf_shared_t	pxgstrf_shared;
	pdgstrf_threadarg_t	*pdgstrf_threadarg;
};

/* SuperLUSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
SuperLUSolver::SuperLUSolver(unsigned nt, integer iMatOrd,
			 std::vector<integer>*const piTmpRow, 
			 std::vector<integer>*const piTmpCol,
			 std::vector<doublereal>*const pdTmpMat,
			 doublereal* pdTmpRhs)
: Aip(piTmpRow),
App(piTmpCol),
Axp(pdTmpMat),
iN(iMatOrd),
iNonZeroes(0),
bFirstSol(true),
sld(0),
nThreads(nt),
thread_data(0)
{
	LinearSolver::pdRhs = pdTmpRhs;
	LinearSolver::pdSol = pdTmpRhs;

	ASSERT(Aip != NULL);
	ASSERT(App != NULL);
	ASSERT(Axp != NULL);
	ASSERT(pdTmpRhs != NULL);
	ASSERT(iN > 0);

	SAFENEW(sld, SuperLUSolverData);

	pthread_mutex_init(&thread_mutex, NULL);
	pthread_cond_init(&thread_cond, NULL);

	/* Set up the dense matrix data structure for B. */
	dCreate_Dense_Matrix(&sld->B, iN, 1, pdTmpRhs, iN, DN, _D, GE);


	SAFENEWARR(thread_data, thread_data_t, nThreads);

	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].pSLUS = this;
		thread_data[t].threadNumber = t;

		if (t == 0) {
			continue;
		}

		sem_init(&thread_data[t].sem, 0, 0);
		if (pthread_create(&thread_data[t].thread, NULL, thread_op, &thread_data[t]) != 0) {
			silent_cerr("SuperLUSolver: pthread_create() failed "
					"for thread " << t
					<< " of " << nThreads << std::endl);
			THROW(ErrGeneric());
		}
	}
}

/* Distruttore */
SuperLUSolver::~SuperLUSolver(void)
{
	/* ------------------------------------------------------------
	 * Clean up and free storage after multithreaded factorization.
	 * ------------------------------------------------------------*/
	pdgstrf_thread_finalize(sld->pdgstrf_threadarg, &sld->pxgstrf_shared, 
			&sld->A, &sld->perm_r[0], &sld->L, &sld->U);

	/* ------------------------------------------------------------
	 * Deallocate storage after factorization.
	 * ------------------------------------------------------------*/
	pdgstrf_finalize(&sld->pdgstrf_options, &sld->AC);
	
	thread_operation = SuperLUSolver::EXIT;
	thread_count = nThreads - 1;

	for (unsigned i = 1; i < nThreads; i++) {
		sem_post(&thread_data[i].sem);
	}

	/* thread cleanup func? */

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);

	pthread_mutex_destroy(&thread_mutex);
	pthread_cond_destroy(&thread_cond);

	/* other cleanup... */
}

void *
SuperLUSolver::thread_op(void *arg)
{
	thread_data_t *td = (thread_data_t *)arg;

	/* deal with signals ... */
	sigset_t newset /* , oldset */ ;
	sigemptyset(&newset);
	sigaddset(&newset, SIGTERM);
	sigaddset(&newset, SIGINT);
	sigaddset(&newset, SIGHUP);
	pthread_sigmask(SIG_BLOCK, &newset, /* &oldset */ NULL);

	bool bKeepGoing(true);

	while (bKeepGoing) {
		sem_wait(&td->sem);

		switch (td->pSLUS->thread_operation) {
		case SuperLUSolver::FACTOR:
			(void)pdgstrf_thread(td->pdgstrf_threadarg);
			break;

		case SuperLUSolver::EXIT:
			bKeepGoing = false;
			break;

		default:
			silent_cerr("SuperLUSolver: unhandled op"
					<< std::endl);
			THROW(ErrGeneric());
		}

		td->pSLUS->EndOfOp();
	}

	pthread_exit(NULL);
}

void
SuperLUSolver::EndOfOp(void)
{
	bool last;

	/* decrement the thread counter */
	pthread_mutex_lock(&thread_mutex);
	thread_count--;
	last = (thread_count == 0);

	/* if last thread, signal to restart */
	if (last) {
		pthread_cond_signal(&thread_cond);
		pthread_mutex_unlock(&thread_mutex);
	}
}

#ifdef DEBUG
void 
SuperLUSolver::IsValid(void) const
{
	ASSERT(Aip != NULL);
	ASSERT(App != NULL);
	ASSERT(Axp != NULL);
	ASSERT(iN > 0); 
}
#endif /* DEBUG */

/* Fattorizza la matrice */
void
SuperLUSolver::Factor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	ASSERT(iNonZeroes > 0);

	yes_no_t	refact = YES,
			usepr = NO;
	doublereal	u = 1.0,
			drop_tol = 0.0;
	void		*work = NULL;
	int		info = 0, lwork = 0;
	int		panel_size = sp_ienv(1),
			relax = sp_ienv(2);

	if (bFirstSol) {
		/* ------------------------------------------------------------
		 * Allocate storage and initialize statistics variables. 
		 * ------------------------------------------------------------*/
		refact = NO;

		/* Set up the sparse matrix data structure for A. */
		dCreate_CompCol_Matrix(&sld->A, iN, iN, iNonZeroes,
				&(*Axp)[0], &(*Aip)[0], &(*App)[0],
				NC, _D, GE);

		/* ------------------------------------------------------------
		 * Get column permutation vector perm_c[], according to permc_spec:
		 * permc_spec = 0: use the natural ordering 
		 * permc_spec = 1: use minimum degree ordering on structure of A'*A
		 * permc_spec = 2: use minimum degree ordering on structure of A'+A
		 * !!! permc_spec = 3: use approximate minimum degree column order !!!
		 * ------------------------------------------------------------*/
		int	permc_spec = 0;

		sld->perm_c.resize(iN);
		sld->perm_r.resize(iN);
		get_perm_c(permc_spec, &sld->A, &sld->perm_c[0]);

		StatAlloc(iN, nThreads, panel_size, relax, &sld->Gstat);
	}

	StatInit(iN, nThreads, &sld->Gstat);

	/* ------------------------------------------------------------
	 * Initialize the option structure pdgstrf_options using the
	 * user-input parameters;
	 * Apply perm_c to the columns of original A to form AC.
	 * ------------------------------------------------------------*/
	pdgstrf_init(nThreads, refact, panel_size, relax,
			u, usepr, drop_tol,
			&sld->perm_c[0], &sld->perm_r[0],
			work, lwork, &sld->A, &sld->AC,
			&sld->pdgstrf_options, &sld->Gstat);

	if (bFirstSol) {
		bFirstSol = false;

		/* --------------------------------------------------------------
		 * Initializes the parallel data structures for pdgstrf_thread().
		 * --------------------------------------------------------------*/
		sld->pdgstrf_threadarg = pdgstrf_thread_init(&sld->AC,
				&sld->L, &sld->U, &sld->pdgstrf_options,
				&sld->pxgstrf_shared, &sld->Gstat, &info);
		for (unsigned t = 0; t < nThreads; t++) {
			thread_data[t].pdgstrf_threadarg = &sld->pdgstrf_threadarg[t];
		}
	}

	thread_operation = SuperLUSolver::FACTOR;
	thread_count = nThreads - 1;

	for (unsigned t = 1; t < nThreads; t++) {
		sem_post(&thread_data[t].sem);
	}

	(void)pdgstrf_thread(thread_data[0].pdgstrf_threadarg);

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);
}

/* Risolve */
void
SuperLUSolver::Solve(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	
	if (bHasBeenReset) {
      		((SuperLUSolver *)this)->Factor();
      		bHasBeenReset = false;
	}

#if 0
	y12solve(&iN, &((*Axp)[0]), &iCurSize, LinearSolver::pdRhs,
			    pdPIVOT, pic,
			    piHA, &iN,
			    iIFLAG, &iIFAIL);
#endif
}

/* Index Form */
void
SuperLUSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<int>& Ar, std::vector<int>& Ac,
		std::vector<int>& Ap) const
{
	if (!bHasBeenReset) {
		return;
	}
	
	iNonZeroes = mh.MakeCompressedColumnForm(Ax, Ar, Ap, 0);
	ASSERT(iNonZeroes > 0);
}

/* SuperLUSolver - end */


/* SuperLUSparseSolutionManager - begin: code */

/* Costruttore */
SuperLUSparseSolutionManager::SuperLUSparseSolutionManager(integer iSize, 
		integer iWorkSpaceSize,
		const doublereal& dPivotFactor,
		unsigned nt)
: iMatSize(iSize), 
iColStart(iSize + 1),
MH(iSize),
pVH(NULL)
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));


   	/* Valore di default */
   	if (iWorkSpaceSize == 0) {
		/*
		 * y12 requires at least 3*numzeros to store factors
		 * for multiple backsubs
		 */
      		iWorkSpaceSize = 3*iSize*iSize;
   	}
	
	integer iPivot;
	if (dPivotFactor == 0.) {
		iPivot = 0;
	} else {
		iPivot = 1;
	}

   	/* Alloca arrays */
	dVec.resize(iMatSize, 0.);
   
   	/* Alloca handlers ecc. */
   	SAFENEWWITHCONSTRUCTOR(pVH,
			       MyVectorHandler,
			       MyVectorHandler(iMatSize, &(dVec[0])));

	iRow.reserve(iWorkSpaceSize);
	dMat.reserve(iWorkSpaceSize);

   	SAFENEWWITHCONSTRUCTOR(SolutionManager::pLS, 
			       SuperLUSolver,
			       SuperLUSolver(nt, iMatSize,
			       		   &iRow, &iColStart,
					   &dMat, &(dVec[0])));
   
	pLS->SetSolutionManager(this);

#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
SuperLUSparseSolutionManager::~SuperLUSparseSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	if (pVH != NULL) {      
      		SAFEDELETE(pVH);
   	}
   
   	/* Dealloca roba, tra cui i thread */
}

#ifdef DEBUG
/* Test di validita' del manager */
void 
SuperLUSparseSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatSize > 0);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pVH));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLS));
#endif /* DEBUG_MEMMANAGER */
   
   	ASSERT((pVH->IsValid(), 1));
   	ASSERT((pLS->IsValid(), 1));
}
#endif /* DEBUG */

/* Inizializza il gestore delle matrici */
void
SuperLUSparseSolutionManager::MatrInit(const doublereal& dResetVal)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	MatrReset(dResetVal);
	pLS->Init();
}

void
SuperLUSparseSolutionManager::MatrReset(const doublereal& d)
{
	MH.Reset(d);
}

void
SuperLUSparseSolutionManager::MakeCompressedColumnForm(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
	pLS->MakeCompactForm(MH, dMat, iRow, iCol, iColStart);
}

/* Risolve il problema */
void
SuperLUSparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
   	MakeCompressedColumnForm();

#if 0
	std::cerr << "### after MakeIndexForm:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < iColStart.size(); i++) {
		std::cerr << i << " " << iColStart[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	unsigned c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << iRow[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == iColStart[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif

   	pLS->Solve();

#if 0
	std::cerr << "### after Solve:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < iColStart.size(); i++) {
		std::cerr << i << " " << iColStart[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << iRow[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == iColStart[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif
}

/* SuperLUSparseSolutionManager - end */

/* SuperLUSparseCCSolutionManager - begin */

template <class CC>
SuperLUSparseCCSolutionManager<CC>::SuperLUSparseCCSolutionManager(integer Dim,
		integer dummy, doublereal dPivot, unsigned nt)
: SuperLUSparseSolutionManager(Dim, dummy, dPivot, true),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
SuperLUSparseCCSolutionManager<CC>::~SuperLUSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MatrReset(const doublereal& d)
{
	if (!CCReady) {
		MH.Reset(d);
	} else {
		Ac->Reset(d);
	}
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(MH, dMat, iRow, iCol, iColStart);

		ASSERT(Ac == 0);

		SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(dMat, iRow, iColStart));

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MatrInitialize(const doublereal& d)
{
	SAFEDELETE(Ac);
	Ac = 0;

	CCReady = false;

	MatrInit();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
SuperLUSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &MH;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class SuperLUSparseCCSolutionManager<CColMatrixHandler<0> >;
template class SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* SuperLUSparseCCSolutionManager - end */

#endif /* USE_SUPERLU_MT */

