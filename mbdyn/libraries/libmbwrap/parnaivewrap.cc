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
 *                          ParNaive C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

// #include "spmh.h"
// #include "spmapmh.h"
// #include "dirccmh.h"
// #include "ccmh.h"

// extern "C" {
// #include <pdsp_defs.h>
// #include <util.h>
// 
// extern void *pdgstrf_thread(void *);
// extern void pdgstrf_finalize(pdgstrf_options_t *, SuperMatrix*);
// }

#include "parnaivewrap.h"


struct ParNaiveSolverData {
// 	SuperMatrix		A,
// 				AC,
// 				L,
// 				U,
// 				B;
// 	Gstat_t			Gstat;
	std::vector<int>	perm_c, /* column permutation vector */
				perm_r; /* row permutations from partial pivoting */
// 	pdgstrf_options_t	pdgstrf_options;
// 	pxgstrf_shared_t	pxgstrf_shared;
// 	pdgstrf_threadarg_t	*pdgstrf_threadarg;
};

/* ParNaiveSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
ParNaiveSolver::ParNaiveSolver(unsigned nt, const integer &size, 
		const doublereal& dMP,
		NaiveMatrixHandler *const a)
: LinearSolver(0),
iSize(size),
dMinPiv(dMP),
piv(size),
A(a),
//bFirstSol(true),
//bRegenerateMatrix(true),
sld(0),
nThreads(nt),
thread_data(0)
{
	ASSERT(iN > 0);

	SAFENEW(sld, ParNaiveSolverData);
	
	/*
	 * This means it's the first run
	 * FIXME: create a dependence on the library's internals
	 */
// 	sld->A.Store = NULL;

	pthread_mutex_init(&thread_mutex, NULL);
	pthread_cond_init(&thread_cond, NULL);

	SAFENEWARR(thread_data, thread_data_t, nThreads);

	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].pSLUS = this;
		thread_data[t].threadNumber = t;

		if (t == 0) {
			continue;
		}

		sem_init(&thread_data[t].sem, 0, 0);
		if (pthread_create(&thread_data[t].thread, NULL, thread_op, &thread_data[t]) != 0) {
			silent_cerr("ParNaiveSolver: pthread_create() failed "
					"for thread " << t
					<< " of " << nThreads << std::endl);
			THROW(ErrGeneric());
		}
	}
}

/* Distruttore */
ParNaiveSolver::~ParNaiveSolver(void)
{
#if 0
	/* ------------------------------------------------------------
	 * Clean up and free storage after multithreaded factorization.
	 * ------------------------------------------------------------*/
	pdgstrf_thread_finalize(sld->pdgstrf_threadarg, &sld->pxgstrf_shared, 
			&sld->A, &sld->perm_r[0], &sld->L, &sld->U);

	/* ------------------------------------------------------------
	 * Deallocate storage after factorization.
	 * ------------------------------------------------------------*/
	pdgstrf_finalize(&sld->pdgstrf_options, &sld->AC);
#endif
	
	thread_operation = ParNaiveSolver::EXIT;
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

	for (unsigned i = 1; i < nThreads; i++) {
		sem_destroy(&thread_data[i].sem);
	}

	pthread_mutex_destroy(&thread_mutex);
	pthread_cond_destroy(&thread_cond);

	/* other cleanup... */
}

void *
ParNaiveSolver::thread_op(void *arg)
{
	thread_data_t *td = (thread_data_t *)arg;

	silent_cout("ParNaiveSolver: thread " << td->threadNumber
			<< " [self=" << pthread_self()
			<< ",pid=" << getpid() << "]"
			<< " starting..." << std::endl);

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
		case ParNaiveSolver::FACTOR:
// 			(void)pdgstrf_thread(td->pdgstrf_threadarg);
			break;

		case ParNaiveSolver::EXIT:
			bKeepGoing = false;
			break;

		default:
			silent_cerr("ParNaiveSolver: unhandled op"
					<< std::endl);
			THROW(ErrGeneric());
		}

		td->pSLUS->EndOfOp();
	}

	pthread_exit(NULL);
}

void
ParNaiveSolver::EndOfOp(void)
{
	bool last;

	/* decrement the thread counter */
	pthread_mutex_lock(&thread_mutex);
	thread_count--;
	last = (thread_count == 0);

	/* if last thread, signal to restart */
	if (last) {
#if 0
		pthread_cond_broadcast(&thread_cond);
#else
		pthread_cond_signal(&thread_cond);
#endif
	}
	pthread_mutex_unlock(&thread_mutex);
}

#ifdef DEBUG
void 
ParNaiveSolver::IsValid(void) const
{
	ASSERT(Aip != NULL);
	ASSERT(App != NULL);
	ASSERT(Axp != NULL);
	ASSERT(iN > 0); 
}
#endif /* DEBUG */

/* Fattorizza la matrice */
void
ParNaiveSolver::Factor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	ASSERT(iNonZeroes > 0);

// 	yes_no_t	refact = YES,
// 			usepr = NO;
// 	doublereal	u = dMinPiv,
// 			drop_tol = 0.0;
// 	void		*work = NULL;
// 	int		info = 0, lwork = 0;
// 	int		panel_size = sp_ienv(1),
// 			relax = sp_ienv(2);

// 	if (bRegenerateMatrix) {
// 		refact = NO;
// 		
// 		/* NOTE: we could use sld->A.Store == NULL */
// 		if (bFirstSol) {
// 			ASSERT(Astore == NULL);
// 
// 			/* ---------------------------------------------------
// 			 * Allocate storage and initialize statistics variables. 
// 			 * ---------------------------------------------------*/
// 			/* Set up the dense matrix data structure for B. */
// 			dCreate_Dense_Matrix(&sld->B, iN, 1,
// 					LinearSolver::pdRhs,
// 					iN, DN, _D, GE);
// 
// 			/* Set up the sparse matrix data structure for A. */
// 			dCreate_CompCol_Matrix(&sld->A, iN, iN, iNonZeroes,
// 					Axp, Aip, App, NC, _D, GE);
// 
// 			StatAlloc(iN, nThreads, panel_size, relax, &sld->Gstat);
// 			StatInit(iN, nThreads, &sld->Gstat);
// 			
// 			sld->perm_c.resize(iN);
// 			sld->perm_r.resize(iN);
// 
// 			bFirstSol = false;	/* never change this again */
// 
// 		} else {
// 			NCformat *Astore = (NCformat *) sld->A.Store;
// 
// 			ASSERT(Astore);
// 
// 			Astore->nnz = iNonZeroes;
// 			Astore->nzval = Axp;
// 			Astore->rowind = Aip;
// 			Astore->colptr = App;
// 		}
// 
// 		/* --------------------------------------------------
// 		 * Get column permutation vector perm_c[], according
// 		 * 	to permc_spec:
// 		 * permc_spec = 0: use the natural ordering 
// 		 * permc_spec = 1: use minimum degree ordering
// 		 * 	on structure of A'*A
// 		 * permc_spec = 2: use minimum degree ordering
// 		 * 	on structure of A'+A
// 		 * !!! permc_spec = 3: use approximate minimum
// 		 * 	degree column order !!!
// 		 * --------------------------------------------------*/
// 
// 		/*
// 		 * According to Umfpack's use of AMD:
// 		 *
// 		 * symmetric matrix:
// 		 *	AMD: A^T + A permutation
// 		 *
// 		 * Non symmetric matrix:
// 		 *	COLAMD: A^T * A permutation
// 		 *
// 		 * so we use permc_spec = 1
// 		 */
// 	
// 		int	permc_spec = 1;
// 		int	*pc = &(sld->perm_c[0]);
// 		get_perm_c(permc_spec, &sld->A, pc);
// 
// 		bRegenerateMatrix = false;
// 	}

	int		*pr = &(sld->perm_r[0]),
			*pc = &(sld->perm_c[0]);

// 	/* ------------------------------------------------------------
// 	 * Initialize the option structure pdgstrf_options using the
// 	 * user-input parameters;
// 	 * Apply perm_c to the columns of original A to form AC.
// 	 * ------------------------------------------------------------*/
// 	pdgstrf_init(nThreads, refact, panel_size, relax,
// 			u, usepr, drop_tol, pc, pr,
// 			work, lwork, &sld->A, &sld->AC,
// 			&sld->pdgstrf_options, &sld->Gstat);
// 
// 
// 	/* --------------------------------------------------------------
// 	 * Initializes the parallel data structures for pdgstrf_thread().
// 	 * --------------------------------------------------------------*/
// 	sld->pdgstrf_threadarg = pdgstrf_thread_init(&sld->AC,
// 			&sld->L, &sld->U, &sld->pdgstrf_options,
// 			&sld->pxgstrf_shared, &sld->Gstat, &info);

// 	if (info != 0) {
// 		silent_cerr("ParNaiveSolver::Factor: pdgstrf_thread_init failed"
// 				<< std::endl);
// 	}
		
	for (unsigned t = 0; t < nThreads; t++) {
// 		thread_data[t].pdgstrf_threadarg = &sld->pdgstrf_threadarg[t];
	}

	thread_operation = ParNaiveSolver::FACTOR;
	thread_count = nThreads - 1;

	for (unsigned t = 1; t < nThreads; t++) {
		sem_post(&thread_data[t].sem);
	}

// 	(void)pdgstrf_thread(thread_data[0].pdgstrf_threadarg);

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);

	/* ------------------------------------------------------------
	 * Clean up and free storage after multithreaded factorization.
	 * ------------------------------------------------------------*/
	 
	 //!!!!!!
	 //TODO:
	 //!!!!!!
	 
// 	pdgstrf_thread_finalize(sld->pdgstrf_threadarg, &sld->pxgstrf_shared, 
// 			&sld->A, pr, &sld->L, &sld->U);
}

/* Risolve */
void
ParNaiveSolver::Solve(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	
	if (bHasBeenReset) {
      		((ParNaiveSolver *)this)->Factor();
      		bHasBeenReset = false;
	}

	/* ------------------------------------------------------------
	 * Solve the system A*X=B, overwriting B with X.
	 * ------------------------------------------------------------*/
// 	trans_t		trans = NOTRANS;
// 	int		info = 0;

	int		*pr = &(sld->perm_r[0]),
			*pc = &(sld->perm_c[0]);

// 	dgstrs(trans, &sld->L, &sld->U, pr, pc,
// 			&sld->B, &sld->Gstat, &info);
}

// /* Index Form */
// void
// ParNaiveSolver::MakeCompactForm(SparseMatrixHandler& mh,
// 		std::vector<doublereal>& Ax,
// 		std::vector<integer>& Ar, std::vector<integer>& Ac,
// 		std::vector<integer>& Ap) const
// {
// 	/* no need to rebuild matrix */
// 	if (!bHasBeenReset) {
// 		return;
// 	}
// 	
// 	iNonZeroes = mh.MakeCompressedColumnForm(Ax, Ar, Ap, 0);
// 	ASSERT(iNonZeroes > 0);
// 
// 	Axp = &Ax[0];
// 	Aip = &Ar[0];
// 	App = &Ap[0];
// 
// 	/* rebuild matrix ... (CC is broken) */
// 	bRegenerateMatrix = true;
// 
// #if 0
// 	Destroy_CompCol_Matrix(&sld->A);
// #endif
// }

/* ParNaiveSolver - end */


/* ParNaiveSparseSolutionManager - begin: code */

/* Costruttore */
ParNaiveSparseSolutionManager::ParNaiveSparseSolutionManager(unsigned nt,
		const integer Dim, const doublereal dMP)
: A(Dim),
VH(Dim),
XH(Dim)
{
   	ASSERT(Dim > 0);
   	ASSERT((dMP >= 0.0) && (dMP <= 1.0));

   	SAFENEWWITHCONSTRUCTOR(SolutionManager::pLS, 
			       ParNaiveSolver,
			       ParNaiveSolver(nt, Dim, dMP, &A));
   
	pLS->ChangeResPoint(VH.pdGetVec());
	pLS->ChangeSolPoint(XH.pdGetVec());
	pLS->SetSolutionManager(this);

#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
ParNaiveSparseSolutionManager::~ParNaiveSparseSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	/* Dealloca roba, tra cui i thread */
}

#ifdef DEBUG
/* Test di validita' del manager */
void 
ParNaiveSparseSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatSize > 0);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(VH.pdGetVec()));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(XH.pdGetVec()));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLS));
#endif /* DEBUG_MEMMANAGER */
   
   	ASSERT((VH.IsValid(), 1));
   	ASSERT((XH.IsValid(), 1));
   	ASSERT((pLS->IsValid(), 1));
}
#endif /* DEBUG */

/* Inizializza il gestore delle matrici */
void
ParNaiveSparseSolutionManager::MatrReset(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	//MH.Reset();
	pLS->Reset();
}

// void
// ParNaiveSparseSolutionManager::MakeCompressedColumnForm(void)
// {
// #ifdef DEBUG
//    	IsValid();
// #endif /* DEBUG */
// 
// 	/* FIXME: move this to the matrix handler! */
// 	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
// }
 
/* Risolve il problema */
void
ParNaiveSparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
//    	MakeCompressedColumnForm();

#if 0
	std::cerr << "### after MakeIndexForm:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < Ap.size(); i++) {
		std::cerr << i << " " << Ap[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	unsigned c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << Ai[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == Ap[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif

   	pLS->Solve();

#if 0
	std::cerr << "### after Solve:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < Ap.size(); i++) {
		std::cerr << i << " " << Ap[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << Ai[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == Ap[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif
}

/* ParNaiveSparseSolutionManager - end */



