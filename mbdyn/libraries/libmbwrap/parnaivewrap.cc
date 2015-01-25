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
/*
 * Copyright (C) 2009
 *
 * Marco Morandini
 *
 */

/*****************************************************************************
 *                                                                           *
 *                          ParNaive C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_NAIVE_MULTITHREAD

/* FIXME: incompatible with RTAI at present */
#ifndef USE_RTAI

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <cstdio>

#include <unistd.h>
#include <signal.h>

#include "parnaivewrap.h"
#include "mthrdslv.h"
#include "task2cpu.h"

#include "pmthrdslv.h"

/* ParNaiveSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
ParNaiveSolver::ParNaiveSolver(unsigned nt, const integer &size, 
		const doublereal& dMP,
		NaiveMatrixHandler *const a)
: LinearSolver(0),
iSize(size),
dMinPiv(dMP),
ppril(0),
pnril(0),
A(a),
nThreads(nt),
thread_data(0)
{
	ASSERT(iN > 0);

	piv.resize(iSize);
	fwd.resize(iSize);
	todo.resize(iSize);
	row_locks.resize(iSize + 2);
	col_locks.resize(iSize, AO_TS_INITIALIZER);

	pthread_mutex_init(&thread_mutex, NULL);
	pthread_cond_init(&thread_cond, NULL);

	SAFENEWARR(ppril, integer *, iSize);
	ppril[0] = 0;
	SAFENEWARR(ppril[0], integer, iSize*iSize);
	for (integer i = 1; i < iSize; i++) {
		ppril[i] = ppril[i - 1] + iSize;
	}
	SAFENEWARR(pnril, integer, iSize);
#ifdef HAVE_MEMSET_H
	memset(pnril, 0, sizeof(integer)*iSize);
#else /* ! HAVE_MEMSET_H */
	for (integer row = 0; row < iSize; row++) {
		pnril[row] = 0;
	}
#endif /* ! HAVE_MEMSET_H */		


	SAFENEWARRNOFILL(thread_data, thread_data_t, nThreads);
	
	(void)mbdyn_task2cpu(nThreads - 1);

	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].pSLUS = this;
		thread_data[t].threadNumber = t;
		thread_data[t].retval = 0;

		sem_init(&thread_data[t].sem, 0, 0);
		if (pthread_create(&thread_data[t].thread, NULL, thread_op, &thread_data[t]) != 0) {
			silent_cerr("ParNaiveSolver: pthread_create() failed "
					"for thread " << t
					<< " of " << nThreads << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

/* Distruttore */
ParNaiveSolver::~ParNaiveSolver(void)
{
	
	thread_operation = ParNaiveSolver::EXIT;
	thread_count = nThreads;

	for (unsigned i = 0; i < nThreads; i++) {
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

	if (ppril) {
		if (ppril[0]) {
			SAFEDELETEARR(ppril[0]);
		}
		SAFEDELETEARR(ppril);
	}
	if (pnril) {
		SAFEDELETEARR(pnril);
	}

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

	(void)mbdyn_task2cpu(td->threadNumber);

	bool bKeepGoing(true);

	while (bKeepGoing) {
		sem_wait(&td->sem);

		switch (td->pSLUS->thread_operation) {
		case ParNaiveSolver::FACTOR:
			td->retval = pnaivfct(td->pSLUS->A->ppdRows,
				td->pSLUS->A->iSize,
				td->pSLUS->A->piNzr,
				td->pSLUS->A->ppiRows,
				td->pSLUS->A->piNzc,
				td->pSLUS->A->ppiCols,
				td->pSLUS->pnril,
				td->pSLUS->ppril, 
				td->pSLUS->A->ppnonzero,
				&td->pSLUS->piv[0],
				&td->pSLUS->todo[0],
				td->pSLUS->dMinPiv,
				&td->pSLUS->row_locks[0],
				&td->pSLUS->col_locks[0],
				td->threadNumber,
				td->pSLUS->nThreads
			);

			if (td->retval) {
				if (td->retval & NAIVE_ENULCOL) {
					silent_cerr("NaiveSolver: NAIVE_ENULCOL("
							<< (td->retval & ~NAIVE_ENULCOL) << ")" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (td->retval & NAIVE_ENOPIV) {
					silent_cerr("NaiveSolver: NAIVE_ENOPIV("
							<< (td->retval & ~NAIVE_ENOPIV) << ")" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* default */
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case ParNaiveSolver::SOLVE:
			pnaivslv(td->pSLUS->A->ppdRows,
				td->pSLUS->A->iSize,
				td->pSLUS->A->piNzc,
				td->pSLUS->A->ppiCols,
				td->pSLUS->LinearSolver::pdRhs,
				&td->pSLUS->piv[0],
				&td->pSLUS->fwd[0],
				td->pSLUS->LinearSolver::pdSol,
				&td->pSLUS->row_locks[0],
				td->threadNumber,
				td->pSLUS->nThreads
			);
			break;

		case ParNaiveSolver::EXIT:
			bKeepGoing = false;
			break;

		default:
			silent_cerr("ParNaiveSolver: unhandled op"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	thread_operation = ParNaiveSolver::FACTOR;
	thread_count = nThreads;

	for (int i = 0; i < iSize; i++) {
			piv[i] = -1;
			todo[i] = -1;
			row_locks[i] = 0;
			pnril[0] = 0;
	}

	/* NOTE: these are for pivot_lock and sync_lock */
	/* FIXME: bottleneck? */
	row_locks[iSize] = 0;
	row_locks[iSize + 1] = 0;

	/* NOTE: no need to reset col_locks because they're
	 * always left equal to zero after use */

	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].retval = 0;
		sem_post(&thread_data[t].sem);
	}

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);

}

/* Risolve */
void
ParNaiveSolver::Solve(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	
	if (bHasBeenReset) {
      		const_cast<ParNaiveSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}

	for (int i = 0; i < iSize + 2; i++) {
		row_locks[i] = 0;
	}

	thread_operation = ParNaiveSolver::SOLVE;
	thread_count = nThreads;
	
	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].retval = 0;
		sem_post(&thread_data[t].sem);
	}


	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);
}

void
ParNaiveSolver::SetMat(NaiveMatrixHandler *const a)
{
	A = a;
}

/* ParNaiveSolver - end */


/* ParNaiveSparseSolutionManager - begin: code */

/* Costruttore */
ParNaiveSparseSolutionManager::ParNaiveSparseSolutionManager(unsigned nt,
		const integer Dim, const doublereal dMP)
: A(0),
VH(Dim),
XH(Dim)
{
   	ASSERT(Dim > 0);
   	ASSERT((dMP >= 0.0) && (dMP <= 1.0));

	SAFENEWWITHCONSTRUCTOR(A, NaiveMatrixHandler, NaiveMatrixHandler(Dim));
   	SAFENEWWITHCONSTRUCTOR(SolutionManager::pLS, 
			       ParNaiveSolver,
			       ParNaiveSolver(nt, Dim, dMP, A));
   
	pLS->pdSetResVec(VH.pdGetVec());
	pLS->pdSetSolVec(XH.pdGetVec());
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
   
   	/* Dealloca roba, ... */
	if (A != 0) {
		SAFEDELETE(A);
		A = 0;
	}

   	/* ... tra cui i thread */
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

	pLS->Reset();
}

 
/* Risolve il problema */
void
ParNaiveSparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */


   	pLS->Solve();

}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
ParNaiveSparseSolutionManager::pMatHdl(void) const
{
	return A;
}

VectorHandler*
ParNaiveSparseSolutionManager::pResHdl(void) const
{
#ifdef DEBUG
	VH.IsValid();
#endif /* DEBUG */
	return &VH;
}

/* Rende disponibile l'handler per la soluzione */
VectorHandler*
ParNaiveSparseSolutionManager::pSolHdl(void) const
{
#ifdef DEBUG
	XH.IsValid();
#endif /* DEBUG */
	return &XH;
}

/* ParNaiveSparseSolutionManager - end */


/* ParNaiveSparsePermSolutionManager - begin */

extern "C" {
#include "colamd.h"
}

ParNaiveSparsePermSolutionManager::ParNaiveSparsePermSolutionManager(
	unsigned nt,
	const integer Dim,
	const doublereal dMP)
: ParNaiveSparseSolutionManager(nt, Dim, dMP),
dMinPiv(dMP),
ePermState(PERM_NO)
{
	perm.resize(Dim, 0);
	invperm.resize(Dim, 0);

	SAFEDELETE(A);
	A = 0;
	SAFENEWWITHCONSTRUCTOR(A, NaivePermMatrixHandler, NaivePermMatrixHandler(Dim, perm, invperm));

	dynamic_cast<ParNaiveSolver *>(pLS)->SetMat(A);

	MatrInitialize();
}

ParNaiveSparsePermSolutionManager::~ParNaiveSparsePermSolutionManager(void) 
{
	NO_OP;
}

void
ParNaiveSparsePermSolutionManager::MatrReset(void)
{
	if (ePermState == PERM_INTERMEDIATE) {
		ePermState = PERM_READY;

		pLS->pdSetResVec(VH.pdGetVec());
		pLS->pdSetSolVec(XH.pdGetVec());

		pLS->SetSolutionManager(this);
	}

	pLS->Reset();
}

void
ParNaiveSparsePermSolutionManager::ComputePermutation(void)
{
	std::vector<integer> Ai;
	A->MakeCCStructure(Ai, invperm);
	doublereal knobs [COLAMD_KNOBS];
	integer stats [COLAMD_STATS];
	integer Alen = mbdyn_colamd_recommended (Ai.size(), A->iGetNumRows(), A->iGetNumCols());
	Ai.resize(Alen);
	mbdyn_colamd_set_defaults(knobs);
	if (!mbdyn_colamd(A->iGetNumRows(), A->iGetNumCols(), Alen,
		&(Ai[0]), &(invperm[0]), knobs, stats)) {
		silent_cerr("colamd permutation failed" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	for (integer i = 0; i < A->iGetNumRows(); i++) {
		perm[invperm[i]] = i;
	}
	ePermState = PERM_INTERMEDIATE;
}

void
ParNaiveSparsePermSolutionManager::BackPerm(void)
{
	for (integer i = 0; i < A->iGetNumCols(); i++) {
		XH(invperm[i] + 1) = VH(i + 1);
	}
}


/* Risolve il sistema  Fattorizzazione + Backward Substitution */
void
ParNaiveSparsePermSolutionManager::Solve(void)
{
	if (ePermState == PERM_NO) {
		ComputePermutation();

	} else {
		pLS->pdSetSolVec(VH.pdGetVec());
	}

	pLS->Solve();

	if (ePermState == PERM_READY) {
		BackPerm();
		pLS->pdSetSolVec(XH.pdGetVec());
	}
}

/* Inizializzatore "speciale" */
void
ParNaiveSparsePermSolutionManager::MatrInitialize()
{
	ePermState = PERM_NO;

	for (integer i = 0; i < A->iGetNumRows(); i++) {
		perm[i] = i;
		invperm[i] = i;
	}
 
	MatrReset();
}
	
/* ParParNaiveSparsePermSolutionManager - end */


#endif /* ! USE_RTAI */

#endif /* USE_NAIVE_MULTITHREAD */
