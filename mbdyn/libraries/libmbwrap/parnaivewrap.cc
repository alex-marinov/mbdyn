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

#ifdef USE_NAIVE_MULTITHREAD

/* FIXME: incompatible with RTAI at present */
#ifndef USE_RTAI

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <stdio.h>

#include <unistd.h>
#include <signal.h>

#include "parnaivewrap.h"

extern "C" {
int pnaivfct(doublereal** a,
	integer neq,
	integer *nzr, integer** ri,
	integer *nzc, integer** ci,
	char ** nz,
	integer *piv,
	integer *todo,
	doublereal minpiv,
	unsigned long *locks,
	int task,
	int NCPU);

void pnaivslv(doublereal** a,
	integer neq,
	integer *nzc,
	integer** ci,
	doublereal *rhs,
	integer *piv,
	doublereal *fwd,
	doublereal *sol,
	unsigned long *locks,
	int task,
	int NCPU);
}

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
nThreads(nt),
thread_data(0)
{
	ASSERT(iN > 0);

	
	piv.resize(iSize);
	fwd.resize(iSize);
	todo.resize(iSize);
	locks.resize(iSize+2);

	pthread_mutex_init(&thread_mutex, NULL);
	pthread_cond_init(&thread_cond, NULL);

	SAFENEWARR(thread_data, thread_data_t, nThreads);
	
	do {
		int fd;
		fd = open("/dev/TASK2CPU",O_RDWR);
		if (fd <= 0) {
			silent_cerr("Error opening /dev/TASK2CPU" << std::endl);
		}
		ioctl(fd, 0, nThreads);
		close(fd);
	} while(false);

	for (unsigned t = 0; t < nThreads; t++) {
		thread_data[t].pSLUS = this;
		thread_data[t].threadNumber = t;

		sem_init(&thread_data[t].sem, 0, 0);
		if (pthread_create(&thread_data[t].thread, NULL, thread_op, &thread_data[t]) != 0) {
			silent_cerr("ParNaiveSolver: pthread_create() failed "
					"for thread " << t
					<< " of " << nThreads << std::endl);
			throw ErrGeneric();
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

	do {
		int fd;
		fd = open("/dev/TASK2CPU",O_RDWR);
		if (fd <= 0) {
			silent_cerr("Error opening /dev/TASK2CPU" << std::endl);
		}
		ioctl(fd, 0, td->threadNumber);
		close(fd);
	} while(false);


	bool bKeepGoing(true);

	while (bKeepGoing) {
		sem_wait(&td->sem);

		switch (td->pSLUS->thread_operation) {
		case ParNaiveSolver::FACTOR:
			pnaivfct(td->pSLUS->A->ppdRows,
				td->pSLUS->A->iSize,
				td->pSLUS->A->piNzr,
				td->pSLUS->A->ppiRows,
				td->pSLUS->A->piNzc,
				td->pSLUS->A->ppiCols,
				td->pSLUS->A->ppnonzero,
				&(td->pSLUS->piv[0]),
				&(td->pSLUS->todo[0]),
				td->pSLUS->dMinPiv,
				&(td->pSLUS->locks[0]),
				td->threadNumber,
				//std::min(td->threadNumber,1U),
				td->pSLUS->nThreads
			);
			break;

		case ParNaiveSolver::SOLVE:
			pnaivslv(td->pSLUS->A->ppdRows,
				td->pSLUS->A->iSize,
				td->pSLUS->A->piNzc,
				td->pSLUS->A->ppiCols,
				td->pSLUS->LinearSolver::pdRhs,
				&(td->pSLUS->piv[0]),
				&(td->pSLUS->fwd[0]),
				td->pSLUS->LinearSolver::pdSol,
				&(td->pSLUS->locks[0]),
				td->threadNumber,
				//std::min(td->threadNumber,1U),
				td->pSLUS->nThreads
			);
			break;

		case ParNaiveSolver::EXIT:
			bKeepGoing = false;
			break;

		default:
			silent_cerr("ParNaiveSolver: unhandled op"
					<< std::endl);
			throw ErrGeneric();
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

	//prepara i dati
	for (int i = 0; i < iSize; i++) {
			piv[i] = -1;
			todo[i] = -1;
	}
	for (int i = 0; i < iSize+2; i++) {
			locks[i] = 0;
	}

	for (unsigned t = 0; t < nThreads; t++) {
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
      		((ParNaiveSolver *)this)->Factor();
      		bHasBeenReset = false;
	} else {
		for (int i = 0; i < iSize+2; i++) {
			locks[i] = 0;
		}
	}

	thread_operation = ParNaiveSolver::SOLVE;
	thread_count = nThreads;

	
	for (unsigned t = 0; t < nThreads; t++) {
		sem_post(&thread_data[t].sem);
	}


	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);
}


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

VectorHandler* ParNaiveSparseSolutionManager::pResHdl(void) const {
#ifdef DEBUG
	VH.IsValid();
#endif /* DEBUG */
	return &VH;
};

	/* Rende disponibile l'handler per la soluzione */
VectorHandler* ParNaiveSparseSolutionManager::pSolHdl(void) const {
#ifdef DEBUG
	XH.IsValid();
#endif /* DEBUG */
	return &XH;
};

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
PermReady(false),
Ap(0)
{
	perm.resize(Dim,0);
	invperm.resize(Dim,0);
}

ParNaiveSparsePermSolutionManager::~ParNaiveSparsePermSolutionManager(void) 
{
	if (Ap) {
		SAFEDELETE(Ap);
	}
}

void
ParNaiveSparsePermSolutionManager::MatrReset(void)
{
	if (Ap) {
		PermReady = true;

		pLS->ChangeResPoint(VH.pdGetVec());
		pLS->ChangeSolPoint(XH.pdGetVec());

		pLS->SetSolutionManager(this);
	}
	pLS->Reset();
}

void
ParNaiveSparsePermSolutionManager::ComputePermutation(void) {
	std::vector<integer> Ai;
	A.MakeCCStructure(Ai, invperm);
	doublereal knobs [COLAMD_KNOBS];
	integer stats [COLAMD_STATS];
	integer Alen = colamd_recommended (Ai.size(), A.iGetNumRows(), A.iGetNumCols());
	Ai.resize(Alen);
	colamd_set_defaults(knobs);
	if (!colamd(A.iGetNumRows(), A.iGetNumCols(), Alen,
		&(Ai[0]), &(invperm[0]), knobs, stats)) {
		silent_cerr("colamd permutation failed" << std::endl);
		throw ErrGeneric();
	}
	for (integer i = 0; i < A.iGetNumRows(); i++) {
		perm[invperm[i]] = i;
	}
}

void
ParNaiveSparsePermSolutionManager::BackPerm(void) {
	for (integer i = 0; i < A.iGetNumCols(); i++) {
		XH.PutCoef(invperm[i] + 1, VH.dGetCoef(i + 1));
	}
}


/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
ParNaiveSparsePermSolutionManager::Solve(void)
{
	if ((!PermReady) && (!Ap)) {
		ComputePermutation();
		if (Ap) {
			SAFEDELETE(Ap);
		}
		SAFENEWWITHCONSTRUCTOR(Ap, NaivePermMatrixHandler, 
			NaivePermMatrixHandler(&A, &(perm[0])));
	} else {
		pLS->ChangeSolPoint(VH.pdGetVec());
	}
	pLS->Solve();
	if (PermReady) {
		BackPerm();
		pLS->ChangeSolPoint(XH.pdGetVec());
	}
}

/* Inizializzatore "speciale" */
void
ParNaiveSparsePermSolutionManager::MatrInitialize()
{
	PermReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
MatrixHandler*
ParNaiveSparsePermSolutionManager::pMatHdl(void) const
{
	if (!PermReady) {
		return &A;
	}

	ASSERT(Ap != 0);
	return Ap;
}



/* ParParNaiveSparsePermSolutionManager - end */


#endif /* ! USE_RTAI */

#endif /* USE_NAIVE_MULTITHREAD */
