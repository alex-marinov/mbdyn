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

/* datamanager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_MULTITHREAD

extern "C" {
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif /* HAVE_SYS_TIMES_H */
#ifdef HAVE_SCHED_H
#include <sched.h>
#endif /* HAVE_SCHED_H */
}

#include "mtdataman.h"
#include "umfpackwrap.h"

/* MultiThreadDataManager - begin */


/*
 * costruttore: inizializza l'oggetto, legge i dati e crea le strutture di
 * gestione di Dof, nodi, elementi e drivers.
 */

MultiThreadDataManager::MultiThreadDataManager(MBDynParser& HP, 
		unsigned OF,
		doublereal dInitialTime,
		const char* sInputFileName, 
		const char* sOutputFileName,
		bool bAbortAfterInput,
		unsigned nt)
:
DataManager(HP, OF, dInitialTime, sInputFileName, sOutputFileName,
		bAbortAfterInput),
nThreads(nt),
CCReady(CC_NO),
thread_data(0),
op(MultiThreadDataManager::UNKNOWN_OP),
thread_count(0),
propagate_ErrMatrixRebuild(sig_atomic_t(false))
{
#if 0	/* no effects ... */
	struct sched_param	sp;
	int			policy = SCHED_FIFO;
	int			rc;

	rc = sched_getparam(0, &sp);
	if (rc != 0) {
		silent_cerr("sched_getparam() failed: " << errno << std::endl);
		THROW(ErrGeneric());
	}

	int pmin = sched_get_priority_min(policy);
	int pmax = sched_get_priority_max(policy);

	silent_cout("current priority is " << sp.sched_priority
			<< " {" << pmin << "," << pmax << "}" << std::endl);

	if (sp.sched_priority > pmax || sp.sched_priority < pmin) {
		sp.sched_priority = pmax;
	}

	rc = sched_setscheduler(0, policy, &sp);
	if (rc != 0) {
		silent_cerr("sched_setscheduler() unable "
				"to set SCHED_FIFO scheduling policy: "
				<< errno
				<< std::endl);
		THROW(ErrGeneric());
	}
#endif
	if (pthread_mutex_init(&thread_mutex, NULL)) {
		silent_cerr("MultiThreadDataManager::MultiThreadDataManager(): "
				"mutex init failed" << std::endl);
		THROW(ErrGeneric());
	}
		
	if (pthread_cond_init(&thread_cond, NULL)) {
		silent_cerr("MultiThreadDataManager::MultiThreadDataManager(): "
				"cond init failed" << std::endl);
		THROW(ErrGeneric());
	}

	ThreadSpawn();
}

MultiThreadDataManager::~MultiThreadDataManager(void)
{
	pthread_mutex_destroy(&thread_mutex);
	pthread_cond_destroy(&thread_cond);

}

clock_t
MultiThreadDataManager::ThreadDestroy(void)
{
	if (thread_data == 0) {
		return 0;
	}

	clock_t cputime = 0;

	op = MultiThreadDataManager::OP_EXIT;
	thread_count = nThreads - 1;

	for (unsigned i = 1; i < nThreads; i++) {
		void *retval = NULL;

		sem_post(&thread_data[i].sem);
		if (pthread_join(thread_data[i].thread, &retval)) {
			silent_cerr("pthread_join() failed on thread " << i
					<< std::endl);
			/* already shutting down ... */
		}

		cputime += thread_data[i].cputime;
	}

	thread_cleanup(&thread_data[0]);

	SAFEDELETEARR(thread_data);

	return cputime;
}

void *
MultiThreadDataManager::thread(void *p)
{
	MultiThreadDataManager::ThreadData *arg
		= (MultiThreadDataManager::ThreadData *)p;

	silent_cout("thread " << arg->threadNumber 
			<< " starting..." << std::endl);
	
	bool bKeepGoing = true;

	/* deal with signals ... */
	sigset_t newset /* , oldset */ ;
	sigemptyset(&newset);
	sigaddset(&newset, SIGTERM);
	sigaddset(&newset, SIGINT);
	sigaddset(&newset, SIGHUP);
	pthread_sigmask(SIG_BLOCK, &newset, /* &oldset */ NULL);

	while (bKeepGoing) {
		/* stop here until told to start */
		/* NOTE: here
		 * - the requested operation must be set;
		 * - the appropriate operation args must be set
		 * - the thread_count must be set to nThreads - 1
		 */
		sem_wait(&arg->sem);

		DEBUGCERR("thread " << arg->threadNumber << ": "
				"op " << arg->pDM->op << std::endl);

		/* select requested operation */
		switch (arg->pDM->op) {
		case MultiThreadDataManager::OP_ASSJAC:
			arg->pJacHdl->Reset(0.);
			try {
				arg->pDM->DataManager::AssJac(*(arg->pJacHdl),
						arg->dCoef,
						arg->ElemIter,
						*arg->pWorkMat);

			} catch (MatrixHandler::ErrRebuildMatrix) {
				silent_cerr("thread " << arg->threadNumber
						<< " caught ErrRebuildMatrix"
						<< std::endl);

				mbdyn_compare_and_swap(arg->pDM->propagate_ErrMatrixRebuild,
						sig_atomic_t(true), sig_atomic_t(false));

			} catch (...) {
				throw;
			}
			break;

		case MultiThreadDataManager::OP_ASSRES:
			arg->pResHdl->Reset(0.);
			arg->pDM->DataManager::AssRes(*(arg->pResHdl),
					arg->dCoef,
					arg->ElemIter,
					*arg->pWorkVec);
			break;

		case MultiThreadDataManager::OP_EXIT:
			/* cleanup */
			thread_cleanup(arg);
			bKeepGoing = false;
			break;

		default:
			silent_cerr("unhandled op" << std::endl);
			THROW(ErrGeneric());
		}

		/* decrease counter and signal if last
		 * (mutex + cond) */
		arg->pDM->EndOfOp();
	}

	/* all threads are joined */
	pthread_exit(NULL);
}

void
MultiThreadDataManager::thread_cleanup(ThreadData *arg)
{
	/* cleanup */
	SAFEDELETE(arg->pWorkMatA);
	SAFEDELETE(arg->pWorkMatB);
	SAFEDELETE(arg->pWorkVec);
	SAFEDELETEARR(arg->piWorkIndex);
	SAFEDELETEARR(arg->pdWorkMat);
	if (arg->threadNumber > 0) {
		if (arg->pJacHdl) {
			SAFEDELETE(arg->pJacHdl);
		}
		SAFEDELETE(arg->pResHdl);
	}
	sem_destroy(&arg->sem);

#ifdef HAVE_SYS_TIMES_H	
	/* Tempo di CPU impiegato */
	struct tms tmsbuf;
	times(&tmsbuf);

	pedantic_cout("Thread " << arg->threadNumber << ":" << std::endl
		<< "\tutime:  " << tmsbuf.tms_utime << std::endl
		<< "\tstime:  " << tmsbuf.tms_stime << std::endl
		<< "\tcutime: " << tmsbuf.tms_cutime << std::endl
		<< "\tcstime: " << tmsbuf.tms_cstime << std::endl);
			
	arg->cputime = tmsbuf.tms_utime + tmsbuf.tms_cutime
		+ tmsbuf.tms_stime + tmsbuf.tms_cstime;
#endif /* HAVE_SYS_TIMES_H */
}

void
MultiThreadDataManager::EndOfOp(void)
{
	bool last;
	
	/* decrement the thread counter */
	pthread_mutex_lock(&thread_mutex);
	thread_count--;
	last = (thread_count == 0);

	/* if last thread, signal to restart */
	if (last) {
		pthread_cond_signal(&thread_cond);
		// pthread_cond_broadcast(&thread_cond);
	}

	pthread_mutex_unlock(&thread_mutex);
}

/* starts the helper threads */
void
MultiThreadDataManager::ThreadSpawn(void)
{
	ASSERT(nThreads > 1);

	SAFENEWARR(thread_data, MultiThreadDataManager::ThreadData, nThreads);
	
	for (unsigned i = 0; i < nThreads; i++) {
		/* callback data */
		thread_data[i].pDM = this;
		sem_init(&thread_data[i].sem, 0, 0);
		thread_data[i].threadNumber = i;
		thread_data[i].ElemIter.Init(ppElems, iTotElem);

		/* thread workspace */
		SAFENEWARR(thread_data[i].piWorkIndex, integer, 2*iWorkIntSize);
		SAFENEWARR(thread_data[i].pdWorkMat, doublereal, 2*iWorkDoubleSize);

		/* SubMatrixHandlers */
		SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkMatA,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iWorkIntSize,
					iWorkDoubleSize,
					thread_data[i].piWorkIndex,
					thread_data[i].pdWorkMat));

		SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkMatB,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iWorkIntSize,
					iWorkDoubleSize,
					thread_data[i].piWorkIndex + iWorkIntSize,
					thread_data[i].pdWorkMat + iWorkDoubleSize));

		thread_data[i].pWorkMat = thread_data[i].pWorkMatA;

		SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkVec,
				MySubVectorHandler,
				MySubVectorHandler(iWorkIntSize,
					thread_data[i].piWorkIndex,
					thread_data[i].pdWorkMat));

		if (i == 0) continue;

		/* set by AssJac when in CC form */
		thread_data[i].pJacHdl = 0;

		SAFENEWWITHCONSTRUCTOR(thread_data[i].pResHdl,
				MyVectorHandler, MyVectorHandler(iTotDofs));

		/* create thread */
		if (pthread_create(&thread_data[i].thread, NULL, thread,
					&thread_data[i]) != 0) {
			silent_cerr("pthread_create() failed "
					"for thread " << i
					<< " of " << nThreads << std::endl);
			THROW(ErrGeneric());
		}
	}
}

void
MultiThreadDataManager::ResetInUse(bool b)
{
	DEBUGCOUT("Entering MultiThreadDataManager::ResetInUse()" << std::endl);

	Elem* pTmpEl = NULL;
	if (ElemIter.bGetFirst(pTmpEl)) {
		do {
			pTmpEl->SetInUse(b);
		} while (ElemIter.bGetNext(pTmpEl));
	}
}

void
MultiThreadDataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
	ASSERT(thread_data != NULL);

	propagate_ErrMatrixRebuild = sig_atomic_t(false);

	CompactSparseMatrixHandler *pMH
		= dynamic_cast<CompactSparseMatrixHandler *>(&JacHdl);

	while (false) {
retry:;
		CCReady = CC_NO;
		for (unsigned i = 1; i < nThreads; i++) {
			SAFEDELETE(thread_data[i].pJacHdl);
			thread_data[i].pJacHdl = 0;
		}
	}

	switch (CCReady) {
	case CC_NO:
		DEBUGCERR("CC_NO => CC_FIRST" << std::endl);

		ASSERT(dynamic_cast<SpMapMatrixHandler *>(&JacHdl));

		DataManager::AssJac(JacHdl, dCoef, ElemIter, *pWorkMat);
		CCReady = CC_FIRST;

		return;

	case CC_FIRST: {

#if 0
		ASSERT(pMH);
#endif
		if (pMH == 0) {
			goto retry;
		}

		DEBUGCERR("CC_FIRST => CC_YES" << std::endl);

		for (unsigned i = 1; i < nThreads; i++) {
			thread_data[i].pJacHdl = pMH->Copy();
		}

		CCReady = CC_YES;

		break;
	}

	case CC_YES:

#if 0
		ASSERT(pMH);
#endif
		if (pMH == 0) {
			goto retry;
		}

		DEBUGCERR("CC_YES" << std::endl);

		break;

	default:
		THROW(ErrGeneric());

	}

	thread_data[0].ElemIter.ResetAccessData();
	op = MultiThreadDataManager::OP_ASSJAC;
	thread_count = nThreads - 1;

	for (unsigned i = 1; i < nThreads; i++) {
		thread_data[i].dCoef = dCoef;
	
		sem_post(&thread_data[i].sem);
	}

	try {
		DataManager::AssJac(JacHdl, dCoef, thread_data[0].ElemIter,
				*thread_data[0].pWorkMat);

	} catch (MatrixHandler::ErrRebuildMatrix) {
		silent_cerr("Thread " << thread_data[0].threadNumber
				<< " caught ErrRebuildMatrix"
				<< std::endl);

		mbdyn_compare_and_swap(propagate_ErrMatrixRebuild,
				sig_atomic_t(true), sig_atomic_t(false));

	} catch (...) {
		throw;
	}

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);

	if (propagate_ErrMatrixRebuild) {
		for (unsigned i = 1; i < nThreads; i++) {
			SAFEDELETE(thread_data[i].pJacHdl);
			thread_data[i].pJacHdl = 0;
		}
		CCReady = CC_NO;

		throw MatrixHandler::ErrRebuildMatrix();
	}

	for (unsigned i = 1; i < nThreads; i++) {
		pMH->AddUnchecked(*thread_data[i].pJacHdl);
	}
}

void
MultiThreadDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
	ASSERT(thread_data != NULL);

	thread_data[0].ElemIter.ResetAccessData();
	op = MultiThreadDataManager::OP_ASSRES;
	thread_count = nThreads - 1;

	for (unsigned i = 1; i < nThreads; i++) {
		thread_data[i].dCoef = dCoef;
	
		sem_post(&thread_data[i].sem);
	}

	DataManager::AssRes(ResHdl, dCoef, thread_data[0].ElemIter,
			*thread_data[0].pWorkVec);

	pthread_mutex_lock(&thread_mutex);
	if (thread_count > 0) {
		pthread_cond_wait(&thread_cond, &thread_mutex);
	}
	pthread_mutex_unlock(&thread_mutex);

	for (unsigned i = 1; i < nThreads; i++) {
		ResHdl += *thread_data[i].pResHdl;
	}
}

clock_t
MultiThreadDataManager::GetCPUTime(void) const
{
	return ((MultiThreadDataManager *)this)->ThreadDestroy();
}

#endif /* USE_MULTITHREAD */

