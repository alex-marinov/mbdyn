/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
ptd(NULL),
op(MultiThreadDataManager::UNKNOWN_OP),
dataman_thread_count(0)
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
	ThreadSpawn();
}

MultiThreadDataManager::~MultiThreadDataManager(void)
{
	NO_OP;
}

clock_t
MultiThreadDataManager::ThreadDestroy(void)
{
	if (ptd == NULL) {
		return 0;
	}

	clock_t cputime = 0;

	op = MultiThreadDataManager::OP_EXIT;
	dataman_thread_count = nThreads - 1;

	for (unsigned i = 1; i < nThreads; i++) {
		void *retval = NULL;

		sem_post(&ptd[i].sem);
		if (pthread_join(ptd[i].thread, &retval)) {
			silent_cerr("pthread_join() failed on thread " << i
					<< std::endl);
			/* already shutting down ... */
		}

		cputime += ptd[i].cputime;
	}

	dataman_thread_cleanup(&ptd[0]);

	SAFEDELETEARR(ptd);

	return cputime;
}

void *
MultiThreadDataManager::dataman_thread(void *p)
{
	MultiThreadDataManager::PerThreadData *arg
		= (MultiThreadDataManager::PerThreadData *)p;
	bool bKeepGoing = true;

	while (bKeepGoing) {
		/* stop here until told to start */
		/* NOTE: here
		 * - the requested operation must be set;
		 * - the appropriate operation args must be set
		 * - the dataman_thread_count must be set to nThreads
		 * - the caller must be ready to wait on dataman_thread_cond
		 */
		sem_wait(&arg->sem);

		// std::cerr << "thread " << arg->threadNumber << ": "
		// 	"op " << arg->pDM->op << std::endl;

		/* select requested operation */
		switch (arg->pDM->op) {
		case MultiThreadDataManager::OP_ASSJAC:
			arg->pJacHdl->Reset(0.);
			arg->pDM->DataManager::AssJac(*(arg->pJacHdl),
					arg->dCoef,
					(VecIter<Elem *> *)&arg->ElemIter,
					*arg->pWorkMat);
			break;

		case MultiThreadDataManager::OP_ASSRES:
			arg->pResHdl->Reset(0.);
			arg->pDM->DataManager::AssRes(*(arg->pResHdl),
					arg->dCoef,
					(VecIter<Elem *> *)&arg->ElemIter,
					*arg->pWorkVec);
			break;

		case MultiThreadDataManager::OP_EXIT:
			/* cleanup */
			dataman_thread_cleanup(arg);
			bKeepGoing = false;
			break;

		default:
			std::cerr << "unhandled op" << std::endl;
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
MultiThreadDataManager::dataman_thread_cleanup(PerThreadData *arg)
{
	/* cleanup */
	SAFEDELETE(arg->pWorkMatA);
	SAFEDELETE(arg->pWorkMatB);
	SAFEDELETE(arg->pWorkVec);
	SAFEDELETEARR(arg->piWorkIndex);
	SAFEDELETEARR(arg->pdWorkMat);
	SAFEDELETE(arg->pJacHdl);
	SAFEDELETE(arg->pResHdl);
	sem_destroy(&arg->sem);

#ifdef HAVE_SYS_TIMES_H	
	/* Tempo di CPU impiegato */
	struct tms tmsbuf;
	times(&tmsbuf);

#if 0
	std::cerr
		<< "utime:  " << tmsbuf.tms_utime << std::endl
		<< "stime:  " << tmsbuf.tms_stime << std::endl
		<< "cutime: " << tmsbuf.tms_cutime << std::endl
		<< "cstime: " << tmsbuf.tms_cstime << std::endl;
#endif
			
	arg->cputime = tmsbuf.tms_utime + tmsbuf.tms_cutime
		+ tmsbuf.tms_stime + tmsbuf.tms_cstime;
#endif /* HAVE_SYS_TIMES_H */
}

void
MultiThreadDataManager::EndOfOp(void)
{
	bool last;
	
	/* decrement the thread counter */
	pthread_mutex_lock(&dataman_thread_mutex);
	dataman_thread_count--;
	last = (dataman_thread_count == 0);
	pthread_mutex_unlock(&dataman_thread_mutex);

	/* if last thread, signal to restart */
	if (last) {
		pthread_cond_signal(&dataman_thread_cond);
	}
}

/* starts the helper threads */
void
MultiThreadDataManager::ThreadSpawn(void)
{
	ASSERT(nThreads > 1);

	SAFENEWARR(ptd, MultiThreadDataManager::PerThreadData, nThreads);
	
	for (unsigned i = 0; i < nThreads; i++) {
		/* callback data */
		ptd[i].pDM = this;
		sem_init(&ptd[i].sem, 0, 0);
		ptd[i].threadNumber = i;
		ptd[i].ElemIter.Init(ppElems, iTotElem);

		/* thread workspace */
		SAFENEWARR(ptd[i].piWorkIndex, integer, 2*iWorkIntSize);
		SAFENEWARR(ptd[i].pdWorkMat, doublereal, 2*iWorkDoubleSize);

		/* SubMatrixHandlers */
		SAFENEWWITHCONSTRUCTOR(ptd[i].pWorkMatA,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iWorkIntSize,
					iWorkDoubleSize,
					ptd[i].piWorkIndex,
					ptd[i].pdWorkMat));

		SAFENEWWITHCONSTRUCTOR(ptd[i].pWorkMatB,
				VariableSubMatrixHandler,
				VariableSubMatrixHandler(iWorkIntSize,
					iWorkDoubleSize,
					ptd[i].piWorkIndex + iWorkIntSize,
					ptd[i].pdWorkMat + iWorkDoubleSize));

		ptd[i].pWorkMat = ptd[i].pWorkMatA;

		SAFENEWWITHCONSTRUCTOR(ptd[i].pWorkVec,
				MySubVectorHandler,
				MySubVectorHandler(iWorkIntSize,
					ptd[i].piWorkIndex, ptd[i].pdWorkMat));

		if (i == 0) continue;

		/* set by AssJac when in CC form */
		ptd[i].pJacHdl = NULL;

		SAFENEWWITHCONSTRUCTOR(ptd[i].pResHdl,
				MyVectorHandler, MyVectorHandler(iTotDofs));

		/* create thread */
		if (pthread_create(&ptd[i].thread, NULL, dataman_thread,
					&ptd[i]) != 0) {
			std::cerr << "pthread_create() failed "
				"for thread " << i << " of " << nThreads 
				<< std::endl;
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
	ASSERT(ptd != NULL);

	switch (CCReady) {
	case CC_NO:
		if (dynamic_cast<SpMapMatrixHandler *>(&JacHdl) == NULL) {
			THROW(ErrGeneric());
		}

		DataManager::AssJac(JacHdl, dCoef, &ElemIter, *pWorkMat);
		CCReady = CC_FIRST;
		// std::cerr << "CC_NO => CC_FIRST" << std::endl;
		return;

	case CC_FIRST: {
		CColMatrixHandler *pMH = dynamic_cast<CColMatrixHandler *>(&JacHdl);
		if (pMH == NULL) {
			THROW(ErrGeneric());
		}

		for (unsigned i = 1; i < nThreads; i++) {
			ptd[i].pJacHdl = pMH->Copy();
		}
		CCReady = CC_YES;
		// std::cerr << "CC_FIRST => CC_YES" << std::endl;
		break;
	}

	case CC_YES:
		// std::cerr << "CC_YES" << std::endl;
		break;

	default:
		THROW(ErrGeneric());

	}

	ResetInUse(false);
	op = MultiThreadDataManager::OP_ASSJAC;
	dataman_thread_count = nThreads - 1;

	pthread_mutex_lock(&dataman_thread_mutex);

	for (unsigned i = 1; i < nThreads; i++) {
		ptd[i].dCoef = dCoef;
	
		sem_post(&ptd[i].sem);
	}

	DataManager::AssJac(JacHdl, dCoef,
			(VecIter<Elem *> *)&ptd[0].ElemIter,
			*ptd[0].pWorkMat);

	pthread_cond_wait(&dataman_thread_cond, &dataman_thread_mutex);
	pthread_mutex_unlock(&dataman_thread_mutex);

	CColMatrixHandler *pMH = dynamic_cast<CColMatrixHandler *>(&JacHdl);
	for (unsigned t = 1; t < nThreads; t++) {
		pMH->AddUnchecked(*ptd[t].pJacHdl);
	}
}

void
MultiThreadDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
	ASSERT(ptd != NULL);

	ResetInUse(false);
	op = MultiThreadDataManager::OP_ASSRES;
	dataman_thread_count = nThreads - 1;

	pthread_mutex_lock(&dataman_thread_mutex);

	for (unsigned i = 1; i < nThreads; i++) {
		ptd[i].dCoef = dCoef;
	
		sem_post(&ptd[i].sem);
	}

	DataManager::AssRes(ResHdl, dCoef,
			(VecIter<Elem *> *)&ptd[0].ElemIter,
			*ptd[0].pWorkVec);

	pthread_cond_wait(&dataman_thread_cond, &dataman_thread_mutex);
	pthread_mutex_unlock(&dataman_thread_mutex);

	for (unsigned t = 1; t < nThreads; t++) {
		ResHdl += *ptd[t].pResHdl;
	}
}

clock_t
MultiThreadDataManager::GetCPUTime(void) const
{
	return ((MultiThreadDataManager *)this)->ThreadDestroy();
}

#endif /* USE_MULTITHREAD */

