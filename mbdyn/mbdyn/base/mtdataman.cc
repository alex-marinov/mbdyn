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

#include "mtdataman.h"

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
ptd(NULL),
op(MultiThreadDataManager::UNKNOWN_OP),
dataman_helper_count(0)
{
	MultiThreadSpawn();
} /* End of MultiThreadDataManager::MultiThreadDataManager() */

MultiThreadDataManager::~MultiThreadDataManager(void)
{
	NO_OP;
}

void *
MultiThreadDataManager::dataman_helper(void *p)
{
	MultiThreadDataManager::PerThreadData *arg
		= (MultiThreadDataManager::PerThreadData *)p;
	bool bKeepGoing = true;

	while (bKeepGoing) {
		/* stop here until told to start */
		/* NOTE: here
		 * - the requested operation must be set;
		 * - the appropriate operation args must be set
		 * - the dataman_helper_count must be set to nThreads
		 * - the caller must be ready to wait on dataman_helper_cond
		 */
		sem_wait(&arg->sem);

		/* select requested operation */
		switch (arg->pDM->op) {
		case MultiThreadDataManager::OP_ASSJAC:
			arg->pDM->DataManager::AssJac(*(arg->pJacHdl),
					arg->dCoef,
					(VecIter<Elem *> *)&arg->ElemIter,
					*(arg->pWorkMat));
			break;

		case MultiThreadDataManager::OP_ASSMATS:
			arg->pDM->DataManager::AssMats(*(arg->pMatA),
					*(arg->pMatB),
					(VecIter<Elem *> *)&arg->ElemIter,
					*(arg->pWorkMatA), *(arg->pWorkMatB));
			break;

		case MultiThreadDataManager::OP_ASSRES:
			arg->pDM->DataManager::AssRes(*(arg->pResHdl),
					arg->dCoef,
					(VecIter<Elem *> *)&arg->ElemIter,
					*(arg->pWorkVec));
			break;

		case MultiThreadDataManager::OP_EXIT:
			/* cleanup */
			SAFEDELETE(arg->pWorkMatA);
			SAFEDELETE(arg->pWorkMatB);
			SAFEDELETE(arg->pWorkVec);
			SAFEDELETEARR(arg->piWorkIndex);
			SAFEDELETEARR(arg->pdWorkMat);
			sem_destroy(&arg->sem);

			bKeepGoing = false;
			break;

		default:
			std::cerr << "unhandled op" << std::endl;
			THROW(ErrGeneric());
		}

		/* decrement the thread counter */
		pthread_mutex_lock(&arg->pDM->dataman_helper_mutex);
		arg->pDM->dataman_helper_count--;
		pthread_mutex_unlock(&arg->pDM->dataman_helper_mutex);

		/* if last thread, signal to restart */
		if (arg->pDM->dataman_helper_count == 0) {
			pthread_cond_signal(&arg->pDM->dataman_helper_cond);
		}
	}

	return NULL;
}

/* starts the helper threads */
void
MultiThreadDataManager::MultiThreadSpawn(void)
{
	ASSERT(nThreads > 1);

	SAFENEWARR(ptd, MultiThreadDataManager::PerThreadData, nThreads);
	
	for (unsigned i = 0; i < nThreads; i++) {
		/* callback data */
		ptd[i].pDM = this;
		ptd[i].threadNumber = i;
		sem_init(&ptd[i].sem, 0, 0);
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
					ptd[i].piWorkIndex+iWorkIntSize,
					ptd[i].pdWorkMat+iWorkDoubleSize));

		ptd[i].pWorkMat = ptd[i].pWorkMatA;

		SAFENEWWITHCONSTRUCTOR(ptd[i].pWorkVec,
				MySubVectorHandler,
				MySubVectorHandler(iWorkIntSize,
					ptd[i].piWorkIndex, ptd[i].pdWorkMat));

		/* create thread */
		if (pthread_create(&ptd[i].thread, NULL, dataman_helper,
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
	if (ptd == NULL) {
		DataManager::AssJac(JacHdl, dCoef);
		return;
	}

	ResetInUse(false);
	op = MultiThreadDataManager::OP_ASSJAC;
	dataman_helper_count = nThreads;

	pthread_mutex_lock(&dataman_helper_mutex);

	for (unsigned i = 0; i < nThreads; i++) {
		ptd[i].pJacHdl = &JacHdl;
		ptd[i].dCoef = dCoef;
	
		sem_post(&ptd[i].sem);
	}

	pthread_cond_wait(&dataman_helper_cond, &dataman_helper_mutex);
	pthread_mutex_unlock(&dataman_helper_mutex);
}

void
MultiThreadDataManager::AssMats(MatrixHandler& MatA, MatrixHandler& MatB)
{
	if (ptd == NULL) {
		DataManager::AssMats(MatA, MatB);
		return;
	}

	ResetInUse(false);
	op = MultiThreadDataManager::OP_ASSMATS;
	dataman_helper_count = nThreads;

	pthread_mutex_lock(&dataman_helper_mutex);

	for (unsigned i = 0; i < nThreads; i++) {
		ptd[i].pMatA = &MatA;
		ptd[i].pMatB = &MatB;
	
		sem_post(&ptd[i].sem);
	}

	pthread_cond_wait(&dataman_helper_cond, &dataman_helper_mutex);
	pthread_mutex_unlock(&dataman_helper_mutex);
}

void
MultiThreadDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
	if (ptd == NULL) {
		DataManager::AssRes(ResHdl, dCoef);
		return;
	}

	ResetInUse(false);
	op = MultiThreadDataManager::OP_ASSRES;
	dataman_helper_count = nThreads;

	pthread_mutex_lock(&dataman_helper_mutex);

	for (unsigned i = 0; i < nThreads; i++) {
		ptd[i].pResHdl = &ResHdl;
		ptd[i].dCoef = dCoef;
	
		sem_post(&ptd[i].sem);
	}

	pthread_cond_wait(&dataman_helper_cond, &dataman_helper_mutex);
	pthread_mutex_unlock(&dataman_helper_mutex);
}

#endif /* USE_MULTITHREAD */
