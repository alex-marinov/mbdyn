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

/* gestore dei dati */

#ifndef MTDATAMAN_H
#define MTDATAMAN_H

#ifdef USE_MULTITHREAD

#include <pthread.h>
#include <semaphore.h>

#include "dataman.h"

/* MultiThreadDataManager - begin */

class MultiThreadDataManager : public DataManager {
protected:
	/* from input file, or auto-detected */
	unsigned int nThreads;

	/* per-thread specific data */
	struct PerThreadData {
		MultiThreadDataManager *pDM;
		unsigned int threadNumber;
		pthread_t thread;
		sem_t sem;

		mutable MT_VecIter<Elem *> ElemIter;
	
		integer	*piWorkIndex;    /* array di lavoro */
		doublereal *pdWorkMat;   /* matrice di lavoro */

		VariableSubMatrixHandler *pWorkMatA;  /* SubMatrix di lavoro */
		VariableSubMatrixHandler *pWorkMatB;
		VariableSubMatrixHandler *pWorkMat;
		MySubVectorHandler *pWorkVec;

		MatrixHandler* pJacHdl;
		VectorHandler* pResHdl;
		MatrixHandler* pMatA;
		MatrixHandler* pMatB;
		doublereal dCoef;
	} *ptd;

	enum DataManagerOp {
		UNKNOWN_OP = -1,

		OP_ASSJAC,
		OP_ASSMATS,
		OP_ASSRES,
		OP_BEFOREPREDICT,
		OP_AFTERPREDICT,
		OP_AFTERCONVERGENCE,

		OP_EXIT,

		LAST_OP
	} op;

	/* will be replaced by barriers ... */
	unsigned dataman_thread_count;

	/* this can be replaced by a barrier ... */
	pthread_mutex_t	dataman_thread_mutex;
	pthread_cond_t	dataman_thread_cond;

	void EndOfOp(void);

	/* thread function */
	static void *dataman_thread(void *arg);
	static void dataman_thread_cleanup(PerThreadData *arg);

	/* starts the helper threads */
	void MultiThreadSpawn(void);

	/* reset InUse flag(s) before multithread execution */
	void ResetInUse(bool b = false);

public:
	/* costruttore - legge i dati e costruisce le relative strutture */
	MultiThreadDataManager(MBDynParser& HP,
			unsigned OF,
			doublereal dInitialTime,
			const char* sInputFileName,
			const char* sOutputFileName,
			bool bAbortAfterInput,
			unsigned nt);

	/* distruttore */
	virtual ~MultiThreadDataManager(void);

	/* Assembla lo jacobiano */
	virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef);

	/* Assembla le matrici per gli autovalori */
	virtual void AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl);

	/* Assembla il residuo */
	virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef);
};

/* MultiThreadDataManager - end */

#endif /* USE_MULTITHREAD */

#endif /* MTDATAMAN_H */

