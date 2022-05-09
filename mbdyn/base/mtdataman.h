/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "ac/pthread.h"		/* includes POSIX semaphores */

#include "dataman.h"
#include "spmh.h"

#ifdef USE_NAIVE_MULTITHREAD
#include "naivemh.h"
#endif

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient_spmh.h"
#endif

class Solver;

/* MultiThreadDataManager - begin */

class MultiThreadDataManager : public DataManager {
protected:
        // nThreads is now in DataManager
        enum {
                ASS_UNKNOWN = -1,

                ASS_CC,			/* use native column-compressed form */
#ifdef USE_NAIVE_MULTITHREAD
                ASS_NAIVE,		/* use native H-P sparse solver */
#endif
#ifdef USE_SPARSE_AUTODIFF
                ASS_GRAD,
                ASS_GRAD_PROD,
#endif
                ASS_DEFAULT,
                ASS_LAST
        } AssMode;

        /* steps of CC computation */
        enum {
                CC_NO,
                CC_FIRST,
                CC_YES
        } CCReady;

        /* per-thread specific data */
        struct ThreadData {
                MultiThreadDataManager *pDM;
                integer threadNumber;
                integer iCPUIndex;
                pthread_t thread;
                sem_t sem;
                std::exception_ptr except;
                mutable MT_VecIter<Elem *> ElemIter;

                VariableSubMatrixHandler *pWorkMatA;	/* Working SubMatrix */
                VariableSubMatrixHandler *pWorkMatB;
                VariableSubMatrixHandler *pWorkMat;	/* same as pWorkMatA */
                MySubVectorHandler *pWorkVec;

                /* for CC assembly */
                CompactSparseMatrixHandler* pJacHdl;
#ifdef USE_NAIVE_MULTITHREAD
                /* for Naive assembly */
                NaiveMatrixHandler** ppNaiveJacHdl;
#endif
#ifdef USE_SPARSE_AUTODIFF
                SpGradientSparseMatrixWrapper oGradJacHdl;
                const VectorHandler* pY;
                VectorHandler* pJacProd;
#endif
                AO_TS_t* lock;
#ifdef MBDYN_X_MT_ASSRES
                VectorHandler* pResHdl;
#endif
                VectorHandler* pAbsResHdl;
                MatrixHandler* pMatA;
                MatrixHandler* pMatB;
                doublereal dCoef;
        } *thread_data;

        enum DataManagerOp {
                OP_UNKNOWN = -1,

                OP_ASSJAC_CC,
#ifdef USE_NAIVE_MULTITHREAD
                OP_ASSJAC_NAIVE,
                OP_SUM_NAIVE,
#endif
#ifdef USE_SPARSE_AUTODIFF
                OP_ASSJAC_GRAD,
                OP_ASSJAC_PROD,
#endif
                /* used only #ifdef MBDYN_X_MT_ASSRES */
                OP_ASSRES,

                /* not used yet */
                OP_ASSMATS,
                OP_BEFOREPREDICT,
                OP_AFTERPREDICT,
                OP_AFTERCONVERGENCE,
                /* end of not used yet */

                OP_EXIT,

                LAST_OP
        } op;

        /* will be replaced by barriers ... */
        unsigned thread_count;

        /* this can be replaced by a barrier ... */
        pthread_mutex_t	thread_mutex;
        pthread_cond_t	thread_cond;

        /* this is used to propagate ErrMatrixRebuild ... */
        AO_TS_t	propagate_ErrMatrixRebuild;

        void EndOfOp(void);

        /* thread function */
        static void *thread(void *arg);
        static void thread_cleanup(ThreadData *arg);

        /* starts the helper threads */
        void ThreadSpawn(void);
        void ThreadDestroy(void);

        /* specialized assembly */
        virtual void CCAssJac(MatrixHandler& JacHdl, doublereal dCoef);
#ifdef USE_NAIVE_MULTITHREAD
        virtual void NaiveAssJac(NaiveMatrixHandler& JacHdl, doublereal dCoef);
        virtual void NaiveAssJacInit(NaiveMatrixHandler& JacHdl, doublereal dCoef);
#endif
#ifdef USE_SPARSE_AUTODIFF
        void GradAssJac(SpGradientSparseMatrixHandler& JacHdl, doublereal dCoef);
        void GradAssJacProd(VectorHandler& JacY, const VectorHandler& Y, doublereal dCoef);
        virtual void AssJac(VectorHandler& JacY, const VectorHandler& Y, doublereal dCoef) override;
#endif
        static void SetAffinity(const ThreadData& oThread);
public:
        /* costruttore - legge i dati e costruisce le relative strutture */
        MultiThreadDataManager(MBDynParser& HP,
                        unsigned OF,
                        Solver* pS,
                        doublereal dInitialTime,
                        const char* sOutputFileName,
                        const char* sInputFileName,
                        bool bAbortAfterInput,
                        unsigned nt);

        /* distruttore */
        virtual ~MultiThreadDataManager(void);

        /* Assembla lo jacobiano */
        virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef);

#ifdef MBDYN_X_MT_ASSRES
        /* Assembla il residuo */
        virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef, VectorHandler*const pAbsResHdl = 0)
                /*throw(ChangedEquationStructure)*/;
#endif /* MBDYN_X_MT_ASSRES */
};

/* MultiThreadDataManager - end */

#endif /* USE_MULTITHREAD */

#endif /* MTDATAMAN_H */
