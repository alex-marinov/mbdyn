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

/* datamanager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_MULTITHREAD

extern "C" {
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif /* HAVE_SYS_TIMES_H */
#ifdef HAVE_SCHED_H
#include <sched.h>
#endif /* HAVE_SCHED_H */


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif
#include <stdio.h>
}

#include <cerrno>

#include "mtdataman.h"
#include "spmapmh.h"
#include "task2cpu.h"

static inline void
do_lock(volatile AO_TS_t *p)
{
        while (mbdyn_test_and_set(p) == AO_TS_SET);
}

static inline void
do_unlock(AO_TS_t *p)
{
        AO_CLEAR(p);
}

#ifdef USE_NAIVE_MULTITHREAD
static void
naivepsad(doublereal **ga, integer **gri,
                integer *gnzr, integer **gci, integer *gnzc, char **gnz,
                doublereal **a, integer **ci, integer *nzc,
                integer from, integer to, AO_TS_t *lock)
{

        for (integer r = from; r < to; r++) {
                integer nc = nzc[r];

                if (nc) {
                        doublereal *pgar = ga[r];
                        doublereal *par  = a[r];

                        for (integer i = 0; i < nc; i++) {
                                integer c = ci[r][i];

                                if (gnz[r][c]) {
                                        pgar[c] += par[c];

                                } else {
                                        pgar[c] = par[c];
                                        gci[r][gnzc[r]] = c;
                                        /* This can only be set to 1 from 0,
                                         * so concurrency is harmless
                                         */
                                        gnz[r][c] = 1;

                                        do_lock(&lock[c]);

                                        gri[c][gnzr[c]] = r;
                                        gnzr[c]++;

                                        do_unlock(&lock[c]);
                                        gnzc[r]++;
                                }
                        }
                }
        }
}
#endif

/* MultiThreadDataManager - begin */


/*
 * costruttore: inizializza l'oggetto, legge i dati e crea le strutture di
 * gestione di Dof, nodi, elementi e drivers.
 */

MultiThreadDataManager::MultiThreadDataManager(MBDynParser& HP,
                unsigned OF,
                Solver* pS,
                doublereal dInitialTime,
                const char* sOutputFileName,
                const char* sInputFileName,
                bool bAbortAfterInput,
                unsigned nThreads)
:
DataManager(HP, OF, pS, dInitialTime, sOutputFileName, sInputFileName, bAbortAfterInput),
AssMode(ASS_UNKNOWN),
CCReady(CC_NO),
thread_data(0),
op(MultiThreadDataManager::OP_UNKNOWN),
thread_count(0),
propagate_ErrMatrixRebuild(AO_TS_INITIALIZER)
{
        DataManager::nThreads = nThreads;

#if 0	/* no effects ... */
        struct sched_param	sp;
        int			policy = SCHED_FIFO;
        int			rc;

        rc = sched_getparam(0, &sp);
        if (rc != 0) {
                silent_cerr("sched_getparam() failed: " << errno << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
#endif

        if (pthread_mutex_init(&thread_mutex, NULL)) {
                silent_cerr("MultiThreadDataManager::MultiThreadDataManager(): "
                                "mutex init failed" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (pthread_cond_init(&thread_cond, NULL)) {
                silent_cerr("MultiThreadDataManager::MultiThreadDataManager(): "
                                "cond init failed" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        ThreadSpawn();
}

MultiThreadDataManager::~MultiThreadDataManager(void)
{
        ThreadDestroy();
        pthread_mutex_destroy(&thread_mutex);
        pthread_cond_destroy(&thread_cond);
}

void MultiThreadDataManager::ThreadDestroy(void)
{
        if (thread_data == 0) {
                return;
        }

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
        }

        if (thread_data[0].lock) {
                SAFEDELETEARR(thread_data[0].lock);
        }
        thread_cleanup(&thread_data[0]);

        SAFEDELETEARR(thread_data);
        thread_data = 0;
}


void *
MultiThreadDataManager::thread(void *p)
{
        MultiThreadDataManager::ThreadData *arg
                = (MultiThreadDataManager::ThreadData *)p;

        silent_cout("MultiThreadDataManager: thread " << arg->threadNumber
                        << " [self=" << pthread_self()
                        << ",pid=" << getpid() << "]"
                        << " starting..." << std::endl);

        bool bKeepGoing = true;

#ifdef HAVE_PTHREAD_SIGMASK
        /* deal with signals ... */
        sigset_t newset /* , oldset */ ;
        sigemptyset(&newset);
        sigaddset(&newset, SIGTERM);
        sigaddset(&newset, SIGINT);
        sigaddset(&newset, SIGHUP);
        pthread_sigmask(SIG_BLOCK, &newset, /* &oldset */ NULL);
#endif

        SetAffinity(*arg);

        while (bKeepGoing) {
             /* stop here until told to start */
             /*
              * NOTE: here
              * - the requested operation must be set;
              * - the appropriate operation args must be set
              * - the thread_count must be set to nThreads - 1
              */
             sem_wait(&arg->sem);

             try {
                  DEBUGCOUT("thread " << arg->threadNumber << ": "
                            "op " << arg->pDM->op << std::endl);

                  /* select requested operation */
                  switch (arg->pDM->op) {
                  case MultiThreadDataManager::OP_ASSJAC_CC:
                       //arg->pJacHdl->Reset();
                       try {
                            arg->pDM->DataManager::AssJac(*arg->pJacHdl,
                                                          arg->dCoef,
                                                          arg->ElemIter,
                                                          *arg->pWorkMat);

                       } catch (MatrixHandler::ErrRebuildMatrix& e) {
                            silent_cerr("thread " << arg->threadNumber
                                        << " caught ErrRebuildMatrix"
                                        << std::endl);

                            mbdyn_test_and_set(&arg->pDM->propagate_ErrMatrixRebuild);

                       } catch (...) {
                            throw;
                       }
                       break;
#ifdef USE_NAIVE_MULTITHREAD
                  case MultiThreadDataManager::OP_ASSJAC_NAIVE:
#if 0
                       arg->ppNaiveJacHdl[arg->threadNumber]->Reset();
#endif
                       /* NOTE: Naive should never throw
                        * ErrRebuildMatrix ... */
                       arg->pDM->DataManager::AssJac(*arg->ppNaiveJacHdl[arg->threadNumber],
                                                     arg->dCoef,
                                                     arg->ElemIter,
                                                     *arg->pWorkMat);
                       break;

                  case MultiThreadDataManager::OP_SUM_NAIVE:
                  {
                       /* FIXME: if the naive matrix is permuted (colamd),
                        * this should not impact the parallel assembly,
                        * because all the matrices refer to the same
                        * permutation vector */
                       NaiveMatrixHandler* to = arg->ppNaiveJacHdl[0];
                       integer nn = to->iGetNumRows();
                       integer iFrom = (nn*(arg->threadNumber))/arg->pDM->nThreads;
                       integer iTo = (nn*(arg->threadNumber + 1))/arg->pDM->nThreads;
                       for (unsigned int matrix = 1; matrix < arg->pDM->nThreads; matrix++) {
                            NaiveMatrixHandler* from = arg->ppNaiveJacHdl[matrix];
                            naivepsad(to->ppdRows,
                                      to->ppiRows, to->piNzr,
                                      to->ppiCols, to->piNzc, to->ppnonzero,
                                      from->ppdRows, from->ppiCols, from->piNzc,
                                      iFrom, iTo, arg->lock);
                       }
                       break;
                  }
#endif
#ifdef USE_SPARSE_AUTODIFF
                  case MultiThreadDataManager::OP_ASSJAC_GRAD:
                  {
                       arg->pDM->DataManager::AssJac(arg->oGradJacHdl,
                                                     arg->dCoef,
                                                     arg->ElemIter,
                                                     *arg->pWorkMat);
                       break;
                  }
#endif
#ifdef MBDYN_X_MT_ASSRES
                  case MultiThreadDataManager::OP_ASSRES:
                       arg->pResHdl->Reset();
                       if (arg->pAbsResHdl) arg->pAbsResHdl->Reset();
                       arg->pDM->DataManager::AssRes(*arg->pResHdl,
                                                     arg->dCoef,
                                                     arg->ElemIter,
                                                     *arg->pWorkVec,
                                                     arg->pAbsResHdl);
                       break;
#endif /* MBDYN_X_MT_ASSRES */

                  case MultiThreadDataManager::OP_EXIT:
                       /* cleanup */
                       thread_cleanup(arg);
                       bKeepGoing = false;
                       break;

                  default:
                       silent_cerr("MultiThreadDataManager: unhandled op"
                                   << std::endl);
                       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                  }

             } catch (...) {
                  arg->except = std::current_exception();
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

        if (arg->threadNumber > 0) {
                if (arg->pJacHdl) {
                        SAFEDELETE(arg->pJacHdl);
                }
#ifdef MBDYN_X_MT_ASSRES
                SAFEDELETE(arg->pResHdl);
#endif
#ifdef USE_NAIVE_MULTITHREAD
                if (arg->ppNaiveJacHdl && arg->ppNaiveJacHdl[arg->threadNumber]) {
                        SAFEDELETE(arg->ppNaiveJacHdl[arg->threadNumber]);
                        arg->ppNaiveJacHdl[arg->threadNumber] = nullptr;
                }
#endif
                if (arg->pAbsResHdl) {
                     SAFEDELETE(arg->pAbsResHdl);
                }
        } else {
#ifdef USE_NAIVE_MULTITHREAD
                if (arg->ppNaiveJacHdl) {
                        // can be nonzero only when in Naive form
                        SAFEDELETEARR(arg->ppNaiveJacHdl);
                        arg->ppNaiveJacHdl = nullptr;
                }
#endif
        }

#ifdef USE_SPARSE_AUTODIFF
        arg->oGradJacHdl.~SpGradientSparseMatrixWrapper();
#endif

        sem_destroy(&arg->sem);

#ifdef HAVE_SYS_TIMES_H
        /* Tempo di CPU impiegato */
        struct tms tmsbuf;
        times(&tmsbuf);

        pedantic_cout("thread " << arg->threadNumber << ":" << std::endl
                << "\tutime:  " << tmsbuf.tms_utime << std::endl
                << "\tstime:  " << tmsbuf.tms_stime << std::endl
                << "\tcutime: " << tmsbuf.tms_cutime << std::endl
                << "\tcstime: " << tmsbuf.tms_cstime << std::endl);
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

        SAFENEWARRNOFILL(thread_data, MultiThreadDataManager::ThreadData, nThreads);

        const Task2CPU& oCPUSet = Task2CPU::GetGlobalState();
        int iCPUIndex = oCPUSet.iGetFirstCPU();
        const unsigned uNumCPUs = oCPUSet.iGetCount();

        for (unsigned i = 0; i < nThreads; i++) {
                /* callback data */
                thread_data[i].pDM = this;
                sem_init(&thread_data[i].sem, 0, 0);
                thread_data[i].threadNumber = i;

                if (uNumCPUs >= nThreads) {
                     thread_data[i].iCPUIndex = iCPUIndex;
                     iCPUIndex = oCPUSet.iGetNextCPU(iCPUIndex);
                } else {
                     thread_data[i].iCPUIndex = -1;
                }

                thread_data[i].ElemIter.Init(&Elems[0], Elems.size());
                thread_data[i].lock = 0;

                /* SubMatrixHandlers */
                thread_data[i].pWorkMatA = 0;
                SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkMatA,
                                VariableSubMatrixHandler,
                                VariableSubMatrixHandler(iMaxWorkNumRowsJac,
                                                         iMaxWorkNumColsJac,
                                                         iMaxWorkNumItemsJac));

                thread_data[i].pWorkMatB = 0;
                SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkMatB,
                                VariableSubMatrixHandler,
                                VariableSubMatrixHandler(iMaxWorkNumRowsJac,
                                                         iMaxWorkNumColsJac,
                                                         iMaxWorkNumItemsJac));

                thread_data[i].pWorkMat = thread_data[i].pWorkMatA;

                thread_data[i].pWorkVec = 0;
                SAFENEWWITHCONSTRUCTOR(thread_data[i].pWorkVec,
                                MySubVectorHandler,
                                MySubVectorHandler(iMaxWorkNumRowsRes));

                /* set by AssJac when in CC form */
                thread_data[i].pJacHdl = 0;
#ifdef USE_NAIVE_MULTITHREAD
                /* set by AssJac when in Naive form */
                thread_data[i].ppNaiveJacHdl = 0;
#endif
#ifdef USE_SPARSE_AUTODIFF
                new(&thread_data[i].oGradJacHdl) SpGradientSparseMatrixWrapper(nullptr);
#endif
#ifdef MBDYN_X_MT_ASSRES
                /* set below */
                thread_data[i].pResHdl = 0;
#endif
                thread_data[i].pAbsResHdl = 0;

                /* to be sure... */
                thread_data[i].pMatA = 0;
                thread_data[i].pMatB = 0;

                if (i == 0) {
                        continue;
                }
#ifdef MBDYN_X_MT_ASSRES
                SAFENEWWITHCONSTRUCTOR(thread_data[i].pResHdl,
                                MyVectorHandler, MyVectorHandler(iTotDofs));

                SAFENEWWITHCONSTRUCTOR(thread_data[i].pAbsResHdl,
                                MyVectorHandler, MyVectorHandler(iTotDofs));
#endif
                /* create thread */
                if (pthread_create(&thread_data[i].thread, NULL, thread,
                                        &thread_data[i]) != 0) {
                        silent_cerr("pthread_create() failed "
                                        "for thread " << i
                                        << " of " << nThreads << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }


        SetAffinity(thread_data[0]);
}

void
MultiThreadDataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
#ifdef USE_SPARSE_AUTODIFF
        SpGradientSparseMatrixHandler* pGradJacHdl = nullptr;
#endif
#ifdef USE_NAIVE_MULTITHREAD
        NaiveMatrixHandler* pNaiveJacHdl = nullptr;
#endif
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
        NodesUpdateJac(dCoef, NodeIter);
#endif
        if (dynamic_cast<CompactSparseMatrixHandler*>(&JacHdl) || dynamic_cast<SpMapMatrixHandler*>(&JacHdl)) {
                AssMode = ASS_CC;
                CCAssJac(JacHdl, dCoef);
#ifdef USE_SPARSE_AUTODIFF
        } else if ((pGradJacHdl = dynamic_cast<SpGradientSparseMatrixHandler*>(&JacHdl))) {
                AssMode = ASS_GRAD;
                GradAssJac(*pGradJacHdl, dCoef);
#endif
#ifdef USE_NAIVE_MULTITHREAD
        } else if ((pNaiveJacHdl = dynamic_cast<NaiveMatrixHandler*>(&JacHdl))) {
                AssMode = ASS_NAIVE;
                NaiveAssJacInit(*pNaiveJacHdl, dCoef);
#endif
        } else {
                AssMode = ASS_DEFAULT;
                // Single threaded assembly needed for eigenanalysis using lapack
                DataManager::AssJac(JacHdl, dCoef);
        }
}

void
MultiThreadDataManager::CCAssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
        ASSERT(thread_data != NULL);

        AO_CLEAR(&propagate_ErrMatrixRebuild);

        auto *pMH = dynamic_cast<CompactSparseMatrixHandler*>(&JacHdl);

        while (false) {
retry:;
                CCReady = CC_NO;
                for (unsigned i = 1; i < nThreads; i++) {
                        if (thread_data[i].pJacHdl) {
                                SAFEDELETE(thread_data[i].pJacHdl);
                                thread_data[i].pJacHdl = 0;
                        }
                }
        }

        switch (CCReady) {
        case CC_NO:
                DEBUGCERR("CC_NO => CC_FIRST" << std::endl);

                ASSERT(dynamic_cast<SpMapMatrixHandler *>(&JacHdl) != 0);

                DataManager::AssJac(JacHdl, dCoef, ElemIter, *pWorkMat);
                CCReady = CC_FIRST;

                return;

        case CC_FIRST:
                if (pMH == 0) {
                        goto retry;
                }

                DEBUGCERR("CC_FIRST => CC_YES" << std::endl);

                for (unsigned i = 1; i < nThreads; i++) {
                        thread_data[i].pJacHdl = pMH->Copy();
                }

                CCReady = CC_YES;

                break;

        case CC_YES:
                if (pMH == 0) {
                        goto retry;
                }

                DEBUGCERR("CC_YES" << std::endl);

                break;

        default:
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);

        }

        thread_data[0].ElemIter.ResetAccessData();
        op = MultiThreadDataManager::OP_ASSJAC_CC;
        thread_count = nThreads - 1;

        for (unsigned i = 0; i < nThreads; ++i) {
             thread_data[i].except = std::exception_ptr{};
        }
        
        for (unsigned i = 1; i < nThreads; i++) {
                thread_data[i].dCoef = dCoef;

                sem_post(&thread_data[i].sem);
        }

        try {
                DataManager::AssJac(JacHdl, dCoef, thread_data[0].ElemIter,
                                    *thread_data[0].pWorkMat);

        } catch (MatrixHandler::ErrRebuildMatrix& e) {
                silent_cerr("thread " << thread_data[0].threadNumber
                                << " caught ErrRebuildMatrix"
                                << std::endl);

                mbdyn_test_and_set(&propagate_ErrMatrixRebuild);
        } catch (...) {
             thread_data[0].except = std::current_exception();
        }

        pthread_mutex_lock(&thread_mutex);

        if (thread_count > 0) {
                pthread_cond_wait(&thread_cond, &thread_mutex);
        }

        pthread_mutex_unlock(&thread_mutex);

        if (propagate_ErrMatrixRebuild == AO_TS_SET) {
                for (unsigned i = 1; i < nThreads; i++) {
                        SAFEDELETE(thread_data[i].pJacHdl);
                        thread_data[i].pJacHdl = 0;
                }

                CCReady = CC_NO;

                throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
        }

        for (unsigned i = 0; i < nThreads; ++i) {
             if (thread_data[i].except) {
                  std::rethrow_exception(thread_data[i].except);
             }
        }
        
        for (unsigned i = 1; i < nThreads; i++) {
                pMH->AddUnchecked(*thread_data[i].pJacHdl);
        }
}

#ifdef USE_NAIVE_MULTITHREAD
void MultiThreadDataManager::NaiveAssJacInit(NaiveMatrixHandler& JacHdl, doublereal dCoef)
{
     while (true) {
          if (thread_data[0].ppNaiveJacHdl) {
               if (&JacHdl == thread_data[0].ppNaiveJacHdl[0]) {
                    /*
                     * NOTE: this test is here to detect whether JacHdl changed
                     * (typically it changes any time Solver::SetupSolmans() is called)
                     * However, just looking at the pointer does not suffice;
                     * we also need to check that "perm" did not change.
                     * For this purpose, we compare the address of "perm"
                     * in JacHdl and in the MatrixHandler of the second thread.
                     */
                    NaivePermMatrixHandler *pNaivePermJacHdl = dynamic_cast<NaivePermMatrixHandler *>(&JacHdl);
                    bool bDoAssemble(false);
                    if (!pNaivePermJacHdl) {
                         bDoAssemble = true;

                    } else {
                         NaivePermMatrixHandler *pNaivePermJacHdl2 = dynamic_cast<NaivePermMatrixHandler *>(thread_data[0].ppNaiveJacHdl[1]);
                         ASSERT(pNaivePermJacHdl2 != 0);
                         if (&pNaivePermJacHdl->GetPerm() == &pNaivePermJacHdl2->GetPerm()) {
                              bDoAssemble = true;
                         }
                    }

                    if (bDoAssemble) {
                         NaiveAssJac(JacHdl, dCoef);
                         break;
                    }
               }

               for (unsigned i = 1; i < nThreads; i++) {
                    if (thread_data[0].ppNaiveJacHdl[i]) {
                         SAFEDELETE(thread_data[0].ppNaiveJacHdl[i]);
                         thread_data[0].ppNaiveJacHdl[i] = 0;
                    }
               }

               SAFEDELETEARR(thread_data[0].lock);
               thread_data[0].lock = 0;
          }
          /* use JacHdl as matrix for the first thread,
           * and create copies for the other threads;
           * each thread sees the array of all the matrices,
           * and uses only its own for element assembly,
           * all for per-thread matrix summation */
          SAFENEWARR(thread_data[0].lock, AO_TS_t, JacHdl.iGetNumRows());

          for (integer i = 0; i < JacHdl.iGetNumRows(); i++) {
               thread_data[0].lock[i] = AO_TS_INITIALIZER;
          }

          if (thread_data[0].ppNaiveJacHdl) {
               SAFEDELETEARR(thread_data[0].ppNaiveJacHdl);
          }

          thread_data[0].ppNaiveJacHdl = 0;
          SAFENEWARR(thread_data[0].ppNaiveJacHdl,
                     NaiveMatrixHandler*, nThreads);
          thread_data[0].ppNaiveJacHdl[0] = &JacHdl;

          NaivePermMatrixHandler *pNaivePermJacHdl = dynamic_cast<NaivePermMatrixHandler *>(&JacHdl);

          for (unsigned i = 1; i < nThreads; i++) {
               thread_data[i].lock = thread_data[0].lock;
               thread_data[i].ppNaiveJacHdl = thread_data[0].ppNaiveJacHdl;
               thread_data[0].ppNaiveJacHdl[i] = 0;

               if (pNaivePermJacHdl) {
                    SAFENEWWITHCONSTRUCTOR(thread_data[0].ppNaiveJacHdl[i],
                                           NaivePermMatrixHandler,
                                           NaivePermMatrixHandler(JacHdl.iGetNumRows(),
                                                                  pNaivePermJacHdl->GetPerm(),
                                                                  pNaivePermJacHdl->GetInvPerm()));

               } else {
                    SAFENEWWITHCONSTRUCTOR(thread_data[0].ppNaiveJacHdl[i],
                                           NaiveMatrixHandler,
                                           NaiveMatrixHandler(JacHdl.iGetNumRows()));
               }
          }
     }
}

void
MultiThreadDataManager::NaiveAssJac(NaiveMatrixHandler& JacHdl, doublereal dCoef)
{
        ASSERT(thread_data != NULL);

        /* Assemble per-thread matrix */
        thread_data[0].ElemIter.ResetAccessData();
        op = MultiThreadDataManager::OP_ASSJAC_NAIVE;
        thread_count = nThreads - 1;

        for (unsigned i = 0; i < nThreads; ++i) {
             thread_data[i].except = std::exception_ptr{};
        }

        for (unsigned i = 1; i < nThreads; i++) {
                thread_data[i].dCoef = dCoef;

                sem_post(&thread_data[i].sem);
        }

        /* FIXME Right now it's already done before calling AssJac;
         * needs be moved here to improve parallel performances... */
#if 0
        thread_data[0].ppNaiveJacHdl[0]->Reset();
#endif
        try {
             DataManager::AssJac(*thread_data[0].ppNaiveJacHdl[0],
                                 dCoef,
                                 thread_data[0].ElemIter,
                                 *thread_data[0].pWorkMat);
        } catch (...) {
             thread_data[0].except = std::current_exception();
        }
        
        pthread_mutex_lock(&thread_mutex);
        if (thread_count > 0) {
                pthread_cond_wait(&thread_cond, &thread_mutex);
        }
        pthread_mutex_unlock(&thread_mutex);

        for (unsigned i = 0; i < nThreads; ++i) {
             if (thread_data[i].except) {
                  std::rethrow_exception(thread_data[i].except);
             }
        }
        
        /* Sum per-thread matrices */
        op = MultiThreadDataManager::OP_SUM_NAIVE;
        thread_count = nThreads - 1;
        for (unsigned i = 1; i < nThreads; i++) {
                sem_post(&thread_data[i].sem);
        }

        NaiveMatrixHandler* to = thread_data[0].ppNaiveJacHdl[0];
        integer nn = to->iGetNumRows();
        integer iFrom = 0;
        integer iTo = nn/nThreads;
        for (unsigned matrix = 1; matrix < nThreads; matrix++) {
                NaiveMatrixHandler* from = thread_data[0].ppNaiveJacHdl[matrix];
                naivepsad(to->ppdRows,
                                to->ppiRows, to->piNzr,
                                to->ppiCols, to->piNzc, to->ppnonzero,
                                from->ppdRows, from->ppiCols, from->piNzc,
                                iFrom, iTo, thread_data[0].lock);
        }

        pthread_mutex_lock(&thread_mutex);
        if (thread_count > 0) {
                pthread_cond_wait(&thread_cond, &thread_mutex);
        }

        pthread_mutex_unlock(&thread_mutex);
}
#endif

#ifdef USE_SPARSE_AUTODIFF
void
MultiThreadDataManager::GradAssJac(SpGradientSparseMatrixHandler& JacHdl, doublereal dCoef)
{
        ASSERT(thread_data != NULL);

        AO_CLEAR(&propagate_ErrMatrixRebuild);

        JacHdl.Reset(); // FIXME: Matrix cannot be reset in parallel by DataManager::AssJac

        thread_data[0].ElemIter.ResetAccessData();
        thread_data[0].oGradJacHdl.SetMatrixHandler(&JacHdl);
        op = MultiThreadDataManager::OP_ASSJAC_GRAD;
        thread_count = nThreads - 1;

        for (unsigned i = 0; i < nThreads; ++i) {
             thread_data[i].except = std::exception_ptr{};
        }
        
        for (unsigned i = 1; i < nThreads; i++) {
                thread_data[i].dCoef = dCoef;
                thread_data[i].oGradJacHdl.SetMatrixHandler(&JacHdl);
                sem_post(&thread_data[i].sem);
        }

        try {
             DataManager::AssJac(thread_data[0].oGradJacHdl, dCoef, thread_data[0].ElemIter,
                                 *thread_data[0].pWorkMat);
        } catch (...) {
             thread_data[0].except = std::current_exception();
        }
        
        pthread_mutex_lock(&thread_mutex);

        if (thread_count > 0) {
                pthread_cond_wait(&thread_cond, &thread_mutex);
        }

        pthread_mutex_unlock(&thread_mutex);

        for (unsigned i = 0; i < nThreads; ++i) {
             if (thread_data[i].except) {
                  std::rethrow_exception(thread_data[i].except);
             }
        }
}
#endif

void MultiThreadDataManager::SetAffinity(const ThreadData& oThread)
{
     if (oThread.iCPUIndex >= 0) {
          Task2CPU oCPUSet;

          pedantic_cerr("Setting affinity of thread " << oThread.threadNumber << " to CPU " << oThread.iCPUIndex << " ...\n");

          oCPUSet.SetCPU(oThread.iCPUIndex);

          if (!oCPUSet.bSetAffinity()) {
               silent_cerr("Failed to set affinity of thread " << oThread.threadNumber
                           << " to CPU " << oThread.iCPUIndex << "\n");
          }
     }
}

#ifdef MBDYN_X_MT_ASSRES
void
MultiThreadDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef, const VectorHandler* pAbsResHdl)
        /*throw(ChangedEquationStructure)*/
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
                        *thread_data[0].pWorkVec, *thread_data[0].pAbsResHdl);

        pthread_mutex_lock(&thread_mutex);
        if (thread_count > 0) {
                pthread_cond_wait(&thread_cond, &thread_mutex);
        }
        pthread_mutex_unlock(&thread_mutex);

        for (unsigned i = 1; i < nThreads; i++) {
                ResHdl += *thread_data[i].pResHdl;
        }
        if (pAbsResHdl) {
                for (unsigned i = 1; i < nThreads; i++) {
                        *thread_data[i].pAbsResHdl.AddAbsValuesTo(ResHdl);
                }
        }
}
#endif /* MBDYN_X_MT_ASSRES */

#endif /* USE_MULTITHREAD */
