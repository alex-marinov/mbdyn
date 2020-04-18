/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2018(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include "mbconfig.h"

#ifdef USE_PASTIX
#include <algorithm>
#include "pastixwrap.h"

PastixSolver::PastixSolver(SolutionManager* pSM, integer iDim, integer iNumIter, integer iNumThreads)
    :LinearSolver(pSM),
     iDim(iDim),
     iNumIter(iNumIter),
     iNumThreads(iNumThreads),
     Axp(0),
     Aip(0),
     App(0),
     pastix_data(0),
     perm(iDim, 0),
     iperm(iDim, 0),
     bDoOrdering(true),
     iNumNonZeros(-1)
{
    std::fill(&iparm[0], &iparm[IPARM_SIZE], 0);

    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK] = API_TASK_INIT;

    std::fill(&dparm[0], &dparm[DPARM_SIZE], 0.);
}

PastixSolver::~PastixSolver()
{
    if (pastix_data) {
        try {
            PastixCall(API_TASK_CLEAN, API_TASK_CLEAN);
        } catch (const std::exception& oErr) {
            silent_cerr("warning: Failed to clean up pastix memory: " << oErr.what() << std::endl);
        }
    }
}

#ifdef DEBUG
void PastixSolver::IsValid(void) const
{
    ASSERT(NO_ERR == CheckMatrix());    
}
#endif /* DEBUG */


pastix_int_t PastixSolver::PastixCall(const pastix_int_t iStartTask, const pastix_int_t iEndTask) const
{
    iparm[IPARM_START_TASK] = iStartTask;
    iparm[IPARM_END_TASK] = iEndTask;

    pastix(&pastix_data,
           MPI_COMM_WORLD,
           iDim,
           App,
           Aip,
           Axp,
           &perm[0],
           &iperm[0],
           pdSol,
           1,
           iparm,
           dparm);

    if (iparm[IPARM_ERROR_NUMBER] != NO_ERR) {
        static const struct {
            pastix_int_t iTask;
            char szTask[18];
        } rgTasks[8] = {
            {API_TASK_INIT,     "API_TASK_INIT"},
            {API_TASK_ORDERING, "API_TASK_ORDERING"},
            {API_TASK_SYMBFACT, "API_TASK_SYMBFACT"},
            {API_TASK_ANALYSE,  "API_TASK_ANALYSE"},
            {API_TASK_NUMFACT,  "API_TASK_NUMFACT"},
            {API_TASK_SOLVE,    "API_TASK_SOLVE"},
            {API_TASK_REFINE,   "API_TASK_REFINE"},
            {API_TASK_CLEAN,    "API_TASK_CLEAN"}
        };

        static const char szUnknownTask[] = "<unknown task>";
        const char* pszStartTask = szUnknownTask;
        const char* pszEndTask = szUnknownTask;

        for (auto i = std::begin(rgTasks); i != std::end(rgTasks); ++i) {
            if (i->iTask == iStartTask) {
                pszStartTask = i->szTask;
            }

            if (i->iTask == iEndTask) {
                pszEndTask = i->szTask;
            }
        }
        
        silent_cerr("PastixCall("
                    << pszStartTask << ", "
                    << pszEndTask
                    << ") failed with status "
                    << iparm[IPARM_ERROR_NUMBER] << std::endl);
        
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    return iparm[IPARM_ERROR_NUMBER];
}

pastix_int_t PastixSolver::CheckMatrix() const
{
    pastix_int_t status = pastix_checkMatrix(MPI_COMM_WORLD,
                                             API_VERBOSE_NOT,
                                             API_SYM_NO,
                                             API_NO,
                                             iDim,
                                             &App,
                                             &Aip,
                                             &Axp,
                                             0,
                                             1);

    if (NO_ERR != status) {
        silent_cerr("pastix_checkMatrix failed with status " << status << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    return status;
}

void PastixSolver::Solve(void) const
{
    // Right hand side will be overwritten by the solution by pastix
    std::copy(&pdRhs[0], &pdRhs[iDim], pdSol);
    
    if (pastix_data == 0) {
        CheckMatrix();
        
        iparm[IPARM_MODIFY_PARAMETER] = API_NO;

        PastixCall(API_TASK_INIT, API_TASK_INIT);

        iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
        iparm[IPARM_SYM] = API_SYM_NO;
        iparm[IPARM_THREAD_NBR] = iNumThreads;
        iparm[IPARM_FACTORIZATION] = API_FACT_LU;
        iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
        iparm[IPARM_ORDERING] = API_ORDER_SCOTCH;
        iparm[IPARM_BINDTHRD] = API_BIND_NO;
        iparm[IPARM_ITERMAX] = iNumIter;
        iparm[IPARM_RHS_MAKING] = API_RHS_B;        
    } else {
        ASSERT(NO_ERR == CheckMatrix());
    }

    pastix_int_t iStartTask;

    if (pastix_data == 0 || bDoOrdering) {
        // First call or nonezero pattern has changed
        iStartTask = API_TASK_ORDERING;
    } else if (bHasBeenReset) {
        // Assume same nonezero pattern (e.g. full Newton Raphson iteration)
        iStartTask = API_TASK_NUMFACT;
    } else {
        // Same matrix (e.g. modified Newton Raphson iteration)
        iStartTask = API_TASK_SOLVE;
    }

    pastix_int_t iEndTask = (iNumIter > 0) ? API_TASK_REFINE : API_TASK_SOLVE;
    
    PastixCall(iStartTask, iEndTask);

    if (iNumIter > 0 && iparm[IPARM_NBITER] >= iNumIter) {
        silent_cerr("Pastix warning: iparm[IPARM_NBITER] = " << iparm[IPARM_NBITER]
                    << " >= iparm[IPARM_ITERMAX] = " << iparm[IPARM_ITERMAX] <<  std::endl
                    << "Pastix warning: dparm[DPARM_RELATIVE_ERROR] = " << dparm[DPARM_RELATIVE_ERROR] << std::endl
                    << "Pastix warning: The flag \"max iterations\" of the linear solver is too small or the condition number of the Jacobian matrix is too high" << std::endl);
    }
    
    bDoOrdering = false;
    bHasBeenReset = false;
}

void PastixSolver::MakeCompactForm(SparseMatrixHandler& mh,
                                   std::vector<doublereal>& Ax,
                                   std::vector<integer>& Ai,
                                   std::vector<integer>& Adummy,
                                   std::vector<integer>& Ap) const
{
    if (!bHasBeenReset) {
        return;
    }

    mh.MakeCompressedColumnForm(Ax, Ai, Ap, 1);
    
    if (iNumNonZeros != mh.Nz()) {
        // Force a new ordering step (e.g. needed if we are using a SpMapMatrixHandler)
        bDoOrdering = true;
    }

    iNumNonZeros = mh.Nz();
    
    Axp = &Ax[0];
    Aip = &Ai[0];
    App = &Ap[0];
}

PastixSolutionManager::PastixSolutionManager(integer iDim, integer iNumThreads, integer iNumIter, const ScaleOpt& scale)
    :x(iDim),
     b(iDim),
     xVH(iDim, &x[0]),
     bVH(iDim, &b[0]),
     scale(scale),
     pMatScale(0),
     A(iDim)
{
    SAFENEWWITHCONSTRUCTOR(pLS,
                           PastixSolver,
                           PastixSolver(this, iDim, iNumIter, iNumThreads));

    pLS->pdSetResVec(&b[0]);
    pLS->pdSetSolVec(&x[0]);
}

PastixSolutionManager::~PastixSolutionManager(void)
{
    if (pMatScale) {
        SAFEDELETE(pMatScale);
        pMatScale = 0;
    }
}

#ifdef DEBUG
void PastixSolutionManager::IsValid(void) const
{
    ASSERT(b.size() == x.size());
    ASSERT(A.iGetNumRows() == A.iGetNumCols());
    ASSERT(b.size() == A.iGetNumRows());
    ASSERT(x.size() == A.iGetNumCols());
    
    pLS->IsValid();
}
#endif /* DEBUG */

void PastixSolutionManager::MatrReset(void)
{
    pLS->Reset();
}

void PastixSolutionManager::MatrInitialize(void)
{
    MatrReset();
}

void PastixSolutionManager::ForceSymmetricGraph(SpMapMatrixHandler& A) const
{
    if (pLS->bReset()) {
        for (auto i = A.begin(); i != A.end(); ++i) {
            // Metis and Scotch ordering require a symmetric graph
            // (e.g. for each matrix element A(i, j) another element A(j, i) must exist)
            A(i->iCol + 1, i->iRow + 1);
        }
    }
}

void PastixSolutionManager::MakeCompressedColumnForm(void)
{
    ForceSymmetricGraph(A);
    ScaleMatrixAndRightHandSide<SpMapMatrixHandler>(A);   
    pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

template <typename MH>
void PastixSolutionManager::ScaleMatrixAndRightHandSide(MH& mh)
{
    if (scale.when != SCALEW_NEVER) {
        MatrixScale<MH>& rMatScale = GetMatrixScale<MH>();

        if (pLS->bReset()) {
            if (!rMatScale.bGetInitialized()
                || scale.when == SolutionManager::SCALEW_ALWAYS) {
                // (re)compute
                rMatScale.ComputeScaleFactors(mh);
            }
            // in any case scale matrix and right-hand-side
            rMatScale.ScaleMatrix(mh);

            if (silent_err) {
                rMatScale.Report(std::cerr);
            }
        }

        rMatScale.ScaleRightHandSide(bVH);
    }
}

template <typename MH>
MatrixScale<MH>& PastixSolutionManager::GetMatrixScale()
{
    if (pMatScale == 0) {
        pMatScale = MatrixScale<MH>::Allocate(scale);
    }

    // Will throw std::bad_cast if the type does not match
    return dynamic_cast<MatrixScale<MH>&>(*pMatScale);
}

void PastixSolutionManager::ScaleSolution(void)
{
    if (scale.when != SCALEW_NEVER) {
        ASSERT(pMatScale != 0);

        pMatScale->ScaleSolution(xVH);
    }
}

void PastixSolutionManager::Solve(void)
{
    MakeCompressedColumnForm();
    
    pLS->Solve();

    ScaleSolution();
}

MatrixHandler* PastixSolutionManager::pMatHdl(void) const
{
    return &A;
}

VectorHandler* PastixSolutionManager::pResHdl(void) const
{
    return &bVH;
}

VectorHandler* PastixSolutionManager::pSolHdl(void) const
{
    return &xVH;
}

template <class CC>
PastixCCSolutionManager<CC>::PastixCCSolutionManager(integer iDim,
                                                     integer iNumThreads,
                                                     integer iNumIter,
                                                     const ScaleOpt& scale)
    : PastixSolutionManager(iDim, iNumThreads, iNumIter, scale),
      CCReady(false),
      Ac(0)
{

}

template <class CC>
PastixCCSolutionManager<CC>::~PastixCCSolutionManager(void) 
{
    if (Ac) {
        SAFEDELETE(Ac);
    }
}

template <class CC>
void
PastixCCSolutionManager<CC>::MatrReset(void)
{
    pLS->Reset();
}

template <class CC>
void
PastixCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
    if (!CCReady) {
        ForceSymmetricGraph(A);
        pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

        if (Ac == 0) {
            SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
        }

        CCReady = true;
    }

    ScaleMatrixAndRightHandSide(*Ac);
}

template <class CC>
void
PastixCCSolutionManager<CC>::MatrInitialize()
{
    CCReady = false;

    if (Ac) {
        // If a DirCColMatrixHandler is in use and matrix scaling is enabled
        // an uncaught exception (MatrixHandler::ErrRebuildMatrix) will be thrown
        // if zero entries in the matrix become nonzero.
        // For that reason we have to reinitialize Ac!
        SAFEDELETE(Ac);
        Ac = 0;
    }

    MatrReset();
    pGetSolver()->Initialize();
}
	
template <class CC>
MatrixHandler*
PastixCCSolutionManager<CC>::pMatHdl(void) const
{
    if (!CCReady) {
        return &A;
    }

    ASSERT(Ac != 0);
    return Ac;
}

template class PastixCCSolutionManager<CColMatrixHandler<1> >;
template class PastixCCSolutionManager<DirCColMatrixHandler<1> >;

#endif
