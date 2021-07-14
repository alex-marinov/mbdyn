/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2021
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
  Copyright (C) 2021(-2021) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include "mbconfig.h"

#ifdef USE_PARDISO

#include <algorithm>

#ifdef USE_OMP_SET_NUM_THREADS
#include <omp.h>
#endif

#include "pardisowrap.h"

PardisoSolver::PardisoSolver(SolutionManager* pSM, integer iDim, integer iNumThreads, integer iNumIter, const SolutionManager::ScaleOpt& scale, integer iVerbose)
     :LinearSolver(pSM),
      pAx(nullptr),
      pAi(nullptr),
      pAp(nullptr),
      iNumNz(-1),
      phase(-1),
      maxfct(1),
      mnum(1),
      mtype(11),
      n(-1),
      nrhs(1),
      msglvl(iVerbose)
{
     std::fill(std::begin(iparm), std::end(iparm), 0);
     std::fill(std::begin(pt), std::end(pt), nullptr);

     iparm[0] = 1; // Use default values.
     iparm[1] = 3; // The parallel (OpenMP) version of the nested dissection algorithm
     iparm[7] = iNumIter; // Iterative refinement step
     iparm[9] = 13; // Pivoting perturbation

     switch (scale.algorithm) {
     case SolutionManager::SCALEA_NONE:
          break;
     default:
          iparm[10] = 1; // Enable scaling
     }

     iparm[12] = 1; // Improved accuracy using (non-) symmetric weighted matching.

#ifdef USE_OMP_SET_NUM_THREADS
     omp_set_num_threads(iNumThreads);
#endif
}

PardisoSolver::~PardisoSolver()
{
     MKL_INT ierror;

     phase = -1;

     pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &ierror);

     ASSERT(ierror == 0);
}

void PardisoSolver::Solve(void) const

{
     const MKL_INT iNumNzA = pAp[n] - pAp[0];
     const bool bForceOrderingStep = iNumNz != iNumNzA;

     MKL_INT ierror;

     if (bHasBeenReset) {
          if (bForceOrderingStep) {
               phase = 11; // Analysis step

               pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

               if (ierror != 0) {
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }

               iNumNz = iNumNzA;
          }

          phase = 22; // Numerical factorization step

          pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

          if (ierror != 0) {
               throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
          }

          bHasBeenReset = false;
     }

     phase = 33; // Solve, iterative refinement step

     pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

     if (ierror) {
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

MKL_INT PardisoSolver::MakeCompactForm(SparseMatrixHandler& mh,
                                       std::vector<doublereal>& Ax,
                                       std::vector<MKL_INT>& Ai,
                                       std::vector<MKL_INT>& Ap) const
{
     MKL_INT iNumNzA = mh.MakeCompressedRowForm(Ax, Ai, Ap, 1);

     n = mh.iGetNumCols();

     pAx = &Ax.front();
     pAi = &Ai.front();
     pAp = &Ap.front();

     return iNumNzA;
}

template <typename MatrixHandlerType>
PardisoSolutionManager<MatrixHandlerType>::PardisoSolutionManager(integer iDim, integer iNumThreads, integer iNumIter, const ScaleOpt& scale, integer iVerbose)
    :x(iDim),
     b(iDim),
     A(iDim, iDim)
{
    SAFENEWWITHCONSTRUCTOR(pLS,
                           PardisoSolver,
                           PardisoSolver(this, iDim, iNumThreads, iNumIter, scale, iVerbose));

    pLS->pdSetResVec(b.pdGetVec());
    pLS->pdSetSolVec(x.pdGetVec());
}

template <typename MatrixHandlerType>
PardisoSolutionManager<MatrixHandlerType>::~PardisoSolutionManager(void)
{

}

#ifdef DEBUG
template <typename MatrixHandlerType>
void PardisoSolutionManager<MatrixHandlerType>::IsValid() const
{
}
#endif

template <typename MatrixHandlerType>
void PardisoSolutionManager<MatrixHandlerType>::MatrReset(void)
{
    pLS->Reset();
}

template <typename MatrixHandlerType>
void PardisoSolutionManager<MatrixHandlerType>::MatrInitialize(void)
{
    MatrReset();
}

template <typename MatrixHandlerType>
void PardisoSolutionManager<MatrixHandlerType>::MakeCompressedRowForm(void)
{
     pGetSolver()->MakeCompactForm(A, Ax, Ai, Ap);
}

template <typename MatrixHandlerType>
void PardisoSolutionManager<MatrixHandlerType>::Solve(void)
{
    MakeCompressedRowForm();

    pLS->Solve();
}

template <typename MatrixHandlerType>
MatrixHandler* PardisoSolutionManager<MatrixHandlerType>::pMatHdl(void) const
{
    return &A;
}

template <typename MatrixHandlerType>
VectorHandler* PardisoSolutionManager<MatrixHandlerType>::pResHdl(void) const
{
    return &b;
}

template <typename MatrixHandlerType>
VectorHandler* PardisoSolutionManager<MatrixHandlerType>::pSolHdl(void) const
{
    return &x;
}

template class PardisoSolutionManager<SpMapMatrixHandler>;

#ifdef USE_SPARSE_AUTODIFF
template class PardisoSolutionManager<SpGradientSparseMatrixHandler>;
#endif

#endif
