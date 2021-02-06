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
#include <omp.h>

#include "task2cpu.h"
#include "linsol.h"
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
      msglvl(iVerbose),
      bDoOrdering(true)
{
     std::fill(std::begin(iparm), std::end(iparm), 0);
     std::fill(std::begin(pt), std::end(pt), nullptr);

     iparm[0] = 1;
     iparm[7] = iNumIter;
     iparm[9] = 13;

     switch (scale.algorithm) {
     case SolutionManager::SCALEA_NONE:
     case SolutionManager::SCALEA_UNDEF:          
          break;
     default:
          iparm[10] = 1;
     }

     iparm[12] = 1;
     
     omp_set_num_threads(iNumThreads);
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
     MKL_INT ierror;

     static_assert(sizeof(MKL_INT) == 4);
     
     if (bHasBeenReset) {
          if (bDoOrdering) {
               phase = 11;
          
               pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

               if (ierror != 0) {
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }

               bDoOrdering = false;
          }

          phase = 22;
          pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

          if (ierror != 0) {
               throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
          }

          bHasBeenReset = false;
     }
     
     phase = 33;

     pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

     if (ierror != 0) {
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

MKL_INT PardisoSolver::MakeCompactForm(SparseMatrixHandler& mh,
                                       std::vector<doublereal>& Ax,
                                       std::vector<MKL_INT>& Ai,
                                       std::vector<MKL_INT>& Ap) const
{
     MKL_INT iNumNzCurr = mh.MakeCompressedRowForm(Ax, Ai, Ap, 1);

     if (iNumNz != iNumNzCurr) {
          bDoOrdering = true;
     }
     
     iNumNz = iNumNzCurr;
     n = mh.iGetNumCols();
     
     pAx = &Ax.front();
     pAi = &Ai.front();
     pAp = &Ap.front();
     
     return iNumNz;
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
