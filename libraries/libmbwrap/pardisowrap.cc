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
  Copyright (C) 2021(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include "mbconfig.h"

#ifdef USE_PARDISO

#include <algorithm>
#include <cstring>
#include <functional>
#include <memory>

#ifdef USE_MKL_SET_NUM_THREADS_LOCAL
#include <mkl/mkl_service.h>
#endif

#include "pardisowrap.h"
#include "cscmhtpl.h"

template <typename MKL_INT_TYPE>
PardisoSolver<MKL_INT_TYPE>::PardisoSolver(SolutionManager* pSM, integer iDim, doublereal dPivot, integer iNumThreads, integer iNumIter, integer iVerbose)
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
      AX(iDim)
{
     std::fill(std::begin(iparm), std::end(iparm), 0);
     std::fill(std::begin(pt), std::end(pt), nullptr);

     iparm[0] = 1; // Use default values.
     iparm[1] = 3; // The parallel (OpenMP) version of the nested dissection algorithm
     iparm[7] = iNumIter; // Iterative refinement step
     iparm[9] = -log10(dPivot); // Pivoting perturbation
     iparm[10] = 1; // Enable builtin scaling
     iparm[12] = 1; // Improved accuracy using (non-) symmetric weighted matching.

#ifdef USE_MKL_SET_NUM_THREADS_LOCAL
     mkl_set_num_threads_local(iNumThreads);
#endif
}

template <typename MKL_INT_TYPE>
PardisoSolver<MKL_INT_TYPE>::~PardisoSolver()
{
     MKL_INT_TYPE ierror;

     phase = -1;

     SolverType::pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &ierror);

     ASSERT(ierror == 0);
}

template <typename MKL_INT_TYPE>
void PardisoSolver<MKL_INT_TYPE>::Solve(void) const

{
    const MKL_INT_TYPE iNumNzA = pAp[n] - pAp[0];

    MKL_INT_TYPE ierror = 0;
    constexpr doublereal dTolBW = sqrt(std::numeric_limits<doublereal>::epsilon());
    doublereal dErrBW = std::numeric_limits<doublereal>::max();
    bool bSymbolicFactor = iNumNzA != iNumNz;
    bool bRepeatFactor = false;
    
    do {
         if (bRepeatFactor) {
              bHasBeenReset = bSymbolicFactor = true;
         }
      
         if (bHasBeenReset) {        
              if (bSymbolicFactor) {
                   phase = 11; // Analysis step

                   SolverType::pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

                   if (ierror != 0) {
                        silent_cerr("Pardiso symbolic factorization failed with status " << ierror << "\n");
                        throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
                   }

                   iNumNz = iNumNzA;
              }
          
              phase = 22; // Numerical factorization step

              SolverType::pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

              if (ierror != 0) {
                   silent_cerr("Pardiso numeric factorization failed with status " << ierror << "\n");
                   throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
              }

              bHasBeenReset = false;
         }

         phase = 33; // Solve, iterative refinement step

         SolverType::pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, pAx, pAp, pAi, nullptr, &nrhs, iparm, &msglvl, pdRhs, pdSol, &ierror);

         const CSCMatrixHandlerTpl<doublereal, MKL_INT_TYPE, 1> A(pAx, pAi, pAp, n, iNumNzA);

         const MyVectorHandler X(n, pdSol);
     
         A.MatTVecMul(AX, X);

         const doublereal* const pAX = AX.pdGetVec();
         doublereal normAXpB = 0., normAXmB = 0.;

         for (MKL_INT_TYPE i = 0; i < n; ++i) {
              normAXpB += std::pow(pAX[i] + pdRhs[i], 2);
              normAXmB += std::pow(pAX[i] - pdRhs[i], 2);
         }

         dErrBW = sqrt(normAXmB / normAXpB);

         if (msglvl) {
              silent_cerr("Pardiso: backward error: " << dErrBW << "\n");
         }

         bRepeatFactor = (ierror != 0 || dErrBW > dTolBW);
    } while (bRepeatFactor && !bSymbolicFactor);

    if (ierror) {
         silent_cerr("Pardiso solve failed with status " << ierror << "\n");
         throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
    }
}

template <typename MKL_INT_TYPE>
MKL_INT_TYPE PardisoSolver<MKL_INT_TYPE>::MakeCompactForm(SparseMatrixHandler& mh,
                                                          std::vector<doublereal>& Ax,
                                                          std::vector<MH_INT_TYPE>& Ai,
                                                          std::vector<MH_INT_TYPE>& Ap) const
{
     MH_INT_TYPE iNumNzA = mh.MakeCompressedRowForm(Ax, Ai, Ap, 1);

     n = mh.iGetNumCols();

     pAx = &Ax.front();

     static_assert(sizeof(MH_INT_TYPE) == sizeof(MKL_INT_TYPE));
     static_assert(std::numeric_limits<MH_INT_TYPE>::is_integer);
     static_assert(std::numeric_limits<MKL_INT_TYPE>::is_integer);
     
     pAi = reinterpret_cast<MKL_INT_TYPE*>(&Ai.front());
     pAp = reinterpret_cast<MKL_INT_TYPE*>(&Ap.front());

     return iNumNzA;
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::PardisoSolutionManager(integer iDim, doublereal dPivot, integer iNumThreads, integer iNumIter, integer iVerbose)
    :x(iDim),
     b(iDim),
     A(iDim, iDim)
{
    SAFENEWWITHCONSTRUCTOR(pLS,
                           PardisoSolver<MKL_INT_TYPE>,
                           PardisoSolver<MKL_INT_TYPE>(this, iDim, dPivot, iNumThreads, iNumIter, iVerbose));

    pLS->pdSetResVec(b.pdGetVec());
    pLS->pdSetSolVec(x.pdGetVec());
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::~PardisoSolutionManager(void)
{

}

#ifdef DEBUG
template <typename MatrixHandlerType, typename MKL_INT_TYPE>
void PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::IsValid() const
{
}
#endif

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
void PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::MatrReset(void)
{
    pLS->Reset();
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
void PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::MatrInitialize(void)
{
    MatrReset();
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
void PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::MakeCompressedRowForm(void)
{
    pGetSolver()->MakeCompactForm(A, Ax, Ai, Ap);
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
void PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::Solve(void)
{
    MakeCompressedRowForm();

    pLS->Solve();
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
MatrixHandler* PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::pMatHdl(void) const
{
    return &A;
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
VectorHandler* PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::pResHdl(void) const
{
    return &b;
}

template <typename MatrixHandlerType, typename MKL_INT_TYPE>
VectorHandler* PardisoSolutionManager<MatrixHandlerType, MKL_INT_TYPE>::pSolHdl(void) const
{
    return &x;
}

template class PardisoSolutionManager<SpMapMatrixHandler, MKL_INT>;
template class PardisoSolutionManager<SpMapMatrixHandler, long long>;
template class PardisoSolutionManager<SpGradientSparseMatrixHandler, MKL_INT>;
template class PardisoSolutionManager<SpGradientSparseMatrixHandler, long long>;
#endif
