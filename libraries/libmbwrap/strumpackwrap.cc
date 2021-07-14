/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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
  Copyright (C) 2020(-2020) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/


#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_STRUMPACK

#include "myassert.h"
#include "linsol.h"
#include "strumpackwrap.h"

StrumpackSolver::StrumpackSolver(SolutionManager* pSolMan,
				 integer iNumThreads,
				 integer iNumIter,
				 unsigned uSolverFlags,
				 integer iVerbose)
     :LinearSolver(pSolMan), iNumNonZeros(-1), bDoOrderingStep(true)
{
     oSolver.options().set_verbose(iVerbose);
     oSolver.options().set_maxit(std::max(1, iNumIter));

     const unsigned uPermFlag = uSolverFlags & LinSol::SOLVER_FLAGS_PERM_MASK;
     
     switch (uPermFlag) {
     case LinSol::SOLVER_FLAGS_ALLOWS_SCOTCH:
	  oSolver.options().set_reordering_method(strumpack::ReorderingStrategy::SCOTCH);
	  break;
     case LinSol::SOLVER_FLAGS_ALLOWS_METIS:
	  oSolver.options().set_reordering_method(strumpack::ReorderingStrategy::METIS);
	  break;
     }

     const unsigned uCompressionFlag = uSolverFlags & LinSol::SOLVER_FLAGS_COMPRESSION_MASK;

     switch (uCompressionFlag) {
     case LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_HSS:
	  oSolver.options().set_compression(strumpack::CompressionType::HSS);
	  break;
     case LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_BLR:
	  oSolver.options().set_compression(strumpack::CompressionType::BLR);
	  break;
     case LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_HODLR:
	  oSolver.options().set_compression(strumpack::CompressionType::HODLR);
	  break;
     }
     
     strumpack::params::num_threads = iNumThreads;
}

StrumpackSolver::~StrumpackSolver(void)
{
}

void StrumpackSolver::Reset(void)
{
     bHasBeenReset = true;
}

void StrumpackSolver::ResetSymbolic()
{
     bDoOrderingStep = true;
}
     
void StrumpackSolver::MakeCompactForm(class SparseMatrixHandler& mh,
				      std::vector<doublereal>& Ax,
				      std::vector<integer>& Ai,
				      std::vector<integer>&,
				      std::vector<integer>& Ap) const
{
     mh.MakeCompressedRowForm(Ax, Ai, Ap);
     
     oSolver.set_csr_matrix(mh.iGetNumCols(), &Ap.front(), &Ai.front(), &Ax.front());

     if (mh.Nz() != iNumNonZeros) {
	  bDoOrderingStep = true;
	  iNumNonZeros = mh.Nz();
     }
}

void StrumpackSolver::Solve(void) const
{
     if (bDoOrderingStep) {
	  oSolver.reorder();
	  bDoOrderingStep = false;
     }
     
     if (bHasBeenReset) {
	  oSolver.factor();
	  bHasBeenReset = false;
     }
     
     oSolver.solve(pdGetResVec(), pdGetSolVec());
}

template <typename MatrixHandlerType>
StrumpackSolutionManager<MatrixHandlerType>::StrumpackSolutionManager(integer Dim,
								      integer iNumThreads,
								      integer iNumIter,
								      const ScaleOpt& scale,
								      unsigned uSolverFlags,
								      integer iVerbose)
     :A(Dim, Dim), x(Dim), b(Dim), pMatScale(nullptr), scale(scale)
{
     SAFENEWWITHCONSTRUCTOR(pLS,
			    StrumpackSolver,
			    StrumpackSolver(this, iNumThreads, iNumIter, uSolverFlags, iVerbose));

     pLS->pdSetResVec(b.pdGetVec());
     pLS->pdSetSolVec(x.pdGetVec());
}

template <typename MatrixHandlerType>
StrumpackSolutionManager<MatrixHandlerType>::~StrumpackSolutionManager(void)
{
     if (pMatScale) {
	  SAFEDELETE(pMatScale);
     }
}

#ifdef DEBUG
template <typename MatrixHandlerType>
void StrumpackSolutionManager<MatrixHandlerType>::IsValid(void) const
{
     pLS->IsValid();
}
#endif /* DEBUG */

template <typename MatrixHandlerType>
void StrumpackSolutionManager<MatrixHandlerType>::MatrReset(void)
{
     pLS->Reset();
}

template <typename MatrixHandlerType>
void StrumpackSolutionManager<MatrixHandlerType>::MatrInitialize()
{
     pGetSolver()->ResetSymbolic();     
}

template <typename MatrixHandlerType>
void StrumpackSolutionManager<MatrixHandlerType>::Solve(void)
{
     std::vector<integer>* pDummy = nullptr;

     ScaleMatrixAndRightHandSide(A);
     
     pLS->MakeCompactForm(A, Ax, Ai, *pDummy, Ap);
     
     pLS->Solve();

     ScaleSolution();
}

template <typename MatrixHandlerType>
MatrixHandler* StrumpackSolutionManager<MatrixHandlerType>::pMatHdl(void) const
{
     return &A;
}

template <typename MatrixHandlerType>
MyVectorHandler* StrumpackSolutionManager<MatrixHandlerType>::pResHdl(void) const
{
     return &b;
}

template <typename MatrixHandlerType>
MyVectorHandler* StrumpackSolutionManager<MatrixHandlerType>::pSolHdl(void) const
{
     return &x;
}

template <typename MatrixHandlerType>
StrumpackSolver* StrumpackSolutionManager<MatrixHandlerType>::pGetSolver() const
{
     return static_cast<StrumpackSolver*>(pLS);
}

template <typename MatrixHandlerType>
template <typename MH>
void StrumpackSolutionManager<MatrixHandlerType>::ScaleMatrixAndRightHandSide(MH& mh)
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

        rMatScale.ScaleRightHandSide(b);
    }
}

template <typename MatrixHandlerType>
template <typename MH>
MatrixScale<MH>& StrumpackSolutionManager<MatrixHandlerType>::GetMatrixScale()
{
    if (pMatScale == nullptr) {
        pMatScale = MatrixScale<MH>::Allocate(scale);
    }

    // Will throw std::bad_cast if the type does not match
    return dynamic_cast<MatrixScale<MH>&>(*pMatScale);
}

template <typename MatrixHandlerType>
void StrumpackSolutionManager<MatrixHandlerType>::ScaleSolution(void)
{
    if (scale.when != SCALEW_NEVER) {
        ASSERT(pMatScale != nullptr);

        pMatScale->ScaleSolution(x);
    }
}

template class StrumpackSolutionManager<SpMapMatrixHandler>;

#ifdef USE_SPARSE_AUTODIFF
template class StrumpackSolutionManager<SpGradientSparseMatrixHandler>;
#endif

#endif
