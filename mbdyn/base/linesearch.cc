/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
  * Portions Copyright (C) 2003-2017
  * Giuseppe Quaranta   <quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare
  */

 /*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2011(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
  */

/*
  References:
  Numerical recipes in C: the art of scientific computing / William H. Press [et al.]. – 2nd ed.
  ISBN 0-521-43108-5
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <unistd.h>
#include "ls.h"
#include "solver.h"
#include "linesearch.h"
#include "solman.h"
#ifdef USE_MPI
#include "mbcomm.h"
#include "schsolman.h"
#endif /* USE_MPI */

#include "dofown.h"
#include "output.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>

#ifdef DEBUG
#undef ASSERT
#define ASSERT(expr) assert(expr)
#define TRACE(expr) silent_cerr(__FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ":" << expr)
#else
#define TRACE(expr) static_cast<void>(0)
#endif

#define TRACE_VAR(varname) TRACE(#varname "=" << (varname) << '\n')
#define TRACE_FLAG(varname, flag) TRACE(#varname "&" #flag "=" << ((varname) & (flag)) << '\n')

LineSearchParameters::LineSearchParameters(void)
     :CommonNonlinearSolverParam(DEFAULT_FLAGS, 0, false),
       dTolX(1e-7),
       dTolMin(1e-8),
       iMaxIterations(200),
       dMaxStep(100.),
       dAlphaFull(1e-4),
       dAlphaModified(0.6),
       dLambdaMin(1e-2),
       dLambdaFactMin(1e-1),
       dDivergenceCheck(1.),
       dMinStepScale(1e-3),
       dTimeStepTol(0.1),
       dUpdateRatio(0.05),
       dMcpTol(DBL_EPSILON),
       dMcpSigma(0.9),
       dMcpRho(1e-8),
       dMcpP(2.1)
{
     NO_OP;
}

LineSearchSolver::LineSearchSolver(DataManager* pDM,
                                   const NonlinearSolverTestOptions& options,
                                   const struct LineSearchParameters& param)
     : NonlinearSolver(options),
       LineSearchParameters(param),
       lambda2(std::numeric_limits<doublereal>::max()),
       f2(std::numeric_limits<doublereal>::max()),
       pRes(0),
       pSol(0),
       pNLP(0),
       pS(0),
       pSM(0),
       pDM(pDM),
       iRebuildJac(0)
{
     TRACE_VAR(dTolX);
     TRACE_VAR(dTolMin);
     TRACE_VAR(iMaxIterations);
     TRACE_VAR(dMaxStep);
     TRACE_VAR(dAlphaFull);
     TRACE_VAR(dAlphaModified);
     TRACE_VAR(dLambdaMin);
     TRACE_VAR(dLambdaFactMin);
     TRACE_VAR(dMinStepScale);
     TRACE_VAR(dDivergenceCheck);
     TRACE_VAR(bHonorJacRequest);
     TRACE_VAR(uFlags);
     TRACE_FLAG(uFlags, ALGORITHM_CUBIC);
     TRACE_FLAG(uFlags, ALGORITHM_FACTOR);
     TRACE_FLAG(uFlags, SCALE_NEWTON_STEP);
     TRACE_FLAG(uFlags, RELATIVE_LAMBDA_MIN);
     TRACE_FLAG(uFlags, ZERO_GRADIENT_CONTINUE);
     TRACE_FLAG(uFlags, DIVERGENCE_CHECK);
     TRACE_FLAG(uFlags, PRINT_CONVERGENCE_INFO);
}

LineSearchSolver::~LineSearchSolver(void)
{
     NO_OP;
}

void LineSearchSolver::Attach(const NonlinearProblem* pNLP, Solver* pS)
{
     ASSERT(pNLP != NULL);
     ASSERT(pS != NULL);

     if (this->pNLP != pNLP) {
          ResetCond();
     }

     pSM = pS->pGetSolutionManager();
     this->pS = pS;
     this->pNLP = pNLP;
     pRes = pSM->pResHdl();
     pAbsRes = pGetResTest()->GetAbsRes();
     pSol = pSM->pSolHdl();
     Size = pRes->iGetSize();

     if (g.iGetSize() != Size) {
          TRACE("Resize temporary vectors ...\n");
          g.Resize(Size);
          p.Resize(Size);
     }
}

doublereal LineSearchSolver::dGetMinNewtonInc(const VectorHandler& dX) const {
     const VectorHandler& X = *pDM->GetpXCurr();
     const VectorHandler& XP = *pDM->GetpXPCurr();
     doublereal dTest = 0.;

     for (integer i = 1; i <= Size; ++i) {
          const doublereal Xi = pDM->GetDofType(i) == DofOrder::ALGEBRAIC ? X(i) : XP(i);
          
          const doublereal dTemp = std::abs(dX(i)) / std::max(std::abs(Xi), 1.);
          
          if (dTemp > dTest) {
               dTest = dTemp;
          }
     }

     return std::min(dTolX / dTest, 1.);
}

doublereal LineSearchSolver::dGetLambdaNext(doublereal lambda, const doublereal dSlope, const doublereal fPrev, const doublereal fCurr) {
     TRACE_VAR(dSlope);
     TRACE_VAR(fPrev);
     TRACE_VAR(fCurr);
     TRACE("dLambdaPrev=" << lambda << "\n");

     if (uFlags & ALGORITHM_CUBIC && std::isfinite(fCurr)) {
          doublereal tmplam;

          TRACE("Calculate new value for lambda ...\n");

          if (lambda == 1.) {
               tmplam = -dSlope / (2 * (fCurr - fPrev - dSlope));
               TRACE_VAR(tmplam);
          } else {
               const doublereal rhs1 = fCurr - fPrev - lambda * dSlope;
               const doublereal rhs2 = f2 - fPrev - lambda2 * dSlope;
               const doublereal a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
               const doublereal b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);

               if (a == 0.) {
                    tmplam = -dSlope / (2. * b);
                    TRACE_VAR(tmplam);
               } else {
                    const doublereal disc = b * b - 3. * a * dSlope;

                    if (disc < 0.) {
                         tmplam = 0.5 * lambda;
                         TRACE_VAR(tmplam);
                    } else if (b <= 0.) {
                         tmplam = (-b + sqrt(disc)) / (3. * a);
                         TRACE_VAR(tmplam);
                    } else {
                         tmplam = -dSlope / (b + sqrt(disc));
                         TRACE_VAR(tmplam);
                    }

                    if (tmplam > 0.5 * lambda) {
                         tmplam = 0.5 * lambda;
                         TRACE_VAR(tmplam);
                    }
               }
          }
          lambda2 = lambda;
          f2 = fCurr;
          TRACE_VAR(tmplam);
          lambda = std::max(tmplam, dLambdaFactMin * lambda);
     } else {
          lambda *= dLambdaFactMin;
     }

     TRACE("dLambda=" << lambda << "\n");

     return lambda;
}

bool LineSearchSolver::bCheckZeroGradient(doublereal fCurr,
                                          doublereal dErr,
                                          doublereal dTol,
                                          integer iIterCnt) const {
     const VectorHandler& X = *pDM->GetpXCurr();
     const VectorHandler& XP = *pDM->GetpXPCurr();

     doublereal test = 0.;
     const doublereal den = std::max(fCurr, 0.5 * Size);

     for (integer i = 1; i <= Size; ++i) {
          const doublereal absX = std::max(std::abs(X(i)),
                                           std::abs(XP(i)));
          const doublereal temp = std::abs(g(i)) * std::max(absX, 1.) / den;

          if (temp > test) {
               test = temp;
          }
     }

     if (test < dTolMin) { // we are at a local minimum of the function f
          if (uFlags & VERBOSE_MODE) {
               silent_cerr("line search warning: Zero gradient detected at time step "
                           << pDM->dGetTime() << " at iteration "
                           << iIterCnt << " (spurious convergence) test=" << test
                           << " < tol=" << dTolMin
                           << "\tErr(n)=" << dErr
                           << " > Tol = " << dTol << '\n');
          }

          if (uFlags & ZERO_GRADIENT_CONTINUE) {
               return true;
          } else {
               throw ZeroGradient(MBDYN_EXCEPT_ARGS);
          }
     } else {
          if (uFlags & VERBOSE_MODE) {
               silent_cerr("line search warning: lambda min"
                           " has been reached at time step " << pDM->dGetTime()
                           << " at iteration " << iIterCnt
                           << " but the gradient is not zero" << '\n');
          }

          if (uFlags & ABORT_AT_LAMBDA_MIN) {
               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }
     }

     return false;
}

bool LineSearchSolver::bCheckDivergence(const doublereal dErrFactor, const doublereal dErrPrev, const doublereal dErr, const integer iIterCnt) const
{
     const bool bDivergence = dErrFactor > dDivergenceCheck;

     if (bDivergence) {
          if (uFlags & DIVERGENCE_CHECK) {
               if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: The residual could not be decreased"
                                " at time step " << pDM->dGetTime()
                                << " at iteration " << iIterCnt << '\n');

                    if (iIterCnt > 1) {
                         silent_cerr("Err(n-1)=" << dErrPrev
                                     << "\tErr(n)=" << dErr
                                     << "\tErr(n)/Err(n-1)=" << dErr/dErrPrev
                                     << "\tErr(n)/Err(1)=" << dErrFactor << '\n');
                    }
               }
               throw ResidualNotDecreased(MBDYN_EXCEPT_ARGS);
          }
     }

     return bDivergence;
}

doublereal LineSearchSolver::dGetLambdaMin(doublereal& dSlope, const bool bRebuildJac, const VectorHandler& p, const integer iIterCnt, doublereal fCurr) const
{
     TRACE_VAR(dSlope);

     doublereal dLambdaMinCurr = -1.;

     if (dSlope >= 0. && fCurr > 0.) {
          if (uFlags & VERBOSE_MODE) {
               silent_cerr("line search warning: slope = " << dSlope << " is not negative at time step "
                           << pDM->dGetTime()
                           << " at iteration " << iIterCnt << '\n'
                           << "This could be a roundoff problem" << '\n');
          }

          if (uFlags & NON_NEGATIVE_SLOPE_CONTINUE) {
               // It seems to be a numerical problem.
               // Line search may not work in this situation.
               // Resort to the ordinary Newton Raphson algorithm
               dSlope = 0.;
               dLambdaMinCurr = 1.;
          } else {
               throw SlopeNotNegative(MBDYN_EXCEPT_ARGS);
          }
     }

     if (bRebuildJac) {
          if (dLambdaMinCurr < 0) {
               // dLambdaMinCurr has to be detected
               if (uFlags & RELATIVE_LAMBDA_MIN) {
                    dLambdaMinCurr = std::max(dLambdaMin, dGetMinNewtonInc(p));
               } else {
                    dLambdaMinCurr = dLambdaMin;
               }
          }
     } else {
          dLambdaMinCurr = 1.;
     }

     return dLambdaMinCurr;
}

doublereal LineSearchSolver::dGetMaxNewtonStep(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr) const
{
     return (uFlags & SCALE_NEWTON_STEP)
          ? dMaxStep * std::max(sqrt(XPrimeCurr.Dot() + XCurr.Dot()), static_cast<doublereal>(Size))
          : std::numeric_limits<doublereal>::max();
}

void LineSearchSolver::ScaleNewtonStep(const doublereal stpmax, VectorHandler& dX, const integer iIterCnt) const
{
     if (uFlags & SCALE_NEWTON_STEP) {
          const doublereal dNormSol = dX.Norm();

          if (dNormSol > stpmax) {
               const doublereal dScale = std::max(dMinStepScale, stpmax / dNormSol);
               if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: "
                                "Newton increment is reduced by factor " << dScale
                                << " at time step " << pDM->dGetTime()
                                << " at iteration " << iIterCnt
                                << " The time step is probably too large" << '\n');
               }
               ASSERT(stpmax >= 0);
               ASSERT(dScale <= 1.);
               ASSERT(dScale > 0.);
               dX *= dScale;
          }
     }
}

void LineSearchSolver::OutputIteration(const doublereal dErr, const integer iIterCnt, const bool bRebuildJac, const CPUTime& oCPU)
{
#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
     {
          if (outputIters() || outputSolverConditionNumber()) {
               if (outputIters()) {
                    silent_cout("\tIteration(" << iIterCnt << ") "
                                << std::setw(12) << dErr
                                << ' ' << (bRebuildJac ? 'J' : ' '));
               }

               if (outputSolverConditionNumber()) {
                    silent_cout(" cond=");
                    doublereal dCond;
                    if (pSM->bGetConditionNumber(dCond)) {
                         silent_cout(dCond);
                         if (outputSolverConditionStat()) {
                              AddCond(dCond);
                              silent_cout(" " << dGetCondMin() << " " << dGetCondMax() << " " << dGetCondAvg());
                         }
                    } else {
                         silent_cout("NA");
                    }
               }

               if (outputCPUTime()) {
                    silent_cout(" CPU:" << oCPU.Residual
                                << '+' << oCPU.Jacobian
                                << '+' << oCPU.LinearSolver);
               }

               silent_cout(std::endl); // Flush stdout only once per iteration
          }
     }
}

void LineSearchSolver::OutputLineSearch(const integer iIterCnt,
                                        const integer iLineSearchIter,
                                        const doublereal fCurr,
                                        const doublereal dErr,
                                        const doublereal dLambda,
                                        const doublereal dSlope) const
{
     if (outputIters() && (uFlags & PRINT_CONVERGENCE_INFO)) {
#ifdef USE_MPI
          if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
          {
               silent_cout("\t\tf(" << iIterCnt << ":"
                           << iLineSearchIter << ")="
                           << std::setw(12) << fCurr
                           << "\tErr=" << std::setw(12) << dErr
                           << "\tlambda=" << std::setw(12) << dLambda
                           <<"\tslope=" << dSlope << '\n');
          }
     }
}

void
LineSearchSolver::Residual(doublereal& f, integer iIterCnt)
{
#ifdef USE_EXTERNAL
     SendExternal();
#endif /* USE_EXTERNAL */

     pRes->Reset();

     if (pAbsRes != 0) {
          pAbsRes->Reset();
     }

     try {
          TRACE("Assemble residual\n");
          pNLP->Residual(pRes, pAbsRes);
     }
     catch (const SolutionDataManager::ChangedEquationStructure&) {
          if (bHonorJacRequest) {
               iRebuildJac = 0;
          }
     }

     if (outputRes()) {
          pS->PrintResidual(*pRes, iIterCnt);
     }

     f = 0.5 * pRes->Dot();
}

void
LineSearchSolver::Jacobian(MatrixHandler* pJac)
{
     SolutionManager *const pSM = pS->pGetSolutionManager();

     const integer iMaxIterRebuild = 10;

     pSM->MatrReset();
     integer iIter = 0;

     do {
          try {
               TRACE("Assemble Jacobian\n");
               pNLP->Jacobian(pJac);
          } catch (const MatrixHandler::ErrRebuildMatrix&) {
               silent_cout("LineSearchSolver: "
                           "rebuilding matrix..."
                           << '\n');

               /* need to rebuild the matrix... */
               pSM->MatrInitialize();
               continue;
          }
          break;
     } while (++iIter < iMaxIterRebuild);

     if (iIter >= iMaxIterRebuild) {
          silent_cerr("Maximum number of iterations exceeded when rebuilding the Jacobian matrix" << '\n');
          throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
     }

     TotJac++;

#ifdef USE_MPI
     if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
     {
          if (outputJac()) {
               silent_cout("Jacobian:" << '\n');

               if (silent_out) {
                    pJac->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
               }
          }

          if (outputMatrixConditionNumber()) {
               silent_cout(" cond=" << pJac->ConditionNumber(GetCondMatNorm()) << '\n');
          }
     }
}

LineSearchFull::LineSearchFull(DataManager* pDM,
                               const NonlinearSolverTestOptions& options,
                               const struct LineSearchParameters& param)
     :LineSearchSolver(pDM, options, param)
{
}

LineSearchFull::~LineSearchFull(void)
{
}

bool
LineSearchFull::LineSearch(const doublereal dMaxStep,
                           const doublereal fPrev,
                           doublereal& fCurr,
                           const integer iIterCnt)
{
     ScaleNewtonStep(dMaxStep, *pSol, iIterCnt);

     p = *pSol; // save the Newton increment

     doublereal dSlope = g.InnerProd(p);

     const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, true, p, iIterCnt, fCurr);

     doublereal dLambda = 1., dLambdaPrev = 0;
     integer iLineSearchIter = 0;
     doublereal dLambdaMax = 1.;

     SetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX, dLambdaMax);

     TRACE_VAR(dLambdaMinCurr);

     do {
          if (iLineSearchIter > 0) {
               TRACE("Start new step from Xold, XPold with lambda = " << dLambda << " ...\n");

               // scale the Newton increment by lambda and restore previous state
               for (integer i = 1; i <= Size; ++i) {
                    (*pSol)(i) = (dLambda - dLambdaPrev) * p(i);
               }
          }

          TRACE("Update the nonlinear problem ... pSol->Norm()=" << pSol->Norm() << '\n');
          SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, iLineSearchIter);
          SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, dLambda);

          pNLP->Update(pSol);

          Residual(fCurr, iIterCnt);

          if (iLineSearchIter == 0) {
               dLambdaMax = std::max(dLambdaMinCurr, GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX));

               ASSERT(dLambdaMax <= 1.);
               ASSERT(dLambdaMax >= 0.);
          }

          ++iLineSearchIter;

          TRACE("New value for f:" << fCurr << '\n');

          doublereal dErr = 0., dErrDiff = 0.;
          bool bResTestFinite = false;

          try {
               MakeResTest(pS, pNLP, *pRes, 0., dErr, dErrDiff);
               bResTestFinite = true;
          } catch (const NonlinearSolver::ErrSimulationDiverged&) {
               if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: residual test failed (residual=" << dErr << ")!"  << std::endl);
               }
          }

          OutputLineSearch(iIterCnt, iLineSearchIter, fCurr / fPrev, dErr, dLambda, dSlope);

          if (bResTestFinite && dLambda <= dLambdaMax) {
               pS->CheckTimeStepLimit(dErr, dErrDiff);

               if (fCurr <= fPrev + dAlphaFull * dLambda * dSlope) {
                    TRACE("Sufficient decrease in f: backtrack\n");
                    return false;
               } else if (dLambda <= dLambdaMinCurr) {
                    TRACE("Checking for convergence: lambda=" << dLambda << " < lambdaMin=" << dLambdaMinCurr << '\n');
                    return true; // check for convergence
               }
          }

          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }

          dLambdaPrev = dLambda; // save lambda; so we can restore the previous state
          dLambda = dGetLambdaNext(dLambdaPrev, dSlope, fPrev, fCurr);

          if (dLambdaMax < dLambda) {
               pedantic_cout("lambda reduced from " << dLambda << " to " << dLambdaMax
                             << " because of element request" << std::endl);
               dLambda = dLambdaMax;
          }

     } while (iLineSearchIter < iMaxIterations);

     if (uFlags & VERBOSE_MODE) {
          silent_cerr("line search warning: Maximum number of line search iterations="
                      << iMaxIterations
                      << " exceeded in line search at time step " << pDM->dGetTime()
                      << " at iteration " << iIterCnt << '\n');
     }

     throw MaxIterations(MBDYN_EXCEPT_ARGS);
}

void
LineSearchFull::Solve(const NonlinearProblem *pNonLinProblem,
                      Solver *pS,
                      const integer iMaxIter,
                      const doublereal& Tol,
                      integer& iIterCnt,
                      doublereal& dErr,
                      const doublereal& SolTol,
                      doublereal& dSolErr)
{
     Attach(pNonLinProblem, pS);
     iIterCnt = 0;
     dSolErr = 0.;
     dErr = 0.;
     doublereal dErrDiff = 0.;
     doublereal fCurr;
     CPUTime oCPU(*this);

     oCPU.Residual.Tic();

     Residual(fCurr, iIterCnt);

     TRACE("\t\tf(0) = " << fCurr << '\n');

     bool bResConverged = MakeResTest(pS, pNLP, *pRes, 1e-2 * Tol, dErr, dErrDiff); // use a more stringent test for the first iteration
     bool bSolConverged = pGetSolTest()->GetType() == NonlinearSolverTest::NONE;
     doublereal dErrPrev = std::numeric_limits<doublereal>::max(); // disable error test for the first iteration
     const doublereal dErr0 = dErr;

     OutputIteration(dErr, iIterCnt, false, oCPU);

     pS->CheckTimeStepLimit(dErr, dErrDiff);

     if (bResConverged && bSolConverged) {
          // Attention: Do not throw ConvergenceOnSolution here!
          // It would prevent that AfterConvergence is called in StepIntegrator::Advance.
          return;
     }

     const doublereal dMaxStep = dGetMaxNewtonStep(*pDM->GetpXCurr(), *pDM->GetpXPCurr());

     TRACE_VAR(dMaxStep);

     try {
          while (true) {
               oCPU.Jacobian.Tic(oCPU.Residual);

               Jacobian(pSM->pMatHdl());

               ASSERT(pSM->pMatHdl()->iGetNumCols() == Size);
               ASSERT(pSM->pMatHdl()->iGetNumCols() == pRes->iGetSize());
               ASSERT(pSM->pMatHdl()->iGetNumRows() == pRes->iGetSize());

               oCPU.Residual.Tic(oCPU.Jacobian);

               // compute gradient g = \nabla f = fjac^T \, fvec = -Jac^T \, pRes
               g.Reset();
               pSM->pMatHdl()->MatTVecDecMul(g, *pRes); // Attention: must be called before Solve() if row scaling is used

               const doublereal fPrev = fCurr;

               oCPU.LinearSolver.Tic(oCPU.Residual);

               pSM->Solve();

               ++iIterCnt;

               if (outputSol()) {
                    pS->PrintSolution(*pSol, iIterCnt);
               }

               oCPU.Residual.Tic(oCPU.LinearSolver);

               const bool bCheck = LineSearch(dMaxStep, fPrev, fCurr, iIterCnt);

               bResConverged = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);

               const doublereal dErrFactor = dErr / dErr0;

               OutputIteration(dErr, iIterCnt, true, oCPU);

               bSolConverged = MakeSolTest(pS, *pSol, SolTol, dSolErr);

               if (outputIters()) {
#ifdef USE_MPI
                    if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
                    {
                         silent_cout("\t\tSolErr " << dSolErr << '\n');
                    }
               }

               if (bResConverged && bSolConverged) {
                    // Attention: Do not throw ConvergenceOnSolution here!
                    // It would prevent that AfterConvergence is called in StepIntegrator::Advance.
                    return;
               }

               if (bCheck) { // lambda <= dLambdaMinCurr: check for gradient zero
                    bCheckZeroGradient(fCurr, dErr, Tol, iIterCnt);
               }

               bCheckDivergence(dErrFactor, dErrPrev, dErr, iIterCnt);

               if (iIterCnt >= std::abs(iMaxIter)) {
                    if (iMaxIter < 0 && dErrFactor < 1.) {
                         return;
                    }
                    if (outputBailout()) {
                         pS->PrintResidual(*pRes, iIterCnt);
                    }
                    throw NoConvergence(MBDYN_EXCEPT_ARGS);
               }

               dErrPrev = dErr;

               // allow to bail out in case of multiple CTRL^C
               if (mbdyn_stop_at_end_of_iteration()) {
                    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
               }
          }
     } catch (const LinearSolver::ErrNoPivot&) {
          throw ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
     }
}

LineSearchModified::LineSearchModified(DataManager* pDM,
                                       const NonlinearSolverTestOptions& options,
                                       const struct LineSearchParameters& param)
     :LineSearchSolver(pDM, options, param),
      dTimePrev(pDM->dGetTime()), dTimeStepPrev(0.)
{
}

LineSearchModified::~LineSearchModified(void)
{

}

void LineSearchModified::Solve(const NonlinearProblem* const pNLP,
                               Solver* const pS,
                               const integer iMaxIter,
                               const doublereal& Tol,
                               integer& iIterCnt,
                               doublereal& dErr,
                               const doublereal& SolTol,
                               doublereal& dSolErr)
{
     if (this->pNLP != pNLP) {
          // Force update of Jacobian matrix in order to avoid problems when scaling an empty matrix
          iRebuildJac = 0;
     }

     Attach(pNLP, pS);
     dSolErr = 0.;
     dErr = 0.;

     if (!bKeepJacAcrossSteps) {
          iRebuildJac = 0;
     }

     const doublereal dTimeCurr = pDM->dGetTime();
     const doublereal dTimeStepCurr = dTimeCurr - dTimePrev;
     const doublereal dTimeStepDiff = (dTimeStepCurr - dTimeStepPrev) / dTimeStepPrev;

     TRACE_VAR(dTimeCurr);
     TRACE_VAR(dTimePrev);
     TRACE_VAR(dTimeStepPrev);
     TRACE_VAR(dTimeStepCurr);
     TRACE_VAR(dTimeStepDiff);
     TRACE_VAR(dTimeStepTol);

     if (std::abs(dTimeStepDiff) > dTimeStepTol) {
          TRACE("Time step changed too much: force rebuild of Jacobian\n");
          iRebuildJac = 0;
     }

     CPUTime oCPU(*this);

     oCPU.Residual.Tic();

     iIterCnt = 0;
     bool bRebuildJac = false;
     bool bResConverged = pGetResTest()->GetType() == NonlinearSolverTest::NONE;
     bool bSolConverged = pGetSolTest()->GetType() == NonlinearSolverTest::NONE;
     bool bUpdateResidual = true;
     bool bDivergence = false;
     doublereal fCurr, dErrDiff, dErr0 = 0;
     const doublereal dMaxStep = dGetMaxNewtonStep(*pDM->GetpXCurr(), *pDM->GetpXPCurr());

     try {
          while (true) {
               if (bUpdateResidual) {
                    Residual(fCurr, iIterCnt);

                    TRACE_VAR(fCurr);

                    bUpdateResidual = false;

                    bResConverged = MakeResTest(pS, pNLP, *pRes, iIterCnt > 0 ? Tol : 1e-2 * Tol, dErr, dErrDiff);

                    OutputIteration(dErr, iIterCnt, bRebuildJac, oCPU);

                    pS->CheckTimeStepLimit(dErr, dErrDiff);

                    if (bResConverged && bSolConverged) {
                         goto exit_success;
                    }

                    if (iIterCnt == 0) {
                         dErr0 = dErr;
                    }
               }

               TRACE_VAR(fCurr);

               const doublereal fPrev = fCurr;
               const doublereal dErrPrev = dErr;

               oCPU.Jacobian.Tic(oCPU.Residual);

               doublereal dAlphaCurr;

               bRebuildJac = iRebuildJac <= 0 || bDivergence;

               if (bRebuildJac) {
                    Jacobian(pSM->pMatHdl());
                    iRebuildJac = bDivergence ? 0 : iIterationsBeforeAssembly;
                    dTimeStepPrev = dTimeStepCurr;
                    dAlphaCurr = dAlphaFull;
               } else {
                    TRACE("Reuse Jacobian, update of Jacobian after " << iRebuildJac << " iterations\n");
                    --iRebuildJac;
                    dAlphaCurr = dAlphaModified;
               }

               // Attention: must be called before Solve() if row scaling is used
               // Attention: Jacobian matrix must not be replaced by LU factors during a previous Solve()
               g.Reset();
               pSM->pMatHdl()->MatTVecDecMul(g, *pRes);

               TRACE("Solving system of linear equations ...\n");

               oCPU.LinearSolver.Tic(oCPU.Jacobian);

               pSM->Solve();

               TRACE("Linear solver completed\n");

               ++iIterCnt;

               if (outputSol()) {
                    pS->PrintSolution(*pSol, iIterCnt);
               }

               bSolConverged = MakeSolTest(pS, *pSol, SolTol, dSolErr);

               if (outputIters()) {
#ifdef USE_MPI
                    if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
                    {
                         silent_cout("\t\tSolErr " << dSolErr << '\n');
                    }
               }

               if (bSolConverged && bResConverged) {
                    pNLP->Update(pSol);

                    goto exit_success;
               }

               oCPU.Residual.Tic(oCPU.LinearSolver);

               doublereal dSlope = g.InnerProd(*pSol);

               TRACE_VAR(dSlope);

               if (dSlope >= 0) {
                    if (!bRebuildJac) {
                         TRACE("Slope is not negative: Jacobian or gradient might be out of date\n");
                         iRebuildJac = 0;
                         bUpdateResidual = true; // Attention: If we are using automatic row scaling, the residual will be destroyed during Solve()
                         continue;
                    } else {
                         TRACE("Jacobian has been updated but slope is still not negative\n");
                         bDivergence = true; // Not sure if strictly required but should be safe
                    }
               }

               ScaleNewtonStep(dMaxStep, *pSol, iIterCnt);

               p = *pSol;

               doublereal dLambda = 1., dLambdaPrev = 0.;
               const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, bRebuildJac, p, iIterCnt, fCurr);

               TRACE_VAR(dLambdaMinCurr);

               integer iLineSearchIter = 0;

               doublereal dLambdaMax = 1.;

               SetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX, dLambdaMax);

               while (true) {
                    TRACE_VAR(dLambda);
                    TRACE_VAR(dLambdaPrev);

                    if (iLineSearchIter > 0) {
                         for (integer i = 1; i <= Size; ++i) {
                              (*pSol)(i) = (dLambda - dLambdaPrev) * p(i);
                         }
                    }

                    SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, iLineSearchIter);
                    SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, dLambda);

                    pNLP->Update(pSol);

                    Residual(fCurr, iIterCnt);

                    if (iLineSearchIter == 0) {
                         dLambdaMax = std::max(dLambdaMinCurr, GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX));

                         ASSERT(dLambdaMax <= 1.);
                         ASSERT(dLambdaMax >= 0.);
                    }

                    ++iLineSearchIter;

                    TRACE_VAR(fCurr);
                    TRACE_VAR(fPrev);
                    TRACE_VAR(dAlphaCurr);
                    TRACE_VAR(dSlope);

                    bool bResTestFinite = false;

                    try {
                         bResConverged = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);
                         bResTestFinite = true;
                    } catch (const NonlinearSolver::ErrSimulationDiverged&) {
                         bResConverged = false;

                         if (uFlags & VERBOSE_MODE) {
                              silent_cerr("line search warning: residual test failed (residual=" << dErr << ")!"  << std::endl);
                         }
                    }

                    OutputLineSearch(iIterCnt, iLineSearchIter, fCurr / fPrev, dErr, dLambda, dSlope);

                    if (bResConverged && bSolConverged && dLambda <= dLambdaMax) {
                         break;
                    }

                    if (!bRebuildJac && sqrt(fCurr / fPrev) > dUpdateRatio) {
                         // Force update of Jacobian because convergence is too slow!
                         iRebuildJac = 0;
                    }

                    TRACE_VAR(fCurr);
                    TRACE_VAR(fPrev);
                    TRACE_VAR(dSlope);
                    TRACE_VAR(dAlphaCurr);
                    TRACE("dAlphaCurr * dLambda * dSlope + fPrev=" << dAlphaCurr * dLambda * dSlope + fPrev << "\n");

                    if (dLambda <= dLambdaMax) {
                         if (bResTestFinite && fCurr < dAlphaCurr * dLambda * dSlope + fPrev) {
                              TRACE("Sufficient decrease in f: backtrack\n");
                              break;
                         } else {
                              bDivergence = true;

                              if (!bRebuildJac) {
                                   if (uFlags & VERBOSE_MODE) {
                                        silent_cerr("line search warning: Divergent solution detected!\n");
                                   }

                                   iRebuildJac = 0;

                                   for (integer i = 1; i <= Size; ++i) {
                                        (*pSol)(i) = -dLambda * p(i);
                                   }

                                   pNLP->Update(pSol);
                                   bUpdateResidual = true;
                                   break;
                              }
                         }

                         if (dLambda <= dLambdaMinCurr) {
                              break;
                         }
                    }

                    if (iLineSearchIter >= iMaxIterations) {
                         throw MaxIterations(MBDYN_EXCEPT_ARGS);
                    }

                    dLambdaPrev = dLambda;
                    dLambda = dGetLambdaNext(dLambdaPrev, dSlope, fPrev, fCurr);

                    if (dLambdaMax < dLambda) {
                         pedantic_cout("lambda reduced from " << dLambda << " to " << dLambdaMax
                                       << " because of element request" << std::endl);
                         dLambda = dLambdaMax;
                    }
               }

               TRACE_VAR(dErr);
               TRACE_VAR(dErrDiff);
               TRACE_VAR(Tol);
               TRACE_VAR(bResConverged);

               OutputIteration(dErr, iIterCnt, bRebuildJac, oCPU);

               if (bResConverged && bSolConverged) {
                    goto exit_success;
               }

               if (dLambda < dLambdaMinCurr) {
                    bCheckZeroGradient(fCurr, dErr, Tol, iIterCnt);
               }

               if (bRebuildJac) {
                    bCheckDivergence(dErr / dErr0, dErrPrev, dErr, iIterCnt);
               }

               pS->CheckTimeStepLimit(dErr, dErrDiff);

               if (iIterCnt >= std::abs(iMaxIter)) {
                    if (outputBailout()) {
                         pS->PrintResidual(*pRes, iIterCnt);
                    }
                    throw NoConvergence(MBDYN_EXCEPT_ARGS);
               }

               // allow to bail out in case of multiple CTRL^C
               if (mbdyn_stop_at_end_of_iteration()) {
                    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
               }
          }
     } catch (const LinearSolver::ErrNoPivot&) {
          throw ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
     }

exit_success:
     dTimePrev = dTimeCurr;
     TRACE_VAR(dTimePrev);
     // Attention: Do not throw ConvergenceOnSolution here!
     // It would prevent that AfterConvergence is called in StepIntegrator::Advance.
}

LineSearchBFGS::LineSearchBFGS(DataManager* pDM,
                               const NonlinearSolverTestOptions& options,
                               const struct LineSearchParameters& param)
     :LineSearchSolver(pDM, options, param)

{
}

LineSearchBFGS::~LineSearchBFGS(void)
{
}

void LineSearchBFGS::Attach(const NonlinearProblem* pNLP, Solver* pS)
{
     LineSearchSolver::Attach(pNLP, pS);

     if (t.iGetSize() != Size) {
          t.Resize(Size);
          s.Resize(Size);
          w.Resize(Size);
          FCurr.Resize(Size);
          FPrev.Resize(Size);
     }
}

void LineSearchBFGS::Solve(const NonlinearProblem *pNLP,
                           Solver *pS,
                           const integer iMaxIter,
                           const doublereal& Tol,
                           integer& iIterCnt,
                           doublereal& dErr,
                           const doublereal& SolTol,
                           doublereal& dSolErr)
{
     if (this->pNLP != pNLP) {
          // Force update of Jacobian matrix in order to avoid problems when scaling an empty matrix
          iRebuildJac = 0;
     }

     Attach(pNLP, pS);
     dSolErr = 0.;
     dErr = 0.;

     if (!bKeepJacAcrossSteps) {
          iRebuildJac = 0;
     }

     auto* const pSM = dynamic_cast<QrSolutionManager*>(this->pSM);

     if (!pSM) {
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     CPUTime oCPU(*this);

     oCPU.Residual.Tic();

     iIterCnt = 0;
     bool bResConverged = pGetResTest()->GetType() == NonlinearSolverTest::NONE;
     bool bSolConverged = pGetSolTest()->GetType() == NonlinearSolverTest::NONE;
     bool bRebuildJac = false;
     bool bUpdateResidual = true;
     bool bDivergence = false;
     doublereal fCurr, dErrDiff, dErr0 = 0.;
     const doublereal dMaxStep = dGetMaxNewtonStep(*pDM->GetpXCurr(), *pDM->GetpXPCurr());

     try {
          while (true) {
               if (bUpdateResidual) {
                    Residual(fCurr, iIterCnt);

                    bUpdateResidual = false;

                    bResConverged = MakeResTest(pS, pNLP, *pRes, iIterCnt > 0 ? Tol : 1e-2 * Tol, dErr, dErrDiff);

                    OutputIteration(dErr, iIterCnt, bRebuildJac, oCPU);

                    pS->CheckTimeStepLimit(dErr, dErrDiff);

                    if (bResConverged && bSolConverged) {
                         FCurr = *pRes;
                         goto exit_success;
                    }

                    if (iIterCnt == 0) {
                         dErr0 = dErr;
                    }
               }

               FCurr = *pRes; // *pRes will be overwritten in pSM->Solve()

               const doublereal dErrPrev = dErr;

               oCPU.Jacobian.Tic(oCPU.Residual);

               bRebuildJac = iRebuildJac <= 0 || bDivergence;

               doublereal dAlphaCurr;

               if (bRebuildJac) {
                    Jacobian(pSM->pMatHdl());

                    oCPU.LinearSolver.Tic(oCPU.Jacobian);

                    pSM->InitQR();

                    oCPU.Jacobian.Tic(oCPU.LinearSolver);

                    iRebuildJac = bDivergence ? 0 : iIterationsBeforeAssembly;
                    dAlphaCurr = dAlphaFull;
               } else {
                    if (iIterCnt > 0) {
                         t.Reset();

                         pSM->MatVecOp(QrSolutionManager::OP_A_MINUS_R_B, t, s); // t = -R * s

                         for (integer i = 1; i <= Size; ++i) {
                              w(i) = FCurr(i) - FPrev(i);
                         }

                         pSM->MatVecOp(QrSolutionManager::OP_A_MINUS_Q_B, w, t); // w = dF - Q * t = dF + Q * R * s

                         bool bSkip = true;

                         const auto EPS = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9);

                         for (integer i = 1; i <= Size; ++i) {
                              if (fabs(w(i)) >= EPS * (fabs(FCurr(i)) + fabs(FPrev(i)))) {
                                   bSkip = false;
                              } else {
                                   w(i) = 0.;
                              }
                         }

                         if (!bSkip) {
                              s *= (-1. / s.Dot());
                              pSM->UpdateQR(w, s);
                         }
                    }

                    --iRebuildJac;
                    dAlphaCurr = dAlphaModified;
               }

               if (bResConverged && bSolConverged) {
                    goto exit_success;
               }

               p.Reset();

               pSM->MatVecOp(QrSolutionManager::OP_A_PLUS_QT_B, p, FCurr); // p = Q^T * Fcurr

               g.Reset();

               pSM->MatVecOp(QrSolutionManager::OP_A_MINUS_RT_B, g, p); // g = -(Q * R)^T * Fcurr = -R^T * p

               FPrev = FCurr;

               const doublereal fPrev = fCurr;

               oCPU.LinearSolver.Tic(oCPU.Jacobian);

               *pSol = p;

               pSM->SolveR(); // R * s = p

               s = *pSol;

               ++iIterCnt;

               if (outputSol()) {
                    pS->PrintSolution(*pSol, iIterCnt);
               }

               bSolConverged = MakeSolTest(pS, *pSol, SolTol, dSolErr);

               if (outputIters()) {
#ifdef USE_MPI
                    if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
                    {
                         silent_cout("\t\tSolErr " << dSolErr << '\n');
                    }
               }

               if (bResConverged && bSolConverged) {
                    // Use our current solution to update the Jacobian for the next step
                    // unless the Jacobian will be rebuild anyway
                    if (!bKeepJacAcrossSteps || iRebuildJac <= 0) {
                         pNLP->Update(pSol);
                         goto exit_success;
                    }
               }

               oCPU.Residual.Tic(oCPU.LinearSolver);

               doublereal dSlope = g.InnerProd(*pSol);

               if (dSlope >= 0 && fCurr) {
                    iRebuildJac = 0;

                    if (!bRebuildJac) {
                         if (uFlags & VERBOSE_MODE) {
                              silent_cerr("line search warning: Slope = "
                                          << dSlope
                                          << " is not negative: Jacobian or gradient might be out of date\n");
                         }

                         *pRes = FCurr;
                         continue;
                    } else {
                         if (uFlags & VERBOSE_MODE) {
                              silent_cerr("line search warning: Jacobian has been updated but slope = "
                                          << dSlope
                                          << " is still not negative\n");
                         }
                         bDivergence = true; // Not sure if strictly required but should be safe
                    }
               }

               ScaleNewtonStep(dMaxStep, *pSol, iIterCnt);

               p = *pSol;

               doublereal dLambda = 1., dLambdaPrev = 0.;
               const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, bRebuildJac, p, iIterCnt, fCurr);

               integer iLineSearchIter = 0;
               doublereal dLambdaMax = 1.;

               SetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX, dLambdaMax);

               while (true) {
                    if (iLineSearchIter > 0) {
                         pSol->ScalarMul(p, dLambda - dLambdaPrev);
                    }

                    SetNonlinearSolverHint(LINESEARCH_ITERATION_CURR, iLineSearchIter);
                    SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, dLambda);

                    pNLP->Update(pSol);

                    Residual(fCurr, iIterCnt);

                    if (iLineSearchIter == 0) {
                         dLambdaMax = std::max(dLambdaMinCurr, GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX));
                    }

                    ++iLineSearchIter;

                    bool bResTestFinite = false;

                    try {
                         bResConverged = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);
                         bResTestFinite = true;
                    } catch (const NonlinearSolver::ErrSimulationDiverged&) {
                         bResConverged = false;
                         bDivergence = true;
                         iRebuildJac = 0;

                         if (uFlags & VERBOSE_MODE) {
                              silent_cerr("line search warning: residual test failed (residual=" << dErr << ")!"  << std::endl);
                         }
                    }

                    OutputLineSearch(iIterCnt, iLineSearchIter, fCurr / fPrev, dErr, dLambda, dSlope);

                    if (dLambda <= dLambdaMax) {
                         if (bResConverged && bSolConverged) {
                              break;
                         }

                         if (bResTestFinite && fCurr < dAlphaCurr * dLambda * dSlope + fPrev) {
                              TRACE("Sufficient decrease in f: backtrack\n");
                              break;
                         }

                         if (dLambda <= dLambdaMinCurr) {
                              bDivergence = dErr / dErrPrev > dDivergenceCheck;

                              if (bDivergence && !bRebuildJac) {
                                   if (uFlags & VERBOSE_MODE) {
                                        silent_cerr("line search warning: Divergent solution detected!\n");
                                   }

                                   iRebuildJac = 0;

                                   pSol->ScalarMul(p, -dLambda);

                                   pNLP->Update(pSol);
                                   bUpdateResidual = true;
                              } else if (!bRebuildJac && sqrt(fCurr / fPrev) > dUpdateRatio) {
                                   if (uFlags & VERBOSE_MODE) {
                                        silent_cerr("line search warning: Jacobian update forced because rate of convergence is too slow!\n");
                                   }
                                   iRebuildJac = 0;
                              }

                              break;
                         }
                    }

                    if (iLineSearchIter >= iMaxIterations) {
                         throw MaxIterations(MBDYN_EXCEPT_ARGS);
                    }

                    dLambdaPrev = dLambda;
                    dLambda = dGetLambdaNext(dLambdaPrev, dSlope, fPrev, fCurr);

                    if (dLambdaMax < dLambda) {
                         pedantic_cout("lambda reduced from " << dLambda << " to " << dLambdaMax
                                       << " because of element request" << std::endl);
                         dLambda = dLambdaMax;
                    }
               }

               s.ScalarMul(p, dLambda);

               OutputIteration(dErr, iIterCnt, bRebuildJac, oCPU);

               if (bResConverged && bSolConverged) {
                    // Use our current solution to update the Jacobian
                    // unless it will be rebuild anyway
                    if (!bKeepJacAcrossSteps || iRebuildJac <= 0) {
                         goto exit_success;
                    }
               } else {
                    if (dLambda < dLambdaMinCurr) {
                         if (bCheckZeroGradient(fCurr, dErr, Tol, iIterCnt)) {
                              iRebuildJac = 0;
                         }
                    }

                    if (bCheckDivergence(dErr / dErr0, dErrPrev, dErr, iIterCnt)) {
                         iRebuildJac = 0;
                    }

                    pS->CheckTimeStepLimit(dErr, dErrDiff);

                    if (iIterCnt >= std::abs(iMaxIter)) {
                         if (outputBailout()) {
                              pS->PrintResidual(*pRes, iIterCnt);
                         }

                         throw NoConvergence(MBDYN_EXCEPT_ARGS);
                    }
               }

               // allow to bail out in case of multiple CTRL^C
               if (mbdyn_stop_at_end_of_iteration()) {
                    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
               }
          }
     } catch (const LinearSolver::ErrNoPivot&) {
          throw ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
     }

exit_success:
     // Attention: Do not throw ConvergenceOnSolution here!
     // It would prevent that AfterConvergence is called in StepIntegrator::Advance.
     ;
}

LineSearchMCP::LineSearchMCP(DataManager* pDM,
                             const NonlinearSolverTestOptions& options,
                             const struct LineSearchParameters& param)
     :LineSearchSolver(pDM, options, param)
{
}

LineSearchMCP::~LineSearchMCP()
{
}

void LineSearchMCP::ComputeH(const VectorHandler& z, const VectorHandler& F, const SpGradientSparseMatrixHandler& nablaFMCP, MatrixHandler& H)
{
     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          workV1(i) = (eEqType == DofOrder::INEQUALITY)
               ? (fabs(z(i)) < dMcpTol) && (fabs(F(i)) < dMcpTol)
               : 0.;
     }

     workV2.Reset();

     nablaFMCP.MatTVecDecMul(workV2, workV1);

     H.Reset();

     sp_grad::SpGradient g;

     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          if (eEqType == DofOrder::EQUALITY) {
               H.SubItem(i, nablaFMCP.GetRow(i));
          } else if (workV1(i)) {
               const doublereal normi = sqrt(workV1(i) * workV1(i) + workV2(i) * workV2(i));
               g.Reset(0., i, workV1(i) / normi - 1.);
               H.AddItem(i, g - nablaFMCP.GetRow(i) * (workV2(i) / normi - 1.));
          } else {
               const doublereal normi = sqrt(z(i) * z(i) + F(i) * F(i));
               g.Reset(0., i, z(i) / normi - 1.);
               H.AddItem(i, g - nablaFMCP.GetRow(i) * (F(i) / normi - 1.));
          }
     }
}

void LineSearchMCP::ComputeH(const VectorHandler& z, const VectorHandler& F, MatrixHandler& H, std::vector<doublereal>& rgRowScale, std::vector<doublereal>& rgColScale)
{
     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          workV1(i) = (eEqType == DofOrder::INEQUALITY)
               ? (fabs(z(i)) < dMcpTol) && (fabs(F(i)) < dMcpTol)
               : 0.;
     }

     workV2.Reset();

     H.MatTVecDecMul(workV2, workV1);

     ASSERT(rgRowScale.size() == static_cast<size_t>(Size));
     ASSERT(rgColScale.size() == 0);

     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          if (eEqType == DofOrder::EQUALITY) {
               rgRowScale[i - 1] = -1.;
          } else if (workV1(i)) {
               const doublereal normi = sqrt(workV1(i) * workV1(i) + workV2(i) * workV2(i));
               rgRowScale[i - 1] = -(workV2(i) / normi - 1.);
          } else {
               const doublereal normi = sqrt(z(i) * z(i) + F(i) * F(i));
               rgRowScale[i - 1] = -(F(i) / normi - 1.);
          }
     }

     H.PacMat();

     H.Scale(rgRowScale, rgColScale);

     sp_grad::SpGradient g;

     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          if (eEqType == DofOrder::EQUALITY) {
               continue;
          } else if (workV1(i)) {
               const doublereal normi = sqrt(workV1(i) * workV1(i) + workV2(i) * workV2(i));
               g.Reset(0., i, workV1(i) / normi - 1.);
               H.AddItem(i, g);
          } else {
               const doublereal normi = sqrt(z(i) * z(i) + F(i) * F(i));
               g.Reset(0., i, z(i) / normi - 1.);
               H.AddItem(i, g);
          }
     }     
}

void LineSearchMCP::ComputeHDesc(const SpGradientSparseMatrixHandler& nablaFMCP, const VectorHandler& z, const VectorHandler& F, MatrixHandler& H) const
{
     H.Reset();

     sp_grad::SpGradient g;

     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          if (eEqType == DofOrder::INEQUALITY && z(i) <= F(i)) {
               g.Reset(0., i, 1.);
               H.AddItem(i, g);
          } else {
               H.SubItem(i, nablaFMCP.GetRow(i));
          }
     }

     H.PacMat();
}

void LineSearchMCP::ComputeRHSDesc(const VectorHandler& z, const VectorHandler& F, VectorHandler& Fmin) const
{
     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          Fmin(i) = (eEqType == DofOrder::EQUALITY)
               ? -F(i)
               : (z(i) <= F(i)) ? -z(i) : -F(i);
     }
}

void LineSearchMCP::ComputeFMCP(VectorHandler& F)
{
     F.Reset();

     if (pAbsRes) {
          pAbsRes->Reset();
     }

     try {
          TRACE("Assemble residual\n");
          pNLP->Residual(&F, pAbsRes);
     } catch (const SolutionDataManager::ChangedEquationStructure&) {
          if (bHonorJacRequest) {
               iRebuildJac = 0;
          }
     }
}

doublereal LineSearchMCP::ComputeFMerit(const VectorHandler& z, const VectorHandler& F, VectorHandler& FMerit) const
{
     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          FMerit(i) = (eEqType == DofOrder::EQUALITY)
               ? F(i)
               : sqrt(z(i) * z(i) + F(i) * F(i)) - (z(i) + F(i));
     }

     return 0.5 * FMerit.Dot();
}

void LineSearchMCP::Update(const VectorHandler& DeltaZ, VectorHandler& z) const
{
     pNLP->Update(&DeltaZ);

     const VectorHandler& XCurr = *pDM->GetpXCurr();

     for (integer i = 1; i <= Size; ++i) {
          DofOrder::Equality eEqType = pDM->GetEqualityType(i);

          z(i) = (eEqType == DofOrder::INEQUALITY)
               ? XCurr(i)
               : 0.;
     }
}

doublereal LineSearchMCP::ApplyInc(VectorHandler& z, VectorHandler& F, VectorHandler& FMerit, VectorHandler& DeltaZ, const VectorHandler& d, const doublereal dInc)
{
     DeltaZ.ScalarMul(d, dInc);
     Update(DeltaZ, z);
     ComputeFMCP(F);
     return ComputeFMerit(z, F, FMerit);
}

doublereal LineSearchMCP::LineSearch(const doublereal dThetaPrev,
                                     const doublereal preRHS,
                                     VectorHandler& z,
                                     VectorHandler& F,
                                     VectorHandler& FMerit,
                                     VectorHandler& DeltaZ,
                                     const integer iIterCnt)
{
     doublereal dLambda = 2., dLambdaPrev = 0.;

     p = DeltaZ;

     integer iLineSearchIter = 0;
     doublereal dThetaCurr = dThetaPrev;
     doublereal dErr, dErrDiff;
     const doublereal dLambdaMinCurr = std::max(dLambdaMin, dGetMinNewtonInc(p));

     while (dLambda >= dLambdaMinCurr) {
          dThetaCurr = ApplyInc(z, F, FMerit, DeltaZ, p, dLambda - dLambdaPrev);

          try {
               MakeResTest(pS, pNLP, FMerit, 0., dErr, dErrDiff);
          } catch (const NonlinearSolver::ErrSimulationDiverged&) {
               if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: residual test failed (residual=" << dErr << ")\n");
               }
          }
          
          OutputLineSearch(iIterCnt, iLineSearchIter, dThetaCurr / dThetaPrev, dErr, dLambda, preRHS);

          if (dThetaCurr < dThetaPrev + dLambda * dAlphaFull * preRHS) {
               return dThetaCurr;
          }

          dLambdaPrev = dLambda;
          dLambda *= 0.5; // FIXME: consider dLambdaFactMin?

          if (++iLineSearchIter > iMaxIterations) {
               break;
          }
     }

     return dThetaCurr;
}

void LineSearchMCP::Attach(const NonlinearProblem* pNLPNew, Solver* pS)
{
     LineSearchSolver::Attach(pNLPNew, pS);

     if (zH.iGetSize() != Size) {
          zH.ResizeReset(Size);
          FH.ResizeReset(Size);
          FMeritH.ResizeReset(Size);
          workV1.ResizeReset(Size);
          workV2.ResizeReset(Size);
          JacThetaFMeritH.ResizeReset(Size);
     }
}


MCPNewtonMinFB::MCPNewtonMinFB(DataManager* pDM,
                               const NonlinearSolverTestOptions& options,
                               const struct LineSearchParameters& param)
     :LineSearchMCP(pDM, options, param),
      nablaFMCPH(0, 0)
{
}

MCPNewtonMinFB::~MCPNewtonMinFB()
{
}

void MCPNewtonMinFB::Solve(const NonlinearProblem *pNLP,
                           Solver *pS,
                           const integer iMaxIter,
                           const doublereal& Tol,
                           integer& iIterCnt,
                           doublereal& dErr,
                           const doublereal& SolTol,
                           doublereal& dSolErr)
{
     iIterCnt = 0;
     doublereal dErrDiff = dSolErr = dErr = 0.;

     Attach(pNLP, pS);

     MatrixHandler& H = *pSM->pMatHdl();
     VectorHandler& Fmin = *pSM->pResHdl();
     VectorHandler& DeltaZ = *pSM->pSolHdl();

     CPUTime oCPU(*this);

     oCPU.Residual.Tic();

     ComputeFMCP(FH);

     DeltaZ.Reset();
     Update(DeltaZ, zH);

     doublereal dThetaPrev = ComputeFMerit(zH, FH, FMeritH);

     if (outputRes()) {
          pS->PrintResidual(FMeritH, iIterCnt);
     }

     const bool bResConverged = MakeResTest(pS, pNLP, FMeritH, Tol, dErr, dErrDiff);     

     OutputIteration(dErr, iIterCnt, false, oCPU);
     
     while (true) {
          ComputeRHSDesc(zH, FH, Fmin);
          
          oCPU.Jacobian.Tic(oCPU.Residual);
          
          Jacobian(&nablaFMCPH);
          
          ComputeHDesc(nablaFMCPH, zH, FH, H);

          oCPU.LinearSolver.Tic(oCPU.Jacobian);

          pSM->Solve(); // solve H * DeltaZ = Fmin for DeltaZ

          if (outputSol()) {
               pS->PrintSolution(DeltaZ, iIterCnt);
          }

          oCPU.Residual.Tic(oCPU.LinearSolver);

          const bool bSolConverged = MakeSolTest(pS, DeltaZ, SolTol, dSolErr);

          doublereal dThetaCurr = ComputeFMerit(zH, FH, FMeritH);
          ComputeH(zH, FH, nablaFMCPH, H);
          H.MatTVecMul(JacThetaFMeritH, FMeritH);
          Update(DeltaZ, zH);

          ComputeFMCP(FH);

          dThetaCurr = ComputeFMerit(zH, FH, FMeritH);

          if (dThetaCurr > dMcpSigma * dThetaPrev) {
               doublereal preRHS = JacThetaFMeritH.InnerProd(DeltaZ);
               const doublereal threshold = -dMcpRho * std::pow(DeltaZ.Norm(), dMcpP);

               p.ScalarMul(DeltaZ, -1.);

               Update(p, zH);

               if (preRHS > threshold) {
                    DeltaZ.ScalarMul(JacThetaFMeritH, -1.);
                    preRHS = -JacThetaFMeritH.Dot();
               }

               dThetaCurr = LineSearch(dThetaPrev, preRHS, zH, FH, FMeritH, DeltaZ, iIterCnt);
          }

          const bool bResConverged = MakeResTest(pS, pNLP, FMeritH, Tol, dErr, dErrDiff);

          ++iIterCnt;

          if (outputRes()) {
               pS->PrintResidual(FMeritH, iIterCnt);
          }

          OutputIteration(dErr, iIterCnt, true, oCPU);

          if (outputIters()) {
#ifdef USE_MPI
               if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
               {
                    silent_cout("\t\tSolErr " << dSolErr << '\n');
               }
          }

          if (bResConverged && bSolConverged) {
               break;
          }

          if (iIterCnt >= std::abs(iMaxIter)) {
               if (outputBailout()) {
                    pS->PrintResidual(*pRes, iIterCnt);
               }

               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }

          // allow to bail out in case of multiple CTRL^C
          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }

          dThetaPrev = dThetaCurr;
     }
}

void MCPNewtonMinFB::Attach(const NonlinearProblem* pNLPNew, Solver* pS)
{
     LineSearchMCP::Attach(pNLPNew, pS);

     if (nablaFMCPH.iGetNumRows() != Size) {
          nablaFMCPH.ResizeReset(Size, Size);
     }
}

MCPNewtonFB::MCPNewtonFB(DataManager* pDM,
                         const NonlinearSolverTestOptions& options,
                         const struct LineSearchParameters& param)
     :LineSearchMCP(pDM, options, param)
{
}

MCPNewtonFB::~MCPNewtonFB()
{
}

void MCPNewtonFB::Solve(const NonlinearProblem *pNLP,
                        Solver *pS,
                        const integer iMaxIter,
                        const doublereal& Tol,
                        integer& iIterCnt,
                        doublereal& dErr,
                        const doublereal& SolTol,
                        doublereal& dSolErr)
{
     iIterCnt = 0;
     doublereal dErrDiff = dSolErr = dErr = 0.;

     Attach(pNLP, pS);

     MatrixHandler& H = *pSM->pMatHdl();
     VectorHandler& Fmin = *pSM->pResHdl();
     VectorHandler& DeltaZ = *pSM->pSolHdl();

     CPUTime oCPU(*this);

     oCPU.Residual.Tic();

     ComputeFMCP(FH);

     DeltaZ.Reset();
     Update(DeltaZ, zH);

     doublereal dThetaPrev = ComputeFMerit(zH, FH, FMeritH);

     if (outputRes()) {
          pS->PrintResidual(FMeritH, iIterCnt);
     }

     const bool bResConverged = MakeResTest(pS, pNLP, FMeritH, Tol, dErr, dErrDiff);     

     OutputIteration(dErr, iIterCnt, false, oCPU);
     
     while (true) {
          oCPU.Jacobian.Tic(oCPU.Residual);
          
          Jacobian(&H);
          
          ComputeH(zH, FH, H, rgRowScale, rgColScale);
          
          H.MatTVecMul(JacThetaFMeritH, FMeritH);

          Fmin.ScalarMul(FMeritH, -1.);

          oCPU.LinearSolver.Tic(oCPU.Jacobian);

          pSM->Solve(); // solve H * DeltaZ = Fmin for DeltaZ

          if (outputSol()) {
               pS->PrintSolution(DeltaZ, iIterCnt);
          }

          const bool bSolConverged = MakeSolTest(pS, DeltaZ, SolTol, dSolErr);

          oCPU.Residual.Tic(oCPU.LinearSolver);

          Update(DeltaZ, zH);

          ComputeFMCP(FH);
          
          doublereal dThetaCurr = ComputeFMerit(zH, FH, FMeritH);

          if (dThetaCurr > dMcpSigma * dThetaPrev) {
               doublereal preRHS = JacThetaFMeritH.InnerProd(DeltaZ);
               const doublereal threshold = -dMcpRho * std::pow(DeltaZ.Norm(), dMcpP);

               p.ScalarMul(DeltaZ, -1.);

               Update(p, zH);

               if (preRHS > threshold) {
                    DeltaZ.ScalarMul(JacThetaFMeritH, -1.);
                    preRHS = -JacThetaFMeritH.Dot();
               }

               dThetaCurr = LineSearch(dThetaPrev, preRHS, zH, FH, FMeritH, DeltaZ, iIterCnt);
          }

          const bool bResConverged = MakeResTest(pS, pNLP, FMeritH, Tol, dErr, dErrDiff);

          ++iIterCnt;

          if (outputRes()) {
               pS->PrintResidual(FMeritH, iIterCnt);
          }

          OutputIteration(dErr, iIterCnt, true, oCPU);

          if (outputIters()) {
#ifdef USE_MPI
               if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
               {
                    silent_cout("\t\tSolErr " << dSolErr << '\n');
               }
          }

          if (bResConverged && bSolConverged) {
               break;
          }

          if (iIterCnt >= std::abs(iMaxIter)) {
               if (outputBailout()) {
                    pS->PrintResidual(*pRes, iIterCnt);
               }

               throw NoConvergence(MBDYN_EXCEPT_ARGS);
          }

          // allow to bail out in case of multiple CTRL^C
          if (mbdyn_stop_at_end_of_iteration()) {
               throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
          }

          dThetaPrev = dThetaCurr;
     }
}

void MCPNewtonFB::Attach(const NonlinearProblem* pNLP, Solver* pS)
{
     LineSearchMCP::Attach(pNLP, pS);

     rgRowScale.resize(Size);
}
