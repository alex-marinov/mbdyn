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
 
 /* 
  * Portions Copyright (C) 2003-2017
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

 /*
 AUTHOR: Reinhard Resch <r.resch@a1.net>
  Copyright (C) 2011(-2019) all rights reserved.

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
#include "clock_time.h"
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
: dTolX(1e-7),
dTolMin(1e-8),
iMaxIterations(200),
dMaxStep(100.),
      dAlphaFull(1e-4),
      dAlphaModified(0.6),
dLambdaMin(1e-2),
dLambdaFactMin(1e-1),
dDivergenceCheck(1.),
dMinStepScale(1e-3),
uFlags(DIVERGENCE_CHECK
	| ALGORITHM_CUBIC
	| RELATIVE_LAMBDA_MIN
	| SCALE_NEWTON_STEP
             | ABORT_AT_LAMBDA_MIN),
      iIterationsBeforeAssembly(0),
      bKeepJac(false),
      dTimeStepTol(0.1),
      dUpdateRatio(0.05)
{
	NO_OP;
}

LineSearchSolver::LineSearchSolver(DataManager* pDM,
	const NonlinearSolverOptions& options,
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
        const doublereal dTemp = std::abs(dX(i)) / std::max(std::max(std::abs(X(i)),
                                                                     std::abs(XP(i))), 1.);
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

doublereal LineSearchSolver::dGetLambdaMin(doublereal& dSlope, const bool bRebuildJac, const VectorHandler& p, const integer iIterCnt) const
{
    TRACE_VAR(dSlope);

    doublereal dLambdaMinCurr = -1.;

    if (dSlope >= 0.) {
        if (uFlags & VERBOSE_MODE) {
            silent_cerr("line search warning: slope = " << dSlope << " >= 0"
                        " at time step " << pDM->dGetTime()
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

	try {
        TRACE("Assemble residual\n");
    	pNLP->Residual(pRes);
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
LineSearchSolver::Jacobian()
{
	SolutionManager *const pSM = pS->pGetSolutionManager();

	const integer iMaxIterRebuild = 10;

	pSM->MatrReset();
	integer iIter = 0;

	do {
		try {
            TRACE("Assemble Jacobian\n");
			pNLP->Jacobian(pSM->pMatHdl());
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
                pSM->pMatHdl()->Print(std::cout, MatrixHandler::MAT_PRINT_TRIPLET);
            }
		}

		if (outputMatrixConditionNumber()) {
            silent_cout(" cond=" << pSM->pMatHdl()->ConditionNumber(GetCondMatNorm()) << '\n');
		}
	}
}

LineSearchFull::LineSearchFull(DataManager* pDM,
                               const NonlinearSolverOptions& options,
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

    const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, true, p, iIterCnt);

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
            dLambdaMax = GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX);

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

        OutputLineSearch(iIterCnt, iLineSearchIter, fCurr, dErr, dLambda, dSlope);

        if (bResTestFinite) {
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

    bool bTest = MakeResTest(pS, pNLP, *pRes, 1e-2 * Tol, dErr, dErrDiff); // use a more stringent test for the first iteration
    doublereal dErrPrev = std::numeric_limits<doublereal>::max(); // disable error test for the first iteration
    const doublereal dErr0 = dErr;

    OutputIteration(dErr, iIterCnt, false, oCPU);

	pS->CheckTimeStepLimit(dErr, dErrDiff);

	if (bTest) {
		return;
	}

    const doublereal dMaxStep = dGetMaxNewtonStep(*pDM->GetpXCurr(), *pDM->GetpXPCurr());

    TRACE_VAR(dMaxStep);
        
    try {   
	while (true) {
            oCPU.Jacobian.Tic(oCPU.Residual);

		Jacobian();

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

		bTest = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);

            const doublereal dErrFactor = dErr / dErr0;

            OutputIteration(dErr, iIterCnt, true, oCPU);
        
            if (bTest) {
                return;
		}

            bTest = MakeSolTest(pS, *pSol, SolTol, dSolErr);

            if (outputIters()) {
#ifdef USE_MPI
		if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
		{
                    silent_cout("\t\tSolErr "
                                << dSolErr << '\n');
                }
				}

            if (bTest) {
                if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: Convergence on solution "
                                "at time step " << pDM->dGetTime()
                                << "at iteration " << iIterCnt
                                << "\tErr(n)=" << dErr);

                    if (iIterCnt >= 1) {
                        silent_cerr("\tErr(n-1)=" << dErrPrev
                                    << "\tErr(n)/Err(n-1)=" << dErr/dErrPrev
                                    << "\tErr(n)/Err(1) " << dErrFactor);
						}

                    silent_cerr('\n');
					}
                throw ConvergenceOnSolution(MBDYN_EXCEPT_ARGS);
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
                                       const NonlinearSolverOptions& options,
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
    
    if (!bKeepJac) {
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
    bool bSolConverged = false;
    bool bUpdateResidual = true;
    bool bResConverged = false;
    bool bDivergence = false;
    doublereal fCurr, dErrDiff, dErr0 = 0;
    const doublereal dMaxStep = dGetMaxNewtonStep(*pDM->GetpXCurr(), *pDM->GetpXPCurr());

    try {
        while (true) {
            if (bUpdateResidual) {
                TRACE_VAR(fCurr);
            
                Residual(fCurr, iIterCnt);
            
                TRACE_VAR(fCurr);
            
                bUpdateResidual = false;
            
                bResConverged = MakeResTest(pS, pNLP, *pRes, iIterCnt > 0 ? Tol : 1e-2 * Tol, dErr, dErrDiff);
    
                OutputIteration(dErr, iIterCnt, bRebuildJac, oCPU);

                pS->CheckTimeStepLimit(dErr, dErrDiff);
            
                if (bResConverged) {
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
                Jacobian();
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
#ifdef DEBUG
            std::cout << "\n\nJac=[\n";
            for (integer i = 1; i <= pSM->pMatHdl()->iGetNumRows(); ++i) {
                for (integer j = 1; j <= pSM->pMatHdl()->iGetNumCols(); ++j) {
                    doublereal d = pSM->pMatHdl()->dGetCoef(i, j);
                    ASSERT(std::isfinite(d));
                    if (d) {
                        std::cout << i << ", " << j << ", " << d << ";" << std::endl;
                    }
                }
            }
            std::cout << "];\n\n";
#endif
            oCPU.LinearSolver.Tic(oCPU.Jacobian);
        
            pSM->Solve();

            TRACE("Linear solver completed\n");
        
            ++iIterCnt;
        
            if (outputSol()) {
                pS->PrintSolution(*pSol, iIterCnt);
            }
        
            bSolConverged = MakeSolTest(pS, *pSol, SolTol, dSolErr);
        
            if (bSolConverged) {
			if (uFlags & VERBOSE_MODE) {
				silent_cerr("line search warning: Convergence on solution "
                    "at time step " << pDM->dGetTime()
					<< "at iteration " << iIterCnt
                                << "\tErr(n)=" << dSolErr);
                }

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
            const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, bRebuildJac, p, iIterCnt);
        
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
                    dLambdaMax = GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX);
            
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

                OutputLineSearch(iIterCnt, iLineSearchIter, fCurr, dErr, dLambda, dSlope);
            
                if (bResConverged) {
                    break;
				}

                if (!bRebuildJac && sqrt(fCurr / fPrev) > dUpdateRatio) {
                    // Force update of Jacobian because convergence is too slow!
                    iRebuildJac = 0;
                }
            
                if (bResTestFinite && fCurr < dAlphaCurr * dSlope + fPrev) {
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

            if (bResConverged) {
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

    if (bSolConverged) {
        TRACE("convergence on solution\n");
        throw ConvergenceOnSolution(MBDYN_EXCEPT_ARGS);
    }
}

LineSearchBFGS::LineSearchBFGS(DataManager* pDM,
                               const NonlinearSolverOptions& options,
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

    if (!bKeepJac) {
        iRebuildJac = 0;
    }

    auto* const pSM = dynamic_cast<QrSolutionManager*>(this->pSM);

    if (!pSM) {
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    CPUTime oCPU(*this);

    oCPU.Residual.Tic();

    iIterCnt = 0;
    bool bRebuildJac = false;
    bool bSolConverged = false;
    bool bUpdateResidual = true;
    bool bResConverged = false;
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

                if (bResConverged) {
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
                Jacobian();

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

                    constexpr auto EPS = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9);

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

            if (bResConverged || bSolConverged) {
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

            if (bSolConverged) {
                if (uFlags & VERBOSE_MODE) {
                    silent_cerr("line search warning: Convergence on solution "
                                "at time step " << pDM->dGetTime()
                                << "at iteration " << iIterCnt
                                << "\tErr(n)=" << dSolErr);
                }

                // Use our current solution to update the Jacobian for the next step
                // unless the Jacobian will be rebuild anyway
                if (!bKeepJac || iRebuildJac <= 0) {                        
                        pNLP->Update(pSol);
                        goto exit_success;
                }
            }

            oCPU.Residual.Tic(oCPU.LinearSolver);

            doublereal dSlope = g.InnerProd(*pSol);

            if (dSlope >= 0) {
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
            const doublereal dLambdaMinCurr = dGetLambdaMin(dSlope, bRebuildJac, p, iIterCnt);

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
                    dLambdaMax = GetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX);
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

                OutputLineSearch(iIterCnt, iLineSearchIter, fCurr, dErr, dLambda, dSlope);

                if (bResConverged) {
                    break;
                }

                if (bResTestFinite && fCurr < dAlphaCurr * dSlope + fPrev) {
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

            if (bResConverged || bSolConverged) {
                    // Use our current solution to update the Jacobian
                    // unless it will be rebuild anyway
                    if (!bKeepJac || iRebuildJac <= 0) {
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
    if (bSolConverged) {
        TRACE("convergence on solution\n");
        throw ConvergenceOnSolution(MBDYN_EXCEPT_ARGS);
    }
}
