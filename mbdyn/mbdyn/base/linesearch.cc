/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
  * Portions Copyright (C) 2003-2014
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

 /*
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2014) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
  */
  
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <unistd.h>

#include "solver.h"
#include "linesearch.h"
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

//#define DEBUG

#ifdef DEBUG
#undef ASSERT
#define ASSERT(expr) assert(expr)
#define TRACE(expr) silent_cerr(__FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ":" << expr)
#else
#define TRACE(expr) static_cast<void>(0)
#endif

#define TRACE_VAR(varname) TRACE(#varname "=" << (varname) << std::endl)
#define TRACE_FLAG(varname, flag) TRACE(#varname "&" #flag "=" << ((varname) & (flag)) << std::endl)

LineSearchParameters::LineSearchParameters(void)
: dTolX(1e-7),
dTolMin(1e-8),
iMaxIterations(200),
dMaxStep(100.),
dAlpha(1e-4),
dLambdaMin(1e-2),
dLambdaFactMin(1e-1),
dDivergenceCheck(1.),
dMinStepScale(1e-3),
uFlags(DIVERGENCE_CHECK
	| ALGORITHM_CUBIC
	| RELATIVE_LAMBDA_MIN
	| SCALE_NEWTON_STEP
	| ABORT_AT_LAMBDA_MIN)
{
	NO_OP;
}

LineSearchSolver::LineSearchSolver(DataManager* pDM,
	const NonlinearSolverOptions& options,
	const struct LineSearchParameters& param)
: NonlinearSolver(options),
LineSearchParameters(param),
pRes(0),
pSol(0),
pNLP(0),
pS(0),
pDM(pDM)
{
	TRACE_VAR(dTolX);
	// TRACE_VAR(dTolF);
	TRACE_VAR(dTolMin);
	TRACE_VAR(iMaxIterations);
	TRACE_VAR(dMaxStep);
	TRACE_VAR(dAlpha);
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

void
LineSearchSolver::Residual(doublereal& f, integer iIterCnt)
{
#ifdef USE_EXTERNAL
	SendExternal();
#endif /* USE_EXTERNAL */

	pRes->Reset();

	try {
    	pNLP->Residual(pRes);
	}
	catch (const SolutionDataManager::ChangedEquationStructure&) {
		// The Jacobian will be rebuilt anyway
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
			pNLP->Jacobian(pSM->pMatHdl());
		} catch (const MatrixHandler::ErrRebuildMatrix&) {
			silent_cout("NewtonRaphsonSolver: "
					"rebuilding matrix..."
					<< std::endl);

			/* need to rebuild the matrix... */
			pSM->MatrInitialize();
			continue;
		}
		break;
	} while (++iIter < iMaxIterRebuild);

	if (iIter >= iMaxIterRebuild) {
		silent_cerr("Maximum number of iterations exceeded when rebuilding the Jacobian matrix" << std::endl);
		throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
	}

	TotJac++;

#ifdef USE_MPI
	if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		if (outputJac()) {
			silent_cout("Jacobian:" << std::endl
					<< *(pSM->pMatHdl()));
		}

		if (outputMatrixConditionNumber()) {
			silent_cout(" cond=" << pSM->pMatHdl()->ConditionNumber(GetCondMatNorm()) << std::endl);
		}
	}
}

void
LineSearchSolver::LineSearch(doublereal stpmax, doublereal fold,
	doublereal& f, bool& check, const integer& iIterCnt)
{
	check = false;

	if (uFlags & SCALE_NEWTON_STEP) {
		const doublereal dNormSol = pSol->Norm();

		if (dNormSol > stpmax) {
			const doublereal dScale = std::max(dMinStepScale, stpmax / dNormSol);
			if (uFlags & VERBOSE_MODE) {
				silent_cerr("line search warning: "
	    				"Newton increment is reduced by factor " << dScale
                        << " at time step " << pDM->dGetTime()
	    				<< " at iteration " << iIterCnt
	    				<< " The time step is probably too large" << std::endl);
			}
			ASSERT(stpmax >= 0);
			ASSERT(dScale <= 1.);
			ASSERT(dScale > 0.);
			*pSol *= dScale;
		}
	}

	p = *pSol; // save the Newton increment

	doublereal slope = g.InnerProd(p);

	TRACE_VAR(slope);

	doublereal dLambdaMinEff = -1.;

	if (slope >= 0.) {
		if (uFlags & VERBOSE_MODE) {
			silent_cerr("line search warning: slope=" << slope << " >= 0"
                        " at time step " << pDM->dGetTime()
						<< " at iteration " << iIterCnt << std::endl
						<< "This could be a roundoff problem" << std::endl);
		}

		if (uFlags & NON_NEGATIVE_SLOPE_CONTINUE) {
			// It seems to be a numerical problem.
			// Line search may not work in this situation.
			// Resort to the ordinary Newton Raphson algorithm
			slope = 0.;
			dLambdaMinEff = 1.;
		} else {
			throw SlopeNotNegative(MBDYN_EXCEPT_ARGS);
		}
	}

	if (dLambdaMinEff < 0) {
		// dLambdaMinEff has to be detected
		if (uFlags & RELATIVE_LAMBDA_MIN) {
			doublereal test = 0.;

			for (integer i = 1; i <= Size; ++i) {
				const doublereal temp = std::abs(p(i)) / std::max(std::max(std::abs((*pDM->GetpXCurr())(i)), std::abs((*pDM->GetpXPCurr())(i))), 1.);
				if (temp > test) {
					test = temp;
				}
			}

			dLambdaMinEff = std::max(dLambdaMin, std::min(dTolX / test, 1.));
		} else {
			dLambdaMinEff = dLambdaMin;
		}
	}

	doublereal lambda = 1.;
	doublereal tmplam;
	doublereal lambda2, f2;

#ifdef DEBUG
	lambda2 = f2 = NAN;
#endif

	integer iIter = 0;

	TRACE_VAR(dLambdaMinEff);

	do {
		// FIXME: dLambdaMax not defined
		// ASSERT(lambda <= dLambdaMax);

		if (iIter > 0) {
			TRACE("Start new step from Xold, XPold with lambda = " << lambda << " ..." << std::endl);
			pNLP->Update(&dXneg); // restore the previous state

			pSol->Reset();
			pSol->ScalarMul(p, lambda); // scale the Newton increment by lambda
		}

		for (integer i = 1; i <= Size; ++i)
			dXneg(i) = -(*pSol)(i); // save the increment; so we can restore the previous state

		TRACE("Update the nonlinear problem ... pSol->Norm()=" << pSol->Norm() << std::endl);

		pNLP->Update(pSol);
		Residual(f, iIter);

		++iIter;

		TRACE("New value for f:" << f << std::endl);

		doublereal dErr = 0., dErrDiff = 0.;
		MakeResTest(pS, pNLP, *pRes, 0., dErr, dErrDiff);

		if (outputIters() && (uFlags & PRINT_CONVERGENCE_INFO)) {
#ifdef USE_MPI
			if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
			{
				silent_cout("\t\tf(" << iIterCnt << ":" << iIter << ")=" << std::setw(12) << f
					<< "\tErr=" << std::setw(12) << dErr
					<< "\tlambda=" << std::setw(12) << lambda
					<<"\tslope=" << slope << std::endl);
			}
		}

		pS->CheckTimeStepLimit(dErr, dErrDiff);

		if (f <= fold + dAlpha * lambda * slope) {
			TRACE("Sufficient decrease in f: backtrack" << std::endl);
			return;
		} else if (lambda <= dLambdaMinEff) {
			TRACE("Checking for convergence: lambda=" << lambda << " < lambdaMin=" << dLambdaMinEff << std::endl);
			check = true; // check for convergence
			return;
		} else {
			if (uFlags & ALGORITHM_CUBIC) {
				TRACE("Calculate new value for alam ..." << std::endl);
				if (lambda == 1.) {
					tmplam = -slope / (2 * (f - fold - slope));
					TRACE_VAR(tmplam);
				} else {
					const doublereal rhs1 = f - fold - lambda * slope;
					const doublereal rhs2 = f2 - fold - lambda2 * slope;
					const doublereal a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
					const doublereal b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);

					if (a == 0.) {
						tmplam = -slope / (2. * b);
						TRACE_VAR(tmplam);
					} else {
						const doublereal disc = b * b - 3. * a * slope;

						if (disc < 0.) {
							tmplam = 0.5 * lambda;
							TRACE_VAR(tmplam);
						} else if (b <= 0.) {
							tmplam = (-b + sqrt(disc)) / (3. * a);
							TRACE_VAR(tmplam);
						} else {
							tmplam = -slope / (b + sqrt(disc));
							TRACE_VAR(tmplam);
						}

						if (tmplam > 0.5 * lambda) {
							tmplam = 0.5 * lambda;
							TRACE_VAR(tmplam);
						}
					}
				}
				lambda2 = lambda;
				f2 = f;
				TRACE_VAR(tmplam);
				lambda = std::max(tmplam, dLambdaFactMin * lambda);
			} else {
				lambda *= dLambdaFactMin;
			}
		}

		if (mbdyn_stop_at_end_of_iteration()) {
			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		}
	} while (iIter < iMaxIterations);

	if (uFlags & VERBOSE_MODE) {
		silent_cerr("line search warning: Maximum number of line search iterations="
			<< iMaxIterations
			<< " exceeded in line search at time step " << pDM->dGetTime() 
            << " at iteration " << iIterCnt << std::endl);
	}

	throw MaxIterations(MBDYN_EXCEPT_ARGS);
}

void
LineSearchSolver::Solve(const NonlinearProblem *pNonLinProblem,
		Solver *pS,
		const integer iMaxIter,
		const doublereal& Tol,
		integer& iIterCnt,
		doublereal& dErr,
		const doublereal& SolTol,
		doublereal& dSolErr)
{
	ASSERT(pS != NULL);
	SolutionManager *const pSM = pS->pGetSolutionManager();
	this->pS = pS;
	pNLP = pNonLinProblem;
	pRes = pSM->pResHdl();
	pSol = pSM->pSolHdl();
	Size = pRes->iGetSize();

	if (g.iGetSize() != Size) {
		TRACE("Resize temporary vectors ..." << std::endl);
		g.Resize(Size);
		p.Resize(Size);
		dXneg.Resize(Size);
	}

	bool check = false;
	iIterCnt = 0;
	dSolErr = 0.;
	dErr = 0.;
	doublereal dErrDiff = 0.;
	doublereal f;
	Residual(f, iIterCnt);

	TRACE("\t\tf(0) = " << f << std::endl);

	bool bTest = MakeResTest(pS, pNLP, *pRes, 1e-2 * Tol, dErr, dErrDiff); // use a more stringent test for the first iteration
	doublereal dPrevErr = std::numeric_limits<doublereal>::max(); // disable error test for the first iteration
	doublereal dErrFactor = 1.;

	if (outputIters()) {
#ifdef USE_MPI
		if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
		{
			silent_cout("\tIteration(" << iIterCnt << ") " << std::setw(12) << dErr);

			if (uFlags & PRINT_CONVERGENCE_INFO) {
				silent_cout(" f=" << std::setw(12) << f);
			}

			silent_cout(std::endl);
		}
	}

	pS->CheckTimeStepLimit(dErr, dErrDiff);

	if (bTest) {
		return;
	}

	const doublereal stpmax = (uFlags & SCALE_NEWTON_STEP)
		? dMaxStep * std::max(sqrt(pDM->GetpXPCurr()->Dot() + pDM->GetpXCurr()->Dot()), static_cast<doublereal>(Size))
		: std::numeric_limits<doublereal>::max();

	TRACE_VAR(stpmax);
    
	while (true) {
		Jacobian();
		ASSERT(pSM->pMatHdl()->iGetNumCols() == Size);
		ASSERT(pSM->pMatHdl()->iGetNumCols() == pRes->iGetSize());
		ASSERT(pSM->pMatHdl()->iGetNumRows() == pRes->iGetSize());
		g.Reset();
		pSM->pMatHdl()->MatTVecDecMul(g, *pRes); // compute gradient g = \nabla f = fjac^T \, fvec = -Jac^T \, pRes
		const doublereal fold = f;
		pSM->Solve();

   		if (outputSol()) {
			pS->PrintSolution(*pSol, iIterCnt);
		}
	
		LineSearch(stpmax, fold, f, check, iIterCnt);
		bTest = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);

		if (iIterCnt > 0) {
			dErrFactor *= dErr / dPrevErr;
		}

		iIterCnt++;

#ifdef USE_MPI
		if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
		{
			if (outputIters() || outputSolverConditionNumber()) {
				if (outputIters()) {
					silent_cout("\tIteration(" << iIterCnt << ") " << std::setw(12) << dErr << " J");
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

				if (uFlags & PRINT_CONVERGENCE_INFO) {
					silent_cout(" f=" << std::setw(12) << f);

					if (iIterCnt > 1) {
						silent_cout(" Err(n-1)=" << std::setw(12) << dPrevErr
								 << " Err(n)/Err(n-1)=" << std::setw(12) << dErr / dPrevErr
								 << " Err(n)/Err(1)=" << std::setw(12) << dErrFactor);
					}
				}

				if (outputIters()
						|| outputSolverConditionNumber()
						|| (uFlags & PRINT_CONVERGENCE_INFO)) {
					silent_cout(std::endl);
				}
			}
		}

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
					<< dSolErr << std::endl);
			}
		}

		if (bTest) {
			if (uFlags & VERBOSE_MODE) {
				silent_cerr("line search warning: Convergence on solution "
                    "at time step " << pDM->dGetTime()
					<< "at iteration " << iIterCnt
					<< "\tErr(n)=" << dErr);

				if (iIterCnt >= 1) {
					silent_cerr("\tErr(n-1)=" << dPrevErr
					<< "\tErr(n)/Err(n-1)=" << dErr/dPrevErr
					<< "\tErr(n)/Err(1) " << dErrFactor);
				}

				silent_cerr(std::endl);
			}
			throw ConvergenceOnSolution(MBDYN_EXCEPT_ARGS);
		}

		if (check) { // lambda <= dLambdaMinEff: check for gradient zero
			doublereal test = 0.;
			const doublereal den = std::max(f, 0.5 * Size);

			for (integer i = 1; i <= Size; ++i) {
				const doublereal absX = std::max(std::abs((*pDM->GetpXCurr())(i)),
						 	 	 	 	 	 	 std::abs((*pDM->GetpXPCurr())(i)));
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
						<< " > Tol = " << Tol << std::endl);
				}

				if (uFlags & ZERO_GRADIENT_CONTINUE) {
					return;
				} else {
					throw ZeroGradient(MBDYN_EXCEPT_ARGS);
				}
			} else {
				if (uFlags & VERBOSE_MODE) {
					silent_cerr("line search warning: lambda min"
						" has been reached at time step " << pDM->dGetTime() 
                        << " at iteration " << iIterCnt
						<< " but the gradient is not zero" << std::endl);
				}

				if (uFlags & ABORT_AT_LAMBDA_MIN) {
					throw NoConvergence(MBDYN_EXCEPT_ARGS);
				}
			}
		} else if (dErr > dPrevErr) { // f <= fold + dAlpha * lambda * slope
			// We should not get here since f was decreased also Err should be decreased
			// If not it could be a numerical problem (e.g. the tolerance could be too small)
			if (uFlags & VERBOSE_MODE) {
				silent_cerr("line search warning: f has been reduced during line search but the residual could not be reduced at time step "
						<< pDM->dGetTime()
						<< " at iteration " << iIterCnt << std::endl);
			}

			// Do not throw an exception here because if we specify for example
			//		tolerance: <<Tol>>, test, minmax;
			// it could happen that the the norm of the residual vector is decreased
			// but the maximum residual is not!

			// In any case we will check for divergence if DIVERGENCE_CHECK is enabled.
		}

		if (dErrFactor > dDivergenceCheck) {
			if (uFlags & DIVERGENCE_CHECK) {
				if (uFlags & VERBOSE_MODE) {
					silent_cerr("line search warning: The residual could not be decreased"
            					" at time step " << pDM->dGetTime() 
                                << " at iteration " << iIterCnt << std::endl);

					if (iIterCnt > 1) {
						silent_cerr("Err(n-1)=" << dPrevErr
									<< "\tErr(n)=" << dErr
									<< "\tErr(n)/Err(n-1)=" << dErr/dPrevErr
									<< "\tErr(n)/Err(1)=" << dErrFactor << std::endl);
					}
				}
				throw ResidualNotDecreased(MBDYN_EXCEPT_ARGS);
			}
		}

		if (iIterCnt >= std::abs(iMaxIter)) {
			if (iMaxIter < 0 && dErrFactor < 1.) {
				return;
			}
			if (outputBailout()) {
				pS->PrintResidual(*pRes, iIterCnt);
			}
			throw NoConvergence(MBDYN_EXCEPT_ARGS);
		}

		dPrevErr = dErr;
	}
}

