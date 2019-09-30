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

#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

#include "nonlin.h"
#include "vh.h"

struct LineSearchParameters
{
	doublereal dTolX;
	doublereal dTolMin;
	integer iMaxIterations;
	doublereal dMaxStep;
    doublereal dAlphaFull;
    doublereal dAlphaModified;
	doublereal dLambdaMin;
	doublereal dLambdaFactMin;
    doublereal dDivergenceCheck;
    doublereal dMinStepScale;
	enum { ZERO_GRADIENT_CONTINUE 	   = 0x001,
		   DIVERGENCE_CHECK		  	   = 0x002,
           ALGORITHM_CUBIC		  	   = 0x004,
           ALGORITHM_FACTOR		  	   = 0x008,
           PRINT_CONVERGENCE_INFO 	   = 0x010,
           SCALE_NEWTON_STEP	  	   = 0x020,
           RELATIVE_LAMBDA_MIN    	   = 0x040,
           ABORT_AT_LAMBDA_MIN    	   = 0x080,
           VERBOSE_MODE			  	   = 0x100,
           NON_NEGATIVE_SLOPE_CONTINUE = 0x200,
           ALGORITHM = ALGORITHM_CUBIC | ALGORITHM_FACTOR };
	unsigned uFlags;
    integer iIterationsBeforeAssembly;
    bool bKeepJac;
    doublereal dTimeStepTol;
    doublereal dUpdateRatio;
	LineSearchParameters();
};

class LineSearchSolver: public NonlinearSolver, protected LineSearchParameters
{
private:
    doublereal lambda2;
    doublereal f2;
    
protected:
	VectorHandler* 	pRes;
	VectorHandler* 	pSol;
	MyVectorHandler p;
	MyVectorHandler g;
	const NonlinearProblem* pNLP;
    Solver* pS;
    SolutionManager* pSM;
    DataManager* const pDM;
    integer iRebuildJac;

    struct CPUTime
    {
        explicit CPUTime(LineSearchSolver& oSolver)
            :Residual(oSolver, CPU_RESIDUAL),
             Jacobian(oSolver, CPU_JACOBIAN),
             LinearSolver(oSolver, CPU_LINEAR_SOLVER)
        {
        }
        CPUStopWatch Residual;
        CPUStopWatch Jacobian;
        CPUStopWatch LinearSolver;
    };

public:
	class SlopeNotNegative: public NoConvergence
	{
	public:
		SlopeNotNegative(MBDYN_EXCEPT_ARGS_DECL)
			:NoConvergence(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{

		}
	};

	class ZeroGradient: public NoConvergence
	{
	public:
		ZeroGradient(MBDYN_EXCEPT_ARGS_DECL)
			:NoConvergence(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{

		}
	};

	class MaxIterations: public NoConvergence
	{
	public:
		MaxIterations(MBDYN_EXCEPT_ARGS_DECL)
			:NoConvergence(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{

		}
	};

    class ResidualNotDecreased: public ErrSimulationDiverged
    {
    public:
        ResidualNotDecreased(MBDYN_EXCEPT_ARGS_DECL)
            :ErrSimulationDiverged(MBDYN_EXCEPT_ARGS_PASSTHRU)
        {

        }
    };

	LineSearchSolver(DataManager* pDM,
                     const NonlinearSolverOptions& options,
                     const struct LineSearchParameters& param);
	~LineSearchSolver(void);
	
protected:
    void Attach(const NonlinearProblem* pNLP, Solver* pS);
    doublereal dGetMinNewtonInc(const VectorHandler& dX) const;
    doublereal dGetLambdaNext(doublereal dLambdaCurr, doublereal dSlope, doublereal fPrev, doublereal fCurr);
    doublereal dGetLambdaMin(doublereal& dSlope, bool bRebuildJac, const VectorHandler& p, integer iIterCnt) const;
    doublereal dGetMaxNewtonStep(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr) const;
    void ScaleNewtonStep(doublereal stpmax, VectorHandler& dX, integer iIterCnt) const;
    bool bCheckZeroGradient(doublereal fCurr,
                            doublereal dErr,
                            doublereal dTol,
                            integer iIterCnt) const;
    bool bCheckDivergence(doublereal dErrFactor,
                          doublereal dErrPrev,
                          doublereal dErr,
                          integer iIterCnt) const;
    void OutputIteration(doublereal dErr, integer iIterCnt, bool bRebuildJac, const CPUTime& oCPU);
    void OutputLineSearch(integer iIterCnt,
                          integer iLineSearchIter,
                          doublereal fCurr,
                          doublereal dErr,
                          doublereal dLambda,
                          doublereal dSlope) const;
    void Residual(doublereal& f, integer iIterCnt);
    void Jacobian();
};

class LineSearchFull: public LineSearchSolver
{    
public:
    LineSearchFull(DataManager* pDM,
                   const NonlinearSolverOptions& options,
                   const struct LineSearchParameters& param);
    ~LineSearchFull(void);
	
    virtual void Solve(const NonlinearProblem *pNLP,
			Solver *pS,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr,
			const doublereal& SolTol,
			doublereal& dSolErr);

private:
   bool LineSearch(doublereal dMaxStep, doublereal fPrev, doublereal& fCurr, integer iIterCnt);
};

class LineSearchModified: public LineSearchSolver
{
private:
    doublereal dTimePrev, dTimeStepPrev;
    
public:
    LineSearchModified(DataManager* pDM,
                       const NonlinearSolverOptions& options,
                       const struct LineSearchParameters& param);
    ~LineSearchModified(void);
	
    virtual void Solve(const NonlinearProblem *pNLP,
                       Solver *pS,
                       const integer iMaxIter,
                       const doublereal& Tol,
                       integer& iIterCnt,
                       doublereal& dErr,
                       const doublereal& SolTol,
                       doublereal& dSolErr);
};

class LineSearchBFGS: public LineSearchSolver
{
private:
        void Attach(const NonlinearProblem* pNLP, Solver* pS);
        MyVectorHandler s, t, w, FCurr, FPrev;

public:
        LineSearchBFGS(DataManager* pDM,
                       const NonlinearSolverOptions& options,
                       const struct LineSearchParameters& param);
        ~LineSearchBFGS(void);

        virtual void Solve(const NonlinearProblem *pNLP,
                           Solver *pS,
                           const integer iMaxIter,
                           const doublereal& Tol,
                           integer& iIterCnt,
                           doublereal& dErr,
                           const doublereal& SolTol,
                           doublereal& dSolErr);
};

#endif /* LINE_SEARCH_H */
