/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
  * Portions Copyright (C) 2003-2012
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

 /*
  * Author: Reinhard Resch
  */

#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

#include "nonlin.h"
#include "vh.h"

struct LineSearchParameters
{
	doublereal dTolX;
	doublereal dTolF;
	doublereal dTolMin;
	integer iMaxIterations;
	doublereal dMaxStep;
	doublereal dAlpha;
	doublereal dLambdaMin;
	doublereal dLambdaFactMin;
    doublereal dDivergenceCheck;
    doublereal dMinStepScale;
	enum { ZERO_GRADIENT_CONTINUE = 0x001,
		   DIVERGENCE_CHECK		  = 0x002,
           ALGORITHM_CUBIC		  = 0x004,
           ALGORITHM_FACTOR		  = 0x008,
           PRINT_CONVERGENCE_INFO = 0x010,
           SCALE_NEWTON_STEP	  = 0x020,
           RELATIVE_LAMBDA_MIN    = 0x040,
           ABORT_AT_LAMBDA_MIN    = 0x080,
           VERBOSE_MODE			  = 0x100,
           ALGORITHM = ALGORITHM_CUBIC | ALGORITHM_FACTOR };
	unsigned uFlags;
	LineSearchParameters();
};

class LineSearchSolver: public NonlinearSolver, private LineSearchParameters
{
	VectorHandler* 	pRes;
	VectorHandler* 	pSol;
	MyVectorHandler p;
	MyVectorHandler g;
	MyVectorHandler dXneg;
	const NonlinearProblem* pNLP;
    Solver* pS;
	const DataManager* const pDM;

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
                     bool bHonorJacobianRequest, 
                     const struct LineSearchParameters& param);
	~LineSearchSolver(void);
	
	void Solve(const NonlinearProblem *pNLP,
			Solver *pS,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr,
			const doublereal& SolTol,
			doublereal& dSolErr);

private:
	void LineSearch(doublereal stpmax, doublereal fold, doublereal& f, bool& check, const integer& iIterCnt);
    void Residual(doublereal& f, integer iIterCnt);
    void Jacobian();
};

#endif /* LINE_SEARCH_H */

