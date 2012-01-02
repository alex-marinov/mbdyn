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
  *
  * Copyright (C) 2003-2012
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

#ifndef NR_H
#define NR_H

#include "nonlin.h"

class NewtonRaphsonSolver : public NonlinearSolver
{
	VectorHandler* 	pRes;
	VectorHandler* 	pSol;
	bool bTrueNewtonRaphson;
	integer IterationBeforeAssembly;
	bool bKeepJac;
	integer iPerformedIterations;
	const NonlinearProblem* pPrevNLP;	

public:
	NewtonRaphsonSolver(const bool bTNR,
			const bool bKJ, 
			const integer IterBfAss,
			bool JacReq = false);
	
	~NewtonRaphsonSolver(void);
	
	void Solve(const NonlinearProblem *pNLP,
			Solver *pS,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr,
			const doublereal& SolTol,
			doublereal& dSolErr);
};

#endif /* NONLIN_H */

