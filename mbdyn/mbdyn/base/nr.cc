/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
  * Copyright (C) 2004
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>

#include "solver.h"
#include "nr.h"  
#ifdef USE_MPI
#include "mbcomm.h"
#include "schsolman.h"
#endif /* USE_MPI */

#include "dofown.h"
#include "umfpackwrap.h"
#include "output.h"
#include "ac/float.h"
#include "ac/math.h"

NewtonRaphsonSolver::NewtonRaphsonSolver(const bool bTNR,
		const bool bKJ, 
		const integer IterBfAss)
: pRes(NULL),
pSol(NULL),
bTrueNewtonRaphson(bTNR),
IterationBeforeAssembly(IterBfAss),
bKeepJac(bKJ),
iPerformedIterations(0),
pPrevNLP(NULL)
{
	NO_OP;
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(void)
{
	NO_OP;
}

void
NewtonRaphsonSolver::Solve(const NonlinearProblem *pNLP,
		Solver *pS,
		const integer iMaxIter,
		const doublereal& Tol,
		integer& iIterCnt,
		doublereal& dErr ,
		const doublereal& SolTol,
		doublereal& dSolErr)
{
	ASSERT(pS != NULL);

	SolutionManager *pSM = pS->pGetSolutionManager();
	
	iIterCnt = 0;
	if (!bKeepJac || (pNLP != pPrevNLP)) {
		iPerformedIterations = 0;
	}
	pPrevNLP = pNLP;
	dSolErr = 0.;

	while (true) {
		pRes = pSM->pResHdl();
		pSol = pSM->pSolHdl();
		Size = pRes->iGetSize();
#ifdef 	USE_EXTERNAL 	
		SendExternal();
#endif /* USE_EXTERNAL */
		
		pRes->Reset();
		bool forceJacobian(false);
		try {
	      		pNLP->Residual(pRes);
		}
		catch (SolutionDataManager::ChangedEquationStructure) {
			forceJacobian = true;
		}
		
      		if (outputRes()) {
	 		silent_cout("Residual (" << iIterCnt 
				<< "):" << std::endl);
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			silent_cout("Dof " << std::setw(8)
					<< iTmpCnt << ": " 
					<< pRes->dGetCoef(iTmpCnt)
					<< std::endl);
			}
      		}

		/* FIXME: if Tol == 0., no convergence on residual
		 * is required, so we could simply don't compute
		 * the test; I'm leaving it in place so it appears
		 * in the output (maybe we could conditionally disable 
		 * it? */

		dErr = MakeResTest(pS, *pRes)*pNLP->TestScale(pResTest);

		if (outputIters()) {
#ifdef USE_MPI
			if (MBDynComm.Get_rank() == 0) {
#endif /* USE_MPI */
				silent_cout("\tIteration " << iIterCnt
					<< " " << dErr);
				if (dErr >= Tol && (bTrueNewtonRaphson || (iPerformedIterations%IterationBeforeAssembly == 0))){
					silent_cout(" J");
				}
				silent_cout(std::endl);
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		
      		if (dErr < Tol) {
	 		return;
      		}
      		
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
		}
          
      		iIterCnt++;

		if (bTrueNewtonRaphson || (iPerformedIterations%IterationBeforeAssembly == 0)
			|| forceJacobian) {
      			pSM->MatrReset();
rebuild_matrix:;
			try {
      				pNLP->Jacobian(pSM->pMatHdl());

			} catch (MatrixHandler::ErrRebuildMatrix) {
				silent_cout("NewtonRaphsonSolver: "
						"rebuilding matrix..."
						<< std::endl);

				/* need to rebuild the matrix... */
      				pSM->MatrInitialize();
				goto rebuild_matrix;

			} catch (...) {
				throw;
			}

			TotJac++;
		}
		
		iPerformedIterations++;
		
		if (outputJac()) {
			silent_cout("Jacobian:" << std::endl
					<< *(pSM->pMatHdl()));
		}		
		
		pSM->Solve();

      		if (outputSol()) {      
	 		silent_cout("Solution:" << std::endl);
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			silent_cout("Dof " << std::setw(8) << iTmpCnt
						<< ": "
						<< pSol->dGetCoef(iTmpCnt)
						<< std::endl);
			}
		}		
		
      		pNLP->Update(pSol);
		
		dSolErr = MakeSolTest(pS, *pSol);
		if (outputIters()) {
#ifdef USE_MPI
			if (MBDynComm.Get_rank() == 0) {
#endif /* USE_MPI */
				silent_cout("\t\tSolErr "
					<< dSolErr << std::endl);
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		if (dSolErr < SolTol) {
			THROW(ConvergenceOnSolution());
		}
	}
}

