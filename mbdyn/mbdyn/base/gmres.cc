/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
  * Copyright (C) 2003-2013
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>
#include <cstdlib>

#include "solver.h"
#include "gmres.h"  

#ifdef USE_MPI
#include "mbcomm.h"
#include "schsolman.h"
#endif /* USE_MPI */

#include "dofown.h"
#include "output.h"

Gmres::Gmres(const Preconditioner::PrecondType PType, 
		const integer iPStep,
		doublereal ITol,
		integer MaxIt,
		doublereal etaMx,
		doublereal T,
		const NonlinearSolverOptions& options)
: MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx, T, options),
v(NULL),
s(MaxLinIt + 1), cs(MaxLinIt + 1), sn(MaxLinIt + 1)
{
	SAFENEWARRNOFILL(v, MyVectorHandler, MaxLinIt + 1); 
}
	
Gmres::~Gmres(void)
{
	for (int i = 0; i < MaxLinIt + 1; i++) {
		v[i].Detach();
	}

	SAFEDELETEARR(v);
}
	
void
Gmres::GeneratePlaneRotation(const doublereal &dx, const doublereal &dy, 
		doublereal &cs, doublereal &sn) const 
{
	if (fabs(dy) < std::numeric_limits<doublereal>::epsilon()) {
		cs = 1.0;
		sn = 0.0;

	} else if (fabs(dy) > fabs(dx)) {
		doublereal temp = dx / dy; 
		sn = 1.0 / sqrt( 1.0 + temp*temp );
		cs = temp * sn;

	} else {
		doublereal temp = dy / dx; 
		cs = 1.0 / sqrt( 1.0 + temp*temp );
		sn = temp * cs;
	}
}

void
Gmres::ApplyPlaneRotation(doublereal &dx, doublereal &dy, 
		const doublereal &cs, const doublereal &sn) const 
{ 
	doublereal temp = cs * dx + sn * dy; 
	dy = -sn * dx + cs * dy;
	dx = temp;
}

void
Gmres::Backsolve(VectorHandler& x, integer sz,
		VectorHandler& s, MyVectorHandler* v) 
{ 
	for (int i = sz+1; i > 0; i--) {
    		s.PutCoef(i, s(i) / H(i, i));
    		for (int j = i - 1; j > 0; j--) {
      			s.DecCoef(j, H(j, i) * s(i));
		}
  	}

  	for (int j = 0; j <= sz; j++) {
    		x.ScalarAddMul(v[j], s(j+1));
	}
}

void
Gmres::Solve(const NonlinearProblem* pNLP,
		Solver* pS,
		const integer iMaxIter,
		const doublereal& Tol,
		integer& iIterCnt,
		doublereal& dErr,
		const doublereal& SolTol,
		doublereal& dSolErr)
{
	ASSERT(pNLP != NULL);
	ASSERT(pS != NULL);

	SolutionManager *pSM = pS->pGetSolutionManager();
	
	iIterCnt = 0;
	dSolErr = 0.;
	doublereal dErrDiff = 0.;

	/* external nonlinear iteration */	
	
	/* riassembla sempre lo jacobiano se l'integratore e' nuovo */
	if (pNLP != pPrevNLP) {
		bBuildMat = true;
	}
	
	pPrevNLP = pNLP;
	
        pRes = pSM->pResHdl();
	Size = pRes->iGetSize();

	doublereal eta = etaMax;
	doublereal rateo = 0.;
	doublereal Fnorm = 1.;
        doublereal resid;
	doublereal norm1, norm2;
	VectorHandler* pr;

	/*
	 * these will be resized (actually allocated)
	 * only the first time they're used, unless
	 * the size of the problem changes
	 *
	 * FIXME: need to review this code.
	 */
	w.Resize(Size);
	vHat.Resize(Size);
	dx.Resize(Size);

	H.Resize(Size, MaxLinIt + 1);

	integer TotalIter = 0;

#ifdef DEBUG_ITERATIVE
	std::cout << " Gmres New Step " << std::endl;
#endif /* DEBUG_ITERATIVE */

	doublereal dOldErr = 0.;
	doublereal dErrFactor = 1.;
	while (true) {

#ifdef 	USE_EXTERNAL 	
		SendExternal();
#endif /* USE_EXTERNAL */
		pRes->Reset();
		try {
	      		pNLP->Residual(pRes);
		}
		catch (SolutionDataManager::ChangedEquationStructure) {
			if (bHonorJacRequest) {
				bBuildMat = true;
			}
		}
		
      		if (outputRes()) {
			pS->PrintResidual(*pRes, iIterCnt);
      		}

		bool bTest = MakeResTest(pS, pNLP, *pRes, Tol, dErr, dErrDiff);
		if (iIterCnt > 0) {
			dErrFactor *= dErr/dOldErr;
		}
		dOldErr = dErr;

#ifdef DEBUG_ITERATIVE
      		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (outputIters()) {
#ifdef USE_MPI
			if (!bParallel || MBDynComm.Get_rank() == 0)
#endif /* USE_MPI */
			{
				silent_cout("\tIteration(" << iIterCnt << ") " << dErr);
				if (bBuildMat && !bTest) {
					silent_cout(" J");
				}
				silent_cout(std::endl);
			}
		}
		
		pS->CheckTimeStepLimit(dErr, dErrDiff);

		if (bTest) {
	 		return;
      		}
      		if (!std::isfinite(dErr)) {
			throw ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
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
          	rateo = dErr*dErr/Fnorm;

#ifdef DEBUG_ITERATIVE
		std::cerr << "rateo " << rateo << std::endl;
#endif /* DEBUG_ITERATIVE */

		Fnorm = dErr*dErr;
		
      		iIterCnt++;

		/* inner iteration to solve the linear system */	
	
		/* Gmres(m) Iterative solver */
		DEBUGCOUT("Using Gmres(m) iterative solver" << std::endl);

        	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */
        	/* N.B. *pRes = -F(0) */ 
		
		pr = pRes;
        	doublereal LocTol = eta * dErr;
 
#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocTol " << LocTol << std::endl;
#endif /* DEBUG_ITERATIVE */
	
		resid = dErr;
		dx.Reset();
		
		if (bBuildMat) {
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

			bBuildMat = false;
			TotalIter = 0;
			TotJac++;
			
#ifdef DEBUG_ITERATIVE
			std::cerr << "Jacobian " << std::endl;
#endif /* DEBUG_ITERATIVE */

		}

		int i = 0;
        	v[0].Resize(Size);
		v[0].ScalarMul(*pr, 1./resid);
		s.Reset();
		cs.Reset();
		sn.Reset();
		s.PutCoef(1, resid);
		while ((i < MaxLinIt)) {

#ifdef DEBUG_ITERATIVE
			std::cout << "v[i]:" << i << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< v[i](iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pPM->Precond(v[i], vHat, pSM); 

#ifdef DEBUG_ITERATIVE
			std::cout << "vHat:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< vHat(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pNLP->EvalProd(Tau, *pr, vHat, w);
			
#if 0
			(pSM->pMatHdl())->MatVecMul(w, vHat);
#endif

#ifdef DEBUG_ITERATIVE
			std::cout << "w:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< w(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			norm1 = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm1: " << norm1 << std::endl; 
#endif /* DEBUG_ITERATIVE */

			for (int k = 0; k <= i; k++) {
        			H(k+1, i+1) = w.InnerProd(v[k]);

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k, i): " << k+1 << " " << i+1 << " "<< H(k+1, i+1) << std::endl; 
#endif /* DEBUG_ITERATIVE */

				w.ScalarAddMul(v[k], -H(k+1, i+1));
      			}
			H(i+2, i+1) = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "w.Norm(): " << w.Norm() << std::endl;
			std::cerr << "H(i+1, i): " << i+2 << " " << i+1 << " " << H(i+2, i+1) << std::endl;
#endif /* DEBUG_ITERATIVE */

			norm2 = H(i+2, i+1); 

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm2: " << norm2 << std::endl;
#endif /* DEBUG_ITERATIVE */

			v[i+1].Resize(Size);
			
			/*  Reorthogonalize? */
			if  (.001*norm2/norm1 < std::numeric_limits<doublereal>::epsilon()) {
    				for (int k = 0;  k <= i; k++) {       
					doublereal hr = v[k].InnerProd(w);
        				H(k+1, i+1) += hr;
        				w.ScalarAddMul(v[k], -hr);

#ifdef DEBUG_ITERATIVE
					std::cerr << "Reorthogonalize: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				}
				H(i+2, i+1) = w.Norm();
			}
			if (fabs(H(i+2, i+1)) > std::numeric_limits<doublereal>::epsilon()) { 
				v[i+1].ScalarMul(w, 1./H(i+2, i+1));

#ifdef DEBUG_ITERATIVE
				std::cout << "v[i+1]:" << i+1 << std::endl;

	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< v[i+1](iTmpCnt) << std::endl;
				}
#endif /* DEBUG_ITERATIVE */
			} else {
				/* happy breakdown !! */

#ifdef DEBUG_ITERATIVE
				std::cerr << "happy breakdown: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				v[i+1].Reset();
			} 
			for (int k = 0; k < i; k++) {
        			ApplyPlaneRotation(H(k+1, i+1), H(k+2, i+1),
						cs(k+1),
						sn(k+1));

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k, i): " << k+1 << " " << i+1 << " " << H(k+1, i+1) << std::endl;
				std::cerr << "H(k+1, i): " << k+2 << " " << i+1 << " " << H(k+2, i+1) << std::endl;
#endif /* DEBUG_ITERATIVE */

			}
							
			GeneratePlaneRotation(H(i+1, i+1), H(i+2, i+1),
					cs(i+1), sn(i+1));

#ifdef DEBUG_ITERATIVE
			std::cerr << "cs(i): " << cs(i+1) << std::endl;
			std::cerr << "sn(i): " << sn(i+1) << std::endl;
#endif /* DEBUG_ITERATIVE */

      			ApplyPlaneRotation(H(i+1, i+1), H(i+2, i+1),
					cs(i+1), sn(i+1));
      			ApplyPlaneRotation(s(i+1), s(i+2), cs(i+1), sn(i+1));
			if ((resid = fabs(s(i+2))) < LocTol) {

#ifdef DEBUG_ITERATIVE
				std::cerr << "resid 1: " << resid  << std::endl;
#endif /* DEBUG_ITERATIVE */

				Backsolve(dx, i, s, v);
				pPM->Precond(dx, dx, pSM); 

#ifdef DEBUG_ITERATIVE
				std::cout << "dx:" << std::endl;
	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< dx(iTmpCnt) << std::endl;
				}
#endif /* DEBUG_ITERATIVE */

        			break;
      			}

#ifdef DEBUG_ITERATIVE
			std::cerr << "resid 2 : " << resid << std::endl;
#endif /* DEBUG_ITERATIVE */

			TotalIter++;
			i++;
		}
		if (i == MaxLinIt) {
			Backsolve(dx, MaxLinIt-1, s, v);
			pPM->Precond(dx, dx, pSM);				
		        silent_cerr("Iterative inner solver didn't converge."
				<< " Continuing..." << std::endl);
		}

		/* se ha impiegato troppi passi riassembla lo jacobiano */
		
		if (TotalIter >= PrecondIter) {
			bBuildMat = true;
		}
		/* calcola il nuovo eta */
		doublereal etaNew = gamma * rateo;
		doublereal etaBis;
		if ((etaBis = gamma*eta*eta) > .1) {
			etaNew = (etaNew > etaBis) ? etaNew : etaBis;
		}
		eta = (etaNew < etaMax) ? etaNew : etaMax;
		/* prevent oversolving */
		etaBis = .5*Tol/dErr;
		eta = (eta > etaBis) ? eta : etaBis;

#ifdef DEBUG_ITERATIVE
		std::cerr << "eta " << eta << std::endl;
#endif /* DEBUG_ITERATIVE */
		
		if (outputSol()) {
			pS->PrintSolution(dx, iIterCnt);
		}		
		
      		pNLP->Update(&dx);
		
		bTest = MakeSolTest(pS, dx, SolTol, dSolErr);
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
			throw ConvergenceOnSolution(MBDYN_EXCEPT_ARGS);
      		}

		// allow to bail out in case of multiple CTRL^C
		if (mbdyn_stop_at_end_of_iteration()) {
			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		}
	}
}

