/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
  * Copyright (C) 2003
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
  
  
  
#include <bicg.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <umfpackwrap.h>
#include <unistd.h>
#include <output.h>

BiCGStab::BiCGStab(const Preconditioner::PrecondType PType, 
		const integer iPStep,
		doublereal ITol,
		integer MaxIt,
		doublereal etaMx) 
: MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx)
{
	NO_OP;
}
	
BiCGStab::~BiCGStab(void)
{
	NO_OP;
}

void
BiCGStab::Solve(const NonlinearProblem* pNLP,
		SolutionManager* pSolMan,
		const integer iMaxIter,
		const doublereal Tol,
		const doublereal SolTol,
		integer& iIterCnt,
		doublereal& dErr
#ifdef MBDYN_X_CONVSOL
		, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
	       )
{
	ASSERT(pNLP != NULL);
	ASSERT(pSM != NULL);
	
	iIterCnt = 0;

#ifdef MBDYN_X_CONVSOL
	dSolErr = 0.;
#endif /* MBDYN_X_CONVSOL  */	

	/* external nonlinear iteration */	
	
	/* riassembla sempre lo jacobiano se l'integratore e' nuovo */
	if (pNLP != pPrevNLP) {
		fBuildMat = true;
	}
	
	pPrevNLP = pNLP;
	
	pSM  = pSolMan;
        pRes = pSM->pResHdl();
	Size = pRes->iGetSize();

	doublereal eta = etaMax;
	doublereal rateo = 0.;
	doublereal Fnorm = 1.;
        doublereal resid;
        doublereal rho_1; 
	doublereal rho_2 = 0.; 
	doublereal alpha;
	doublereal beta;
	doublereal omega;
	VectorHandler* pr;
	MyVectorHandler rHat(Size), p(Size), pHat(Size), s(Size), sHat(Size), t(Size), v(Size);
	MyVectorHandler dx(Size); 
	integer TotalIter = 0;

#ifdef DEBUG_ITERATIVE
	std::cout << " BiCGStab New Step " <<std::endl;
#endif /* DEBUG_ITERATIVE */

	while (true) {

#ifdef 	USE_EXTERNAL 	
		SendExternal();
#endif /* USE_EXTERNAL */
		
		pRes->Reset(0.);
      		pNLP->Residual(pRes);
		
      		if (foutRes) {
	 		std::cout << "Residual:" << std::endl;
	 		std::cout << iIterCnt <<std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< pRes->dGetCoef(iTmpCnt) << std::endl;
			}
      		}

#ifdef __HACK_SCALE_RES__
		dErr = pNLP->TestScale(pScale)*MakeTest(*pRes);		
#else /* ! __HACK_SCALE_RES__ */
		dErr = pNLP->TestScale()*MakeTest(*pRes);		
#endif /* ! __HACK_SCALE_RES__ */

#ifdef DEBUG_ITERATIVE
		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (dErr < Tol) {
	 		return;
      		}
      		if (!isfinite(dErr)) {
			THROW(ErrSimulationDiverged());
		}
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
		}
          	rateo = dErr*dErr/Fnorm;
		Fnorm = dErr*dErr;
		
      		iIterCnt++;

		/* inner iteration to solve the linear system */	
	
		/* BiCGSTAB Iterative solver */
		DEBUGCOUT("Using BiCGStab iterative solver" << std::endl);

        	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */
        	/* N.B. *pRes = -F(0) */ 
		
		pr = pRes;
        	doublereal LocTol = eta * dErr;
        	rHat = *pr;

#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocTol " << LocTol << std::endl;
#endif /* DEBUG_ITERATIVE */
		
		rho_1 = dErr*dErr;   /*rhat.InnerProd(r); */
		resid = dErr;
		v.Reset(0.);
		t.Reset(0.);
		p.Reset(0.);
		dx.Reset(0.);
		
		if (fBuildMat) {
			pSM->MatrInit(0.);
			pNLP->Jacobian(pSM->pMatHdl());
			fBuildMat = false;
			TotalIter = 0;
			TotJac++;

#ifdef DEBUG_ITERATIVE			
			std::cerr << "Jacobian " << std::endl;
#endif /* DEBUG_ITERATIVE */

		}

#ifdef DEBUG_ITERATIVE
		std::cerr << "rho_1 " << rho_1 << std::endl;
#endif /* DEBUG_ITERATIVE */

		int It = 0;
        	while ((resid > LocTol) && (It++ < MaxLinIt)) {
			if (It == 1) {
				p = *pr;
			} else {
				rho_1 = rHat.InnerProd(*pr);

#ifdef DEBUG_ITERATIVE
				std::cerr << "rho_1 " << rho_1 << std::endl;
#endif /* DEBUG_ITERATIVE */

	               		if (fabs(rho_1) < DBL_EPSILON) {
                        		std::cout << "Bi-CGStab Iterative Solver breakdown" 
						<<  " rho_1 = 0 " << std::endl;
					break;
				}
				beta = (rho_1/rho_2) * (alpha/omega);

#ifdef DEBUG_ITERATIVE
				std::cerr << "beta " << beta << std::endl;
#endif /* DEBUG_ITERATIVE */

				p.ScalarAddMul(*pr, p.ScalarAddMul(v, -omega), beta);
			}
			/* right preconditioning */
			pPM->Precond(p, pHat, pSM);
			pNLP->EvalProd(Tau, rHat, pHat, v);
# if 0			
			(pSM->pMatHdl())->MatVecMul(v,pHat);
#endif			
#ifdef DEBUG_ITERATIVE
			std::cout << "v:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< v.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			alpha = rHat.InnerProd(v);
			alpha = rho_1 / alpha;

#ifdef DEBUG_ITERATIVE
			std::cerr << "alpha " << alpha << std::endl;
#endif /* DEBUG_ITERATIVE */

			s.ScalarAddMul(*pr, v, -alpha);

#ifdef DEBUG_ITERATIVE
			std::cerr << "s.Norm() " << s.Norm() << std::endl;
#endif /* DEBUG_ITERATIVE */

			if ((resid = s.Norm()) < LocTol) {
				dx.ScalarAddMul(pHat, alpha);
				TotalIter++;
				break;
			}
			pPM->Precond(s, sHat, pSM);
			pNLP->EvalProd(Tau, rHat, sHat, t);
#if 0
			(pSM->pMatHdl())->MatVecMul(t,sHat);
#endif
			omega = t.Norm();
			omega = t.InnerProd(s) / omega;

#ifdef DEBUG_ITERATIVE
			std::cerr << "omega " << omega << std::endl;
#endif /* DEBUG_ITERATIVE */

			dx.ScalarAddMul(pHat, alpha);
			dx.ScalarAddMul(sHat, omega);
			pr->ScalarAddMul(s, t, -omega);
			rho_2 = rho_1;
			resid = pr->Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "resid " << resid << std::endl;
#endif /* DEBUG_ITERATIVE */

			TotalIter++;
                	if (fabs(omega) < DBL_EPSILON) {
                        	std::cout << "Bi-CGStab Iterative Solver breakdown" 
					<<  " omega = 0 " << std::endl;
				break;
			}
			if ( It == MaxLinIt) {
                        	std::cerr << "Iterative inner solver didn't converge."
					<< " Continuing..." << std::endl;
			}
		}
		/* se ha impiegato troppi passi riassembla lo jacobiano */
		
		if (TotalIter >= PrecondIter) {
			fBuildMat = true;
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
		
		if (foutSol) {      
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
					<< dx.dGetCoef(iTmpCnt) << std::endl;
			}
		}		
		
		
		if (foutIters) {
#ifdef USE_MPI
			if (dynamic_cast<SchurSolutionManager*> (pSM) && (MBDynComm.Get_rank() == 0)) {
#endif /* USE_MPI */
				std::cerr << "\tIteration " << iIterCnt
					<< " " << dErr << " J"
					<< std::endl;
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		
      		pNLP->Update(&dx);

		
#ifdef MBDYN_X_CONVSOL
		if (SolTol > 0.) {
			dSolErr = MakeTest(dx);
        		if (dSolErr < SolTol) {
				THROW(ConvergenceOnSolution());
			}
      		}
#endif /* MBDYN_X_CONVSOL */
	}
}

