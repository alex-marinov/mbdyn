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
  
  
  
#include <nonlin.h>  
#ifdef USE_MPI
#include <mpi++.h>
extern MPI::Intracomm MBDynComm;
#include<schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <umfpackwrap.h>
#include <unistd.h>
#include <output.h>

doublereal NewtonRaphsonSolver::MakeTest(const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("NewtonRaphsonSolver::MakeTest");

      
   	doublereal dRes = 0.;
	ASSERT(pSM != NULL);
	
#ifdef USE_MPI
#warning "NonlinSolver MakeTest parallel broken !! "	
#if 0
   	Dof CurrDof; 
	SchurSolutionManager *pSSM;
	if ((pSSM = dynamic_cast<SchurSolutionManager*> (pSM)) != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
			/* else if ALGEBRAIC: non aggiunge nulla */
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

#ifdef __HACK_RES_TEST__
		for (int iCnt = 0; iCnt < iMI; iCnt++) {
			DCount = pMI[iCnt];
			CurrDof = pDofs[DCount-1];
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				doublereal d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
#endif /* __HACK_RES_TEST__ */
			/* else if ALGEBRAIC: non aggiunge nulla */
		}
		
		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(dRes, dXPr);
		
	} else {
#endif		
#endif /* USE_MPI */
			ASSERT(Vec.iGetSize() == Size);
#ifdef __HACK_SCALE_RES__
			ASSERT(pScale != NULL);
			ASSERT(pScale->iGetSize == Size);
#endif /* __HACK_SCALE_RES__ */
 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = Vec.dGetCoef(iCntp1);
			doublereal d2 = d*d;

#ifdef __HACK_SCALE_RES__
			doublereal ds = pScale->dGetCoef(iCntp1);
			doublereal ds2 = ds*ds;
			d2 *= ds2;
#endif /* __HACK_SCALE_RES__ */

			dRes += d2;

#ifdef __HACK_RES_TEST__
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(iCntp1);
				d2 = d*d;

#ifdef __HACK_SCALE_RES__
				d2 *= ds2;         

#endif /* __HACK_SCALE_RES__ */

				dXPr += d2;
			}
			/* else if ALGEBRAIC: non aggiunge nulla */ 
#endif /* __HACK_RES_TEST__ */
		}
#ifdef USE_MPI
#if 0 		
	}
#endif
#endif /* USE_MPI */

#ifdef __HACK_RES_TEST__
	dRes /= (1.+dXPr);
#endif /* __HACK_RES_TEST__ */

   	return sqrt(dRes);
}

NewtonRaphsonSolver::NewtonRaphsonSolver(const flag fTNR, 
		const integer IterBfAss)
: pSM(NULL),
pRes(NULL),
pSol(NULL),
pJac(NULL),
fTrueNewtonRaphson(fTNR),
IterationBeforeAssembly(IterBfAss) { };


void
NewtonRaphsonSolver::Solve(const NonlinearProblem* pNLP,
		SolutionManager* pSolMan,
		const integer iMaxIter,
		const doublereal Toll,
		const doublereal SolToll,
		integer& iIterCnt,
		doublereal& dErr
#ifdef MBDYN_X_CONVSOL
		, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
		)
{
	ASSERT(pNLP != NULL);
	ASSERT(pSM != NULL);
		
	pSM  = pSolMan;
	pJac = pSM->pMatHdl();
        pRes = pSM->pResHdl();
        pSol = pSM->pSolHdl();
	Size = pRes->iGetSize();
	
	iIterCnt = 0;
	integer iPerformedIterations = 0;
#ifdef MBDYN_X_CONVSOL
	dSolErr = 0.;
#endif /* MBDYN_X_CONVSOL  */	

	while (1) {

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

		dErr = MakeTest(*pRes);

      		if (dErr < Toll) {
	 		return;
      		}
      		
		if (!isfinite(dErr)) {
			THROW(ErrSimulationDiverged());
		}
		
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
		}
          
      		iIterCnt++;

		if (fTrueNewtonRaphson || (iPerformedIterations%IterationBeforeAssembly == 0)) {
      			pSM->MatrInit(0.);
      			pNLP->Jacobian(pJac);
			TotJac++;
		}
		
		iPerformedIterations++;
		
		if (foutJac) {
#ifdef USE_UMFPACK
			if (dynamic_cast<UmfpackSparseLUSolutionManager*>(pSM) == 0) {
#endif /* USE_UMFPACK */
			std::cout << "Warning, Jacobian output "
					"avaliable only "
					"with umfpack solver"
					<< std::endl;
#ifdef USE_UMFPACK
			} else {
			 	std::cout << "Jacobian:" << std::endl
					<< *(pSM->pMatHdl());
		 	}
#endif /* USE_UMFPACK */
      		}		
		
		pSM->Solve();

      		if (foutSol) {      
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
					<< pSol->dGetCoef(iTmpCnt) << std::endl;
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
		
      		pNLP->Update(pSol);

		
#ifdef MBDYN_X_CONVSOL
		if (SolToll > 0.) {
			dSolErr = MakeTest(*pSol);
        		if (dSolErr < dSolToll) {
				THROW(ConvergenceOnSolution());
			}
      		}
#endif /* MBDYN_X_CONVSOL */
	}
	
	
}




MatrixFreeSolver::MatrixFreeSolver(const Preconditioner::PrecondType PType, 
	const integer iPStep,
	doublereal ITol,
	integer MaxIt,
	doublereal etaMx) 
:pSM(NULL),
pPM(NULL),
pRes(NULL),
IterToll(ITol),
MaxLinIt(MaxIt),
Tau(defaultTau),
gamma(defaultGamma),
etaMax(etaMx),
PrecondIter(iPStep),
fBuildMat(true),
pPrevNLP(NULL)
{
	
	switch(PType) {
		
	case Preconditioner::FULLJACOBIAN:
			
		pPM  = new FullJacobianPr();				
		break;
	
	default:
		std::cerr << "Unknown Preconditioner type; Aborting " << std::endl;
		THROW(ErrGeneric()); 
	}
}


doublereal MatrixFreeSolver::MakeTest(const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("MatrixFreeSolver::MakeTest");
   
      
   	doublereal dRes = 0.;
	ASSERT(pSM != NULL);
	
#ifdef USE_MPI
#warning "NonlinSolver MakeTest parallel broken !! "	
#if 0
   	Dof CurrDof;
	SchurSolutionManager *pSSM;
	if ((pSSM = dynamic_cast<SchurSolutionManager*> (pSM)) != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
			/* else if ALGEBRAIC: non aggiunge nulla */
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

#ifdef __HACK_RES_TEST__				
		for (int iCnt = 0; iCnt < iMI; iCnt++) {
			DCount = pMI[iCnt];
			CurrDof = pDofs[DCount-1];
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				doublereal d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
#endif /* __HACK_RES_TEST__ */
			/* else if ALGEBRAIC: non aggiunge nulla */
		}
		
		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(dRes, dXPr);
		
	} else {
#endif		
#endif /* USE_MPI */
			ASSERT(Vec.iGetSize() == Size);
#ifdef __HACK_SCALE_RES__
			ASSERT(pScale != NULL);
			ASSERT(pScale->iGetSize == Size);
#endif /* __HACK_SCALE_RES__ */
 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = Vec.dGetCoef(iCntp1);
			doublereal d2 = d*d;

#ifdef __HACK_SCALE_RES__
			doublereal ds = pScale->dGetCoef(iCntp1);
			doublereal ds2 = ds*ds;
			d2 *= ds2;
#endif /* __HACK_SCALE_RES__ */

			dRes += d2;

#ifdef __HACK_RES_TEST__
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(iCntp1);
				d2 = d*d;

#ifdef __HACK_SCALE_RES__
				d2 *= ds2;         

#endif /* __HACK_SCALE_RES__ */

				dXPr += d2;
			}
			/* else if ALGEBRAIC: non aggiunge nulla */ 
#endif /* __HACK_RES_TEST__ */
		
		}
		
#ifdef USE_MPI
#if 0 		
	}
#endif
#endif /* USE_MPI */

#ifdef __HACK_RES_TEST__
	dRes /= (1.+dXPr);
#endif /* __HACK_RES_TEST__ */


   	return sqrt(dRes);
}



void BiCGStab::Solve(const NonlinearProblem* pNLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal Toll,
			const doublereal SolToll,
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
	std::cout << " New Step " <<std::endl; 			
#endif /* DEBUG_ITERATIVE */

	while (1) {

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

		dErr = MakeTest(*pRes);		

#ifdef DEBUG_ITERATIVE
		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (dErr < Toll) {
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
        	doublereal LocToll = eta * dErr;
        	rHat = *pr;

#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocToll " << LocToll << std::endl;
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
        	while ((resid > LocToll) && (It++ < MaxLinIt)) {
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

			if ((resid = s.Norm()) < LocToll) {
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
		etaBis = .5*Toll/dErr;
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
		if (SolToll > 0.) {
			dSolErr = MakeTest(dx);
        		if (dSolErr < dSolToll) {
				THROW(ConvergenceOnSolution());
			}
      		}
#endif /* MBDYN_X_CONVSOL */
	}
	
				
			
}

void Gmres::Solve(const NonlinearProblem* pNLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal Toll,
			const doublereal SolToll,
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
	doublereal norm1, norm2;
	VectorHandler* pr;
	MyVectorHandler s(MaxLinIt+1), cs(MaxLinIt+1), sn(MaxLinIt+1), w(Size);
	MyVectorHandler* v = new MyVectorHandler[MaxLinIt+1]; 
	MyVectorHandler vHat(Size); 
	MyVectorHandler dx(Size); 
	UpHessMatrix H(MaxLinIt+1);
	integer TotalIter = 0;

#ifdef DEBUG_ITERATIVE
	std::cout << " New Step " <<std::endl; 			
#endif /* DEBUG_ITERATIVE */

	while (1) {

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

		dErr = MakeTest(*pRes);		

#ifdef DEBUG_ITERATIVE
      		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (dErr < Toll) {
	 		return;
      		}
      		if (!isfinite(dErr)) {
			THROW(ErrSimulationDiverged());
		}
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
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
        	doublereal LocToll = eta * dErr;
 
#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocToll " << LocToll << std::endl;
#endif /* DEBUG_ITERATIVE */
	
		resid = dErr;
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

		int i = 0;
        	v[0].Resize(Size);
		v[0].ScalarMul(*pr, 1./resid);
		s.Reset(0.);
		cs.Reset(0.);
		sn.Reset(0.);
		s.fPutCoef(1, resid);
		while ((i < MaxLinIt)) {

#ifdef DEBUG_ITERATIVE
			std::cout << "v[i]:" << i << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< v[i].dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pPM->Precond(v[i], vHat, pSM); 

#ifdef DEBUG_ITERATIVE
			std::cout << "vHat:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< vHat.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pNLP->EvalProd(Tau, *pr, vHat, w);

#ifdef DEBUG_ITERATIVE
			std::cout << "w:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< w.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			norm1 = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm1: " << norm1 << std::endl; 
#endif /* DEBUG_ITERATIVE */

			for (int k = 0; k <= i; k++) {
        			H(k, i) = w.InnerProd(v[k]);

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k,i): " << k << " " << i << " "<< H(k,i) << std::endl; 
#endif /* DEBUG_ITERATIVE */

				w.ScalarAddMul(v[k], -H(k, i));
      			}
			H(i+1, i) = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "w.Norm(): " << w.Norm() << std::endl;
			std::cerr << "H(i+1, i): " << i+1 << " " << i << " " << H(i+1, i) << std::endl;
#endif /* DEBUG_ITERATIVE */

			norm2 = H(i+1, i); 

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm2: " << norm2 << std::endl;
#endif /* DEBUG_ITERATIVE */

			v[i+1].Resize(Size);
			
			/*  Reorthogonalize? */
			if  (.001*norm2/norm1 < DBL_EPSILON) {
    				for (int k = 0;  k <= i; k++) {       
					doublereal hr = v[k].InnerProd(w);
        				H(k,i) = H(k,i)+hr;
        				w.ScalarAddMul(v[k], -hr);

#ifdef DEBUG_ITERATIVE
					std::cerr << "Reorthogonalize: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				}
				H(i+1, i) = w.Norm();
			}
			if (fabs(H(i+1, i)) > DBL_EPSILON) { 
				v[i+1].ScalarMul(w, 1./H(i+1, i));

#ifdef DEBUG_ITERATIVE
				std::cout << "v[i+1]:" << i+1 << std::endl;
#endif /* DEBUG_ITERATIVE */

	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< v[i+1].dGetCoef(iTmpCnt) << std::endl;
				}
			} else {
				/* happy breakdown !! */

#ifdef DEBUG_ITERATIVE
				std::cerr << "happy breakdown: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				v[i+1].Reset(0.);
			} 
			for (int k = 0; k < i; k++) {
        			ApplyPlaneRotation(H(k,i), H(k+1,i), cs.dGetCoef(k+1), sn.dGetCoef(k+1));

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k, i): " << k << " " << i << " " << H(k, i) << std::endl;
				std::cerr << "H(k+1, i): " << k+1 << " " << i << " " << H(k+1, i) << std::endl;
#endif /* DEBUG_ITERATIVE */

			}
							
			GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i+1), sn(i+1));

#ifdef DEBUG_ITERATIVE
			std::cerr << "cs(i): " << cs.dGetCoef(i+1) << std::endl;
			std::cerr << "sn(i): " << sn.dGetCoef(i+1) << std::endl;
#endif /* DEBUG_ITERATIVE */

      			ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i+1), sn(i+1));
      			ApplyPlaneRotation(s(i+1), s(i+2), cs(i+1), sn(i+1));
				
			if ((resid = fabs(s.dGetCoef(i+2))) < LocToll) {

#ifdef DEBUG_ITERATIVE
				std::cerr << "resid 1: " << resid  << std::endl;
#endif /* DEBUG_ITERATIVE */

				Backsolve(dx, i, H, s, v);
				pPM->Precond(dx, dx, pSM); 

#ifdef DEBUG_ITERATIVE
				std::cout << "dx:" << std::endl;
	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< dx.dGetCoef(iTmpCnt) << std::endl;
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
			Backsolve(dx, MaxLinIt, H, s, v);
			pPM->Precond(dx, dx, pSM);				
		        std::cerr << "Iterative inner solver didn't converge."
				<< " Continuing..." << std::endl;
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
		etaBis = .5*Toll/dErr;
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
		if (SolToll > 0.) {
			dSolErr = MakeTest(dx);
        		if (dSolErr < dSolToll) {
				THROW(ConvergenceOnSolution());
			}
      		}
#endif /* MBDYN_X_CONVSOL */
	}
	delete [] v;		
			
}

