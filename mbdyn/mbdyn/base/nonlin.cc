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
#include <ac/float.h>

doublereal NewtonRaphsonSolver::MakeTest(const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("NewtonRaphsonSolver::MakeTest");
   
//   	Dof CurrDof;
      
   	doublereal dRes = 0.;
	ASSERT(pSM != NULL);
	
#ifdef USE_MPI
#warning "NonlinSolver MakeTest parallel broken !! "	
#if 0
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
		
/*		for (int iCnt = 0; iCnt < iMI; iCnt++) {
			DCount = pMI[iCnt];
			CurrDof = pDofs[DCount-1];
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				doublereal d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
*/
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

/*			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(iCntp1);
				d2 = d*d;

			#ifdef __HACK_SCALE_RES__
				d2 *= ds2;         
*/
			//#endif /* __HACK_SCALE_RES__ */
/*
				dXPr += d2;
			}
*/
			/* else if ALGEBRAIC: non aggiunge nulla */ 
		}
#ifdef USE_MPI
#if 0 		
	}
#endif
#endif /* USE_MPI */

//   	dRes /= (1.+dXPr);
   	if (!isfinite(dRes)) {      
      		std::cerr << "The simulation diverged; aborting ..." 
			<< std::endl;       
      		THROW(ErrSimulationDiverged());
   	}

   	return sqrt(dRes);
}







NewtonRaphsonSolver::
	NewtonRaphsonSolver(const flag fTNR, 
			const integer IterBfAss)
:pSM(NULL),
pRes(NULL),
pSol(NULL),
pJac(NULL),
fTrueNewtonRaphson(fTNR),
IterationBeforeAssembly(IterBfAss) { };


void NewtonRaphsonSolver::
	Solve(const NonlinearProblem* pNLP,
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




BiCGMatrixFreeSolver::BiCGMatrixFreeSolver(const Preconditioner::PrecondType PType, 
	const integer iPStep,
	doublereal ITol,
	integer MaxIt) 
:pSM(NULL),
pPM(NULL),
pRes(NULL),
IterToll(ITol),
MaxLinIt(MaxIt),
Tau(defaultTau),
gamma(defaultGamma),
etaMax(defaultEtaMax),
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

void BiCGMatrixFreeSolver::Solve(const NonlinearProblem* pNLP,
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
	
	pSM  = pSolMan;
        pRes = pSM->pResHdl();
	Size = pRes->iGetSize();

	doublereal eta = 0.9;
	doublereal rateo = 0.;
	doublereal Fnorm = 1.;
        doublereal resid;
        doublereal rho_1; 
	doublereal rho_2; 
	doublereal alpha;
	doublereal beta;
	doublereal omega;
	VectorHandler* pr;
	MyVectorHandler rHat(Size), p(Size), s(Size), t(Size), v(Size);
	MyVectorHandler dx(Size); 
	
	if (fBuildMat) {
		pSM->MatrInit(0.);
		pNLP->Jacobian(pSM->pMatHdl());
		fBuildMat = false;
	}
		
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
          	rateo = dErr*dErr/Fnorm;
		Fnorm = dErr*dErr;
		
      		iIterCnt++;

		/* inner iteration to solve the linear system */	
	
		/* BiCGSTAB Iterative solver */
		DEBUGCOUT("Using BiCGStab iterative solver" << std::endl);
		rho_1 = 1.; 
		alpha = 1.;
		omega = 1.;

        	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */
        	/* N.B. *pRes = -F(0) */ 
		
		pr = pRes;
        	doublereal LocToll = eta * dErr;
        	rHat = *pr;

		rho_2 = dErr*dErr;   /*rhat.InnerProd(r); */
		resid = dErr;
		v.Reset(0.);
		t.Reset(0.);
		p.Reset(0.);
		dx.Reset(0.);
		
		int It = 0;
        	while ((resid > LocToll) && (It++ < MaxLinIt)) {
                	if (fabs(omega) < DBL_EPSILON) {
                	//if (omega == 0.) {
                        	std::cout << "Bi-CGStab Iterative Solver breakdown" 
					<<  " omega = 0 " << std::endl;
				THROW(ErrGeneric());
			}
			if (It == 1) {
				p = *pr;
			} else {
				beta = (rho_2/rho_1) * (alpha/omega);
				p.ScalarAddMul(*pr, p.ScalarAddMul(v, -omega), beta);
			}

			/* preconditioning */
			pPM->Precond(p, pSM);

			pNLP->EvalProd(Tau, rHat, p, v);
			alpha = rHat.InnerProd(v);
			if (fabs(alpha) < DBL_EPSILON) {
			//if ( alpha == 0)  {
				std::cout << "Bi-CGSTAB Iterative Solver breakdown"
					<< " r_0 * v = 0 " << std::endl;
				THROW(ErrGeneric());
			}
			alpha = rho_2 / alpha;
			s.ScalarAddMul(*pr, v, -alpha);
			pPM->Precond(s, pSM);
			pNLP->EvalProd(Tau, rHat, s, t);
			omega = sqrt(t.InnerProd(t));
			if (fabs(omega) < DBL_EPSILON) {
			//if ( omega == 0)  {
				std::cout << "Bi-CGSTAB Iterative Solver breakdown"
					<< "  t(Vec) = 0 " << std::endl;
				THROW(ErrGeneric());
			}
			omega = t.InnerProd(s) / omega;
			rho_1 = rho_2;
			rho_2 = rHat.InnerProd(t);
			rho_2 = - omega * rho_2;
			dx.ScalarAddMul(p, alpha);
			dx.ScalarAddMul(s, omega);
			pr->ScalarAddMul(s, t, -omega);
			resid = pr->Norm();	
			
			if ( It == MaxLinIt) {
                        	std::cerr << "Iterative inner solver didn't converge."
					<< " Continuing..." << std::endl;
			}
		}
		
		/* se ha impiegato troppi passi riassembla lo jacobiano */
		if (It > PrecondIter) {
			fBuildMat = true;
		}
		/* calcolo il nuovo eta */
		
		doublereal etaNew = gamma * rateo;
		doublereal etaBis;
		if ((etaBis = gamma*eta*eta) > .1) {
			etaNew = (etaNew > etaBis) ? etaNew : etaBis;
		}
		eta = (etaNew < etaMax) ? etaNew : etaMax;
		
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

doublereal BiCGMatrixFreeSolver::MakeTest(const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("BiCGMatrixFreeSolver::MakeTest");
   
//   	Dof CurrDof;
      
   	doublereal dRes = 0.;
	ASSERT(pSM != NULL);
	
#ifdef USE_MPI
#warning "NonlinSolver MakeTest parallel broken !! "	
#if 0
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
		
/*		for (int iCnt = 0; iCnt < iMI; iCnt++) {
			DCount = pMI[iCnt];
			CurrDof = pDofs[DCount-1];
			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				doublereal d = XP.dGetCoef(DCount);
				dXPr += d*d;
			}
*/
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

/*			if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
				d = XP.dGetCoef(iCntp1);
				d2 = d*d;

			#ifdef __HACK_SCALE_RES__
				d2 *= ds2;         
*/
			//#endif /* __HACK_SCALE_RES__ */
/*
				dXPr += d2;
			}
*/
			/* else if ALGEBRAIC: non aggiunge nulla */ 
		}
#ifdef USE_MPI
#if 0 		
	}
#endif
#endif /* USE_MPI */

//   	dRes /= (1.+dXPr);
   	if (!isfinite(dRes)) {      
      		std::cerr << "The simulation diverged; aborting ..." 
			<< std::endl;       
      		THROW(ErrSimulationDiverged());
   	}

   	return sqrt(dRes);
}
