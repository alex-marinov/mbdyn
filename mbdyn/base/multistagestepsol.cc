/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2021
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
 * Author: Huimin Zhang <huimin.zhang@polimi.it> 2021
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "schurdataman.h"
#include "external.h"
#include "ls.h"
#include "solver.h"
#include "invsolver.h"
#include "stepsol.h"
#include "multistagestepsol.h"
#include "stepsol.hc"

/* Stage2Integrator - begin */
Stage2Integrator::Stage2Integrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, 2, bmod_res_test),
pXPrev(NULL),
pXPrimePrev(NULL),
pXInte(NULL),
pXPrimeInte(NULL)
{	
	NO_OP;
}

Stage2Integrator::~Stage2Integrator(void)
{
	NO_OP;
}


void Stage2Integrator::PredictDofforStage1(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXnm1 = pXPrev->operator()(DCount);
		doublereal dXPnm1 = pXPrimePrev->operator()(DCount);
		doublereal dXPn = dPredDerforStage1(dXnm1, dXPnm1);
		doublereal dXn = dPredStateforStage1(dXnm1, dXPn, dXPnm1);

		pXPrimeCurr->PutCoef(DCount, dXPn);
		pXCurr->PutCoef(DCount, dXn);

	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXnm1 = pXPrev->operator()(DCount);
		doublereal dXInm1 =
			pXPrimePrev->operator()(DCount);

		doublereal dXn = dPredDerAlgforStage1(dXInm1, dXnm1);
		doublereal dXIn = dPredStateAlgforStage1(dXInm1, dXn, dXnm1);

		pXCurr->PutCoef(DCount, dXn);
		pXPrimeCurr->PutCoef(DCount, dXIn);

	} else {
		silent_cerr("Stage2Integrator::"
			"PredictDofforStage1(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	//std::cout<<"PredictforStage1, DCount = "<< DCount <<"; Order = "<< Order << "; pXPrev = "<<pXPrev->operator()(DCount)<<"; pXCurr = "<<pXCurr->operator()(DCount)<<std::endl;
	//std::cout<<"PredictforStage1, DCount = "<< DCount <<"; Order = "<< Order << "; pXPrimePrev = "<<pXPrimePrev->operator()(DCount)<<"; pXPrimeCurr = "<<pXPrimeCurr->operator()(DCount)<<std::endl;
	

}

void
Stage2Integrator::PredictforStage1(void)
{
   	DEBUGCOUTFNAME("Stage2Integrator::PredictforStage1");
   	ASSERT(pDM != NULL);
	UpdateLoop(this, &Stage2Integrator::PredictDofforStage1);
}

void Stage2Integrator::PredictDofforStage2(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXnm1 = pXInte->operator()(DCount);
		doublereal dXnm2 = pXPrev->operator()(DCount);
		doublereal dXPnm1 = pXPrimeInte->operator()(DCount);
		doublereal dXPnm2 = pXPrimePrev->operator()(DCount);
		doublereal dXPn = dPredDerforStage2(dXnm1, dXnm2, dXPnm1, dXPnm2);
		doublereal dXn = dPredStateforStage2(dXnm1, dXnm2, dXPn, dXPnm1, dXPnm2);

		pXPrimeCurr->PutCoef(DCount, dXPn);
		pXCurr->PutCoef(DCount, dXn);
	

	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXnm1 = pXInte->operator()(DCount);
        doublereal dXnm2 = pXPrev->operator()(DCount);
		doublereal dXInm1 =
			pXPrimeInte->operator()(DCount);
		doublereal dXInm2 =
			pXPrimePrev->operator()(DCount);

		doublereal dXn = dPredDerAlgforStage2(dXInm1, dXInm2, dXnm1, dXnm2);
		doublereal dXIn = dPredStateAlgforStage2(dXInm1, dXInm2, dXn, dXnm1, dXnm2);

		pXCurr->PutCoef(DCount, dXn);
		pXPrimeCurr->PutCoef(DCount, dXIn);

	} else {
		silent_cerr("Stage2Integrator::"
			"PredictDofforStage2(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	//std::cout<<"PredictforStage2, DCount = "<< DCount <<"; Order = "<< Order << "; pXPrev = "<<pXPrev->operator()(DCount) << "; pXInte = "<<pXInte->operator()(DCount)<<"; pXCurr = "<<pXCurr->operator()(DCount)<<std::endl;
	//std::cout<<"PredictforStage2, DCount = "<< DCount <<"; Order = "<< Order << "; pXPrimePrev = "<<pXPrimePrev->operator()(DCount) << "; pXPrimeInte = "<<pXPrimeInte->operator()(DCount)<<"; pXPrimeCurr = "<<pXPrimeCurr->operator()(DCount)<<std::endl;
	
}

void
Stage2Integrator::PredictforStage2(void)
{
   	DEBUGCOUTFNAME("Stage2Integrator::PredictforStage2");
   	ASSERT(pDM != NULL);
	UpdateLoop(this, &Stage2Integrator::PredictDofforStage2);
}

doublereal
Stage2Integrator::Advance(Solver* pS,
		const doublereal TStep,
		const doublereal dAph, const StepChange StType,
		std::deque<MyVectorHandler*>& qX,
	 	std::deque<MyVectorHandler*>& qXPrime,
		MyVectorHandler*const pX,
		MyVectorHandler*const pXPrime,
		integer& EffIter,
		doublereal& Err,
		doublereal& SolErr)
{
	ASSERT(pDM != NULL);
	pXCurr  = pX;
	pXPrev  = qX[0];
	 
	pXPrimeCurr  = pXPrime;
	pXPrimePrev  = qXPrime[0];

    /*First-stage*/
    SetCoef(TStep, dAph, StType);
	PredictforStage1();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();

#ifdef DEBUG
	integer iNumDofs = pDM->iGetNumDofs();
	if (outputPred) {
		std::cout << "After prediction, time=" << pDM->dGetTime() << std::endl;
		std::cout << "Dof:      |    XCurr  ";
		for (unsigned idx = 0; idx < qX.size(); idx++) {
			std::cout << "|  XPrev[" << idx << "] ";
		}
		std::cout << "|   XPrime  ";
		for (unsigned idx = 0; idx < qXPrime.size(); idx++) {
			std::cout << "| XPPrev[" << idx << "] ";
		}
		std::cout << "|" << std::endl;
		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
    			std::cout << std::setw(8) << iTmpCnt << ": ";
			std::cout << std::setw(12) << pX->operator()(iTmpCnt);
			for (unsigned int ivec = 0; ivec < qX.size(); ivec++) {
				std::cout << std::setw(12)
					<< (qX[ivec])->operator()(iTmpCnt);
			}
			std::cout << std::setw(12) << pXPrime->operator()(iTmpCnt);
			for (unsigned int ivec = 0; ivec < qXPrime.size(); ivec++) {
				std::cout << std::setw(12)
					<< (qXPrime[ivec])->operator()(iTmpCnt);
			}
			std::cout << " " << pDM->DataManager::GetDofDescription(iTmpCnt) << std::endl;
 		}
	}
#endif /* DEBUG */

	Err = 0.;
	pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol,
    			EffIter, Err, dSolTol, SolErr);
	integer EffIter1 = EffIter;
	doublereal Err1 = Err;
	/* if it gets here, it surely converged */
	pDM->AfterConvergence();

    /*Second-stage*/
	pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *qX[0], *qXPrime[0]);
	qX.push_front(qX.back());
	qX.pop_back();
	qXPrime.push_front(qXPrime.back());
	qXPrime.pop_back();

	/* copy from pX, pXPrime to qx[0], qxPrime[0] */
	for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
		qX[0]->PutCoef(i, pXCurr->operator()(i));
		qXPrime[0]->PutCoef(i, pXPrimeCurr->operator()(i));
	}

	pXInte = qX[0];
	pXPrev = qX[1];

	pXPrimeInte = qXPrime[0];
	pXPrimePrev = qXPrime[1]; 

    SetCoef(TStep, dAph, StType);
	PredictforStage2();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();

#ifdef DEBUG
	integer iNumDofs = pDM->iGetNumDofs();
	if (outputPred) {
		std::cout << "After prediction, time=" << pDM->dGetTime() << std::endl;
		std::cout << "Dof:      |    XCurr  ";
		for (unsigned idx = 0; idx < qX.size(); idx++) {
			std::cout << "|  XPrev[" << idx << "] ";
		}
		std::cout << "|   XPrime  ";
		for (unsigned idx = 0; idx < qXPrime.size(); idx++) {
			std::cout << "| XPPrev[" << idx << "] ";
		}
		std::cout << "|" << std::endl;
		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
    			std::cout << std::setw(8) << iTmpCnt << ": ";
			std::cout << std::setw(12) << pX->operator()(iTmpCnt);
			for (unsigned int ivec = 0; ivec < qX.size(); ivec++) {
				std::cout << std::setw(12)
					<< (qX[ivec])->operator()(iTmpCnt);
			}
			std::cout << std::setw(12) << pXPrime->operator()(iTmpCnt);
			for (unsigned int ivec = 0; ivec < qXPrime.size(); ivec++) {
				std::cout << std::setw(12)
					<< (qXPrime[ivec])->operator()(iTmpCnt);
			}
			std::cout << " " << pDM->DataManager::GetDofDescription(iTmpCnt) << std::endl;
 		}
	}
#endif /* DEBUG */

	Err = 0.;
	pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol,
    			EffIter, Err, dSolTol, SolErr);
	EffIter1+=EffIter;
	EffIter=EffIter1;
	Err1+=Err;
	Err=Err1;
	/* if it gets here, it surely converged */
	pDM->AfterConvergence();


	return Err;
}

/* Stage2Integrator - end */

/* TunableBatheSolver - begin */

TunableBatheSolver::TunableBatheSolver(const doublereal Tl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test)
: Stage2Integrator(iMaxIt, Tl, dSolTl, bmod_res_test),
Rho(pRho), AlgebraicRho(pAlgRho), iStage(1)
{
	ASSERT(pRho != NULL);
	ASSERT(pAlgRho != NULL);
}

TunableBatheSolver::~TunableBatheSolver(void)
{
	NO_OP;
}

void
TunableBatheSolver::SetDriveHandler(const DriveHandler* pDH)
{
	Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
TunableBatheSolver::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange sType)
{	
	if (iStage == 1)
	{
		SetCoefforStage1(dT, dAlpha, sType);
		//std::cout<<"Time = "<<pDM->dGetTime()<<std::endl;
		//std::cout<<"iStep = "<<pDM->pGetDrvHdl()->iGetStep()<<std::endl;
		//std::cout<<"iStage = "<<iStage<<std::endl;
		iStage = iStage + 1;
	} else
	if (iStage == 2)
	{
		SetCoefforStage2(dT, dAlpha, sType);
		//std::cout<<"Time = "<<pDM->dGetTime()<<std::endl;
		//std::cout<<"iStep = "<<pDM->pGetDrvHdl()->iGetStep()<<std::endl;
		//std::cout<<"iStage = "<<iStage<<std::endl;
		iStage = iStage - 1;
	}else {
		silent_cerr("TunableBatheSolver::"
			"SetCoef(): "
			"unknown stage number "
			<< iStage << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
TunableBatheSolver::SetCoefforStage1(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	doublereal dRho = Rho.dGet();
	doublereal dAlgebraicRho = AlgebraicRho.dGet();

	if (dRho == 1.)
	{
		gamma = 1./2.;
	} else 
	{
		gamma = (2. - sqrt(2.+ 2.* dRho))/(1.- dRho);
	}

	ASSERT(pDM != NULL);
	pDM->SetTime(pDM->dGetTime()-(1.-gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
	a[0][DIFFERENTIAL] = 1.;
	a[1][DIFFERENTIAL] = 0.;
	b[0][DIFFERENTIAL] = gamma*dT/2.;
	b[1][DIFFERENTIAL] = gamma*dT/2.;
	b[2][DIFFERENTIAL] = 0.;


	DEBUGCOUT("PredictforStage1()" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << dRho << std::endl 
			<< "a1    = " << a[0][DIFFERENTIAL] << std::endl
			<< "a2    = " << a[1][DIFFERENTIAL] << std::endl
			<< "b0    = " << b[0][DIFFERENTIAL] << std::endl
			<< "b1    = " << b[1][DIFFERENTIAL] << std::endl
			<< "b2    = " << b[2][DIFFERENTIAL] << std::endl);

	a[0][ALGEBRAIC] = a[0][DIFFERENTIAL];
	a[1][ALGEBRAIC] = a[1][DIFFERENTIAL];
	b[0][ALGEBRAIC] = b[0][DIFFERENTIAL];
	b[1][ALGEBRAIC] = b[1][DIFFERENTIAL];
	b[2][ALGEBRAIC] = b[2][DIFFERENTIAL];


	DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "Asymptotic rho =" << dAlgebraicRho << std::endl 
			<< "a1    = " << a[0][ALGEBRAIC] << std::endl
			<< "a2    = " << a[1][ALGEBRAIC] << std::endl
			<< "b0    = " << b[0][ALGEBRAIC] << std::endl
			<< "b1    = " << b[1][ALGEBRAIC] << std::endl
			<< "b2    = " << b[2][ALGEBRAIC] << std::endl);


	db0Differential = b[0][DIFFERENTIAL];
	db0Algebraic = b[0][ALGEBRAIC];

	//std::cout<<"CoefforStage1= "<<a[0][DIFFERENTIAL]<<", "<<a[1][DIFFERENTIAL]<<", "<<b[0][DIFFERENTIAL]<<", "<<b[1][DIFFERENTIAL]<<", "<<b[2][DIFFERENTIAL]<<std::endl;
}

void
TunableBatheSolver::SetCoefforStage2(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	doublereal dRho = Rho.dGet();
	doublereal dAlgebraicRho = AlgebraicRho.dGet();

	ASSERT(pDM != NULL);
	pDM->SetTime(pDM->dGetTime()+(1.-gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
	/*mp[0] = -12.*dAlpha*(1.+dAlpha)/dT;
	mp[1] = -mp[0];
	np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
	np[1] = dAlpha*(2.+3.*dAlpha);*/

	a[0][DIFFERENTIAL] = 0.;
	a[1][DIFFERENTIAL] = 1.;
	b[0][DIFFERENTIAL] = (1. - gamma)*dT/(gamma*dRho- gamma + 2.);
	b[1][DIFFERENTIAL] = (1. + dRho)*dT/(2.*(gamma*dRho- gamma + 2.));
	b[2][DIFFERENTIAL] = (2.*gamma*dRho - dRho + 1.)*dT/(2.*(gamma*dRho- gamma + 2.));


	DEBUGCOUT("PredictforStage2()" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << dRho << std::endl 
			<< "a1    = " << a[0][DIFFERENTIAL] << std::endl
			<< "a2    = " << a[1][DIFFERENTIAL] << std::endl
			<< "b0    = " << b[0][DIFFERENTIAL] << std::endl
			<< "b1    = " << b[1][DIFFERENTIAL] << std::endl
			<< "b2    = " << b[2][DIFFERENTIAL] << std::endl);

	a[0][ALGEBRAIC] = a[0][DIFFERENTIAL];
	a[1][ALGEBRAIC] = a[1][DIFFERENTIAL];
	b[0][ALGEBRAIC] = b[0][DIFFERENTIAL];
	b[1][ALGEBRAIC] = b[1][DIFFERENTIAL];
	b[2][ALGEBRAIC] = b[2][DIFFERENTIAL];


	DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "Asymptotic rho =" << dAlgebraicRho << std::endl 
			<< "a1    = " << a[0][ALGEBRAIC] << std::endl
			<< "a2    = " << a[1][ALGEBRAIC] << std::endl
			<< "b0    = " << b[0][ALGEBRAIC] << std::endl
			<< "b1    = " << b[1][ALGEBRAIC] << std::endl
			<< "b2    = " << b[2][ALGEBRAIC] << std::endl);


	db0Differential = b[0][DIFFERENTIAL];
	db0Algebraic = b[0][ALGEBRAIC];
	//std::cout<<"PredictCoefforStage2= "<<mp[0]<<", "<<mp[1]<<", "<<np[0]<<", "<<np[1]<<std::endl;
	//std::cout<<"CoefforStage2= "<<a[0][DIFFERENTIAL]<<", "<<a[1][DIFFERENTIAL]<<", "<<b[0][DIFFERENTIAL]<<", "<<b[1][DIFFERENTIAL]<<", "<<b[2][DIFFERENTIAL]<<std::endl;

}

doublereal
TunableBatheSolver::dPredDerforStage1(const doublereal& dXm1,
		const doublereal& dXPm1) const
{
	return dXPm1;
}

doublereal
TunableBatheSolver::dPredStateforStage1(const doublereal& dXm1,
		const doublereal& dXP,
		const doublereal& dXPm1) const
{
	return a[0][DIFFERENTIAL]*dXm1+b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1;
}

doublereal
TunableBatheSolver::dPredDerAlgforStage1(const doublereal& dXm1,
		const doublereal& dXPm1) const
{
	return dXPm1;
}

doublereal
TunableBatheSolver::dPredStateAlgforStage1(const doublereal& dXm1,
		const doublereal& dXP,
		const doublereal& dXPm1) const
{
	return a[0][ALGEBRAIC]*dXm1+b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1;
}

doublereal
TunableBatheSolver::dPredDerforStage2(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	/*return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;*/
	return dXPm1;
}

doublereal
TunableBatheSolver::dPredStateforStage2(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2+
	b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1+b[2][DIFFERENTIAL]*dXPm2;
}

doublereal
TunableBatheSolver::dPredDerAlgforStage2(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	/*return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;*/
	return dXPm1;
}

doublereal
TunableBatheSolver::dPredStateAlgforStage2(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	/*return  -a[1][ALGEBRAIC]*dXm1+
	b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2;*/
	return a[0][ALGEBRAIC]*dXm1+a[1][ALGEBRAIC]*dXm2+
	b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2;
}

/* TunableBatheSolver - end */
