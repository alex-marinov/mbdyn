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
  * Copyright (C) 2008
  * Marco Morandini	<morandini@aero.polimi.it>
  *
  * third order integrator; brain-damaged code, mainly due to some
  *                         brain-damaged design decision in mbdyn.
  *                         This will have to change, but will require
  *                         a substantial effort.
  */


#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include "thirdorderstepsol.h"
#include "schurdataman.h" 

#include "stepsol.hc"

ThirdOrderIntegrator::ThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const bool bmod_res_test)
: ImplicitStepIntegrator(iMaxIt, dT, dSolutionTol, 1, 2, bmod_res_test),
pXPrev(0),
pXPrimePrev(0),
bAdvanceCalledFirstTime(true)
{
	pJacxi_xp = &Jacxi_xp;
	pJacxi_x  = &Jacxi_x;
	pJac_xp   = &Jac_xp;
	pJac_x    = &Jac_x;
};

ThirdOrderIntegrator::~ThirdOrderIntegrator(){
	NO_OP;
};

doublereal
ThirdOrderIntegrator::Advance(Solver* pS, 
		const doublereal TStep, 
		const doublereal dAlph, 
		const StepChange StType,
		std::deque<MyVectorHandler*>& qX,
		std::deque<MyVectorHandler*>& qXPrime,
		MyVectorHandler* const pX,
		MyVectorHandler* const pXPrime,
		integer& EffIter,
		doublereal& Err,
		doublereal& SolErr)
{
	if (bAdvanceCalledFirstTime) {
		integer n = pDM->iGetNumDofs();
		Restmp.Resize(n);
		EqIsAlgebraic.resize(n);
		EqIsDifferential.resize(n);
	   	Dof CurrDof;
		DofIterator.bGetFirst(CurrDof);
		for (int iCntm1 = 0; iCntm1 < n;
			iCntm1++, DofIterator.bGetNext(CurrDof)) {
			EqIsAlgebraic[iCntm1] = (
				CurrDof.EqOrder==DofOrder::ALGEBRAIC);
			EqIsDifferential[iCntm1] = (!EqIsAlgebraic[iCntm1]);
		}
		DofIterator.bGetFirst(CurrDof);
		Jacxi_xp.Resize(n, n);
		Jacxi_x.Resize(n, n);
		Jac_xp.Resize(n, n);
		Jac_x.Resize(n, n);
		bAdvanceCalledFirstTime = false;
	}
	pXCurr  = pX;
	pXPrev  = qX[0];

	pXPrimeCurr  = pXPrime;
	pXPrimePrev  = qXPrime[0];
	
	SetCoef(TStep, dAlph, StType);	

	Predict();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	
	Err = 0.;        
	pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol, 
			EffIter, Err, dSolTol, SolErr);
	
	return Err;
};


void ThirdOrderIntegrator::PredictDof_for_AfterPredict(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol) const {
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXnm1 = pXPrev->operator()(DCount);
 		doublereal dXPnm1 = 
			pXPrimePrev->operator()(DCount);
			
 		doublereal dXn = dXnm1+dXPnm1*dT;
		doublereal dXPn = dXPnm1;
	
 		pXPrimeCurr->PutCoef(DCount, dXPn);
 		pXCurr->PutCoef(DCount, dXn);
		
	} else if (Order == DofOrder::ALGEBRAIC) {
 		doublereal dXnm1 = pXPrev->operator()(DCount);
 		//doublereal dXInm1 = pXPrimePrev->operator()(DCount);
 		doublereal dXn = dXnm1;
		doublereal dXIn = dXnm1*dT;
		
 		pXCurr->PutCoef(DCount, dXn);
 		pXPrimeCurr->PutCoef(DCount, dXIn);

	} else {
 		silent_cerr("ThirdOrderIntegrator::PredictDof_for_AfterPredict:"
			<< "unknown order for dof " 
			<< DCount<< std::endl);
 		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}
void ThirdOrderIntegrator::RealPredictDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol) const {
	integer iNumDofs = pDM->iGetNumDofs();
	//simple copy of predicted state
	pXPrimeCurr->PutCoef(DCount+iNumDofs,
		pXPrimeCurr->operator()(DCount));
 	pXCurr->PutCoef(DCount+iNumDofs, pXCurr->operator()(DCount));
	if (Order == DofOrder::DIFFERENTIAL) {
		//doublereal dXPnm1 = 
		//	pXPrimePrev->operator()(DCount);
		
		/* tempo theta*/
// 		doublereal dXn = dXPnm1*(theta-1.)*dT;
// 		pXCurr->IncCoef(DCount+iNumDofs, dXn);
// 		pXCurr->PutCoef(DCount+iNumDofs,
// 			m0*pXCurr->operator()(DCount)+m1*pXPrev->operator()(DCount)
// 			+dT*(n0*pXPrimeCurr->operator()(DCount)+
// 				n1*pXPrimePrev->operator()(DCount)));
		pXCurr->IncCoef(DCount+iNumDofs,
			pXPrimePrev->operator()(DCount)*theta*dT);
      	} else if (Order == DofOrder::ALGEBRAIC) {
		//doublereal dXnm1 = pXPrev->operator()(DCount);
		
		/* tempo theta*/
// 		doublereal dXIn = dXnm1*(theta-1.)*dT;
// 		pXPrimeCurr->IncCoef(DCount+iNumDofs, dXIn);
// 		pXPrimeCurr->PutCoef(DCount+iNumDofs,
// 			m0*pXPrimeCurr->operator()(DCount)+m1*pXPrimePrev->operator()(DCount)
// 			+dT*(n0*pXCurr->operator()(DCount)+
// 				n1*pXPrev->operator()(DCount)));
// 		pXCurr->PutCoef(DCount+iNumDofs,
// 			pXPrev->operator()(DCount));
		pXPrimeCurr->IncCoef(DCount+iNumDofs,
			pXPrev->operator()(DCount)*theta*dT);
	} else {
		silent_cerr("ThirdOrderIntegrator::RealPredictDof: "
			<< "unknown order for dof " 
			<< DCount<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
ThirdOrderIntegrator::Predict(void)
{
   	DEBUGCOUTFNAME("ThirdOrderIntegrator::Predict");
   	ASSERT(pDM != NULL);
   	Dof CurrDof;

#ifdef USE_SCHUR
	SchurDataManager* pSDM;
	if ((pSDM = dynamic_cast<SchurDataManager*> (pDM)) != 0) {
		silent_cerr("Warning: ThirdOrderIntegrator currently is "
			<< "untested with the parallel solver" << std::endl);
	}
#endif // USE_SCHUR

	DofIterator.bGetFirst(CurrDof);
	
   	/* 
	 * Linear combination of previous step state and derivative 
	 * etc....
	*/
	/*
	 * Predict for AfterPredict
	*/
	UpdateLoop(this, &ThirdOrderIntegrator::PredictDof_for_AfterPredict);
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();

	/*
	 * Vero Predict
	 */
	UpdateLoop(this, &ThirdOrderIntegrator::RealPredictDof);
	return;
};

void ThirdOrderIntegrator::Residual(VectorHandler* pRes) const
{
   	DEBUGCOUTFNAME("ThirdOrderIntegrator::Residual");
	ASSERT(pDM != NULL);
	
	integer iNumDofs = pDM->iGetNumDofs();

	MyVectorHandler state, stateder, res;
	
	/* theta*dT */
	state.Attach(iNumDofs, pXCurr->pdGetVec()+iNumDofs);
	stateder.Attach(iNumDofs, pXPrimeCurr->pdGetVec()+iNumDofs);
	res.Attach(iNumDofs, pRes->pdGetVec()+iNumDofs);
	pDM->SetTime(pDM->dGetTime() + theta*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
	pDM->AssRes(res, 1.);
	
	/* dT */
	pDM->SetTime(pDM->dGetTime() - theta*dT);
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	pDM->Update();
	pDM->AssRes(*pRes, 1.);

	return;
};

void ThirdOrderIntegrator::Jacobian(MatrixHandler* pJac) const
{
   	DEBUGCOUTFNAME("ThirdOrderIntegrator::Jacobian");
	ASSERT(pDM != NULL);
	
	integer iNumDofs = pDM->iGetNumDofs();

	MyVectorHandler state, stateder;
	pJac->Reset();
#if 0
	pJacxi_x->Reset();
	pJacxi_xp->Reset();
	pJac_x->Reset();
	pJac_xp->Reset();
#endif
	/* theta*dT */
	state.Attach(iNumDofs, pXCurr->pdGetVec()+iNumDofs);
	stateder.Attach(iNumDofs, pXPrimeCurr->pdGetVec()+iNumDofs);
	pDM->SetTime(pDM->dGetTime() + theta*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
#warning This trick does not work with Coulomb friction (and other elements too)
	/*
	 The reason is that 
	 	a) Colomb friction residual (and other elements) can throw
		   to signal the need of a new jacobian
		b) Have internal states that depend on the order of residual
		   calls
	 However, this is needed for all the elements (almost all)
	 that compute something during residual and use it 
	 later during Jacobain computation
	*/
	pDM->AssRes(Restmp, 1.);
	pDM->AssJac(*pJacxi_x, 1.);
	pDM->AssJac(*pJacxi_xp, 0.);
	
	/* dT */
	pDM->SetTime(pDM->dGetTime() - theta*dT);
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	pDM->Update();
	pDM->AssRes(Restmp, 1.);
	pDM->AssJac(*pJac_x, 1.);
	pDM->AssJac(*pJac_xp, 0.);
	
	
	/* Attenzione: a differenza di quanto riportato a p. 16,
	 * "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications"
	 * qui il tempo finale e' in cima, il tempo theta in basso
	 */ 
	 
	/* 2,2 */
	// doublereal J22_x = (1.+3.*rho)/(6.*rho*(1.+rho))*dT;
	Jacxi_x.MulAndSumWithShift(*pJac, jx22, iNumDofs, iNumDofs);
	Jacxi_xp.FakeThirdOrderMulAndSumWithShift(*pJac, EqIsDifferential, 1. - jx22, iNumDofs, iNumDofs);
	
	/* 2,1 */
	// doublereal J21_x = -1./(6.*rho*(1.+rho)*(1.+rho))*dT;
	Jacxi_x.MulAndSumWithShift(*pJac, jx21, iNumDofs, 0);
	Jacxi_xp.FakeThirdOrderMulAndSumWithShift(*pJac, EqIsDifferential, -jx21, iNumDofs, 0);
	
	/* 1,2 */
	// doublereal J12_x = (1.+rho)*(1.+rho)/(6.*rho)*dT;
	Jac_x.MulAndSumWithShift(*pJac, jx12, 0, iNumDofs);
	Jac_xp.FakeThirdOrderMulAndSumWithShift(*pJac, EqIsDifferential, -jx12, 0, iNumDofs);
	
	/* 1,1 */
	// doublereal J11_x = (2.*rho-1.)/(6.*rho)*dT;
	Jac_x.MulAndSumWithShift(*pJac, jx11, 0, 0);
	Jac_xp.FakeThirdOrderMulAndSumWithShift(*pJac, EqIsDifferential, 1. -jx11, 0, 0);
	
	return;
};

void ThirdOrderIntegrator::UpdateDof(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const {
	integer iNumDofs = pDM->iGetNumDofs();
	doublereal dxp = pSol->operator()(DCount);
	doublereal dxp_xi = pSol->operator()(DCount+iNumDofs);
	if (Order == DofOrder::DIFFERENTIAL) {
		
 		pXPrimeCurr->IncCoef(DCount, dxp);
 		pXPrimeCurr->IncCoef(DCount+iNumDofs, dxp_xi);
		
 		pXCurr->IncCoef(DCount, dT*(w1*dxp_xi+w0*dxp));
 		pXCurr->IncCoef(DCount+iNumDofs, 
			dT*(m0*w1*dxp_xi+(m0*w0+n0)*dxp));
	
	} else if (Order == DofOrder::ALGEBRAIC) {
 		pXCurr->IncCoef(DCount, dxp);
 		pXCurr->IncCoef(DCount+iNumDofs, dxp_xi);
		
 		pXPrimeCurr->IncCoef(DCount, dT*(w1*dxp_xi+w0*dxp));
 		pXPrimeCurr->IncCoef(DCount+iNumDofs, 
			dT*(m0*w1*dxp_xi+(m0*w0+n0)*dxp));
	} else {
 		silent_cerr("unknown order for dof " 
			<< DCount<< std::endl);
 		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
};

void
ThirdOrderIntegrator::Update(const VectorHandler* pSol) const
{
  	DEBUGCOUTFNAME("ThirdOrderIntegrator::Predict");
  	ASSERT(pDM != NULL);

#ifdef USE_SCHUR
	SchurDataManager* pSDM;
	if ((pSDM = dynamic_cast<SchurDataManager*> (pDM)) != 0) {
		silent_cerr("Warning: ThirdOrderIntegrator is untested "
			<< "with the parallel solver" << std::endl);
	}
#endif // USE_SCHUR

	UpdateLoop(this,&ThirdOrderIntegrator::UpdateDof,pSol);	
	pDM->Update();
	return;
};

// /* scale factor for tests */
// doublereal
// ThirdOrderIntegrator::TestScale(const NonlinearSolverTest *pResTest) const
// {
// #ifdef __HACK_RES_TEST__
// 
// #ifdef USE_MPI
// #warning "StepNIntegrator TestScale parallel broken !! "	
// #endif /* USE_MPI */
// 
//    	Dof CurrDof;
// 	doublereal dXPr = 0.;
// 
// 	DofIterator.bGetFirst(CurrDof); 
// 
//    	for (int iCntp1 = 1; iCntp1 <= pXPrimeCurr->iGetSize(); 
// 			iCntp1++, DofIterator.bGetNext(CurrDof)) {
// 
// 		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
// 			doublereal d = pXPrimeCurr->operator()(iCntp1);
// 			doublereal d2 = d*d;
// 
// 			doublereal ds = pResTest->dScaleCoef(iCntp1);
// 			doublereal ds2 = ds*ds;
// 			d2 *= ds2;
// 
// 			dXPr += d2;
// 		}
// 		/* else if ALGEBRAIC: non aggiunge nulla */
// 	}
// 
//    	return 1./(1.+dXPr);
// 
// #else /* ! __HACK_RES_TEST__ */
// 	return 1.;
// #endif /* ! __HACK_RES_TEST__ */
// }


TunableThirdOrderIntegrator::TunableThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const bool bmod_res_test)
: ThirdOrderIntegrator(dT, dSolutionTol, iMaxIt, bmod_res_test),
Rho(pRho)
{
	NO_OP;
}

TunableThirdOrderIntegrator::~TunableThirdOrderIntegrator()
{
	NO_OP;
}

void
TunableThirdOrderIntegrator::SetCoef(doublereal dt,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	dT = dt;

	/* from "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications", pp 8-9
	 */
	rho = Rho.dGet();
	theta = -rho/(1.+rho);
	w0 = (1.+3*theta)/(6.*theta);
	w1 = -1./(6.*theta*(1.+theta));
	w2 = (2.+3.*theta)/(6.*(1.+theta));
	/* from "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications", pp 3
	 */
	m0 = 1.-theta*theta*(3.+2.*theta);
	m1 = theta*theta*(3.+2.*theta);
	n0 = theta*(1.+theta)*(1.+theta);
	n1 = theta*theta*(1.+theta);

	/* Attenzione: a differenza di quanto riportato a p. 16,
	 * "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications"
	 * qui il tempo finale e' in cima, il tempo theta in basso
	 */ 
	jx22 = (1.+3.*rho)/(6.*rho*(1.+rho))*dT;
	jx21 = -1./(6.*rho*std::pow(1.+rho,2.))*dT;
	jx12 = std::pow(1.+rho,2.)/(6.*rho)*dT;
	jx11 = (2.*rho-1.)/(6.*rho)*dT;

	DEBUGCOUT("Tunable Third Order Integrator coefficients:" << std::endl
		<< "\t  rho: " << rho << std::endl
		<< "\ttheta: " << theta << std::endl
		<< "\t w0: " << w0 << std::endl
		<< "\t w1: " << w1 << std::endl
		<< "\t w2: " << w2 << std::endl
		<< "\t   m0: " << m0 << std::endl
		<< "\t   m1: " << m1 << std::endl
		<< "\t   n0: " << n0 << std::endl
		<< "\t   n1: " << n1 << std::endl);
}

void
TunableThirdOrderIntegrator::SetDriveHandler(const DriveHandler* pDH)
{
	Rho.pGetDriveCaller()->SetDrvHdl(pDH);
}

AdHocThirdOrderIntegrator::AdHocThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const bool bmod_res_test)
: ThirdOrderIntegrator(dT, dSolutionTol, iMaxIt, bmod_res_test)
{
	NO_OP;
}

AdHocThirdOrderIntegrator::~AdHocThirdOrderIntegrator()
{
	NO_OP;
}

void
AdHocThirdOrderIntegrator::SetCoef(doublereal dt,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	dT = dt;

	/* from "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications", pp 8-9
	 */
	theta = -2./3.;
	w0 = (2.*theta + 1.)/(2.*theta);
	w1 = -1./(2.*theta);
	w2 = 0.;
	/* from "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications", pp 3
	 */
	m0 = 1.-theta*theta;
	m1 = theta*theta;
	n0 = theta*(1.+theta);
	n1 = 0.;

	/* Attenzione: a differenza di quanto riportato a p. 16,
	 * "Unconditionally stable multistep integration of ordinary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications"
	 * qui il tempo finale e' in cima, il tempo theta in basso
	 */ 
	jx22 = 5./12.*dT;
	jx21 = -1./12.*dT;
	jx12 = 3./4.*dT;
	jx11 = 1./4.*dT;

	DEBUGCOUT("Ad Hoc Third Order Integrator coefficients:" << std::endl
		<< "\t  rho: " << rho << std::endl
		<< "\ttheta: " << theta << std::endl
		<< "\t w0: " << w0 << std::endl
		<< "\t w1: " << w1 << std::endl
		<< "\t w2: " << w2 << std::endl
		<< "\t   m0: " << m0 << std::endl
		<< "\t   m1: " << m1 << std::endl
		<< "\t   n0: " << n0 << std::endl
		<< "\t   n1: " << n1 << std::endl);
}

