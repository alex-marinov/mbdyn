/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
  * Marco Morandini	<morandini@aero.polimi.it>
  *
  * third order integrator; orrible code.
  */


#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */


#include "thirdorderstepsol.h"


ThirdOrderIntegrator::ThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const DriveCaller* pRho)
: ImplicitStepIntegrator(iMaxIt, dT, dSolutionTol, 1, 2),
pXPrev(0),
pXPrimePrev(0),
Rho(pRho),
bAdvanceCalledFirstTime(true)
{
	pJacxi_xp = &Jacxi_xp;
	pJacxi_x  = &Jacxi_x;
	pJac_xp   = &Jac_xp;
	pJac_x    = &Jac_x;
}


void ThirdOrderIntegrator::SetCoef(doublereal dt,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	dT = dt;
	doublereal rho = Rho.dGet();
	theta = -rho/(1.+rho);
	w[0] = (2.+3.*theta)/(6*(1+theta));
	w[1] = -1./(6*theta*(1+theta));
	w[2] = (1.+3*theta)/(6.*theta);
	jx[0][0] = (1+3.*rho)/(6*rho*(1+rho))*dT;
	jx[0][1] = -1./(6*rho*std::pow(1+rho,2.))*dT;
	jx[1][0] = std::pow(1+rho,2.)/(6.*rho)*dT;
	jx[1][1] = (2.*rho-1.)/(6.*rho)*dT;
	jxp[0][0] = 1.;
	jxp[0][1] = 0.;
	jxp[1][0] = 0.;
	jxp[1][1] = 1.;
	
};

doublereal ThirdOrderIntegrator::Advance(const doublereal TStep, 
			const doublereal dAlph, 
			const StepChange StType,
			SolutionManager* pSM,
			NonlinearSolver* pNLS, 
			std::deque<MyVectorHandler*>& qX,
	 		std::deque<MyVectorHandler*>& qXPrime,
			MyVectorHandler*const pX,
 			MyVectorHandler*const pXPrime,
			integer& EffIter
#ifdef MBDYN_X_CONVSOL
			, doublereal& SolErr
#endif /* MBDYN_X_CONVSOL */
			){
	doublereal dErr = 0.;
	if (bAdvanceCalledFirstTime) {
		integer n = pDM->iGetNumDofs();
// 		Res1.Resize(n);
// 		Res2.Resize(n);
		Jacxi_xp.Resize(n);
		Jacxi_x.Resize(n);
		Jac_xp.Resize(n);
		Jac_x.Resize(n);
		bAdvanceCalledFirstTime = false;
	}
	pXCurr  = pX;
	pXPrev  = qX[0];

	pXPrimeCurr  = pXPrime;
	pXPrimePrev  = qXPrime[0];
	Predict();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	pDM->AfterPredict();
	
	return dErr;
}

void ThirdOrderIntegrator::Predict(void) {
   	DEBUGCOUTFNAME("ThirdOrderIntegrator::Predict");
   	ASSERT(pDM != NULL);
   	Dof CurrDof;
	
	DofIterator.fGetFirst(CurrDof);
	integer iNumDofs = pDM->iGetNumDofs();
   	/* 
	 * Combinazione lineare di stato e derivata 
	 * al passo precedente ecc. 
	 */
	 /*
	  * Predict per AfterPredict
	  */
	for (int iCntp1 = 1; iCntp1 <= iNumDofs;
		iCntp1++, DofIterator.fGetNext(CurrDof)) {
		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
			doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXPnm1 = 
				pXPrimePrev->dGetCoef(iCntp1);

	 		doublereal dXn = dXnm1+dXPnm1*dT;
			doublereal dXPn = dXPnm1;
		
	 		pXPrimeCurr->fPutCoef(iCntp1, dXPn);
	 		pXCurr->fPutCoef(iCntp1, dXn);
		
      		} else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXInm1 = pXPrimePrev->dGetCoef(iCntp1);
	 		doublereal dXn = dXnm1;
			doublereal dXIn = dXnm1*dT;
		
	 		pXCurr->fPutCoef(iCntp1, dXn);
	 		pXPrimeCurr->fPutCoef(iCntp1, dXIn);

		} else {
	 		std::cerr << "unknown order for dof " 
				<< iCntp1<< std::endl;
	 		THROW(ErrGeneric());
		}
   	}
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();
	/*
	 * Vero Predict
	 */
	for (int iCntp1 = 1; iCntp1 <= iNumDofs;
		iCntp1++, DofIterator.fGetNext(CurrDof)) {
		//simple copy of predicted state
 		pXPrimeCurr->fPutCoef(iCntp1+iNumDofs,
			pXPrimeCurr->dGetCoef(iCntp1));
 		pXCurr->fPutCoef(iCntp1+iNumDofs, pXCurr->dGetCoef(iCntp1));
		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
	 		doublereal dXPnm1 = 
				pXPrimePrev->dGetCoef(iCntp1);
			/* tempo theta*/
	 		doublereal dXn = dXPnm1*(theta-1.)*dT;
	 		pXCurr->fIncCoef(iCntp1+iNumDofs, dXn);
      		} else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
			/* tempo theta*/
			doublereal dXIn = dXnm1*(theta-1.)*dT;
	 		pXPrimeCurr->fIncCoef(iCntp1+iNumDofs, dXIn);

		} else {
	 		std::cerr << "unknown order for dof " 
				<< iCntp1<< std::endl;
	 		THROW(ErrGeneric());
		}
   	}
}

void ThirdOrderIntegrator::Residual(VectorHandler* pRes) const
{
	ASSERT(pDM != NULL);
	
	integer iNumDofs = pDM->iGetNumDofs();

	MyVectorHandler state, stateder, res;
	
	/* theta*dT */
	state.Attach(iNumDofs,pXCurr->pdGetVec()+iNumDofs);
	stateder.Attach(iNumDofs,pXPrimeCurr->pdGetVec()+iNumDofs);
	res.Attach(iNumDofs,pRes->pdGetVec()+iNumDofs);
	pDM->SetTime(pDM->dGetTime()-(1.-theta)*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
	pDM->AssRes(res, 1.);
	
	/* dT */
	state.Attach(iNumDofs,pXCurr->pdGetVec());
	stateder.Attach(iNumDofs,pXPrimeCurr->pdGetVec());
	res.Attach(iNumDofs,pRes->pdGetVec());
	pDM->SetTime(pDM->dGetTime()+(1.-theta)*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
	pDM->AssRes(res, 1.);
	
}

void ThirdOrderIntegrator::Jacobian(MatrixHandler* pJac) const
{
	ASSERT(pDM != NULL);
	
	integer iNumDofs = pDM->iGetNumDofs();

	MyVectorHandler state, stateder;

	/* theta*dT */
	state.Attach(iNumDofs,pXCurr->pdGetVec()+iNumDofs);
	stateder.Attach(iNumDofs,pXPrimeCurr->pdGetVec()+iNumDofs);
	pDM->SetTime(pDM->dGetTime()-(1.-theta)*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
	pDM->AssJac(*pJacxi_x, 1.);
	pDM->AssJac(*pJacxi_xp, 0.);
	
	/* dT */
	state.Attach(iNumDofs,pXCurr->pdGetVec());
	stateder.Attach(iNumDofs,pXPrimeCurr->pdGetVec());
	pDM->SetTime(pDM->dGetTime()+(1.-theta)*dT);
	pDM->LinkToSolution(state, stateder);
	pDM->Update();
	pDM->AssJac(*pJac_x, 1.);
	pDM->AssJac(*pJac_xp, 0.);
	
	
	/* Attenzione: a differenza di quanto riportato a p. 16,
	 * "Unconditionally stable multistep integration of oridnary
	 * differential and differential-algebraic equations with
	 * controlled algorithmic dissipation for multibody dynamic
	 * applications"
	 * qui il tempo finale e' in cima, il tempo theta in basso
	 */ 
	 
	/* 2,2 */
	doublereal J22_x = (1.+3.*rho)/(6.*rho*(1.+rho))*dT;
	Jacxi_x.MulAndSumWithShift(*pJac,J22_x,iNumDofs,iNumDofs);
	Jacxi_xp.MulAndSumWithShift(*pJac,1.-J22_x,iNumDofs,iNumDofs);
	
	/* 2,1 */
	doublereal J21_x = -1./(6.*rho*(1.+rho)*(1.+rho))*dT;
	Jacxi_x.MulAndSumWithShift(*pJac,J21_x,iNumDofs,0);
	Jacxi_xp.MulAndSumWithShift(*pJac,-J21_x,iNumDofs,0);
	
	/* 1,2 */
	doublereal J12_x = (1.+rho)*(1.+rho)/(6.*rho)*dT;
	Jac_x.MulAndSumWithShift(*pJac,J12_x,0,iNumDofs);
	Jac_xp.MulAndSumWithShift(*pJac,-J12_x,0,iNumDofs);
	
	/* 1,1 */
	doublereal J11_x = (2.*rho-1.)/(6.*rho)*dT;
	Jac_x.MulAndSumWithShift(*pJac,J11_x,0,0);
	Jac_xp.MulAndSumWithShift(*pJac,1.-J11_x,0,0);
}


