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
	NO_OP;
}


void ThirdOrderIntegrator::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
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
	if (bAdvanceCalledFirstTime) {
		integer n = pDM->iGetNumDofs();
		Res1.Resize(n);
		Res2.Resize(n);
		Jac11.Resize(n);
		Jac12.Resize(n);
		Jac21.Resize(n);
		Jac22.Resize(n);
		bAdvanceCalledFirstTime = false;
	}
	pXCurr  = pX;
	pXPrev  = qX[0];

	pXPrimeCurr  = pXPrime;
	pXPrimePrev  = qXPrime[0];
	Predict();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	pDM->AfterPredict();

}

void Predict(void) {
	DofIterator.fGetFirst(CurrDof);
	integer iNumDofs = pDM->iGetNumDofs();
   	/* 
	 * Combinazione lineare di stato e derivata 
	 * al passo precedente ecc. 
	 */
	for (int iCntp1 = 1; iCntp1 <= iNumDofs;
		iCntp1++, DofIterator.fGetNext(CurrDof)) {
		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
			doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXPnm1 = 
				pXPrimePrev->dGetCoef(iCntp1);

	 		doublereal dXPn = dPredDer(dXnm1, dXPnm1);
			doublereal dXn = dPredState(dXnm1, dXPn, dXPnm1);
		
	 		pXPrimeCurr->fPutCoef(iCntp1, dXPn);
	 		pXCurr->fPutCoef(iCntp1, dXn);
		
      		} else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXInm1 = pXPrimePrev->dGetCoef(iCntp1);
	 		doublereal dXn = dPredDerAlg(dXInm1, dXnm1);
			doublereal dXIn = dPredStateAlg(dXInm1, dXn, dXnm1);
		
	 		pXCurr->fPutCoef(iCntp1, dXn);
	 		pXPrimeCurr->fPutCoef(iCntp1, dXIn);

		} else {
	 		std::cerr << "unknown order for dof " 
				<< iCntp1<< std::endl;
	 		THROW(ErrGeneric());
		}
   	}
	
}
