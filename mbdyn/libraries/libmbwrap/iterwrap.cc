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

 /*****************************************************************************
 *                                                                           *
 *                      ITERATIVE SOLUTION MANAGER                           *
 *                                                                           *
 *****************************************************************************/

/* Giuseppe Quaranta <quaranta@aero.polimi.it> */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <iterwrap.h>

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

template <class S> IterativeSolutionManager::IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pdm,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 S* psm,
				 integer WorkSpaceSize,
				 const doublereal dPivot) 
:
Size(iDim),
pSM(psm),
pDM(pdm),
pP(NULL),
pA(NULL),
pxVH(NULL),
pbVH(NULL),
pXCurr(pX),
pXPrimeCurr(pXP),
Toll(Tollerance),
MaxIt(max),
Tau(0.)
{

	DEBUGCOUT("Entering IterativeSolutionManager" << std::endl);
  	ASSERT(iDim > 0);
          
	SAFENEWWITHCONSTRUCTOR(pSM,
				S,
				S(iDim,
				WorkSpaceSize,
				dPivot));

  	pA   = pSM->pMatHdl();
	pxVH = pSM->pResHdl();
	pbVH = pSM->pSolHdl();
	
	SAFENEWWITHCONSTRUCTOR(pP, Preconditioner,
				Preconditioner(iDim, pA, pSM));

};


IterativeSolutionManager::~IterativeSolutionManager(void){

	if (pSM != 0) {
		SAFEDELETE(pSM);
	}
	if (pP != 0) {
		SAFEDELETE(pP);
	}
};

void IterativeSolutionManager::MatrInit(const doublereal& d = 0.) {

	pSM->MatrInit(d);
	return;
};


void IterativeSolutionManager::ComputeTau(void) {
	Tau = 1.e-6;
};

void IterativeSolutionManager::Solve(doublereal dCoef) {

	ComputeTau();
	
	/* BiCGSTAB Iterative solver */
	//std::cout << "USING BICGSTAB SOLVER" << std::endl; 
	doublereal resid;
	doublereal rho_1, rho_2 = 1., alpha, beta, omega; 
  	MyVectorHandler r(Size), p(Size), phat(Size), s(Size), shat(Size), t(Size), v(Size);
        
	/* Vectors to compute the Matrix-Free Appox. of A * x */
	MyVectorHandler XTau(Size), XPTau(Size), TauRes(Size);
  	
	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */  
	r = *pbVH;
	doublereal normb = pbVH->Norm();
        MyVectorHandler rtilde = r;

  	if (normb == 0.0) { normb = 1;}
  
  	if ((resid = r.Norm() / normb) <= Toll) {
    		DEBUGCOUT("Iterative solver Converged in 0 steps. "
			<< "Residual: " << resid << std::endl);
    		return;
  	}

  	for (int i = 1; i <= MaxIt; i++) {
    		rho_1 = rtilde.InnerProd(r);
    		if (rho_1 == 0) {
    			DEBUGCOUT("Iterative Solver Converged in "
				<< i << " Steps, with Residual: " 
				<< r.Norm() / normb  << std::endl);
   			pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);         
      			return;
    		}
    		if (i == 1) {
			p = r;
		} else {
      			beta = (rho_1/rho_2) * (alpha/omega);
			p.ScalarAddMul(r, p.ScalarAddMul(v, -omega), beta);
    		}
    		pP->solve(p,phat);
		
		/*                                                                 *
		 * J(XCurr) * phat = (Res(XCurr + Tau * phat) - Res(XCurr)) / Tau  *
		 * Res(XCurr) = b                                                  *
		 */
		     
		XTau.ScalarAddMul((*pXCurr), phat, Tau);     
		/* ???????????????????????????????????????????? */
		XPTau.ScalarMul(XTau, 1/dCoef); 
   		pDM->LinkToSolution(XTau, XPTau);         
		pDM->AssRes(v, dCoef);
		v -= *pbVH;
		t.ScalarMul(t, 1/Tau);
    		alpha = rho_1 / rtilde.InnerProd(v);
		s.ScalarAddMul(r, v, -alpha);
    		if ((resid = s.Norm()/normb) < Toll) {
   			pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);         
      			pxVH->ScalarAddMul(phat, alpha);
    			DEBUGCOUT("Iterative Solver Converged in "
				<< i << " Steps, with Residual: " 
				<< resid  << std::endl);
      			return;
    		}
    		pP->solve(s,shat);
		
		/*                                                                 *
		 * J(XCurr) * shat = (Res(XCurr + Tau * shat) - Res(XCurr)) / Tau  *
		 * Res(XCurr) = b                                                  *
		 */
		 
		XTau.ScalarAddMul((*pXCurr), shat, Tau);     
		/* ???????????????????????????????????????????? */
		XPTau.ScalarMul(XTau, 1/dCoef);
   		pDM->LinkToSolution(XTau, XPTau);         
		pDM->AssRes(t, dCoef);
		t -= *pbVH;
		t.ScalarMul(t, 1/Tau);
    		omega = t.InnerProd(s) / t.Dot();
      		pxVH->ScalarAddMul(phat, alpha);
      		pxVH->ScalarAddMul(shat, omega);
		r.ScalarAddMul(s, t, -omega); 

    		rho_2 = rho_1;
    		if ((resid = r.Norm() / normb) < Toll) {
   			pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);         
    			DEBUGCOUT("Iterative Solver Converged in "
				<< i << " Steps, with Residual: " 
				<< resid << std::endl);
      			return;
    		}
    		if (omega == 0) {
   			pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
			DEBUGCOUT("Iterative Solver Converged in "
				<< i << " Steps, with Residual: " 
				<< r.Norm() / normb  << std::endl);
      			return;
    		}
  	}

   	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);         
  	std::cerr << "Iterative solver didn't converge." << std::endl;
	THROW(ErrGeneric());
  	return;
	  
};

#ifdef USE_Y12
template IterativeSolutionManager::IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pdm,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 Y12SparseLUSolutionManager* psm,
				 integer WorkSpaceSize,
				 doublereal dPivot);
#endif /* USE_Y12 */ 

#ifdef USE_MESCHACH
template IterativeSolutionManager::IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pdm,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 MeschachSparseLUSolutionManager* psm,
				 integer WorkSpaceSize,
				 doublereal dPivot);
#endif /* USE_MESCHACH */ 

#ifdef USE_HARWELL
template IterativeSolutionManager::IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pdm,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 HarwellSparseLUSolutionManager* psm,
				 integer WorkSpaceSize,
				 doublereal dPivot);
#endif /* USE_HARWELL */ 

#ifdef USE_UMFPACK3
template IterativeSolutionManager::IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pdm,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 Umfpack3SparseLUSolutionManager * psm,
				 integer WorkSpaceSize,
				 doublereal dPivot);
#endif /* USE_UMFPACK3 */ 
