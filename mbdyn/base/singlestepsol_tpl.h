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
 
#ifndef SINGLESTEPSOL_TPL_H
#define SINGLESTEPSOL_TPL_H

#include <unistd.h>
#include <cfloat>
#include <cmath>
#include <deque>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "dataman.h"
#include "dofown.h"
#include "drive.h"
#include "nonlinpb.h"
#include "nonlin.h"

/* tplSingleStepIntegrator - begin */ 

template <unsigned N>
class tplSingleStepIntegrator :   
	public StepNIntegrator
{
public:
	// helper indexes
	enum IDX_XP {
		IDX_XP0 = 0,
		IDX_XPm1 = 1,
		IDX_XPI1 = 2,   // first intermediate variable 
		IDX_XPI2 = 3,   // second intermediate varibale
		IDX_XPI3 = 4,   // ...
		// add as needed...
	};

    enum IDX_B {
		IDX_B0 = 0,
		IDX_Bm1 = 1,
		IDX_BI1 = 2,
		IDX_BI2 = 3,
		IDX_BI3 = 4
		// add as needed
	};

    enum IDX_C {
        IDX_C1 = 0,
        IDX_C2 = 1,
        IDX_C3 = 2
        // add as needed
    };

    enum IDX_D {
        IDX_D1 = 0,
        IDX_D2 = 1,
        IDX_D3 = 2
        // add as needed
    };

    enum IDX_E {
        IDX_E1 = 0,
        IDX_E2 = 1,
        IDX_E3 = 2
        // add as needed
    };

protected:
	const VectorHandler *m_pXPrev;
	const VectorHandler *m_pXPrimePrev;
    VectorHandler *m_pXPrimeIntePrev[N - 1];    // store intermediate variable of last step
    VectorHandler *pXPrimeInteCurr[N - 1];    // store intermediate variable of current step

	doublereal m_b[N + 1][2];
	doublereal m_c[N - 1];
    doublereal m_d[N - 1];
    doublereal m_e[N - 1];
    // add as needed

 
public:
	tplSingleStepIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test);

	virtual ~tplSingleStepIntegrator(void);

	virtual doublereal
	Advance(Solver* pS, 
		const doublereal TStep, 
		const doublereal dAlph, 
		const StepChange StType,
		std::deque<VectorHandler*>& qX,
 		std::deque<VectorHandler*>& qXPrime,
		MyVectorHandler*const pX,
 		MyVectorHandler*const pXPrime,
		integer& EffIter,
		doublereal& Err,
		doublereal& SolErr);

protected:
	void PredictDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;

	virtual void Predict(void);

   	virtual doublereal 
     	dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[N + 1]) const = 0;

   	virtual doublereal 
     	dPredState(const doublereal dXm1,
		const doublereal dXP0mI[N + 1]) const = 0;   

   	virtual doublereal 
     	dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[N + 1]) const = 0;

   	virtual doublereal 
     	dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[N + 1]) const = 0;

	void UpdateInteDof(const int DCount);

    virtual void UpdateInte(void);

    virtual doublereal 
     	dUpdateInte(unsigned uNumber, 
		 const doublereal dXPmI, const doublereal dXPIm, const doublereal dXPmIm) const = 0;

    virtual void SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange NewStep) = 0;
};

template <unsigned N>
tplSingleStepIntegrator<N>::tplSingleStepIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, N, bmod_res_test)
{
	for (unsigned i = 0; i < N - 1; i++) {
		m_pXPrimeIntePrev[i] = 0;
        pXPrimeInteCurr[i] = 0;
	}
}

template <unsigned N>
tplSingleStepIntegrator<N>::~tplSingleStepIntegrator(void)
{
	for (unsigned i = 0; i < N - 1; i++) {
		if (m_pXPrimeIntePrev[i] != 0) {
			SAFEDELETE(m_pXPrimeIntePrev[i]);
			m_pXPrimeIntePrev[i] = 0;
		}

		if (pXPrimeInteCurr[i] != 0) {
			SAFEDELETE(pXPrimeInteCurr[i]);
			pXPrimeInteCurr[i] = 0;
		}
	}
}

template <unsigned N>
void
tplSingleStepIntegrator<N>::PredictDof(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXm1;
		doublereal dXP0mI[N + 1];

        dXm1 = m_pXPrev->operator()(DCount);
        dXP0mI[1] = m_pXPrimePrev->operator()(DCount);
		for (unsigned i = 0; i < N - 1; i++) {
			dXP0mI[i + 2] = m_pXPrimeIntePrev[i]->operator()(DCount);
		}

		dXP0mI[0] = dPredDer(dXm1, dXP0mI);
		doublereal dXn = dPredState(dXm1, dXP0mI);

		pXPrimeCurr->PutCoef(DCount, dXP0mI[0]);
		pXCurr->PutCoef(DCount, dXn);

		//std::cout << "Predict: DCount = " << DCount << "; dXm1 = " << dXm1 << "; dXPm1 = " << dXP0mI[1] << std::endl;
		//std::cout << "dXPI1m = " <<  dXP0mI[2] << "; dXPI2m = " <<  dXP0mI[3] << "; dXPI3m = " <<  dXP0mI[4] << std::endl;
		//std::cout << "dX = " << dXn << "; dXP = " << dXP0mI[0] << std::endl;
		
	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXIm1;
		doublereal dX0mI[N + 1];

        dX0mI[1] = m_pXPrev->operator()(DCount);
        dXIm1 = m_pXPrimePrev->operator()(DCount);
		for (unsigned i = 0; i < N - 1; i++) {
			dX0mI[i + 2] = m_pXPrimeIntePrev[i]->operator()(DCount);
		}

		dX0mI[0] = dPredDerAlg(dXIm1, dX0mI);
		doublereal dXIn = dPredStateAlg(dXIm1, dX0mI);

		pXCurr->PutCoef(DCount, dX0mI[0]);
		pXPrimeCurr->PutCoef(DCount, dXIn);
		

	} else {
		silent_cerr("tplSingleStepIntegrator<" << N << ">::PredictDof(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

template <unsigned N>
void
tplSingleStepIntegrator<N>::Predict(void)
{
   	DEBUGCOUTFNAME("tplSingleStepIntegrator::Predict");
   	ASSERT(pDM != 0);
	UpdateLoop(this, &tplSingleStepIntegrator<N>::PredictDof);
}

template <unsigned N>
void
tplSingleStepIntegrator<N>::UpdateInteDof(const int DCount) 
{
    doublereal dXPmI;
	doublereal dXPIm;
    doublereal dXPmIm;
    for (unsigned i = 0; i < N - 1; i++) {
        dXPmI = m_pXPrimeIntePrev[i]->operator()(DCount);
        if (i == 0) {
            dXPIm = pXPrimeCurr->operator()(DCount);
            dXPmIm = m_pXPrimePrev->operator()(DCount);
        }else 
        {
            dXPIm = pXPrimeInteCurr[i - 1]->operator()(DCount);
            dXPmIm = m_pXPrimeIntePrev[i - 1]->operator()(DCount);
        }
        doublereal dXPInte = dUpdateInte(i + 1, dXPmI, dXPIm, dXPmIm);
        pXPrimeInteCurr[i]->PutCoef(DCount, dXPInte);
    }
	//std::cout << "UpdataInte: DCount = " << DCount << "; dXPI1 = " <<  pXPrimeInteCurr[0]->operator()(DCount) << "; dXPI2 = " << pXPrimeInteCurr[1]->operator()(DCount) <<  "; dXPI3 = " << pXPrimeInteCurr[2]->operator()(DCount) <<std::endl;
}

template <unsigned N>
void
tplSingleStepIntegrator<N>::UpdateInte(void)
{
   	ASSERT(pDM != 0);
	for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
		UpdateInteDof(i);
	}
}

template <unsigned N>
doublereal
tplSingleStepIntegrator<N>::Advance(Solver* pS, 
		const doublereal TStep, 
		const doublereal dAlph, 
		const StepChange StType,
		std::deque<VectorHandler*>& qX,
 		std::deque<VectorHandler*>& qXPrime,
		MyVectorHandler*const pX,
 		MyVectorHandler*const pXPrime,
		integer& EffIter,
		doublereal& Err,
		doublereal& SolErr)
{
	ASSERT(pDM != NULL);

	pXCurr  = pX;
	pXPrimeCurr  = pXPrime;

    m_pXPrev = qX[0];
    m_pXPrimePrev = qXPrime[0];

	if (m_pXPrimeIntePrev[0] == 0){
        for (unsigned S = 0; S < N - 1; S++) {
            SAFENEWWITHCONSTRUCTOR(m_pXPrimeIntePrev[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
                for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
                    m_pXPrimeIntePrev[S]->PutCoef(i, m_pXPrimePrev->operator()(i));
                }
            SAFENEWWITHCONSTRUCTOR(pXPrimeInteCurr[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
        }
    }else
	{
		for (unsigned S = 0; S < N - 1; S++) {
        	for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
            	m_pXPrimeIntePrev[S]->PutCoef(i,qXPrime[S + 1]->operator()(i));
       	 	}
    	}
	}


	/* predizione */
	SetCoef(TStep, dAlph, StType);
	Predict();
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

	/* if it gets here, it surely converged */
	pDM->AfterConvergence();

    UpdateInte();

	for (unsigned S = 0; S < N - 1; S++) {
        	for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
				qX[S]->PutCoef(i,pXCurr->operator()(i));
            	qXPrime[S]->PutCoef(i,pXPrimeInteCurr[S]->operator()(i));
       	 	}
    }

	return Err;
}

#endif /* SINGLESTEPSOL_TPL_H */