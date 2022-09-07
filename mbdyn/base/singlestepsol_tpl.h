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

	enum IDX_N
	{
		IDX_Real = 0, 
		IDX_Imag = 1
		// allow complex coefficients
	};

protected:
	const VectorHandler *m_pXPrev;
	const VectorHandler *m_pXPrimePrev[2];		// need XP of step k-1 and k-2 to update intermediate varibales XPI of step k-1
	VectorHandler *m_pXPrimeIntePrevReal[N - 1];    // store real part of intermediate variables of last step
	VectorHandler *m_pXPrimeIntePrevImag[N - 1];	// store imaginary part of intermediate variables of last step

	doublereal m_b[N + 1][2][2];
	doublereal m_c[N - 1][2];
    doublereal m_d[N - 1][2];
    doublereal m_e[N - 1][2];
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
		const doublereal dXP0mI[N + 1][2]) const = 0;

   	virtual doublereal 
     	dPredState(const doublereal dXm1,
		const doublereal dXP0mI[N + 1][2]) const = 0;   

   	virtual doublereal 
     	dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[N + 1][2]) const = 0;

   	virtual doublereal 
     	dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[N + 1][2]) const = 0;

	void UpdateInteDof(const int DCount);

    virtual void UpdateInte(void);

    virtual doublereal 
     	dUpdateInteReal(unsigned uNumber, 
		 const doublereal dXPmI[2], const doublereal dXPIm[2], const doublereal dXPmIm[2]) const = 0;

	virtual doublereal 
     	dUpdateInteImag(unsigned uNumber, 
		 const doublereal dXPmI[2], const doublereal dXPIm[2], const doublereal dXPmIm[2]) const = 0;

    virtual void SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange NewStep) = 0;
};

template <unsigned N>
tplSingleStepIntegrator<N>::tplSingleStepIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, 2, bmod_res_test)
{
	for (unsigned S = 0; S < N - 1; S++) {
		m_pXPrimeIntePrevReal[S] = 0;
		m_pXPrimeIntePrevImag[S] = 0;
	}
}

template <unsigned N>
tplSingleStepIntegrator<N>::~tplSingleStepIntegrator(void)
{
	for (unsigned S = 0; S < N - 1; S++) {
		if (m_pXPrimeIntePrevReal[S] != 0) {
			SAFEDELETE(m_pXPrimeIntePrevReal[S]);
			m_pXPrimeIntePrevReal[S] = 0;
		}
		if (m_pXPrimeIntePrevImag[S] != 0) {
			SAFEDELETE(m_pXPrimeIntePrevImag[S]);
			m_pXPrimeIntePrevImag[S] = 0;
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
		doublereal dXP0mI[N + 1][2];

        dXm1 = m_pXPrev->operator()(DCount);
        dXP0mI[1][0] = m_pXPrimePrev[0]->operator()(DCount);
		for (unsigned i = 0; i < N - 1; i++)
		{
			dXP0mI[i + 2][0] = m_pXPrimeIntePrevReal[i]->operator()(DCount);
			dXP0mI[i + 2][1] = m_pXPrimeIntePrevImag[i]->operator()(DCount);
		}

		dXP0mI[0][0] = dPredDer(dXm1, dXP0mI);
		doublereal dXn = dPredState(dXm1, dXP0mI);
		
		pXPrimeCurr->PutCoef(DCount, dXP0mI[0][0]);
		pXCurr->PutCoef(DCount, dXn);
		
		
	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXIm1;
		doublereal dX0mI[N + 1][2];

        dX0mI[1][0] = m_pXPrev->operator()(DCount);
        dXIm1 = m_pXPrimePrev[0]->operator()(DCount);
		for (unsigned i = 0; i < N - 1; i++) {
			dX0mI[i + 2][0] = m_pXPrimeIntePrevReal[i]->operator()(DCount);
			dX0mI[i + 2][1] = m_pXPrimeIntePrevImag[i]->operator()(DCount);
		}

		dX0mI[0][0] = dPredDerAlg(dXIm1, dX0mI);
		doublereal dXIn = dPredStateAlg(dXIm1, dX0mI);

		pXCurr->PutCoef(DCount, dX0mI[0][0]);
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
	doublereal dXP1[N - 1][2];
	doublereal dXP2[N - 1][2];
	doublereal dXP3[N - 1][2];
	for (unsigned i = 0; i < N - 1; i++){
		dXP1[i][0] = m_pXPrimeIntePrevReal[i]->operator()(DCount);
		dXP1[i][1] = m_pXPrimeIntePrevImag[i]->operator()(DCount);
		if (i == 0)
		{
			dXP3[i][0] = m_pXPrimePrev[1]->operator()(DCount);
			dXP3[i][1] = 0.;
		}
		else
		{
			dXP3[i][0] = m_pXPrimeIntePrevReal[i - 1]->operator()(DCount);
			dXP3[i][1] = m_pXPrimeIntePrevImag[i - 1]->operator()(DCount);
		}
	}
	for (unsigned i = 0; i < N - 1; i++)
	{
		if (i == 0)
		{
			dXP2[i][0] = m_pXPrimePrev[0]->operator()(DCount);
			dXP2[i][1] = 0.;
		}
		else
		{
			dXP2[i][0] = m_pXPrimeIntePrevReal[i - 1]->operator()(DCount);
			dXP2[i][1] = m_pXPrimeIntePrevImag[i - 1]->operator()(DCount);
		}
		doublereal dXPInteReal = dUpdateInteReal(i + 1, dXP1[i], dXP2[i], dXP3[i]);
		m_pXPrimeIntePrevReal[i]->PutCoef(DCount, dXPInteReal);
		doublereal dXPInteImag = dUpdateInteImag(i + 1, dXP1[i], dXP2[i], dXP3[i]);
		m_pXPrimeIntePrevImag[i]->PutCoef(DCount, dXPInteImag);
	}	
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
	m_pXPrimePrev[0] = qXPrime[0];

	SetCoef(TStep, dAlph, StType);

	if (m_pXPrimeIntePrevReal[0] == 0){
        for (unsigned S = 0; S < N - 1; S++) {
            SAFENEWWITHCONSTRUCTOR(m_pXPrimeIntePrevReal[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
                for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
                    m_pXPrimeIntePrevReal[S]->PutCoef(i, m_pXPrimePrev[0]->operator()(i));
                }
		}
		for (unsigned S = 0; S < N - 1; S++) {
            SAFENEWWITHCONSTRUCTOR(m_pXPrimeIntePrevImag[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
                for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
                    m_pXPrimeIntePrevImag[S]->PutCoef(i, 0.);
                }
		}
		qX.pop_back();	// 	To distinguish ssn, which needs BeforePredict() to achieve different function
	}
	else
	{
		m_pXPrimePrev[1] = qXPrime[1];
		UpdateInte();	// Update intermediate variables of step k-1 using XP of step k-1 and k-2
	}

	/* predizione */
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

	return Err;
}

#endif /* SINGLESTEPSOL_TPL_H */
