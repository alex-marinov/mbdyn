/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
  * Copyright (C) 2003-2017
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * Classe che gestisce l'integrazione di un passo temporale
  *
  * Classe di base virtuale pura: NonlinearProblem 
  *
  *	metodi:
  * 		- Residual, Jacobian e Update che vengono richiesti
  * 		  da nonlinear solver
  *
  * Classe derivata: StepIntegrator
  *	 contiene:
  *		- un puntatore al data manager
  *		- le flag per il tipo di output da implementare
  *		- il numero di iterazioni
  *		- le tolleranze per la verifica della convergenza
  * 	
  *	metodi:
  *		- Predict e After Predict	
  *
  */
 
#ifndef STEPSOL_TPL_H
#define STEPSOL_TPL_H

#include "stepsol.h"
#include "stepsol.hc"


/* Base class for multistep integrators using templates */ 
template <unsigned N>
class tplStepNIntegrator :   
	public StepNIntegrator
{
public:
	// helper indexes
	enum IDX_X {
		IDX_Xm1 = 0,
		IDX_Xm2 = 1,
		IDX_Xm3 = 2,
		IDX_Xm4 = 3,
		IDX_Xm5 = 4
		// add as needed
	};

	enum IDX_XP {
		IDX_XP0 = 0,
		IDX_XPm1 = 1,
		IDX_XPm2 = 2,
		IDX_XPm3 = 3,
		IDX_XPm4 = 4,
		IDX_XPm5 = 5
		// add as needed...
	};

	// usage:
	// 	dXm1mN[IDX_Xm1] -> X in the range from minus 1 to minus N
	// 	dXP0mN[IDX_XPm1] -> XP in the range from 0 to minus N

protected:
	const VectorHandler *m_pXPrev[N];
	const VectorHandler *m_pXPrimePrev[N];

	// access as 
	// 	m_a[IDX_A1-N][DIFFERENTIAL|ALGEBRAIC]
	// 	m_b[IDX_B0-N][DIFFERENTIAL|ALGEBRAIC]
	doublereal m_a[N][2];
	doublereal m_b[N + 1][2];

	// at most N coefficients; often much less are used!
	doublereal m_mp[N];
	doublereal m_np[N];
 
public:
	tplStepNIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const bool bmod_res_test);

	virtual ~tplStepNIntegrator(void);

	virtual doublereal
	Advance(Solver* pS, 
			const doublereal TStep, 
			const doublereal dAlph, 
			const StepChange StType,
			std::deque<MyVectorHandler*>& qX,
	 		std::deque<MyVectorHandler*>& qXPrime,
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

	// dXmN: n-1, n-2, ...
	// dXPmN: n, n-1, n-2, ...
   	virtual doublereal 
     	dPredDer(const doublereal dXm1mN[N],
			const doublereal dXP0mN[N + 1]) const = 0;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
   	virtual doublereal 
     	dPredState(const doublereal dXm1mN[N],
			const doublereal dXP0mN[N + 1]) const = 0;   

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
   	virtual doublereal 
     	dPredDerAlg(const doublereal dXm1mN[N],
			const doublereal dXP0mN[N + 1])  const = 0;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
   	virtual doublereal 
     	dPredStateAlg(const doublereal dXm1mN[N],
			const doublereal dXP0mN[N + 1]) const = 0;

	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
};

template <unsigned N>
tplStepNIntegrator<N>::tplStepNIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, N, bmod_res_test)
{
	for (unsigned i = 0; i < N; i++) {
		m_pXPrev[i] = 0;
		m_pXPrimePrev[i] = 0;
	}
}

template <unsigned N>
tplStepNIntegrator<N>::~tplStepNIntegrator(void)
{
	NO_OP;
}

template <unsigned N>
void
tplStepNIntegrator<N>::PredictDof(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXm1mN[N];
		doublereal dXP0mN[N + 1];

		for (unsigned i = 0; i < N; i++) {
			dXm1mN[i] = m_pXPrev[i]->operator()(DCount);
			dXP0mN[i + 1] = m_pXPrimePrev[i]->operator()(DCount);
		}

		dXP0mN[0] = dPredDer(dXm1mN, dXP0mN);
		doublereal dXn = dPredState(dXm1mN, dXP0mN);

		pXPrimeCurr->PutCoef(DCount, dXP0mN[0]);
		pXCurr->PutCoef(DCount, dXn);

	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXIm1mN[N];
		doublereal dX0mN[N + 1];

		for (unsigned i = 0; i < N; i++) {
			dXIm1mN[i] = m_pXPrimePrev[i]->operator()(DCount);
			dX0mN[i + 1] = m_pXPrev[i]->operator()(DCount);
		}

		dX0mN[0] = dPredDerAlg(dXIm1mN, dX0mN);
		doublereal dXIn = dPredStateAlg(dXIm1mN, dX0mN);

		pXCurr->PutCoef(DCount, dX0mN[0]);
		pXPrimeCurr->PutCoef(DCount, dXIn);

	} else {
		silent_cerr("tplStepNIntegrator<" << 2 << ">::PredictDof(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

template <unsigned N>
void
tplStepNIntegrator<N>::Predict(void)
{
   	DEBUGCOUTFNAME("tplStepNIntegrator::Predict");
   	ASSERT(pDM != 0);
	UpdateLoop(this, &tplStepNIntegrator<N>::PredictDof);
}

template <unsigned N>
doublereal
tplStepNIntegrator<N>::Advance(Solver* pS,
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
	pXPrimeCurr  = pXPrime;

	for (unsigned i = 0; i < N; i++) {
		m_pXPrev[i] = qX[i];
		m_pXPrimePrev[i]  = qXPrime[i];
	}

	/* predizione */
	SetCoef(TStep, dAph, StType);
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

#endif /* STEPSOL_TPL_H */
