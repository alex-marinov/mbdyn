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

#ifndef MULTISTAGESTEPSOL_TPL_H
#define MULTISTAGESTEPSOL_TPL_H

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

/* tplStageNIntegrator - begin */

template <unsigned N>
class tplStageNIntegrator :
	public StepNIntegrator
{
public:
	// k-1     s1      s2      s{S-1}       k
	//  o-------o-------o- ... -o-------o
	// dX[S] = {k-1, s1, s2, ..., s{S-1}}
	// dXP[S+1] = {k, k-1, s1, s2, ..., s{S-1}}
	// S: number of stages

	enum IDX_X {
		IDX_Xm1 = 0,	// previous step
		IDX_Xs1 = 1,	// two stages ends here
		IDX_Xs2 = 2,	// three stages ends here
		IDX_Xs3 = 3,	// ...
		IDX_Xs4 = 4,
		IDX_Xs5 = 5	    // six stages ends here
		// add as needed
	};

	enum IDX_XP {
		IDX_XP0 = 0,	// new step
		IDX_XPm1 = 1,	// previous step
		IDX_XPs1 = 2,	// two stages ends here
		IDX_XPs2 = 3,	// three stages ends here
		IDX_XPs3 = 4,	// ...
		IDX_XPs4 = 5,
		IDX_XPs5 = 6	// six stages ends here
		// add as needed
	};

	enum IDX_A {
		IDX_Am1 = 0,
		IDX_As1 = 1,
		IDX_As2 = 2,
		IDX_As3 = 3,
		IDX_As4 = 4,
		IDX_As5 = 5
		// add as needed
	};

	enum IDX_B {
		IDX_B0 = 0,
		IDX_Bm1 = 1,
		IDX_Bs1 = 2,
		IDX_Bs2 = 3,
		IDX_Bs3 = 4,
		IDX_Bs4 = 5,
		IDX_Bs5 = 6
		// add as needed
	};

protected:
	VectorHandler *m_pX[N];		// state vectors
	VectorHandler *m_pXP[N];	// state derivative vectors

	std::deque<VectorHandler*> m_qXPr, m_qXPPr;	// queues for prediction

	doublereal m_a[N][2];
	doublereal m_b[N + 1][2];

	doublereal m_mp[N];
	doublereal m_np[N];

public:
	tplStageNIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test);

	virtual ~tplStageNIntegrator(void);

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
	template <unsigned S>
	void PredictDofForStageS(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;

	void PredictForStageS(unsigned S);

	virtual doublereal
     	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[N],
		const doublereal dXP0mN[N + 1]) const = 0;

   	virtual doublereal
     	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[N],
		const doublereal dXP0mN[N + 1]) const = 0;

   	virtual doublereal
     	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[N],
		const doublereal dXP0mN[N + 1]) const = 0;

   	virtual doublereal
     	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[N],
		const doublereal dXP0mN[N + 1]) const = 0;

	virtual void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep) = 0;

	// Not needed; use specialized function SetCoefForStageS()
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
};


template <unsigned N>
tplStageNIntegrator<N>::tplStageNIntegrator(const integer MaxIt,
	const doublereal dT,
	const doublereal dSolutionTol,
	const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, 1, bmod_res_test)
{
	for (unsigned i = 0; i < N; i++) {
		m_pX[i] = 0;
		m_pXP[i] = 0;
	}
}

template <unsigned N>
tplStageNIntegrator<N>::~tplStageNIntegrator(void)
{
	m_pX[0] = 0;
	m_pXP[0] = 0;

	// NOTE: start from 1 because tplStageNIntegrator does not own the vector in 0
	for (unsigned i = 1; i < N; i++) {
		if (m_pX[i] != 0) {
			SAFEDELETE(m_pX[i]);
			m_pX[i] = 0;
		}

		if (m_pXP[i] != 0) {
			SAFEDELETE(m_pXP[i]);
			m_pXP[i] = 0;
		}
	}
}

template <unsigned N> template <unsigned S>
void
tplStageNIntegrator<N>::PredictDofForStageS(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXm1mN[N];
		doublereal dXP0mN[N + 1];

		for (unsigned i = 0; i < S; i++) {
			dXm1mN[i] = m_pX[i]->operator()(DCount);
			dXP0mN[i + 1] = m_pXP[i]->operator()(DCount);
		}

		doublereal dXP = dPredDerForStageS(S, dXm1mN, dXP0mN);
		if (S == N) {
			dXP0mN[0] = dXP;
		} else {
			dXP0mN[S + 1] = dXP;
		}
		doublereal dX = dPredStateForStageS(S, dXm1mN, dXP0mN);

		pXPrimeCurr->PutCoef(DCount, dXP);
		pXCurr->PutCoef(DCount, dX);

		//std::cout << "Stage =" << S << std::endl;
		//std::cout << "Predict: DCount =" << DCount << "; Xm=" << dXm1mN[0] << "; XPm=" << dXP0mN[1] << std::endl;
		//std::cout << "Xs1=" << dXm1mN[1] << "; XPs1=" << dXP0mN[2] << std::endl;
		//std::cout << "Xs2=" << dXm1mN[2] << "; XPs2=" << dXP0mN[3] << std::endl;
		//std::cout << "Xs3=" << dXm1mN[3] << "; XPs3=" << dXP0mN[4] << std::endl;
		//std::cout << "Xs4=" << dXm1mN[4] << "; XPs4=" << dXP0mN[5] << std::endl;
		//std::cout << "Xs5=" << dXm1mN[5] << "; XPs5=" << dXP0mN[0] << std::endl;
		//std::cout << "dX = " << dX << "; dXP = " << dXP << std::endl;
	}
	else if (Order == DofOrder::ALGEBRAIC)
	{
		doublereal dXIm1mN[N];
		doublereal dX0mN[N + 1];

		for (unsigned i = 0; i < S; i++) {
			dXIm1mN[i] = m_pXP[i]->operator()(DCount);
			dX0mN[i + 1] = m_pX[i]->operator()(DCount);
		}

		doublereal dX = dPredDerAlgForStageS(S, dXIm1mN, dX0mN);
		if (S == N) {
			dX0mN[0] = dX;
		} else {
			dX0mN[S + 1] = dX;
		}
		doublereal dXI = dPredStateAlgForStageS(S, dXIm1mN, dX0mN);

		pXCurr->PutCoef(DCount, dX);
		pXPrimeCurr->PutCoef(DCount, dXI);

	} else {
		silent_cerr("tplStageNIntegrator::"
			"PredictDofForStage" << S << "(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

template <unsigned N>
void
tplStageNIntegrator<N>::PredictForStageS(unsigned S)
{
   	DEBUGCOUTFNAME("tplStageNIntegrator::PredictForStageS");
   	ASSERT(pDM != 0);

	switch (S) {
	case 1:
		UpdateLoop(this, &tplStageNIntegrator<N>::PredictDofForStageS<1>);
		break;

	case 2:
		UpdateLoop(this, &tplStageNIntegrator<N>::PredictDofForStageS<2>);
		break;

	case 3:
		UpdateLoop(this, &tplStageNIntegrator<N>::PredictDofForStageS<3>);
		break;

	case 4:
		UpdateLoop(this, &tplStageNIntegrator<N>::PredictDofForStageS<4>);
		break;

	case 5:
		UpdateLoop(this, &tplStageNIntegrator<N>::PredictDofForStageS<5>);
		break;

	// add more as needed

	default:
		silent_cerr("Multistage integrators with " << S << "stages not supported" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

template <unsigned N>
doublereal
tplStageNIntegrator<N>::Advance(Solver* pS,
	const doublereal TStep,
	const doublereal dAph, const StepChange StType,
	std::deque<VectorHandler*>& qX,
	std::deque<VectorHandler*>& qXPrime,
	MyVectorHandler*const pX,
	MyVectorHandler*const pXPrime,
	integer& EffIter,
	doublereal& Err,
	doublereal& SolErr)
{
	ASSERT(pDM != 0);

	pXCurr = pX;
	pXPrimeCurr = pXPrime;

	// the solution at the previous step we take from qX, qXPrime
	m_pX[0] = qX[0];
	m_pXP[0] = qXPrime[0];

	ASSERT(m_qXPr.empty());
	ASSERT(m_qXPPr.empty());

	m_qXPr.push_front(qX[0]);
	m_qXPPr.push_front(qXPrime[0]);

	// first time?
	if (m_pX[1] == 0) {
		// we own memory for internal stage(s) (from 1 on; see destructor)
		for (unsigned S = 1; S < N; S++) {
			SAFENEWWITHCONSTRUCTOR(m_pX[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
			SAFENEWWITHCONSTRUCTOR(m_pXP[S],
				MyVectorHandler,
				MyVectorHandler(pDM->iGetNumDofs()));
		}
	}

	// First-stage
	SetCoefForStageS(1, TStep, dAph, StType);
	PredictForStageS(1);
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();

#ifdef DEBUG
	integer iNumDofs = pDM->iGetNumDofs();
	if (outputPred) {
		std::cout << "After prediction, stage " << 1 << " of " << N << " time=" << pDM->dGetTime() << std::endl;
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
	EffIter = 0;
	integer tmpEffIter = EffIter;
	doublereal tmpErr = Err;

	pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol,
    		tmpEffIter, tmpErr, dSolTol, SolErr);

	EffIter += tmpEffIter;
	Err += tmpErr;

	// if it gets here, it surely converged
	pDM->AfterConvergence();

	// Second-stage (and subsequent)
	for (unsigned S = 2; S <= N; S++) {
		// pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *m_pX[S - 2], *m_pXP[S - 2]);

		m_qXPr.push_front(m_pX[S - 2]);
		m_qXPPr.push_front(m_pXP[S - 2]);
		pDM->BeforePredict(*pXCurr, *pXPrimeCurr, m_qXPr, m_qXPPr);

		// copy from pX, pXPrime to m_pX, m_pXP
		for (integer i = 1; i <= pDM->iGetNumDofs(); i++) {
			m_pX[S - 1]->PutCoef(i, pXCurr->operator()(i));
			m_pXP[S - 1]->PutCoef(i, pXPrimeCurr->operator()(i));
		}

		SetCoefForStageS(S, TStep, dAph, StType);
		PredictForStageS(S);
		pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
     	 	pDM->AfterPredict();

#ifdef DEBUG
		integer iNumDofs = pDM->iGetNumDofs();
		if (outputPred) {
			std::cout << "After prediction, stage " << S << " of " << N << " time=" << pDM->dGetTime() << std::endl;
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

		tmpErr = 0.;
		tmpEffIter = 0;
		pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol,
    			tmpEffIter, tmpErr, dSolTol, SolErr);
		EffIter += tmpEffIter;
		Err += tmpErr;

		// if it gets here, it surely converged
		pDM->AfterConvergence();
	}

	while (!m_qXPr.empty()) {
		m_qXPr.pop_back();
		m_qXPPr.pop_back();
	}

	ASSERT(m_qXPr.empty());
	ASSERT(m_qXPPr.empty());

	return Err;
}

// Not needed; use specialized function SetCoefForStageS()
template <unsigned N>
void
tplStageNIntegrator<N>::SetCoef(doublereal dT, 
	doublereal dAlpha,
	enum StepChange NewStep)
{
	NO_OP;
}

/* tplStageNIntegrator - end */

#endif // MULTISTAGESTEPSOL_TPL_H
