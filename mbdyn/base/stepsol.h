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
 
#ifndef STEPSOL_H
#define STEPSOL_H

#include <unistd.h>
#include <array>
#include <cfloat>
#include <cmath>
#include <deque>

/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "dataman.h"
#include "dofown.h"
#include "drive.h"
#include "nonlinpb.h"
#include "nonlin.h"

/* Needed for callback declaration; defined in <mbdyn/base/solver.h> */
class Solver;
class InverseSolver;
 
class StepIntegrator
{

public:
	class ErrGeneric: public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	
	enum { DIFFERENTIAL = 0, ALGEBRAIC = 1 };
	enum StepChange { NEWSTEP, REPEATSTEP };   

protected:	
	DataManager* pDM;
	const DataManager::DofVecType *pDofs; 	/* Dof array, passed by DM */

	bool outputPred;
	
	integer MaxIters;
	doublereal dTol, dSolTol;
	integer steps;
	integer unkstates;
	
	template<class T>
	void UpdateLoop(
		const T* const t,
		void (T::* pUpd)(const int DCount,
			const DofOrder::Order Order,
			const VectorHandler* const pSol) const,
		const VectorHandler* const pSol = 0
	) const;
public:
	StepIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp,
			const integer sts);
	
	virtual ~StepIntegrator(void);

	void SetDataManager(DataManager* pDatMan);
		
	virtual integer GetIntegratorNumPreviousStates(void) const;
	
	virtual integer GetIntegratorNumUnknownStates(void) const;
	
	virtual integer GetIntegratorMaxIters(void) const;
		
	virtual doublereal GetIntegratorDTol(void) const;
	
	virtual doublereal GetIntegratorDSolTol(void) const;

        virtual doublereal dGetCoef(unsigned int iDof) const=0;
     
	virtual void OutputTypes(const bool fpred);
	
	virtual void SetDriveHandler(const DriveHandler* pDH);

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
			doublereal& SolErr) = 0;
};


class ImplicitStepIntegrator:
	public StepIntegrator,
	public NonlinearProblem
{
private:
	/* needed by EvalProd */
	mutable MyVectorHandler XTau;
	mutable MyVectorHandler SavedState;
	mutable MyVectorHandler SavedDerState;
	mutable bool bEvalProdCalledFirstTime;

protected:
	VectorHandler *pXCurr;
	VectorHandler *pXPrimeCurr; 
	bool bModResTest;
public:
	ImplicitStepIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp,
			const integer sts,
			const bool bmod_res_test);
	virtual ~ImplicitStepIntegrator(void);
	
	virtual void
	EvalProd(doublereal Tau, const VectorHandler& f0,
			const VectorHandler& w, VectorHandler& z) const;

	/* scale factor for tests */
	virtual doublereal TestScale(const NonlinearSolverTest *pTest, doublereal& dCoef) const;

};

class DerivativeSolver: 
	public ImplicitStepIntegrator{
private:
	doublereal dCoef;
	const int iMaxIterCoef;
	const doublereal dFactorCoef;

protected:
	void UpdateDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;
public:
	DerivativeSolver(const doublereal Tl, 
			const doublereal dSolTl, 
			const doublereal dC,
			const integer iMaxIt,
			const bool bmod_res_test,
			const integer iMaxIterCoef,
			const doublereal dFactorCoef);

	~DerivativeSolver(void);
	
 	doublereal
	Advance(Solver* pS, 
			const doublereal TStep, 
			const doublereal /* dAph */, 
			const StepChange /* StType */,
			std::deque<VectorHandler*>& qX,
 			std::deque<VectorHandler*>& qXPrime,
			MyVectorHandler*const pX,
 			MyVectorHandler*const pXPrime,
			integer& EffIter,
			doublereal& Err,
			doublereal& SolErr);

 	void Residual(VectorHandler* pRes, VectorHandler* pAbsRes=0) const;

	void Jacobian(MatrixHandler* pJac) const;

	void Jacobian(VectorHandler* pJac, const VectorHandler* pY) const;
	
	void Update(const VectorHandler* pSol) const;

        virtual doublereal dGetCoef(unsigned int iDof) const override;
     
	/* scale factor for tests */
	virtual doublereal TestScale(const NonlinearSolverTest *pTest, doublereal& dAlgebraicEqu) const;
};


/* Base class for integrators of arbitrary order */ 
// FIXME: could probably be merged into the template class tplStepNIntegrator?
class StepNIntegrator :   
	public ImplicitStepIntegrator
{
public:
	doublereal db0Differential;
	doublereal db0Algebraic;

	enum IDX_A {
		IDX_A1 = 0,
		IDX_A2 = 1,
		IDX_A3 = 2,
		IDX_A4 = 3,
		IDX_A5 = 4
		// add as needed
	};

	enum IDX_B {
		IDX_B0 = 0,
		IDX_B1 = 1,
		IDX_B2 = 2,
		IDX_B3 = 3,
		IDX_B4 = 4,
		IDX_B5 = 5
		// add as needed
	};

protected:
	void UpdateDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;
public:
	StepNIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp,
			const bool bmod_res_test);

	virtual ~StepNIntegrator(void);

	virtual void Residual(VectorHandler* pRes, VectorHandler* pAbsRes=0) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
        virtual void Jacobian(VectorHandler* pJac, const VectorHandler* pY) const override;
     
	virtual void Update(const VectorHandler* pSol) const;

        virtual doublereal dGetCoef(unsigned int iDof) const override;
     
	virtual doublereal TestScale(const NonlinearSolverTest *pTest, doublereal& dAlgebraicEqu) const;

protected:
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
};


/* InverseDynamics - Begin*/

class InverseDynamicsStepSolver: 
	public StepIntegrator,
	public NonlinearProblem
{
private:
	/* needed by EvalProd */
	mutable MyVectorHandler XTau;
	mutable MyVectorHandler SavedState;

	/* FIXME: Need this? */
	mutable MyVectorHandler SavedDerState;
	mutable bool bEvalProdCalledFirstTime;

	InverseDynamics::Order iOrder;
	mutable bool m_bJacobian;

protected:
	VectorHandler *pXCurr;
	VectorHandler *pXPrimeCurr; 
	VectorHandler *pXPrimePrimeCurr; 
	VectorHandler *pLambdaCurr;
	bool bModResTest;

public:
	InverseDynamicsStepSolver(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp,
			const integer sts,
			const bool bmod_res_test);

	~InverseDynamicsStepSolver(void);

	virtual void
	EvalProd(doublereal Tau, const VectorHandler& f0,
			const VectorHandler& w, VectorHandler& z) const;

	/* scale factor for tests */
	virtual doublereal TestScale(const NonlinearSolverTest *pTest, doublereal& dCoef) const;

	/* Needed for compatibility with class StepIntegrator */
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
			doublereal& SolErr)
	{
		silent_cerr("InverseDynamicsStepSolver::Advance()");
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};
	
	/* Real Advancer */
	virtual doublereal
	Advance(InverseSolver* pS, 
			const doublereal TStep, 
			const StepChange StType,
			MyVectorHandler*const pX,
 			MyVectorHandler*const pXPrime,
 			MyVectorHandler*const pXPrimePrime,
 			MyVectorHandler*const pLambda,
			integer& EffIter,
			doublereal& Err,
			doublereal& SolErr);
 	
	void Residual(VectorHandler* pRes, VectorHandler* pAbsRes=0) const ;

	void Jacobian(MatrixHandler* pJac) const ;

        virtual void Jacobian(VectorHandler* pJac, const VectorHandler* pY) const override;
	
	void Update(const VectorHandler* pSol) const  ;

        virtual doublereal dGetCoef(unsigned int iDof) const override;
     
	void SetOrder(InverseDynamics::Order iOrder);

	InverseDynamics::Order GetOrder(void) const;

	bool bJacobian(void) const;
};

/* InverseDynamics - End*/

#endif /* STEPSOL_H */
