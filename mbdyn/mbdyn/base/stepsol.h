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
#include <ac/float.h>
#include <ac/math.h>


/* per il debugging */
#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <solman.h>
#include <deque>  
#include <dataman.h> 
#include "dofown.h"
#include "drive.h"
#include <nonlinpb.h>
#include <nonlin.h>
 
class StepIntegrator
{

public:
	class ErrGeneric{};
	
	enum { DIFFERENTIAL = 0, ALGEBRAIC = 1 };
	enum StepChange { NEWSTEP, REPEATSTEP };   

protected:	
	DataManager* pDM;
	VecIter<Dof> DofIterator; 	/* Iteratore per la struttura dei Dof,
					 * passato da DM */

	bool outputPred;
	
	integer MaxIters;
	doublereal dTol, dSolTol;
	integer steps;
	integer unkstates;

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
	
	virtual void OutputTypes(const bool fpred);
	
	virtual void SetDriveHandler(const DriveHandler* pDH) = 0;

	virtual doublereal
	Advance(const doublereal TStep, 
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
			) = 0;
};


class ImplicitStepIntegrator:
	public StepIntegrator,
	public NonlinearProblem
{
private:
	/* needed by EvalProd */
	mutable MyVectorHandler XTau;
	mutable bool bEvalProdCalledFirstTime;
protected:
	VectorHandler *pXCurr;
	VectorHandler *pXPrimeCurr; 
public:
	ImplicitStepIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp,
			const integer sts) :
	StepIntegrator(MaxIt,dT,dSolutionTol,stp,sts),
	bEvalProdCalledFirstTime(true),
	pXCurr(0),
	pXPrimeCurr(0)
	{
		NO_OP;
	};
	virtual ~ImplicitStepIntegrator(){
		NO_OP;
	};
	
	virtual void
	EvalProd(doublereal Tau, const VectorHandler& f0,
			const VectorHandler& w, VectorHandler& z) const;

};

class DerivativeSolver: 
	public ImplicitStepIntegrator{
private:
	doublereal dCoef;

public:
	DerivativeSolver(const doublereal Tl, 
			const doublereal dSolTl, 
			const doublereal dC,
			const integer iMaxIt);

	~DerivativeSolver(void);
	
	void SetDriveHandler(const DriveHandler* /* pDH */ );
	
 	doublereal
	Advance(const doublereal TStep, 
			const doublereal /* dAph */, 
			const StepChange /* StType */,
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
			);

 	void Residual(VectorHandler* pRes) const;

	void Jacobian(MatrixHandler* pJac) const;
	
	void Update(const VectorHandler* pSol) const;

	/* scale factor for tests */
#ifdef __HACK_SCALE_RES__
	virtual doublereal TestScale(const VectorHandler *pScale) const;
#else /* ! __HACK_SCALE_RES__ */
	virtual doublereal TestScale(void) const;
#endif /* ! __HACK_SCALE_RES__ */

 
};


/* classe di base per gli integratori di ordine qualsiasi */ 
class StepNIntegrator :   
	public ImplicitStepIntegrator
{
protected:
	doublereal db0Differential;
	doublereal db0Algebraic;

public:
	StepNIntegrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const integer stp);

	virtual ~StepNIntegrator(void);

	virtual void Residual(VectorHandler* pRes) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
	virtual void Update(const VectorHandler* pSol) const;

	/* scale factor for tests */
#ifdef __HACK_SCALE_RES__
	virtual doublereal TestScale(const VectorHandler *pScale) const;
#else /* ! __HACK_SCALE_RES__ */
	virtual doublereal TestScale(void) const;
#endif /* ! __HACK_SCALE_RES__ */

protected:
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
};

/* classe di base per gli integratori del second'ordine */ 
class Step1Integrator :   
	public StepNIntegrator
{
protected:
	VectorHandler *pXPrev;
	VectorHandler *pXPrimePrev; 

public:
	Step1Integrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol);

	virtual ~Step1Integrator(void);

	virtual doublereal
	Advance(const doublereal TStep, 
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
			);

protected:
	virtual void Predict(void);

	// Overridden by dedicated inline functions
	virtual doublereal 
     	dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXPm1,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
   
   	// Overridden by dedicated inline functions
   	virtual doublereal 
     	dPredictState(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
 
   	virtual doublereal 
     	dPredDer(const doublereal& dXm1,
	      const doublereal& dXPm1) const = 0;
   
   	virtual doublereal 
     	dPredState(const doublereal& dXm1,
		const doublereal& dXP,
		const doublereal& dXPm1) const = 0;   

   	virtual doublereal 
     	dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1)  const = 0;
		    
   	virtual doublereal 
     	dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1) const = 0;
};


class CrankNicholsonSolver: 
	public Step1Integrator
{
public:
	CrankNicholsonSolver(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt);

	~CrankNicholsonSolver(void);

	void SetDriveHandler(const DriveHandler* pDH);
   
protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
  
   	doublereal 
     	dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXPm1,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
	doublereal 
	dPredictState(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
	// Nota: usa predizione lineare per le derivate (massimo ordine possibile)
	doublereal 
	dPredDer(const doublereal& dXm1,
	      const doublereal& dXPm1) const;
   
	doublereal 
	dPredState(const doublereal& dXm1,
		const doublereal& dXP,
		const doublereal& dXPm1) const;
   
	doublereal 
	dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1) const;
   
	doublereal 
	dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1) const;
};

/* classe di base per gli integratori del second'ordine */ 
class Step2Integrator :   
	public StepNIntegrator
{
protected:
	VectorHandler *pXPrev, *pXPrev2;
	VectorHandler *pXPrimePrev, *pXPrimePrev2; 

public:
	Step2Integrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol);

	virtual ~Step2Integrator(void);

	virtual doublereal Advance(const doublereal TStep, 
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
			);

protected:
	virtual void Predict(void);

	// Overridden by dedicated inline functions
	virtual doublereal 
     	dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
   
   	// Overridden by dedicated inline functions
   	virtual doublereal 
     	dPredictState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const = 0;
 
   	virtual doublereal 
     	dPredDer(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const = 0;
   
   	virtual doublereal 
     	dPredState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const = 0;   
   	virtual doublereal 
     	dPredDerAlg(const doublereal& dXm1,
			const doublereal& dXPm1,
			const doublereal& dXPm2)  const = 0;
		    
   	virtual doublereal 
     	dPredStateAlg(const doublereal& dXm1,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const = 0;

	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
};

/* NostroMetodo - begin */

class MultistepSolver: 
	public Step2Integrator
{
protected:
	DriveOwner Rho;
	DriveOwner AlgebraicRho;
   
	doublereal a[2][2];
	doublereal b[3][2];

	doublereal mp[2];
	doublereal np[2];
   
public:
	MultistepSolver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho);

	~MultistepSolver(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal 
	dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;

	doublereal 
	dPredictState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;
   
	// Nota: usa predizione cubica per le derivate (massimo ordine possibile)
	doublereal 
	dPredDer(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;
   
	doublereal 
	dPredState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredDerAlg(const doublereal& dXm1,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredStateAlg(const doublereal& dXm1,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;
};

/* NostroMetodo - end */


/* Hope - begin */

class HopeSolver : 
	public Step2Integrator 
{
protected:
	DriveOwner Rho;
	DriveOwner AlgebraicRho;
   
	flag fStep;
   
	doublereal a[2][2];
	doublereal b[2][2];
   
	doublereal mp[2];
	doublereal np[2];
   
public:
	HopeSolver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho);

	~HopeSolver(void);

protected:
	void SetCoef(doublereal dT,
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal 
	dPredictDerivative(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;

	doublereal 
	dPredictState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const;      

	// Nota: usa predizione cubica per le derivate (massimo ordine possibile)
	doublereal 
	dPredDer(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredDerAlg(const doublereal& dXm1,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredStateAlg(const doublereal& dXm1,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;
};

/* Hope - end */

#endif /* STEPSOL_H */
