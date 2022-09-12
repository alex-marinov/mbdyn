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
 
#ifndef STEPSOL_IMPL_H
#define STEPSOL_IMPL_H

#include "stepsol_tpl.h"

class CrankNicolsonIntegrator: 
	public tplStepNIntegrator<1>
{
public:
	CrankNicolsonIntegrator(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt,
			const bool bmod_res_test);

	~CrankNicolsonIntegrator(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	/* Note: uses linear prediction for derivatives 
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[1],
	      const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[1],
		 const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[1],
		   const doublereal dXP0mN[2]) const;
};


class ImplicitEulerIntegrator: 
	public tplStepNIntegrator<1>
{
public:
	ImplicitEulerIntegrator(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt,
			const bool bmod_res_test);

	~ImplicitEulerIntegrator(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	/* Note: uses linear prediction for derivatives 
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[1],
	      const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[1],
		 const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[1],
		   const doublereal dXP0mN[2]) const;
};


/* 2-step multistep (nostro metodo) - begin */

class Multistep2Solver: 
	public tplStepNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

public:
	Multistep2Solver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho,
			const bool bmod_res_test);

	~Multistep2Solver(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	/* Note: uses cubic prediction for derivatives
	 * (highest possible order) */

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredDer(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredState(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[2],
			const doublereal dXP1mN[2]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;
};

/* 2-step multistep (nostro metodo) - end */


/* Hope - begin */

class HopeSolver : 
	public tplStepNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;
   
	bool m_bStep;
   
public:
	HopeSolver(const doublereal Tl, 
		const doublereal dSolTol, 
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~HopeSolver(void);

protected:
	void SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	/* Note: uses cubic prediction for derivatives
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredState(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredDerAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredStateAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;
};

/* Hope - end */

// The hybrid step integrator allows elements to choose different integration methods for each degree of freedom.
// Simulation entities using this integrator must overwrite SimulationEntity::GetStepIntegrator.
// In that case it is required to call DataManager::dGetStepIntegratorCoef in order to obtain the correct value for dCoef.

class HybridStepIntegrator: public ImplicitStepIntegrator
{
public:
     HybridStepIntegrator(const SolverBase::StepIntegratorType eDefaultIntegrator,
                          const doublereal dTol,
                          const doublereal dSolutionTol,
                          const integer iMaxIterations,
                          const DriveCaller* pRho,
                          const DriveCaller* pAlgRho,
                          const bool bModResTest);

     virtual ~HybridStepIntegrator();

     virtual void
     SetDataManager(DataManager* pDataMan) override;

     virtual void
     SetDriveHandler(const DriveHandler* pDH) override;

     virtual doublereal
     dGetCoef(unsigned int iDof) const override;

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
             doublereal& SolErr) override;

     virtual void
     Update(const VectorHandler* pSol) const override;

     virtual void
     Residual(VectorHandler* pRes, VectorHandler* pAbsRes=0) const override;

     virtual void
     Jacobian(MatrixHandler* pJac) const override;

     virtual void
     Jacobian(VectorHandler* pJac, const VectorHandler* pY) const override;

private:
     void
     SetSolution(std::deque<VectorHandler*>& qX,
                 std::deque<VectorHandler*>& qXPrime,
                 MyVectorHandler* pX,
                 MyVectorHandler* pXPrime);

     struct IntegratorItem {
          IntegratorItem(SolverBase::StepIntegratorType eType,
                         tplStepNIntegratorBase* pInteg)
               :eType(eType),
                pInteg(pInteg) {
          }

          SolverBase::StepIntegratorType eType;
          std::unique_ptr<tplStepNIntegratorBase> pInteg;
     };

     void
     Predict();

     void
     SetCoef(doublereal dT, doublereal dAlpha, StepChange NewStep);

     tplStepNIntegratorBase*
     pAllocateStepIntegrator(SolverBase::StepIntegratorType eStepIntegrator);

     typedef void (tplStepNIntegratorBase::*UpdateFunctionType)(const int DCount,
                                                                const DofOrder::Order Order,
                                                                const VectorHandler* const pSol) const;

     inline void
     UpdateLoop(UpdateFunctionType pfnUpdateFunc, const VectorHandler* const pSol = nullptr) const;

     std::vector<IntegratorItem> rgIntegItems;
     std::array<tplStepNIntegratorBase*, SolverBase::INT_COUNT> rgIntegPtr;
     tplStepNIntegratorBase* pDefaultInteg;
     const SolverBase::StepIntegratorType eDefaultIntegrator;
     DriveOwner m_Rho;
     DriveOwner m_AlgebraicRho;
};

#endif /* STEPSOL_IMPL_H */
