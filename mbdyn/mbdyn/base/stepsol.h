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

public:
	StepIntegrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const integer stp) 
	: pDM(NULL),
	DofIterator(),
	outputPred(false),
	MaxIters(MaxIt),
	dTol(dT),
	dSolTol(dSolutionTol),
	steps(stp)
	{ };
	
	void SetDataManager(DataManager* pDatMan) {
		pDM = pDatMan;
		DofIterator = pDM->GetDofIterator();
		return;
	};
	
	
	virtual ~StepIntegrator(void) { };
		
	virtual integer GetIntegratorStepSize(void) const{
		return steps+1;
	}; 
	
	virtual void OutputTypes(const bool fpred)
	{
		outputPred  = fpred;
	};
	
	
	virtual void SetDriveHandler(const DriveHandler* pDH) = 0;

	virtual doublereal Advance(const doublereal TStep, 
				const doublereal dAlph, 
				const StepChange StType,
				NonlinearSolver* pNLS, 
				std::deque<MyVectorHandler*>& qX,
	 			std::deque<MyVectorHandler*>& qXPrime,
				integer& EffIter
#ifdef MBDYN_X_CONVSOL
				, doublereal& SolErr
#endif /* MBDYN_X_CONVSOL */
				) = 0;
};


class DerivativeSolver: 
	public StepIntegrator,
	public NonlinearProblem
{
private:
	doublereal dCoef;
	VectorHandler *pXCurr;
	VectorHandler *pXPrimeCurr; 

public:
	DerivativeSolver(const doublereal Tl, 
			const doublereal dSolTl, 
			const doublereal dC,
			const integer iMaxIt) 
	:StepIntegrator(iMaxIt, Tl, dSolTl, 1),
	dCoef(dC) { };
	
	~DerivativeSolver(void){ };			

	
	void SetDriveHandler(const DriveHandler* /* pDH */ ) {
      		NO_OP;
   	};
 
		
 	doublereal Advance(const doublereal TStep, 
				const doublereal /* dAph */, 
				const StepChange /* StType */,
				NonlinearSolver* pNLS, 
				std::deque<MyVectorHandler*>& qX,
	 			std::deque<MyVectorHandler*>& qXPrime,
				integer& EffIter
#ifdef MBDYN_X_CONVSOL
				, doublereal& SolErr
#endif /* MBDYN_X_CONVSOL */
			);
 	void Residual(VectorHandler* pRes) const;

	void Jacobian(MatrixHandler* pJac) const;
	
	void Update(const VectorHandler* pSol) const;

 
};


/* classe di base per gli integratori del second'ordine */ 
class Step2Integrator :   
	public StepIntegrator,
	public NonlinearProblem
{
protected:
	doublereal db0Differential;
	doublereal db0Algebraic;
	VectorHandler  *pXCurr, *pXPrev, *pXPrev2;
	VectorHandler  *pXPrimeCurr, *pXPrimePrev, *pXPrimePrev2; 

public:
	Step2Integrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol) 
	:StepIntegrator(MaxIt, dT, dSolutionTol, 2),
	db0Differential(0.),
	db0Algebraic(0.),
	pXCurr(NULL),
	pXPrev(NULL),
	pXPrev2(NULL),
	pXPrimeCurr(NULL),
	pXPrimePrev(NULL),
	pXPrimePrev2(NULL)
	{ };
	

	virtual ~Step2Integrator(void) { };



	virtual doublereal Advance(const doublereal TStep, 
				const doublereal dAlph, 
				const StepChange StType,
				NonlinearSolver* pNLS, 
				std::deque<MyVectorHandler*>& qX,
	 			std::deque<MyVectorHandler*>& qXPrime,
				integer& EffIter
#ifdef MBDYN_X_CONVSOL
				, doublereal& SolErr
#endif /* MBDYN_X_CONVSOL */
				);

	virtual void Residual(VectorHandler* pRes) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
	virtual void Update(const VectorHandler* pSol) const;


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




class CrankNicholsonSolver: 
	public Step2Integrator
{

public:
	
	CrankNicholsonSolver(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt);


	~CrankNicholsonSolver(void) {};


	void SetDriveHandler(const DriveHandler* /* pDH */ ) {
      		NO_OP;
   	};

   
protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange /* NewStep */);
  
   	doublereal 
     	dPredictDerivative(const doublereal& /* dXm1 */,
			const doublereal& /* dXm2 */,
			const doublereal& dXPm1,
			const doublereal& /* dXPm2 */,
			DofOrder::Order o = DofOrder::DIFFERENTIAL) const {
   				return dXPm1;
			};
   
	doublereal 
	dPredictState(const doublereal& dXm1,
		   const doublereal& /* dXm2 */,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& /* dXPm2 */,
		   DofOrder::Order o = DofOrder::DIFFERENTIAL) const {
		   	if (o == DofOrder::ALGEBRAIC) {
				return db0Differential*(dXP+dXPm1);     
			} /* else if (o == DofOrder::DIFFERENTIAL) */   
			return dXm1+db0Differential*(dXP+dXPm1);
		};
   
	// Nota: usa predizione lineare per le derivate (massimo ordine possibile)
	doublereal 
	dPredDer(const doublereal& /* dXm1 */ ,
	      const doublereal& /* dXm2 */ ,
	      const doublereal& dXPm1,
	      const doublereal& /* dXPm2 */ ) const {
		return dXPm1;
	      };
   
	doublereal 
	dPredState(const doublereal& dXm1,
		const doublereal& /* dXm2 */ ,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& /* dXPm2 */ ) const {
		   return dXm1+db0Differential*(dXP+dXPm1);
		};
   
	doublereal 
	dPredDerAlg(const doublereal& /* dXm1 */ ,
		 const doublereal& dXPm1,
		 const doublereal& /* dXPm2 */ ) const {
		    return dXPm1;
		 };
   
	doublereal 
	dPredStateAlg(const doublereal& /* dXm1 */ ,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& /* dXPm2 */ ) const {
		      return db0Differential*(dXP+dXPm1);
		   };     

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
   
   ~MultistepSolver(void) { 
      NO_OP; 
   };

protected:
   
  void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
   
  void SetDriveHandler(const DriveHandler* pDH) {
      Rho.pGetDriveCaller()->SetDrvHdl(pDH);
      AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
   };
   
   
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
	      const doublereal& dXPm2) const {
		 return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
	      };
   
   doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const {
		   return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
		     +b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1+b[2][DIFFERENTIAL]*dXPm2;
		};
   
  doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const {
		    return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		 };
   
   doublereal 
     dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& dXPm2) const {
		      return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2
			-a[1][ALGEBRAIC]*dXm1;
		   };
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
   
   ~HopeSolver(void) { 
      NO_OP; 
   };

protected:
   
   void SetCoef(doublereal dT,
			doublereal dAlpha,
			enum StepChange NewStep);
   
   void SetDriveHandler(const DriveHandler* pDH) {
      Rho.pGetDriveCaller()->SetDrvHdl(pDH);
      AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
   };

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
	      const doublereal& dXPm2) const {
		 return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
	      };
   
   doublereal 
     dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& /* dXPm2 */ ) const {
		   if (fStep) {
		      return dXm1+b[0][ALGEBRAIC]*(dXP+dXPm1);
		   } else {
		      return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2
			+b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1;
		   }
		};
   
   doublereal 
     dPredDerAlg(const doublereal& dXm1,
		 const doublereal& dXPm1,
		 const doublereal& dXPm2) const {
		    return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		 };
   
   doublereal 
     dPredStateAlg(const doublereal& dXm1,
		   const doublereal& dXP,
		   const doublereal& dXPm1,
		   const doublereal& /* dXPm2 */ ) const {
		      if (fStep) {
			 return b[0][ALGEBRAIC]*(dXP+dXPm1);
		      } else {
			 return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1
						   -a[1][ALGEBRAIC]*dXm1;
		      }
		   };
};

/* Hope - end */


#endif /* STEPSOL_H */
