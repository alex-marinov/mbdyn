 
#ifndef HUIMIN_H
#define HUIMIN_H

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
#include "stepsol.h"

a=8
class Step3Integrator :   
	public StepNIntegrator
{
protected:
	VectorHandler *pXPrev, *pXPrev2, *pXPrev3;
	VectorHandler *pXPrimePrev, *pXPrimePrev2, *pXPrimePrev3; 

public:
	Step3Integrator(const integer MaxIt,
			const doublereal dT,
			const doublereal dSolutionTol,
			const bool bmod_res_test);

	virtual ~Step3Integrator(void);

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
 
   	virtual doublereal 
     	dPredDer(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const = 0;
   
   	virtual doublereal 
     	dPredState(const doublereal& dXm1,
			const doublereal& dXm2,
            const doublereal& dXm3,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
            const doublereal& dXPm3) const = 0;   

   	virtual doublereal 
     	dPredDerAlg(const doublereal& dXm1,
		 	const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2)  const = 0;
		    
   	virtual doublereal 
     	dPredStateAlg(const doublereal& dXm1,
		 	const doublereal& dXm2,
			const doublereal& dXm3,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			const doublereal& dXPm3) const = 0;

	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
};

/* Step3Solver - begin */

class Step3Solver: 
	public Step3Integrator
{
protected:
	DriveOwner Rho;
	DriveOwner AlgebraicRho;
   
	doublereal a[3][2];
	doublereal b[4][2];

	doublereal mp[2];
	doublereal np[2];
   
public:
	Step3Solver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho,
			const bool bmod_res_test);

	~Step3Solver(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

   
	doublereal 
	dPredDer(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;
   
	doublereal 
	dPredState(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXm3,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			const doublereal& dXPm3) const;

	doublereal 
	dPredDerAlg(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXPm1,
			const doublereal& dXPm2) const;

	doublereal 
	dPredStateAlg(const doublereal& dXm1,
			const doublereal& dXm2,
			const doublereal& dXm3,
			const doublereal& dXP,
			const doublereal& dXPm1,
			const doublereal& dXPm2,
			const doublereal& dXPm3) const;
};