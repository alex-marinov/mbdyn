#ifndef THIRD_ORDER_STEPSOL_H
#define THIRD_ORDER_STEPSOL_H

#include "stepsol.h"

class ThirdOrderIntegrator :  
	public ImplicitStepIntegrator
{
private:
	VectorHandler *pXCurr;
	VectorHandler *pXPrimeCurr; 
	VectorHandler *pXPrev;
	VectorHandler *pXPrimePrev; 
	DriveOwner Rho;
	
	doublereal rho;
	doublereal theta;
	doublereal w[3];
	doublereal jx[2][2];
	doublereal jxp[2][2];

public:
	ThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const DriveCaller* pRho);

	virtual ~ThirdOrderIntegrator(void);

	virtual void Residual(VectorHandler* pRes) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
	virtual void Update(const VectorHandler* pSol) const;

	/* scale factor for tests */
#ifdef __HACK_SCALE_RES__
	virtual doublereal TestScale(const VectorHandler *pScale) const;
#else /* ! __HACK_SCALE_RES__ */
	virtual doublereal TestScale(void) const;
#endif /* ! __HACK_SCALE_RES__ */

	virtual void
	EvalProd(doublereal Tau, const VectorHandler& f0,
			const VectorHandler& w, VectorHandler& z) const;

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
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
			
	virtual void Predict(void);

};

#endif /* THIRD_ORDER_STEPSOL_H */
