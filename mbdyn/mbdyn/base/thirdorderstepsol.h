#ifndef THIRD_ORDER_STEPSOL_H

#define THIRD_ORDER_STEPSOL_H

#include <solver.h>
#include <stepsol.h>
#include <spmapmh.h>

class ThirdOrderIntegrator :  
	public ImplicitStepIntegrator
{
protected:
	VectorHandler *pXPrev;
	VectorHandler *pXPrimePrev; 

	doublereal dT;
	doublereal rho;
	doublereal theta;
	doublereal w0, w1, w2;
	doublereal jx11, jx12, jx21, jx22;
	doublereal m0, m1, n0, n1;
	
	bool bAdvanceCalledFirstTime;
	mutable MyVectorHandler Restmp;
	std::vector<bool> EqIsAlgebraic, EqIsDifferential;
	SpMapMatrixHandler Jacxi_xp, Jacxi_x, Jac_xp, Jac_x;
	MatrixHandler *pJacxi_xp, *pJacxi_x, *pJac_xp, *pJac_x;

	void PredictDof_for_AfterPredict(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;
	void RealPredictDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;
	void UpdateDof(const int DCount,
		const DofOrder::Order Order,
		const VectorHandler* const pSol = 0) const;
public:
	ThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const bool bmod_res_test);

	virtual ~ThirdOrderIntegrator(void);

	virtual void Residual(VectorHandler* pRes) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
	virtual void Update(const VectorHandler* pSol) const;

	/* scale factor for tests */
//	virtual doublereal TestScale(const NonlinearSolverTest *pTest) const;

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
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep) = 0;
	virtual void Predict(void);

};

class TunableThirdOrderIntegrator :  
	public ThirdOrderIntegrator
{
private:
	DriveOwner Rho;

public:
	TunableThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const bool bmod_res_test);

	virtual ~TunableThirdOrderIntegrator(void);
	void SetDriveHandler(const DriveHandler* pDH);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
};

class AdHocThirdOrderIntegrator :  
	public ThirdOrderIntegrator
{
public:
	AdHocThirdOrderIntegrator(const doublereal dT, 
			const doublereal dSolutionTol, 
			const integer iMaxIt,
			const bool bmod_res_test);

	virtual ~AdHocThirdOrderIntegrator(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
};

#endif /* THIRD_ORDER_STEPSOL_H */
