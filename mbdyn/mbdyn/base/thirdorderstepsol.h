#ifndef THIRD_ORDER_STEPSOL_H

#define THIRD_ORDER_STEPSOL_H

#include "stepsol.h"
#include "spmapmh.h"

class ThirdOrderIntegrator :  
	public ImplicitStepIntegrator
{
private:
	VectorHandler *pXPrev;
	VectorHandler *pXPrimePrev; 
	DriveOwner Rho;
	
	doublereal dT;
	doublereal rho;
	doublereal theta;
	doublereal w[3];
	doublereal jx[2][2];
	doublereal jxp[2][2];
	doublereal m0, m1, n0, n1;
	
	bool bAdvanceCalledFirstTime;
	MyVectorHandler Res1, Res2;
	std::vector<bool> EqIsAlgebraic, EqIsDifferential;
	SpMapMatrixHandler Jacxi_xp, Jacxi_x, Jac_xp, Jac_x;
	MatrixHandler *pJacxi_xp, *pJacxi_x, *pJac_xp, *pJac_x;

protected:
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
			const DriveCaller* pRho,
			const bool bmod_res_test);

	virtual ~ThirdOrderIntegrator(void);

	virtual void Residual(VectorHandler* pRes) const;

	virtual void Jacobian(MatrixHandler* pJac) const;
	
	virtual void Update(const VectorHandler* pSol) const;

	virtual void SetDriveHandler(const DriveHandler* pDH);

	/* scale factor for tests */
//	virtual doublereal TestScale(const NonlinearSolverTest *pTest) const;

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
			integer& EffIter,
			doublereal& Err,
			doublereal& SolErr);

protected:
	virtual void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);
			
	virtual void Predict(void);

};

#endif /* THIRD_ORDER_STEPSOL_H */
