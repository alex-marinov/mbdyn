/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
