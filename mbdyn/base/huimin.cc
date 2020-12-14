/* Step3Integrator - begin */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "schurdataman.h"
#include "external.h"
#include "ls.h"
#include "solver.h"
#include "invsolver.h"
#include "stepsol.h"
#include "huimin.h"

Step3Integrator::Step3Integrator(const integer MaxIt,
		const doublereal dT,
		const doublereal dSolutionTol,
		const bool bmod_res_test)
: StepNIntegrator(MaxIt, dT, dSolutionTol, 3, bmod_res_test),
pXPrev(NULL),
pXPrev2(NULL),
pXPrev3(NULL),
pXPrimePrev(NULL),
pXPrimePrev2(NULL),
pXPrimePrev3(NULL)
{
	NO_OP;
}

Step3Integrator::~Step3Integrator(void)
{
	NO_OP;
}


void Step3Integrator::PredictDof(const int DCount,
	const DofOrder::Order Order,
	const VectorHandler* const pSol) const
{
	if (Order == DofOrder::DIFFERENTIAL) {
		doublereal dXnm1 = pXPrev->operator()(DCount);
		doublereal dXnm2 = pXPrev2->operator()(DCount);
        doublereal dXnm3 = pXPrev3->operator()(DCount);
		doublereal dXPnm1 = pXPrimePrev->operator()(DCount);
		doublereal dXPnm2 =
			pXPrimePrev2->operator()(DCount);
        doublereal dXPnm3 =
			pXPrimePrev3->operator()(DCount);
		doublereal dXPn = dPredDer(dXnm1, dXnm2,
				dXPnm1, dXPnm2);
		doublereal dXn = dPredState(dXnm1, dXnm2, dXnm3
				dXPn, dXPnm1, dXPnm2, dXPnm3);

		pXPrimeCurr->PutCoef(DCount, dXPn);
		pXCurr->PutCoef(DCount, dXn);

	} else if (Order == DofOrder::ALGEBRAIC) {
		doublereal dXnm1 = pXPrev->operator()(DCount);
		doublereal dXnm2 = pXPrev2->operator()(DCount);
		doublereal dXnm3 = pXPrev3->operator()(DCount);
		doublereal dXInm1 =
			pXPrimePrev->operator()(DCount);
		doublereal dXInm2 =
			pXPrimePrev2->operator()(DCount);
		doublereal dXInm3 =
			pXPrimePrev3->operator()(DCount);

		
		doublereal dXn = dPredDerAlg(dXInm1, dXInm2,
				dXnm1, dXnm2);
		doublereal dXIn = dPredStateAlg(dXInm1, dXInm2, dXInm3,
				dXn, dXnm1, dXnm2, dXnm3);

		pXCurr->PutCoef(DCount, dXn);
		pXPrimeCurr->PutCoef(DCount, dXIn);

	} else {
		silent_cerr("Step3Integrator::"
			"PredictDof(): "
			"unknown order for local dof "
			<< DCount << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
Step3Integrator::Predict(void)
{
   	DEBUGCOUTFNAME("Step3Integrator::Predict");
   	ASSERT(pDM != NULL);
	UpdateLoop(this, &Step3Integrator::PredictDof);
}

doublereal
Step3Integrator::Advance(Solver* pS,
		const doublereal TStep,
		const doublereal dAph, const StepChange StType,
		std::deque<MyVectorHandler*>& qX,
	 	std::deque<MyVectorHandler*>& qXPrime,
		MyVectorHandler*const pX,
		MyVectorHandler*const pXPrime,
		integer& EffIter,
		doublereal& Err,
		doublereal& SolErr)
{
	ASSERT(pDM != NULL);
	pXCurr  = pX;
	pXPrev  = qX[0];
	pXPrev2 = qX[1];
    pXPrev3 = qx[2];

	pXPrimeCurr  = pXPrime;
	pXPrimePrev  = qXPrime[0];
	pXPrimePrev2 = qXPrime[1];
    pXPrimePrev3 = qXPrime[2];

	/* prediction */
	SetCoef(TStep, dAph, StType);
	Predict();
	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      	pDM->AfterPredict();

#ifdef DEBUG
	integer iNumDofs = pDM->iGetNumDofs();
	if (outputPred) {
		std::cout << "After prediction, time=" << pDM->dGetTime() << std::endl;
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
	pS->pGetNonlinearSolver()->Solve(this, pS, MaxIters, dTol,
    			EffIter, Err, dSolTol, SolErr);

	/* if it gets here, it surely converged */
	pDM->AfterConvergence();

	return Err;
}

/* Step3Integrator - end */

/* Step3Solver - begin */

Step3Solver::Step3Solver(const doublereal Tl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test)
:Step3Integrator(iMaxIt, Tl, dSolTl, bmod_res_test),
Rho(pRho), AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != NULL);
	ASSERT(pAlgRho != NULL);
}

Step3Solver::~Step3Solver(void)
{
	NO_OP;
}

void
Step3Solver::SetDriveHandler(const DriveHandler* pDH)
{
	Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Step3Solver::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	doublereal dRho = Rho.dGet();
	doublereal dAlgebraicRho = AlgebraicRho.dGet();


	mp[0] = -6.*dAlpha*dAlpha*(1.+dAlpha)/dT;
	mp[1] = -mp[0];
	np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
	np[1] = dAlpha*(2.+3.*dAlpha);

	a[0][DIFFERENTIAL] = 3.*(2.*dRho*dRho-9.*dRho+5.)/(dRho*dRho-5.*dRho+10.);
	a[1][DIFFERENTIAL] = -3.*(5.*dRho*dRho-9.*dRho+2.)/(dRho*dRho-5.*dRho+10.);
	a[2][DIFFERENTIAL] = (10.*dRho*dRho-5.*dRho+1.)/(dRho*dRho-5.*dRho+10.);
	b[0][DIFFERENTIAL] = dT*(6./((1.+dRho)*(dRho*dRho-5.*dRho+10.)));
	b[1][DIFFERENTIAL] = dT*(18.*dRho/((1.+dRho)*(dRho*dRho-5.*dRho+10.)));
	b[2][DIFFERENTIAL] = dT*(18.*dRho*dRho/((1.+dRho)*(dRho*dRho-5.*dRho+10.)));
	b[3][DIFFERENTIAL] = dT*(6.*dRho*dRho*dRho/((1.+dRho)*(dRho*dRho-5.*dRho+10.)));

	DEBUGCOUT("Predict()" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << dRho << std::endl 
			<< "a1    = " << a[0][DIFFERENTIAL] << std::endl
			<< "a2    = " << a[1][DIFFERENTIAL] << std::endl
			<< "a3    = " << a[2][DIFFERENTIAL] << std::endl
			<< "b0    = " << b[0][DIFFERENTIAL] << std::endl
			<< "b1    = " << b[1][DIFFERENTIAL] << std::endl
			<< "b2    = " << b[2][DIFFERENTIAL] << std::endl
			<< "b3    = " << b[3][DIFFERENTIAL] << std::endl);

	if (dAlgebraicRho != dRho) {

		a[0][ALGEBRAIC] = 3.*(2.*dAlgebraicRho*dAlgebraicRho-9.*dAlgebraicRho+5.)/(ddAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.);
		a[1][ALGEBRAIC] = -3.*(5.*dAlgebraicRho*dAlgebraicRho-9.*dAlgebraicRho+2.)/(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.);
		a[2][ALGEBRAIC] = (10.*dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+1.)/(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.);
		b[0][ALGEBRAIC] = dT*(6./((1.+dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.)));
		b[1][ALGEBRAIC] = dT*(18.*dAlgebraicRho/((1.+dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.)));
		b[2][ALGEBRAIC] = dT*(18.*dAlgebraicRho*dAlgebraicRho/((1.+dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.)));
		b[3][ALGEBRAIC] = dT*(6.*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho/((1.+dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho-5.*dAlgebraicRho+10.)));

	} else {
		
		a[0][ALGEBRAIC] = a[0][DIFFERENTIAL];
		a[1][ALGEBRAIC] = a[1][DIFFERENTIAL];
		a[2][ALGEBRAIC] = a[2][DIFFERENTIAL];
		b[0][ALGEBRAIC] = b[0][DIFFERENTIAL];
		b[1][ALGEBRAIC] = b[1][DIFFERENTIAL];
		b[2][ALGEBRAIC] = b[2][DIFFERENTIAL];
		b[3][ALGEBRAIC] = b[3][DIFFERENTIAL];

	}

	DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "Asymptotic rho =" << dAlgebraicRho << std::endl 
			<< "a1    = " << a[0][ALGEBRAIC] << std::endl
			<< "a2    = " << a[1][ALGEBRAIC] << std::endl
			<< "a3    = " << a[2][ALGEBRAIC] << std::endl
			<< "b0    = " << b[0][ALGEBRAIC] << std::endl
			<< "b1    = " << b[1][ALGEBRAIC] << std::endl
			<< "b2    = " << b[2][ALGEBRAIC] << std::endl
			<< "b3    = " << b[3][ALGEBRAIC] << std::endl);

	db0Differential = b[0][DIFFERENTIAL];
	db0Algebraic = b[0][ALGEBRAIC];
}



doublereal
Step3Solver::dPredDer(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
}

doublereal
Step3Solver::dPredState(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXm3,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2,
		const doublereal& dXPm3) const
{
	return a[0][DIFFERENTIAL]*dXm1+a[1][DIFFERENTIAL]*dXm2+a[2][DIFFERENTIAL]*dXm3
		+b[0][DIFFERENTIAL]*dXP+b[1][DIFFERENTIAL]*dXPm1+b[2][DIFFERENTIAL]*dXPm2+b[3][DIFFERENTIAL]*dXPm3;
}

doublereal
Step3Solver::dPredDerAlg(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXPm1,
		const doublereal& dXPm2) const
{
	return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
}

doublereal
Step3Solver::dPredStateAlg(const doublereal& dXm1,
		const doublereal& dXm2,
		const doublereal& dXm3,
		const doublereal& dXP,
		const doublereal& dXPm1,
		const doublereal& dXPm2,
		const doublereal& dXPm3) const
{
	return a[0][ALGEBRAIC]*dXm1+a[1][ALGEBRAIC]*dXm2+a[2][ALGEBRAIC]*dXm3
		+b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2+b[3][ALGEBRAIC]*dXPm3;
}

/* Step3Solver - end */