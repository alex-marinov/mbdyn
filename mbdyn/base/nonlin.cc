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
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "solver.h"
#include "nonlin.h"  
#ifdef USE_MPI
#include "mbcomm.h"
#include "schsolman.h"
#endif /* USE_MPI */
#include "dofown.h"
#include "output.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <unistd.h>

/* NonlinearSolverTest - begin */

NonlinearSolverTest::~NonlinearSolverTest(void)
{
	NO_OP;
}

doublereal
NonlinearSolverTest::MakeTest(Solver *pS, const integer& Size,
		const VectorHandler& Vec, bool bResidual,
		doublereal dScaleAlgEqu, doublereal* pTestDiff)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNorm::MakeTest");

   	doublereal dTest = 0.;

   	if (pTestDiff) {
   		*pTestDiff = 0.;
   	}

#ifdef USE_SCHUR
	/* Only residual test is parallelized; the master node
	 * always knows the entire solution */
	ASSERT(pS != NULL);
	SchurSolutionManager *pSSM = dynamic_cast<SchurSolutionManager *>(pS->pGetSolutionManager());
	if (pSSM) {
		SchurDataManager *pSDM = 
			dynamic_cast<SchurDataManager *>(pS->pGetDataManager());
		ASSERT(pSDM);

		integer iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
		integer *pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);

#if 0
		silent_cout("NonlinearSolverTest::MakeTest("
				<< MBDynComm.Get_rank() << ") "
				"iNumLocDofs=" << iNumLocDofs << std::endl);
#endif
	
		if (bResidual) {
			/*
			 * Chiama la routine di comunicazione per la trasmissione 
			 * del residuo delle interfacce
			 */
			pSSM->StartExchIntRes();

			/* calcola il test per i dofs locali */
			for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
				TestOne(dTest, Vec, pLocDofs[iCnt]);
			}

			/* collect contributions from other nodes,
			 * plus that of the interface; merge them according
			 * to the NonlinearSolverTest type */
			pSSM->ComplExchIntRes(dTest, this);

		} else {
			/*
			 * Chiama la routine di comunicazione per la trasmissione 
			 * del residuo delle interfacce
			 */
			pSSM->StartExchIntSol();

			/* calcola il test per i dofs locali */
			for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
				TestOne(dTest, Vec, pLocDofs[iCnt]);
			}

			/* collect contributions from other nodes,
			 * plus that of the interface; merge them according
			 * to the NonlinearSolverTest type */
			pSSM->ComplExchIntSol(dTest, this);
		}

	} else
#endif // USE_SCHUR
	{
		ASSERT(Vec.iGetSize() == Size);
		const DataManager* const pDM = pS->pGetDataManager();

 	  	for (int iCntp1 = 1; iCntp1 <= Size; iCntp1++) {
 	  		const DofOrder::Order order = pDM->GetEqType(iCntp1);
 	  		const doublereal dCoef = order == DofOrder::DIFFERENTIAL ? 1. : dScaleAlgEqu;

			TestOne(dTest, Vec, iCntp1, dCoef);

			if (pTestDiff && order == DofOrder::DIFFERENTIAL) {
				TestOne(*pTestDiff, Vec, iCntp1, dCoef);
			}
		}
	}

	if (pTestDiff) {
		*pTestDiff = TestPost(*pTestDiff);
	}

	return TestPost(dTest);
}

doublereal
NonlinearSolverTest::TestPost(const doublereal& dRes) const
{
	if (!std::isfinite(dRes)) {
		throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
	}

	return dRes;
}

const doublereal dOne = 1.;

const doublereal&
NonlinearSolverTest::dScaleCoef(const integer& iIndex) const
{
	return ::dOne;
}

/* NonlinearSolverTest - end */

/* NonlinearSolverTestNone - begin */

NonlinearSolverTest::Type NonlinearSolverTestNone::GetType() const
{
     return NONE;
}

doublereal
NonlinearSolverTestNone::MakeTest(Solver *pS, integer Size, 
		const VectorHandler& Vec, bool bResidual,
		doublereal* pTestDiff)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNone::MakeTest");

   	if (pTestDiff) {
   		*pTestDiff = 0.;
   	}

	return 0.;
}

void
NonlinearSolverTestNone::TestOne(doublereal& dRes, 
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	dRes = 0.;
}

void
NonlinearSolverTestNone::TestMerge(doublereal& dResCurr, 
		const doublereal& dResNew) const
{
	dResCurr = 0.;
}

/* NonlinearSolverTestNone - end */

/* NonlinearSolverTestNorm - begin */

NonlinearSolverTest::Type NonlinearSolverTestNorm::GetType() const
{
     return NORM;
}

void
NonlinearSolverTestNorm::TestOne(doublereal& dRes, 
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestNorm::TestMerge(doublereal& dResCurr, 
		const doublereal& dResNew) const
{
	dResCurr += dResNew;
}

doublereal
NonlinearSolverTestNorm::TestPost(const doublereal& dRes) const
{
	/* va qui perche' non posso fare sqrt() su !isfinite() */
	if (!std::isfinite(dRes)) {      
		throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
	}

   	return sqrt(dRes);
}

/* NonlinearSolverTestNorm - end */

/* NonlinearSolverTestRelNorm - begin */

VectorHandler*
NonlinearSolverTestRelNorm::GetAbsRes() {
	return &AbsRes;
}

NonlinearSolverTest::Type NonlinearSolverTestRelNorm::GetType() const
{
     return RELNORM;
}

doublereal
NonlinearSolverTestRelNorm::MakeTest(Solver *pS, const integer &Size, 
		const VectorHandler& Vec, bool bResidual, doublereal dScaleAlgEqu,
		doublereal* pTestDiff)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestRelNorm::MakeTest");

   	/* get norm for absolute residual vector */
	doublereal abs_res_test = NonlinearSolverTest::MakeTest(pS, Size, AbsRes, bResidual, dScaleAlgEqu, pTestDiff);

	/* get norm for residual vector */
	doublereal res_test = NonlinearSolverTest::MakeTest(pS, Size, Vec, bResidual, dScaleAlgEqu, pTestDiff);

	if ( abs_res_test == 0)
		return 0.;

	return res_test/abs_res_test;
}


void
NonlinearSolverTestRelNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestRelNorm::TestMerge(doublereal& dResCurr,
		const doublereal& dResNew) const
{
	dResCurr += dResNew;
}

doublereal
NonlinearSolverTestRelNorm::TestPost(const doublereal& dRes) const
{
	/* va qui perche' non posso fare sqrt() su !isfinite() */
	if (!std::isfinite(dRes)) {
		throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
	}

	return sqrt(dRes);
}

/* NonlinearSolverTestNorm - end */

/* NonlinearSolverTestSepNorm */

VectorHandler*
NonlinearSolverTestSepNorm::GetAbsRes() {
	return &AbsRes;
}

std::map<OutputHandler::Dimensions, std::set<integer>>* 
NonlinearSolverTestSepNorm::GetDimMap() { 
	return &MapOfDimensionIndices; 
};

NonlinearSolverTest::Type NonlinearSolverTestSepNorm::GetType() const
{
     return SEPNORM;
}

doublereal
NonlinearSolverTestSepNorm::MakeTest(Solver *pS, const integer &Size, 
		const VectorHandler& Vec, bool bResidual, doublereal dScaleAlgEqu,
		doublereal* pTestDiff)
{

	std::vector<doublereal> testsVector;
	std::vector<doublereal> testDiffsVector;
	std::vector<doublereal> abs_dTestVector;
	std::vector<doublereal> dTestVector;

   	for (auto it = MapOfDimensionIndices.begin(); it != MapOfDimensionIndices.end(); ++it) {

		doublereal dTest = 0.;
		doublereal abs_dTest = 0.;

		doublereal pTestDiff_temp = 0.;
		doublereal abs_pTestDiff_temp = 0.;


		const DataManager* const pDM = pS->pGetDataManager();
		for (auto i = (*it).second.begin(); i != (*it).second.end(); ++i) {

			const DofOrder::Order order = pDM->GetEqType(*i);
 	  		const doublereal dCoef = order == DofOrder::DIFFERENTIAL ? 1. : dScaleAlgEqu;

			TestOne(dTest, Vec, *i, dCoef);
			TestOne(abs_dTest, AbsRes, *i, dCoef);

			if (pTestDiff && order == DofOrder::DIFFERENTIAL) {
				TestOne(pTestDiff_temp, Vec, *i, dCoef);
				TestOne(abs_pTestDiff_temp, AbsRes, *i, dCoef);
			}

		}

		if (pTestDiff) {
			pTestDiff_temp = TestPost(pTestDiff_temp);
			abs_pTestDiff_temp = TestPost(abs_pTestDiff_temp);

			testDiffsVector.push_back(pTestDiff_temp/abs_pTestDiff_temp);
		}

		dTest = TestPost(dTest);
		dTestVector.push_back(dTest);
		abs_dTest = TestPost(abs_dTest);
		abs_dTestVector.push_back(abs_dTest);
	}

	doublereal abs_dTest_max = *max_element(abs_dTestVector.begin(), abs_dTestVector.end());

	doublereal eps1, eps2;
	eps1 = 1E-1;
	eps2 = 1E-5;

	for ( int i = 0; i < dTestVector.size(); i++) {
		doublereal dTest = dTestVector[i];
		doublereal abs_dTest = abs_dTestVector[i];

		if ((abs_dTest != 0) && !((dTest/abs_dTest > eps1) && (abs_dTest < eps2 * abs_dTest_max))) {
			testsVector.push_back(dTest/abs_dTest);
		} else {
			testsVector.push_back(0.);
		}
	}

	*pTestDiff = *max_element(testDiffsVector.begin(), testDiffsVector.end());

	/* returning the maximum error */
	return *max_element(testsVector.begin(), testsVector.end());;
}

void
NonlinearSolverTestSepNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestSepNorm::TestMerge(doublereal& dResCurr,
		const doublereal& dResNew) const
{
	dResCurr += dResNew;
}

doublereal
NonlinearSolverTestSepNorm::TestPost(const doublereal& dRes) const
{
	/* va qui perche' non posso fare sqrt() su !isfinite() */
	if (!std::isfinite(dRes)) {
		throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
	}

	return sqrt(dRes);
}

/* NonlinearSolverTestSepNorm - end */

/* NonlinearSolverTestMinMax */

NonlinearSolverTest::Type NonlinearSolverTestMinMax::GetType() const
{
     return MINMAX;
}

void
NonlinearSolverTestMinMax::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = fabs(Vec(iIndex)) * dCoef;

	if (!(d < dRes)) { // this will work also if d equal nan
		dRes = d;
	}
}

void
NonlinearSolverTestMinMax::TestMerge(doublereal& dResCurr,
		const doublereal& dResNew) const
{
	ASSERT(dResCurr >= 0.);
	ASSERT(dResNew >= 0.);

	if (dResNew > dResCurr) {
		dResCurr = dResNew;
	}
}

/* NonlinearSolverTestMinMax - end */

/* NonlinearSolverTestScale - begin */

NonlinearSolverTestScale::NonlinearSolverTestScale(const VectorHandler* pScl)
: pScale(pScl)
{
	NO_OP;
}

NonlinearSolverTestScale::~NonlinearSolverTestScale(void)
{
	NO_OP;
}

void
NonlinearSolverTestScale::SetScale(const VectorHandler* pScl)
{
	pScale = pScl;
}

const doublereal&
NonlinearSolverTestScale::dScaleCoef(const integer& iIndex) const
{
	return pScale->operator()(iIndex);
}

/* NonlinearSolverTestScale - end */

/* NonlinearSolverTestScaleNorm - begin */

NonlinearSolverTest::Type NonlinearSolverTestScaleNorm::GetType() const
{
     return NORM;
}

void
NonlinearSolverTestScaleNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * (*pScale)(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestScaleNorm::TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const
{
	NonlinearSolverTestNorm::TestMerge(dResCurr, dResNew);
}

const doublereal&
NonlinearSolverTestScaleNorm::dScaleCoef(const integer& iIndex) const
{
	return NonlinearSolverTestScale::dScaleCoef(iIndex);
}

/* NonlinearSolverTestScaleNorm - end */

/* NonlinearSolverTestScaleRelNorm - begin */

void
NonlinearSolverTestScaleRelNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * (*pScale)(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestScaleRelNorm::TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const
{
	NonlinearSolverTestNorm::TestMerge(dResCurr, dResNew);
}

const doublereal&
NonlinearSolverTestScaleRelNorm::dScaleCoef(const integer& iIndex) const
{
	return NonlinearSolverTestScale::dScaleCoef(iIndex);
}

/* NonlinearSolverTestScaleNorm - end */

/* NonlinearSolverTestScaleSepNorm - begin */

void
NonlinearSolverTestScaleSepNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = Vec(iIndex) * (*pScale)(iIndex) * dCoef;

	dRes += d*d;
}

void
NonlinearSolverTestScaleSepNorm::TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const
{
	NonlinearSolverTestSepNorm::TestMerge(dResCurr, dResNew);
}

const doublereal&
NonlinearSolverTestScaleSepNorm::dScaleCoef(const integer& iIndex) const
{
	return NonlinearSolverTestScale::dScaleCoef(iIndex);
}

/* NonlinearSolverTestScaleSepNorm - end */

/* NonlinearSolverTestScaleMinMax - begin */

NonlinearSolverTest::Type NonlinearSolverTestScaleMinMax::GetType() const
{
     return MINMAX;
}

void
NonlinearSolverTestScaleMinMax::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex, doublereal dCoef) const
{
	doublereal d = fabs(Vec(iIndex) * (*pScale)(iIndex)) * dCoef;

	if (d > dRes) {
		dRes = d;
	}
}

void
NonlinearSolverTestScaleMinMax::TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const
{
	NonlinearSolverTestMinMax::TestMerge(dResCurr, dResNew);
}

const doublereal&
NonlinearSolverTestScaleMinMax::dScaleCoef(const integer& iIndex) const
{
	return NonlinearSolverTestScale::dScaleCoef(iIndex);
}

/* NonlinearSolverTestScaleMinMax - end */

bool
NonlinearSolverTestRange::bIsValid(const integer& iIndex) const
{
	ASSERT(m_iFirstIndex > 0);
	ASSERT(m_iLastIndex > m_iFirstIndex);
	return (iIndex >= m_iFirstIndex && iIndex <= m_iLastIndex);
}

NonlinearSolverTestRange::NonlinearSolverTestRange(NonlinearSolverTest *pTest, integer iFirstIndex, integer iLastIndex)
: m_iFirstIndex(iFirstIndex), m_iLastIndex(iLastIndex), m_pTest(pTest)
{
	NO_OP;
}

NonlinearSolverTestRange::~NonlinearSolverTestRange(void)
{
	delete m_pTest;
}

NonlinearSolverTest::Type NonlinearSolverTestRange::GetType() const
{
     return m_pTest->GetType();
}

#if 0
doublereal
NonlinearSolverTestRange::MakeTest(Solver *pS, const integer& Size,
	const VectorHandler& Vec, bool bResidual)
{
	// return m_pTest->MakeTest(pS, Size, Vec, bResidual);
	return NonlinearSolverTest::MakeTest(pS, Size, Vec, bResidual);
}
#endif

void
NonlinearSolverTestRange::TestOne(doublereal& dRes, const VectorHandler& Vec,
	const integer& iIndex, doublereal dCoef) const
{
	if (bIsValid(iIndex)) {
		m_pTest->TestOne(dRes, Vec, iIndex, dCoef);
	}
}

void
NonlinearSolverTestRange::TestMerge(doublereal& dResCurr, const doublereal& dResNew) const
{
	m_pTest->TestMerge(dResCurr, dResNew);
}

doublereal
NonlinearSolverTestRange::TestPost(const doublereal& dRes) const
{
	return m_pTest->TestPost(dRes);
}

const doublereal&
NonlinearSolverTestRange::dScaleCoef(const integer& iIndex) const
{
	if (bIsValid(iIndex)) {
		return m_pTest->dScaleCoef(iIndex);
	}

	return ::Zero1;
}

void
NonlinearSolverTestRange::SetRange(integer iFirstIndex, integer iLastIndex)
{
	if (iFirstIndex <= 0) {
		silent_cerr("NonlinearSolverTestRange: invalid iFirstIndex=" << iFirstIndex << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iLastIndex <= iFirstIndex) {
		silent_cerr("NonlinearSolverTestRange: invalid iLastIndex=" << iLastIndex << " (iFirstIndex=" << iFirstIndex << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_iFirstIndex = iFirstIndex;
	m_iLastIndex = iLastIndex;
}

NonlinearSolverOptions::NonlinearSolverOptions(bool bHonorJacRequest,
											   enum ScaleFlags eScaleFlags,
											   doublereal dScaleAlgebraic)
:bHonorJacRequest(bHonorJacRequest),
 eScaleFlags(eScaleFlags),
 dScaleAlgebraic(dScaleAlgebraic)
{

}

/* NonlinearSolver - begin */

NonlinearSolver::NonlinearSolver(const NonlinearSolverOptions& options)
: NonlinearSolverOptions(options),
Size(0),
TotJac(0),
#ifdef USE_MPI
bParallel(MPI::Is_initialized()),
#endif /* USE_MPI */
pResTest(0),
pSolTest(0),
iNumCond(0),
dMaxCond(0.),
dMinCond(std::numeric_limits<doublereal>::max()),
dSumCond(0.)
#ifdef USE_EXTERNAL
, ExtStepType(External::ERROR)  
#endif /* USE_EXTERNAL */
{
	std::memset(dTimeCPU, 0, sizeof(dTimeCPU));
        std::memset(&oSolverHints, 0, sizeof(oSolverHints));

        SetNonlinearSolverHint(LINESEARCH_LAMBDA_MAX, 1.0);
        SetNonlinearSolverHint(LINESEARCH_LAMBDA_CURR, 1.0);
}

void
NonlinearSolver::SetTest(NonlinearSolverTest *pr, NonlinearSolverTest *ps)
{
	pResTest = pr;
	pSolTest = ps;
}  

NonlinearSolver::~NonlinearSolver(void)
{
	SAFEDELETE(pResTest);
	SAFEDELETE(pSolTest);
}

std::ostream& NonlinearSolver::PrintSolverTime(std::ostream& os) const
{
     using namespace std;
     using namespace std::chrono;
     typedef duration<float, std::ratio<1, 1> > FloatSec;
	     
     auto prec = os.precision();
     auto fmt = os.flags();
     os.setf(ios::fixed);
     os.precision(3);

     os << "nonlinear solver time:\n";
     os << "\toverall CPU time spent in AssRes:\t" << FloatSec(dTimeCPU[CPU_RESIDUAL]).count() << "s\n";
     os << "\toverall CPU time spent in AssJac:\t" << FloatSec(dTimeCPU[CPU_JACOBIAN]).count() << "s\n";
     os << "\toverall CPU time spent in Solve:\t" << FloatSec(dTimeCPU[CPU_LINEAR_SOLVER]).count() << "s\n";

     nanoseconds total(0);

     for (size_t i = 0; i < CPU_LAST_TYPE; ++i) {
	  total += dTimeCPU[i];
     }

     os << "sum of CPU time spent in nonlinear solver:\t" << FloatSec(total).count() << "s\n";
	  
     os.precision(prec);
     os.flags(fmt);
     return os;
}

integer
NonlinearSolver::TotalAssembledJacobian(void)
{
	return TotJac;
}

bool
NonlinearSolver::MakeResTest(Solver *pS,
	const NonlinearProblem *pNLP,
	const VectorHandler& Vec,
	const doublereal& dTol,
	doublereal& dTest,
	doublereal& dTestDiff)
{
	doublereal dScaleAlgEqu;
	const doublereal dTestScale = pNLP->TestScale(pResTest, dScaleAlgEqu);

	if (eScaleFlags == SCALE_ALGEBRAIC_EQUATIONS_NO) {
		// f = c / dCoef
		dScaleAlgEqu = 1.;
	} else {
		// f = c * dScaleAlgebraic
		dScaleAlgEqu *= dScaleAlgebraic;
	}

	dTest = pResTest->MakeTest(pS, Size, Vec, true, dScaleAlgEqu, &dTestDiff) * dTestScale;
	return ((dTest <= dTol) && pS->pGetDataManager()->IsConverged()); // operator <= will work also for NonlinearSolverTestNone
}

bool
NonlinearSolver::MakeSolTest(Solver *pS,
	const VectorHandler& Vec,
	const doublereal& dTol,
	doublereal& dTest)
{
	dTest = pSolTest->MakeTest(pS, Size, Vec);
	return ((dTest <= dTol) && pS->pGetDataManager()->IsConverged()); // operator <= will work also for NonlinearSolverTestNone
}

#ifdef USE_EXTERNAL

void
NonlinearSolver::SetExternal(const External::ExtMessage Ty)
{
	ExtStepType = Ty;
}

void
NonlinearSolver::SendExternal(void)
{
	switch (ExtStepType) {
	case (External::EMPTY):
		External::SendNull();
		break;

	case (External::REGULAR):
		External::SendRegular();
		break;

	case (External::CLOSE):
		External::SendClose();
		break;

	case (External::ERROR):
	default:
		External::SendError();
	}
}

#endif /* USE_EXTERNAL */

/* NonlinearSolver - end */

