/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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
  * Copyright (C) 2003-2005
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <solver.h>
#include <nonlin.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <unistd.h>
#include <output.h>

/* NonlinearSolverTest - begin */

NonlinearSolverTest::~NonlinearSolverTest(void)
{
	NO_OP;
}

doublereal
NonlinearSolverTest::MakeTest(Solver *pS, const integer& Size,
		const VectorHandler& Vec, bool bResidual)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNorm::MakeTest");

   	doublereal dTest = 0.;

#ifdef USE_MPI
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
#endif /* USE_MPI */
	{
		ASSERT(Vec.iGetSize() == Size);

 	  	for (int iCntp1 = 1; iCntp1 <= Size; iCntp1++) {
			TestOne(dTest, Vec, iCntp1);
		}
	}

	return TestPost(dTest);
}

doublereal
NonlinearSolverTest::TestPost(const doublereal& dRes) const
{
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

doublereal
NonlinearSolverTestNone::MakeTest(Solver *pS, integer Size, 
		const VectorHandler& Vec, bool bResidual)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNone::MakeTest");

	return 0.;
}

void
NonlinearSolverTestNone::TestOne(doublereal& dRes, 
		const VectorHandler& Vec, const integer& iIndex) const
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

void
NonlinearSolverTestNorm::TestOne(doublereal& dRes, 
		const VectorHandler& Vec, const integer& iIndex) const
{
	doublereal d = Vec(iIndex);

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
	if (!isfinite(dRes)) {      
		throw NonlinearSolver::ErrSimulationDiverged();
	}

   	return sqrt(dRes);
}

/* NonlinearSolverTestNorm - end */

/* NonlinearSolverTestMinMax */

void
NonlinearSolverTestMinMax::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex) const
{
	doublereal d = fabs(Vec(iIndex));

	if (d > dRes) {
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
	return pScale->dGetCoef(iIndex);
}

/* NonlinearSolverTestScale - end */

/* NonlinearSolverTestScaleNorm - begin */

void
NonlinearSolverTestScaleNorm::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex) const
{
	doublereal d = Vec(iIndex) * (*pScale)(iIndex);

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

/* NonlinearSolverTestScaleMinMax - begin */

void
NonlinearSolverTestScaleMinMax::TestOne(doublereal& dRes,
		const VectorHandler& Vec, const integer& iIndex) const
{
	doublereal d = fabs(Vec(iIndex) * (*pScale)(iIndex));

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

/* NonlinearSolver - begin */

NonlinearSolver::NonlinearSolver(bool JacReq)
: Size(0),
TotJac(0),
bHonorJacRequest(JacReq),
#ifdef USE_MPI
bParallel(MPI::Is_initialized()),
#endif /* USE_MPI */
pResTest(0),
pSolTest(0)
#ifdef USE_EXTERNAL
, ExtStepType(External::ERROR)  
#endif /* USE_EXTERNAL */
{
	NO_OP;
}

void
NonlinearSolver::SetTest(NonlinearSolverTest *pr, NonlinearSolverTest *ps)
{
	pResTest = pr;
	pSolTest = ps;
}  

NonlinearSolver::~NonlinearSolver(void)
{
	NO_OP;
}

integer
NonlinearSolver::TotalAssembledJacobian(void)
{
	return TotJac;
}

doublereal
NonlinearSolver::MakeResTest(Solver *pS, const VectorHandler& Vec)
{
	return pResTest->MakeTest(pS, Size, Vec, true);
}

doublereal
NonlinearSolver::MakeSolTest(Solver *pS, const VectorHandler& Vec)
{
	return pSolTest->MakeTest(pS, Size, Vec);
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

