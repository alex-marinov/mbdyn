/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
 
#include <nonlin.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <unistd.h>
#include <output.h>

NonlinearSolverTest::~NonlinearSolverTest(void)
{
	NO_OP;
}

const doublereal dOne = 1.;

const doublereal&
NonlinearSolverTest::dScaleCoef(const integer& iIndex) const
{
	return ::dOne;
}

doublereal
NonlinearSolverTestNone::MakeTest(SolutionManager *pSM,
		integer Size, const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNone::MakeTest");

	return 0.;
}

doublereal
NonlinearSolverTestNorm::MakeTest(SolutionManager *pSM,
		integer Size, const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestNorm::MakeTest");

   	doublereal dRes = 0.;
	
#ifdef USE_MPI
#warning "NonlinearSolverTestNorm::MakeTest parallel broken !! "	

#if 0
	ASSERT(pSM != NULL);

	SchurSolutionManager *pSSM = dynamic_cast<SchurSolutionManager *>(pSM);
	if (pSSM != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(Vec);
		
	} else {
#endif

#endif /* USE_MPI */
		ASSERT(Vec.iGetSize() == Size);

 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = Vec.dGetCoef(iCntp1);

			dRes += d*d;
		}
#ifdef USE_MPI
#if 0 
	}
#endif
#endif /* USE_MPI */

	/* va qui perche' non posso fare sqrt() su !isfinite() */
	if (!isfinite(dRes)) {      
		THROW(NonlinearSolver::ErrSimulationDiverged());
	}

   	return sqrt(dRes);
}

doublereal
NonlinearSolverTestMinMax::MakeTest(SolutionManager *pSM,
		integer Size, const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("NonlinearSolverTestMinMax::MakeTest");

   	doublereal dRes = 0.;
	
#ifdef USE_MPI
#warning "NonlinearSolverTestMinMax::MakeTest parallel broken !! "	
#if 0
	ASSERT(pSM != NULL);
	SchurSolutionManager *pSSM;
	if ((pSSM = dynamic_cast<SchurSolutionManager*> (pSM)) != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(dRes, dXPr);
		
	} else {
#endif

#endif /* USE_MPI */
		ASSERT(Vec.iGetSize() == Size);

 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = fabs(Vec.dGetCoef(iCntp1));

			if (d > dRes) {
				dRes = d;
			}
		}
#ifdef USE_MPI
#if 0 
	}
#endif
#endif /* USE_MPI */

	/* va qui perche' non posso fare sqrt() su !isfinite() */
	if (!isfinite(dRes)) {      
		THROW(NonlinearSolver::ErrSimulationDiverged());
	}

   	return dRes;
}

NonlinearSolverTestScale::NonlinearSolverTestScale(const VectorHandler* pScl)
: pScale(pScl)
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

doublereal
NonlinearSolverTestScaleNorm::MakeTest(SolutionManager *pSM,
		integer Size, const VectorHandler& Vec)
{
	DEBUGCOUTFNAME("NonlinearSolverTestScaleNorm::MakeTest");

   	doublereal dRes = 0.;
	
#ifdef USE_MPI
#warning "NonlinearSolverTestScaleNorm::MakeTest parallel broken !! "	
#if 0
	ASSERT(pSM != NULL);
	SchurSolutionManager *pSSM;
	if ((pSSM = dynamic_cast<SchurSolutionManager*> (pSM)) != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(dRes, dXPr);
		
	} else {
#endif

#endif /* USE_MPI */
		ASSERT(Vec.iGetSize() == Size);
		ASSERT(pScale != NULL);
		ASSERT(pScale->iGetSize() == Size);

 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = Vec.dGetCoef(iCntp1);
			doublereal d2 = d*d;

			doublereal ds = pScale->dGetCoef(iCntp1);
			doublereal ds2 = ds*ds;
			d2 *= ds2;

			dRes += d2;
		}
#ifdef USE_MPI
#if 0 
	}
#endif
#endif /* USE_MPI */

	/* FIXME: sicuri che va qui? */
	if (!isfinite(dRes)) {      
		THROW(NonlinearSolver::ErrSimulationDiverged());
	}

   	return sqrt(dRes);
}

doublereal
NonlinearSolverTestScaleMinMax::MakeTest(SolutionManager *pSM,
		integer Size, const VectorHandler& Vec)
{
	DEBUGCOUTFNAME("NonlinearSolverTestScaleMinMax::MakeTest");

   	doublereal dRes = 0.;
	
#ifdef USE_MPI
#warning "NonlinearSolverTestScaleMinMax::MakeTest parallel broken !! "	
#if 0
	ASSERT(pSM != NULL);
	SchurSolutionManager *pSSM;
	if ((pSSM = dynamic_cast<SchurSolutionManager*> (pSM)) != 0) {
		
		/*
		 * Chiama la routine di comunicazione per la trasmissione 
		 * del residuo delle interfacce
		 */
		pSSM->StartExchInt();

		/* calcola il test per i dofs locali */
		int DCount = 0;
		for (int iCnt = 0; iCnt < iNumLocDofs; iCnt++) {
			DCount = pLocDofs[iCnt];
			CurrDof = pDofs[DCount-1];
			doublereal d = Res.dGetCoef(DCount);
			dRes += d*d;
		}

		integer iMI = pSDM->HowManyDofs(SchurDataManager::MYINTERNAL);
		integer *pMI = pSDM->GetDofsList(SchurDataManager::MYINTERNAL);

		/* verifica completamento trasmissioni */
		pSSM->ComplExchInt(dRes, dXPr);
		
	} else {
#endif

#endif /* USE_MPI */
		ASSERT(Vec.iGetSize() == Size);
		ASSERT(pScale != NULL);
		ASSERT(pScale->iGetSize() == Size);

 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = fabs(Vec.dGetCoef(iCntp1))*pScale->dGetCoef(iCntp1);
			if (d > dRes) {
				dRes = d;
			}
		}
#ifdef USE_MPI
#if 0 
	}
#endif
#endif /* USE_MPI */

	/* FIXME: sicuri che va qui? */
	if (!isfinite(dRes)) {      
		THROW(NonlinearSolver::ErrSimulationDiverged());
	}

   	return dRes;
}

NonlinearSolver::NonlinearSolver(void)
: Size(0),
TotJac(0),
pResTest(0)
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
NonlinearSolver::MakeResTest(SolutionManager *pSM, const VectorHandler& Vec)
{
	return pResTest->MakeTest(pSM, Size, Vec);
}

doublereal
NonlinearSolver::MakeSolTest(SolutionManager *pSM, const VectorHandler& Vec)
{
	return pSolTest->MakeTest(pSM, Size, Vec);
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

