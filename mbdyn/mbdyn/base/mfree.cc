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

#include <precond_.h>
#include <mfree.h>
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <umfpackwrap.h>
#include <unistd.h>
#include <output.h>

const doublereal defaultTau = 1.e-7;
const doublereal defaultGamma = 0.9;

MatrixFreeSolver::MatrixFreeSolver(
		const Preconditioner::PrecondType PType, 
		const integer iPStep,
		doublereal ITol,
		integer MaxIt,
		doublereal etaMx) 
: pSM(NULL),
pPM(NULL),
pRes(NULL),
IterTol(ITol),
MaxLinIt(MaxIt),
Tau(defaultTau),
gamma(defaultGamma),
etaMax(etaMx),
PrecondIter(iPStep),
fBuildMat(true),
pPrevNLP(NULL)
{
	
	switch(PType) {
	case Preconditioner::FULLJACOBIAN:
		SAFENEW(pPM, FullJacobianPr);
		break;
	
	default:
		std::cerr << "Unknown Preconditioner type; aborting"
			<< std::endl;
		THROW(ErrGeneric()); 
	}
}

MatrixFreeSolver::~MatrixFreeSolver(void)
{
	NO_OP;
}

doublereal
MatrixFreeSolver::MakeTest(const VectorHandler& Vec)
{
   	DEBUGCOUTFNAME("MatrixFreeSolver::MakeTest");
   
   	doublereal dRes = 0.;
	ASSERT(pSM != NULL);
	
#ifdef USE_MPI
#warning "NonlinSolver MakeTest parallel broken !! "	
#if 0
   	Dof CurrDof;
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
#ifdef __HACK_SCALE_RES__
			ASSERT(pScale != NULL);
			ASSERT(pScale->iGetSize() == Size);
#endif /* __HACK_SCALE_RES__ */
 	  	for (int iCntp1 = 1; iCntp1 <= Size; 
				iCntp1++) {
			doublereal d = Vec.dGetCoef(iCntp1);
			doublereal d2 = d*d;

#ifdef __HACK_SCALE_RES__
			doublereal ds = pScale->dGetCoef(iCntp1);
			doublereal ds2 = ds*ds;
			d2 *= ds2;
#endif /* __HACK_SCALE_RES__ */

			dRes += d2;
		}
		
#ifdef USE_MPI
#if 0 		
	}
#endif
#endif /* USE_MPI */

	/* FIXME: sicuri che va qui? */
	if (!isfinite(dRes)) {      
		THROW(ErrSimulationDiverged());
	}

   	return sqrt(dRes);
}


