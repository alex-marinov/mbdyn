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
 * Copyright 1999-2003 Giuseppe Quaranta <giuquaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * This copyright statement applies to the MPI related code, which was
 * merged from files schur.h/schur.cc
 */

/* 
 *
 * Copyright (C) 2003
 * Giuseppe Quaranta	<quaranta@aero.polimi.it>
 *
 */

/* metodo per la soluzione del modello */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include <ac/float.h>
#include <ac/math.h>

#include <solver.h>
#include <nr.h>
#include <bicg.h>
#include <gmres.h>
#include <solman.h>
#include <vector>

#if defined(HAVE_SIGNAL) && defined(HAVE_SIGNAL_H)
#include <signal.h>
#endif /* HAVE_SIGNAL && HAVE_SIGNAL_H */
  
#ifdef USE_MPI
#include <mbcomm.h>
#ifdef USE_EXTERNAL
#include <external.h>
#endif /* USE_EXTERNAL */
#endif /* USE_MPI */


#if defined(USE_RTAI) 
#include <mbrtai_utils.h>
#if defined(HAVE_SYS_MMAN_H)
#include <sys/mman.h>
#endif /* HAVE_SYS_MMAN_H */
#endif /* USE_RTAI */

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

#ifdef HAVE_SIGNAL
static volatile sig_atomic_t mbdyn_keep_going = 1;
static __sighandler_t mbdyn_sh_term = SIG_DFL;
static __sighandler_t mbdyn_sh_int = SIG_DFL;
static __sighandler_t mbdyn_sh_hup = SIG_DFL;

static void
modify_final_time_handler(int signum)
{
   	::mbdyn_keep_going = 0;
   	switch (signum) {
    	case SIGTERM:
      		signal(signum, ::mbdyn_sh_term);
      		break;
	
    	case SIGINT:
      		signal(signum, ::mbdyn_sh_int);
      		break;
	
    	case SIGHUP:
      		signal(signum, ::mbdyn_sh_hup);
      		break;
   	}
}
#endif /* HAVE_SIGNAL */

#ifdef USE_RTAI
static int
reserve_stack(unsigned long size)
{
	int buf[size];

#ifdef HAVE_MEMSET
	memset(buf, '\0', size*sizeof(int));
#else /* !HAVE_MEMSET */
	for (unsigned long i = 0; i < size; i++) {
		buf[i] = 0;
	}
#endif /* !HAVE_MEMSET */

#ifdef HAVE_MLOCKALL
	return mlockall(MCL_CURRENT | MCL_FUTURE);
#else /* !HAVE_MLOCKALL */
	return 0;
#endif /* !HAVE_MLOCKALL */
}
#endif /* USE_RTAI */


/* Parametri locali */
const integer iDefaultFictitiousStepsNumber = 0;
const doublereal dDefaultFictitiousStepsRatio = 1.e-3;
const integer iDefaultIterationsBeforeAssembly = 2;
const integer iDefaultIterativeSolversMaxSteps = 100;
const integer iDefaultPreconditionerSteps = 20;
const doublereal dDefaultTol = 1.e-6;
const doublereal defaultIterativeEtaMax = 0.9;


/* Costruttore: esegue la simulazione */
Solver::Solver(MBDynParser& HPar,
		const char* sInFName,
		const char* sOutFName,
		flag fPar)
:
CurrStrategy(NOCHANGE),
sInputFileName(NULL),
sOutputFileName(NULL),
HP(HPar),
pStrategyChangeDrive(NULL),
#ifdef __HACK_EIG__
fEigenAnalysis(0),
dEigParam(1.),
fOutputModes(0),
dUpperFreq(FLT_MAX),
dLowerFreq(0.),
#endif /* __HACK_EIG__ */
#ifdef USE_RTAI
bRT(false),
bRTAllowNonRoot(false),
RTMode(MBRTAI_UNKNOWN),
bRTHard(false),
lRTPeriod(-1),
RTSemPtr(NULL),
RTStackSize(1024),
#endif /* USE_RTAI */
#ifdef __HACK_POD__
fPOD(0),
iPODStep(0),
iPODFrames(0),
#endif /*__HACK_POD__*/
iNumPreviousVectors(3),
pdWorkSpace(NULL),
qX(),
qXPrime(),
pDM(NULL),
iNumDofs(0),
dTime(0.),
dInitialTime(0.), 
dFinalTime(0.),
dRefTimeStep(0.),
dInitialTimeStep(1.), 
dMinimumTimeStep(1.),
dMaxTimeStep(1.),
iFictitiousStepsNumber(iDefaultFictitiousStepsNumber),
dFictitiousStepsRatio(dDefaultFictitiousStepsRatio),
fAbortAfterInput(0),
fAbortAfterAssembly(0),
fAbortAfterDerivatives(0),
fAbortAfterFictitiousSteps(0),
iStepsAfterReduction(0),
iStepsAfterRaise(0),
iWeightedPerformedIters(0),
fLastChance(0),
pDerivativeSteps(NULL),
pFirstRegularStep(NULL),
pRegularSteps(NULL),
pFictitiousSteps(NULL),
CurrSolver(LinSol::defaultSolver),
iWorkSpaceSize(0),
dPivotFactor(-1.),
fTrueNewtonRaphson(1),
iIterationsBeforeAssembly(0),
fMatrixFree(0),
MFSolverType(MatrixFreeSolver::BICGSTAB),
dIterTol(dDefaultTol),
PcType(Preconditioner::FULLJACOBIAN),
iPrecondSteps(iDefaultPreconditionerSteps),
iIterativeMaxSteps(iDefaultPreconditionerSteps),
dIterertiveEtaMax(defaultIterativeEtaMax),
#ifdef USE_MPI
fParallel(fPar),
pSDM(NULL),
CurrIntSolver(LinSol::defaultSolver),
iNumLocDofs(0),
iNumIntDofs(0),
pLocDofs(NULL),
pIntDofs(NULL),
pDofs(NULL),
iIWorkSpaceSize(0),
dIPivotFactor(-1.),
pLocalSM(NULL),
pSSM(NULL),
#endif /* USE_MPI */
pSM(NULL),
pCurrSM(NULL),
pNLS(NULL)
{
	DEBUGCOUTFNAME("Solver::Solver");

	ASSERT(sInFName != NULL);
	
	SAFESTRDUP(sInputFileName, sInFName);

	if (sOutFName != NULL) {
		SAFESTRDUP(sOutputFileName, sOutFName);
	}
	
   	/* Legge i dati relativi al metodo di integrazione */
   	ReadData(HP);

#if USE_RTAI
	if (bRT) {
		/* FIXME: if using RTAI, clear out output */
		SetOutputFlags(OUTPUT_NONE);
	}
#endif /* USE_RTAI */
}


void Solver::Run(void)
{
   	DEBUGCOUTFNAME("Solver::Run");

#ifdef USE_MPI
	int MyRank = 0;
	if (fParallel) {
		
		/*
		 * E' necessario poter determinare in questa routine
		 * quale e' il master in modo da far calcolare la soluzione
		 * solo su di esso
		 */
		MyRank = MBDynComm.Get_rank();
		/* chiama il gestore dei dati generali della simulazione */

		/* 
		 * I file di output vengono stampati localmente 
		 * da ogni processo aggiungendo al termine 
		 * dell'OutputFileName il rank del processo 
		 */
		int iRankLength = 3;	/* should be configurable? */
		char* sNewOutName = NULL;
		const char* sOutName = NULL;

		if (sOutputFileName == NULL) {
			int iOutLen = strlen(sInputFileName);
			SAFENEWARR(sNewOutName, char, iOutLen+1+iRankLength+1);
			sOutName = sInputFileName;
		} else {
			int iOutLen = strlen(sOutputFileName);
			SAFENEWARR(sNewOutName, char, iOutLen+1+iRankLength+1);
			sOutName = sOutputFileName;
		}
		sprintf(sNewOutName,"%s.%.*d", sOutName, iRankLength, MyRank);

		DEBUGLCOUT(MYDEBUG_MEM, "creating parallel DataManager" 
				<< std::endl);
		
		SAFENEWWITHCONSTRUCTOR(pSDM,
			SchurDataManager,
			SchurDataManager(HP,
				dInitialTime, 
				sInputFileName,
				sNewOutName,
				fAbortAfterInput));
		pDM = pSDM;

	} else {
#endif /* USE_MPI */

		/* chiama il gestore dei dati generali della simulazione */
		DEBUGLCOUT(MYDEBUG_MEM, "creating DataManager" << std::endl);
		SAFENEWWITHCONSTRUCTOR(pDM,
 				DataManager,
 				DataManager(HP, 
					dInitialTime, 
					sInputFileName,
					sOutputFileName,
					fAbortAfterInput));
#ifdef USE_MPI
	}
#endif /* USE_MPI */
   
   	/* Si fa dare l'std::ostream al file di output per il log */
   	std::ostream& Out = pDM->GetOutFile();

   	if (fAbortAfterInput) {
      		/* Esce */
		pDM->Output(true);
      		Out << "End of Input; no simulation or assembly is required."
			<< std::endl;
      		return;

   	} else if (fAbortAfterAssembly) {
      		/* Fa l'output dell'assemblaggio iniziale e poi esce */
      		pDM->Output(true);
      		Out << "End of Initial Assembly; no simulation is required."
			<< std::endl;
      		return;
   	}

	/* Qui crea le partizioni: principale fra i processi, se parallelo  */
#ifdef USE_MPI
	if (fParallel) {
		pSDM->CreatePartition();
	}
#endif /* USE_MPI */

   	/* Si fa dare il DriveHandler e linka i drivers di rho ecc. */
   	const DriveHandler* pDH = pDM->pGetDrvHdl();
   	pRegularSteps->SetDriveHandler(pDH);
   	pFictitiousSteps->SetDriveHandler(pDH);
   
   	/* Costruisce i vettori della soluzione ai vari passi */
   	DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << std::endl);

#ifdef USE_MPI
	if (fParallel) {
		iNumDofs = pSDM->HowManyDofs(SchurDataManager::TOTAL);
		pDofs = pSDM->pGetDofsList();
		
		iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
		pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);
		iNumIntDofs = pSDM->HowManyDofs(SchurDataManager::INTERNAL);
		pIntDofs = pSDM->GetDofsList(SchurDataManager::INTERNAL);

	} else {
#endif /* USE_MPI */
   		iNumDofs = pDM->iGetNumDofs();
#ifdef USE_MPI
	}
#endif /* USE_MPI */
	
   	ASSERT(iNumDofs > 0);        
	
	integer iNSteps = pRegularSteps->GetIntegratorStepSize();
	integer iFSteps = pFictitiousSteps->GetIntegratorStepSize();
	iNumPreviousVectors = (iNSteps < iFSteps) ? iFSteps : iNSteps;
	
	SAFENEWARR(pdWorkSpace,doublereal, 2*iNumPreviousVectors*iNumDofs);
	
	MyVectorHandler* pX = NULL;
	MyVectorHandler* pXPrime = NULL;
	for (int ivec = 0; ivec < iNumPreviousVectors; ivec++) {  
   		SAFENEWWITHCONSTRUCTOR(pX,
			       	MyVectorHandler,
			       	MyVectorHandler(iNumDofs, pdWorkSpace+ivec*iNumDofs));
		qX.push_back(pX);
   		SAFENEWWITHCONSTRUCTOR(pXPrime,
			       	MyVectorHandler,
			       	MyVectorHandler(iNumDofs, pdWorkSpace+(iNumPreviousVectors+ivec)*iNumDofs));
		qXPrime.push_back(pXPrime);
		pX = NULL;
		pXPrime = NULL;
	}

	/* Resetta i vettori */
   	for (int ivec = 0; ivec < iNumPreviousVectors; ivec++) {
		qX[ivec]->Reset(0.);
		qXPrime[ivec]->Reset(0.);
	}

#ifdef __HACK_POD__
	std::ofstream PodOut;
	if (fPOD) {
		char *PODFileName = NULL;
		
		if (sOutputFileName == NULL) {
			SAFESTRDUP(PODFileName, "MBDyn.POD");
		} else {
			size_t l = strlen(sOutputFileName);
			SAFENEWARR(PODFileName, char, l+sizeof(".POD"));

			memcpy(PODFileName, sOutputFileName, l);
			memcpy(PODFileName+l, ".POD", sizeof(".POD"));
		}
	
		PodOut.open(PODFileName);
		if (!PodOut) {
			std::cerr << "unable to open file \"" << PODFileName
				<< "\"" << std::endl;
			THROW(ErrGeneric());
		}
		SAFEDELETEARR(PODFileName);

#ifdef __HACK_POD_BINARY__
		/* matrix size is coded at the beginning */
		PodOut.write((char *)&(pod.iFrames), sizeof(unsigned long));
		PodOut.write((char *)&iNumDofs, sizeof(unsigned long));
#endif /* __HACK_POD_BINARY__ */
	}
#endif /* __HACK_POD__ */


   	/* Subito collega il DataManager alla soluzione corrente */
   	pDM->LinkToSolution(*(qX[0]), *(qXPrime[0]));         


   	/* costruisce il SolutionManager */
   	DEBUGLCOUT(MYDEBUG_MEM, "creating SolutionManager, size = "
		   << iNumDofs << std::endl);

	
	integer iLWS = iWorkSpaceSize;
	integer iNLD = iNumDofs;
#ifdef USE_MPI
	if (fParallel) {
		iLWS = iWorkSpaceSize*iNumLocDofs/(iNumDofs*iNumDofs);
		iNLD = iNumLocDofs;
	}
#endif /* USE_MPI */

   	switch (CurrSolver) {
     	case LinSol::Y12_SOLVER: 
#ifdef USE_Y12
      		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			Y12SparseLUSolutionManager,
			Y12SparseLUSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			MeschachSparseLUSolutionManager,
			MeschachSparseLUSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
		break;
#else /* !USE_MESCHACH */
		std::cerr << "Configure with --with-meschach "
			"to enable Meschach solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

 	case LinSol::HARWELL_SOLVER:
#ifdef USE_HARWELL
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			HarwellSparseLUSolutionManager,
			HarwellSparseLUSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_HARWELL */
      		std::cerr << "Configure with --with-harwell "
			"to enable Harwell solver" << std::endl;
		THROW(ErrGeneric());
#endif /* !USE_HARWELL */

	case LinSol::UMFPACK_SOLVER:
#ifdef USE_UMFPACK
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			UmfpackSparseLUSolutionManager,
			UmfpackSparseLUSolutionManager(iNLD, 
				0, dPivotFactor));
      		break;
#else /* !USE_UMFPACK */
      		std::cerr << "Configure with --with-umfpack "
			"to enable Umfpack solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_UMFPACK */

	case LinSol::EMPTY_SOLVER:
		break;
				
   	default:
		ASSERT(0);
		THROW(ErrGeneric());

	}

	/*
	 * This is the LOCAL solver if instantiating a parallel
	 * integrator; otherwise it is the MAIN solver
	 */
#ifdef USE_MPI
	if (fParallel) {
		pLocalSM = pCurrSM;

		/* Crea il solutore di Schur globale */
		switch (CurrIntSolver) {
		case LinSol::Y12_SOLVER:
#ifdef USE_Y12
			SAFENEWWITHCONSTRUCTOR(pSSM,
				SchurSolutionManager,
				SchurSolutionManager(iNumDofs, pLocDofs,
					iNumLocDofs,
					pIntDofs, iNumIntDofs,
					pLocalSM,
					(Y12SparseLUSolutionManager*)0,
					iIWorkSpaceSize,
					dIPivotFactor == -1. ? 1. : dIPivotFactor));
			break;
#else /* !USE_Y12 */
			std::cerr << "Configure with --with-y12 "
				"to enable Y12 solver" << std::endl;
			THROW(ErrGeneric());
#endif /* !USE_Y12 */

		case LinSol::HARWELL_SOLVER:
			std::cerr << "Harwell solver cannot be used "
				"as interface solver"
#ifdef USE_MESCHACH
				"; switching to Meschach" 
#endif /* USE_MESCHACH */
				<< std::endl;
#ifndef USE_MESCHACH
			THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

		case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
			SAFENEWWITHCONSTRUCTOR(pSSM,
				SchurSolutionManager,
				SchurSolutionManager(iNumDofs, pLocDofs,
					iNumLocDofs,
					pIntDofs, iNumIntDofs,
					pLocalSM,
					(MeschachSparseLUSolutionManager*)0,
					iIWorkSpaceSize,
					dIPivotFactor == -1. ? 1. : dIPivotFactor));
			break;
#else /* !USE_MESCHACH */
			std::cerr << "Configure with --with-meschach "
				"to enable Meschach solver" << std::endl;
			THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

		case LinSol::UMFPACK_SOLVER:
#ifdef USE_UMFPACK
			SAFENEWWITHCONSTRUCTOR(pSSM,
				SchurSolutionManager,
				SchurSolutionManager(iNumDofs, pLocDofs,
					iNumLocDofs,
					pIntDofs, iNumIntDofs,
					pLocalSM,
					(UmfpackSparseLUSolutionManager*)0,
					0, 
					dIPivotFactor));
			break;
#else /* !USE_UMFPACK */
			std::cerr << "Configure with --with-umfpack "
				"to enable Umfpack solver" << std::endl;
			THROW(ErrGeneric());
#endif /* !USE_UMFPACK */
		case LinSol::EMPTY_SOLVER:
			break;
		}

		pSM = pSSM;
		
	} else {
#endif /* USE_MPI */
		pSM = pCurrSM;
#ifdef USE_MPI
	}
#endif /* USE_MPI */
	
	/*
	 * FIXME: at present there MUST be a pSM
	 * (even for matrix-free nonlinear solvers)
	 */
	if (pSM == NULL) {
		std::cerr << "No linear solver defined" << std::endl;
		THROW(ErrGeneric());
	}
	
	/* a questo punto si costruisce il nonlinear solver passandogli 
	   il solution manager */
	if (fMatrixFree) {
		if (MFSolverType == MatrixFreeSolver::BICGSTAB) {
			SAFENEWWITHCONSTRUCTOR(pNLS,
					BiCGStab,
					BiCGStab(PcType, 
						iPrecondSteps,
						dIterTol, 
						iIterativeMaxSteps,
						dIterertiveEtaMax));
		} else {
			SAFENEWWITHCONSTRUCTOR(pNLS,
					Gmres,
					Gmres(PcType, 
						iPrecondSteps,
						dIterTol, 
						iIterativeMaxSteps,
						dIterertiveEtaMax));
		}			  
	} else {
		SAFENEWWITHCONSTRUCTOR(pNLS,
				NewtonRaphsonSolver,
				NewtonRaphsonSolver(fTrueNewtonRaphson,
					iIterationsBeforeAssembly));  
	}

#ifdef __HACK_SCALE_RES__
	MyVectorHandler Scale(iNumDofs);
	VectorHandler *pScale = &Scale;
	pDM->SetScale(*pScale);
	pNLS->SetScale(pScale);
#endif /* __HACK_SCALE_RES__ */

   	/*
	 * Dell'assemblaggio iniziale dei vincoli se ne occupa il DataManager 
	 * in quanto e' lui il responsabile dei dati della simulazione,
	 * e quindi anche della loro coerenza. Inoltre e' lui a sapere
	 * quali equazioni sono di vincolo o meno.
	 */
   
   	/*
	 * Dialoga con il DataManager per dargli il tempo iniziale 
	 * e per farsi inizializzare i vettori di soluzione e derivata */
	dTime = dInitialTime;
	pDM->SetTime(dTime);
	pDM->SetValue(*(qX[0]), *(qXPrime[0]));
	
#ifdef __HACK_EIG__  
   	if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
	 	Eig();
	 	OneEig.fDone = flag(1);
   	}
#endif /* __HACK_EIG__ */
   
   	integer iTotIter = 0;
	integer iStIter = 0;
   	doublereal dTotErr = 0.;
	doublereal dTest = DBL_MAX;
#ifdef MBDYN_X_CONVSOL
   	doublereal dSolTest = DBL_MAX;
	bool bSolConv = false;
#endif /* MBDYN_X_CONVSOL */	/* calcolo delle derivate */
   	DEBUGLCOUT(MYDEBUG_DERIVATIVES, "derivatives solution step"
			<< std::endl);

#ifdef HAVE_SIGNAL
	/*
	 * FIXME: don't do this if compiling with USE_RTAI
	 * Re FIXME: use sigaction() ...
	 */
   	::mbdyn_sh_term = signal(SIGTERM, modify_final_time_handler);
   	::mbdyn_sh_int = signal(SIGINT, modify_final_time_handler);
   	::mbdyn_sh_hup = signal(SIGHUP, modify_final_time_handler);
   
   	if (::mbdyn_sh_term == SIG_IGN) {
      		signal (SIGTERM, SIG_IGN);
   	}
   	if (::mbdyn_sh_int == SIG_IGN) {
      		signal (SIGINT, SIG_IGN);
   	}
   	if (::mbdyn_sh_hup == SIG_IGN) {
      		signal (SIGHUP, SIG_IGN);
   	}
#endif /* HAVE_SIGNAL */

	/* settaggio degli output Types */
	unsigned OF = OutputFlags;
	if ( DEBUG_LEVEL_MATCH(MYDEBUG_RESIDUAL) ) {
		OF |= OUTPUT_RES;
	}
	if ( DEBUG_LEVEL_MATCH(MYDEBUG_JAC) ) {
		OF |= OUTPUT_JAC;
	}
	if ( DEBUG_LEVEL_MATCH(MYDEBUG_SOL) ) {
		OF |= OUTPUT_SOL;
	}
	pNLS->SetOutputFlags(OF);
	
	pDerivativeSteps->SetDataManager(pDM);
	pFirstFictitiousStep->SetDataManager(pDM);
	pFictitiousSteps->SetDataManager(pDM);
	pFirstRegularStep->SetDataManager(pDM);
	pRegularSteps->SetDataManager(pDM);
	
	pDerivativeSteps->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));
	
	pFirstFictitiousStep->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));
	
	pFictitiousSteps->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));
	
	pFirstRegularStep->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));

	pRegularSteps->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));

#ifdef USE_EXTERNAL
	pNLS->SetExternal(External::EMPTY);
#endif /* USE_EXTERNAL */

#ifdef USE_RTAI
	if (bRT) {
		/* Init RTAI; if init'ed, it will be shut down at exit */
		if (mbdyn_rt_task_init("MBDTSK", 1, 0, 0, &mbdyn_rtai_task)) {
			std::cerr << "unable to init RTAI task" << std::endl;
			THROW(ErrGeneric());
		}
	}
#endif /* USE_RTAI */
	
	try {
		
		dTest = pDerivativeSteps->Advance(0., 1.,
				StepIntegrator::NEWSTEP,
			 	pSM, pNLS, qX, qXPrime, iStIter
#ifdef MBDYN_X_CONVSOL
				, dSolTest
#endif /* MBDYN_X_CONVSOL */
				);
	}
	catch (NonlinearSolver::NoConvergence) {
		std::cerr << std::endl
			<< "Initial derivatives calculation " << iStIter 
			<< " does not converge;" << std::endl
			<< "aborting ..." << std::endl;	 
	 	pDM->Output(true);
	 	THROW(ErrMaxIterations());

	}
	catch (NonlinearSolver::ErrSimulationDiverged) {
		/*
		 * Mettere qui eventuali azioni speciali 
		 * da intraprendere in caso di errore ...
		 */
		THROW(SimulationDiverged());
	}
	catch (NonlinearSolver::ConvergenceOnSolution) {
#ifdef MBDYN_X_CONVSOL
		bSolConv = true;
#endif /* MBDYN_X_CONVSOL */
	}
	catch (...) {
		THROW(ErrGeneric());
	}
	dTotErr  += dTest;
	iTotIter += iStIter;
	  	
	if (outputMsg()) {	
   		Out << "Derivatives solution step at time " << dInitialTime
     			<< " performed in " << iStIter
     			<< " iterations with " << dTest
     			<< " error" << std::endl;
	}
      
   	DEBUGCOUT("Derivatives solution step has been performed successfully"
		  " in " << iStIter << " iterations" << std::endl);
   
   	if (fAbortAfterDerivatives) {
      		/*
		 * Fa l'output della soluzione delle derivate iniziali ed esce
		 */
      		pDM->Output(true);
      		Out << "End of derivatives; no simulation is required."
			<< std::endl;
      		return;
#ifdef HAVE_SIGNAL
   	} else if (!::mbdyn_keep_going) {
      		/*
		 * Fa l'output della soluzione delle derivate iniziali ed esce
		 */
      		pDM->Output(true);
      		Out << "Interrupted during derivatives computation." << std::endl;
      		return;
#endif /* HAVE_SIGNAL */
   	}

	/* Dati comuni a passi fittizi e normali */
   	integer iStep = 1;
   	doublereal dCurrTimeStep = 0.;

   	if (iFictitiousStepsNumber > 0) {

       		/* passi fittizi */
      
      		/*
		 * inizio integrazione: primo passo a predizione lineare
		 * con sottopassi di correzione delle accelerazioni
		 * e delle reazioni vincolari
		 */
      		pDM->BeforePredict(*(qX[0]), *(qXPrime[0]),
				   *(qX[1]), *(qXPrime[1]));
      		Flip();

      		dRefTimeStep = dInitialTimeStep*dFictitiousStepsRatio;
      		dCurrTimeStep = dRefTimeStep;
      		pDM->SetTime(dTime+dCurrTimeStep);    
      	
      		DEBUGLCOUT(MYDEBUG_FSTEPS, "Current time step: "
			   << dCurrTimeStep << std::endl);
		
	 	ASSERT(pFirstFictitiousStep != NULL);
		
		try {
	 		dTest = pFirstFictitiousStep->Advance(dRefTimeStep, 1.,
				StepIntegrator::NEWSTEP,
				pSM, pNLS, qX, qXPrime, iStIter
#ifdef MBDYN_X_CONVSOL
				, dSolTest
#endif /* MBDYN_X_CONVSOL */
				);					
		}
		catch (NonlinearSolver::NoConvergence) {
			std::cerr << std::endl
				<< " first dummy step;" << std::endl
				<< " does not converge;" << std::endl
				<< "time step dt = " << dCurrTimeStep 
				<< " cannot be reduced further;"
				<< std::endl
				<< "aborting ..." << std::endl;
	 		pDM->Output(true);
	 		THROW(ErrMaxIterations());
		}
		catch (NonlinearSolver::ErrSimulationDiverged) {
			/*
		 	 * Mettere qui eventuali azioni speciali 
		 	 * da intraprendere in caso di errore ...
		 	 */
			THROW(SimulationDiverged());
		}
		catch (NonlinearSolver::ConvergenceOnSolution) {
#ifdef MBDYN_X_CONVSOL
			bSolConv = true;
#endif /* MBDYN_X_CONVSOL */
		}
		catch (...) {
			THROW(ErrGeneric());
		}
		pDM->AfterConvergence();
      
      		dRefTimeStep = dCurrTimeStep;
      		dTime += dRefTimeStep;
      
      		dTotErr += dTest;
      		iTotIter += iStIter;

#ifdef HAVE_SIGNAL
      		if (!::mbdyn_keep_going) {
	 		/*
			 * Fa l'output della soluzione delle derivate iniziali
			 * ed esce
			 */
#ifdef DEBUG_FICTITIOUS
	   		pDM->Output(true);
#endif /* DEBUG_FICTITIOUS */
	 		Out << "Interrupted during first dummy step." << std::endl;
	 		return;
      		}
#endif /* HAVE_SIGNAL */
      
#ifdef DEBUG_FICTITIOUS
      		pDM->Output(true);
#endif /* DEBUG_FICTITIOUS */
            
       		/* Passi fittizi successivi */
	
      		for (int iSubStep = 2;
		     iSubStep <= iFictitiousStepsNumber;
		     iSubStep++) {
      			pDM->BeforePredict(*(qX[0]), *(qXPrime[0]),
				   	*(qX[1]), *(qXPrime[1]));
	 		Flip();
	 
	 		DEBUGLCOUT(MYDEBUG_FSTEPS, "Fictitious step "
				   << iSubStep 
				   << "; current time step: " << dCurrTimeStep
				   << std::endl);
	 
	 		ASSERT(pFictitiousSteps!= NULL);
			try {
	 			pDM->SetTime(dTime+dCurrTimeStep);
	 			dTest = pFictitiousSteps->Advance(dRefTimeStep,
						dCurrTimeStep/dRefTimeStep,
					 	StepIntegrator::NEWSTEP, 
				 		pSM, pNLS, qX, qXPrime, iStIter
#ifdef MBDYN_X_CONVSOL
						, dSolTest
#endif /* MBDYN_X_CONVSOL */
						);						
			}
			catch (NonlinearSolver::NoConvergence) {
				std::cerr << std::endl
					<< "Dummy step: " << iSubStep << std::endl
					<< " does not converge;" << std::endl
					<< "time step dt = " << dCurrTimeStep 
					<< " cannot be reduced further;"
					<< std::endl
					<< "aborting ..." << std::endl;
	 			pDM->Output(true);
	 			THROW(ErrMaxIterations());
			}

			catch (NonlinearSolver::ErrSimulationDiverged) {
				/*
		 		 * Mettere qui eventuali azioni speciali 
		 		 * da intraprendere in caso di errore ...
		 		 */
				THROW(SimulationDiverged());
			}			
			catch (NonlinearSolver::ConvergenceOnSolution) {
#ifdef MBDYN_X_CONVSOL
				bSolConv = true;
#endif /* MBDYN_X_CONVSOL */
			}
			catch (...) {
				THROW(ErrGeneric());
			}
			pDM->AfterConvergence();
      
      			dTotErr += dTest;
      			iTotIter += iStIter;
	 
#ifdef DEBUG
	 		if (DEBUG_LEVEL(MYDEBUG_FSTEPS)) {
	    			Out << "Step " << iStep 
					<< " time " << dTime+dCurrTimeStep
					<< " step " << dCurrTimeStep
					<< " iterations " << iStIter
					<< " error " << dTest << std::endl;
	 		}
#endif /* DEBUG */

	 		DEBUGLCOUT(MYDEBUG_FSTEPS, "Substep " << iSubStep 
				   << " of step " << iStep 
				   << " has been completed successfully in "
				   << iStIter << " iterations" << std::endl);
				   
#ifdef HAVE_SIGNAL
	 		if (!::mbdyn_keep_going) {
				/* */
#ifdef DEBUG_FICTITIOUS
	    			pDM->Output();
#endif /* DEBUG_FICTITIOUS */
	    			Out << "Interrupted during dummy steps."
					<< std::endl;
				return;
			}
#endif /* HAVE_SIGNAL */

			dTime += dRefTimeStep;	  
      		}
		if (outputMsg()) {	
      			Out << "Initial solution after dummy steps "
				"at time " << dTime
				<< " performed in " << iStIter
				<< " iterations with " << dTest 
				<< " error" << std::endl;
		}
			
      		DEBUGLCOUT(MYDEBUG_FSTEPS, 
			   "Fictitious steps have been completed successfully"
			   " in " << iStIter << " iterations" << std::endl);
   	} /* Fine dei passi fittizi */

   	/* Output delle "condizioni iniziali" */
   	pDM->Output();

	   
        if (outputMsg()) {	
  	 	Out << "Step " << 0
     			<< " time " << dTime+dCurrTimeStep
     			<< " step " << dCurrTimeStep
     			<< " iterations " << iStIter
     			<< " error " << dTest << std::endl;
	}
   
   	if (fAbortAfterFictitiousSteps) {
      		Out << "End of dummy steps; no simulation is required."
			<< std::endl;
		return;
#ifdef HAVE_SIGNAL
   	} else if (!::mbdyn_keep_going) {
      		/* Fa l'output della soluzione ed esce */
      		Out << "Interrupted during dummy steps." << std::endl;
      		return;
#endif /* HAVE_SIGNAL */
   	}

	/* primo passo regolare */
#ifdef USE_EXTERNAL
	pNLS->SetExternal(External::REGULAR);
#endif /* USE_EXTERNAL */
	
   	iStep = 1; /* Resetto di nuovo iStep */
      
   	DEBUGCOUT("Step " << iStep << " has been completed successfully in "
		  << iStIter << " iterations" << std::endl);

   			
   	DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   
      	pDM->BeforePredict(*(qX[0]), *(qXPrime[0]),
				*(qX[1]), *(qXPrime[1]));
	
	Flip();
	dRefTimeStep = dInitialTimeStep;   
   	dCurrTimeStep = dRefTimeStep;
	 	
	ASSERT(pFirstRegularStep!= NULL);
	StepIntegrator::StepChange CurrStep 
			= StepIntegrator::NEWSTEP;
   
IfFirstStepIsToBeRepeated:
	try {   	
		pDM->SetTime(dTime+dCurrTimeStep);
		dTest = pFirstRegularStep->Advance(dRefTimeStep,
				dCurrTimeStep/dRefTimeStep,
			 	CurrStep, pSM, pNLS, qX, qXPrime, iStIter
#ifdef MBDYN_X_CONVSOL				
				, dSolTest
#endif /* MBDYN_X_CONVSOL */				
				);
	}
	catch (NonlinearSolver::NoConvergence) {
		if (dCurrTimeStep > dMinimumTimeStep) {
			/* Riduce il passo */
			CurrStep = StepIntegrator::REPEATSTEP;
			dCurrTimeStep = NewTimeStep(dCurrTimeStep, 
						iStIter, 
						CurrStep);
			DEBUGCOUT("Changing time step during"
				" first step after "
				<< iStIter << " iterations"
				<< std::endl);
	    		goto IfFirstStepIsToBeRepeated;
	 	}

	    	std::cerr << std::endl
			<< "Maximum iterations number "
			<< iStIter
			<< " has been reached during"
			" first step (time = "
			<< dTime << ");" << std::endl
			<< "time step dt = " << dCurrTimeStep 
			<< " cannot be reduced further;"
			<< std::endl
			<< "aborting ..." << std::endl;
	    	pDM->Output(true);

		THROW(Solver::ErrMaxIterations());
      	}
	catch (NonlinearSolver::ErrSimulationDiverged) {
		/*
		 * Mettere qui eventuali azioni speciali 
		 * da intraprendere in caso di errore ...
		 */

		THROW(SimulationDiverged());
	}
	catch (NonlinearSolver::ConvergenceOnSolution) {
#ifdef MBDYN_X_CONVSOL
		/* FIXME: handle the case of #undef MBDYN_X_CONVSOL */
		bSolConv = true;
#endif /* MBDYN_X_CONVSOL */
	}
	catch (...) {
		THROW(ErrGeneric());
	}

	pDM->AfterConvergence();
   	pDM->Output();
      
#ifdef HAVE_SIGNAL
   	if (!::mbdyn_keep_going) {
      		/* Fa l'output della soluzione al primo passo ed esce */
      		Out << "Interrupted during first step." << std::endl;
      		return;
   	}
#endif /* HAVE_SIGNAL */

	if (outputMsg()) {
      		Out << "Step " << iStep
			<< " time " << dTime+dCurrTimeStep
			<< " step " << dCurrTimeStep
			<< " iterations " << iStIter
			<< " error " << dTest
#ifdef MBDYN_X_CONVSOL
			<< " " << dSolTest
			<< " " << bSolConv
#endif /* MBDYN_X_CONVSOL */
			<< std::endl;
	}

#ifdef MBDYN_X_CONVSOL
	bSolConv = false;
#endif /* MBDYN_X_CONVSOL */

   	dRefTimeStep = dCurrTimeStep;
   	dTime += dRefTimeStep;
   
   	dTotErr += dTest;
   	iTotIter += iStIter;
   
#ifdef __HACK_EIG__  
   	if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
	 	Eig();
		OneEig.fDone = flag(1);
      	}
#endif /* __HACK_EIG__ */

#ifdef USE_RTAI

	if (bRT) {
		/* Need timer */
		if (!mbdyn_rt_is_hard_timer_running() ){
			/* FIXME: ??? */
			mbdyn_rt_set_oneshot_mode();
			mbdyn_start_rt_timer(mbdyn_nano2count(1000000));
		}
			
		if (bRTAllowNonRoot) {
			mbdyn_rt_allow_nonroot_hrt();
		}

		/*
		 * MBDyn can work in two ways:
		 * - internaal timer
		 * - scheduled by an external signal
		 * only the first case is currently implemented
		 */
		if (RTWaitPeriod()) {
			long long t = mbdyn_rt_get_time();
			int r;

			/* Timer should be init'ed */
			ASSERT(t);

			DEBUGCOUT("Task: " << mbdyn_rtai_task 
				<< "; time: " << t 
				<< "; period: " << lRTPeriod << std::endl);
			r = mbdyn_rt_task_make_periodic(mbdyn_rtai_task,
					t, lRTPeriod);

			if (r) {
				std::cerr << "rt_task_make_periodic() failed ("
					<< r << ")" << std::endl;
				THROW(ErrGeneric());
			}
		} else {
			int r;

			/* FIXME: check args 
			 * name should be configurable?
			 * initial value 0: non-blocking
			 */
			r = mbdyn_rt_sem_init("MBDSEM", 0, &RTSemPtr);
			if (r) {
				std::cerr << "rt_sem_init() failed ("
					<< r << ")" << std::endl;
				THROW(ErrGeneric());
			}
		}

		if (bRTHard) {
			/*
			 * FIXME: make hard real time here
			 */
			silent_cout("hard real time not supported yet"
					<< std::endl);
		}

		/* FIXME: should check whether RTStackSize is correclty set? */
		reserve_stack(RTStackSize);
	}
#endif /* USE_RTAI */

    	/* Altri passi regolari */ 
	ASSERT(pRegularSteps != NULL);
	
      	while (1) {
		
		StepIntegrator::StepChange CurrStep 
				= StepIntegrator::NEWSTEP;
	
      		if (dTime >= dFinalTime) {
	 		std::cout << "End of simulation at time "
				<< dTime << " after " 
				<< iStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobians: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl;
			return;
#ifdef HAVE_SIGNAL
      		} else if (!::mbdyn_keep_going) {
	 		std::cout << "Interrupted!" << std::endl
	   			<< "Simulation ended at time "
				<< dTime << " after " 
				<< iStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobians: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl;
	 		return;
#endif /* HAVE_SIGNAL */
      		}
 	
      		iStep++;
      		pDM->BeforePredict(*(qX[0]), *(qXPrime[0]),
				*(qX[1]), *(qXPrime[1]));
	
		Flip();

#ifdef USE_RTAI
		if (bRT) {
			if (RTWaitPeriod()) {
				mbdyn_rt_task_wait_period();
			} else if (RTSemaphore()) {
				/* FIXME: semaphore must be configurable */
				mbdyn_rt_sem_wait(RTSemPtr);
			}
		}
		
#endif /* USE_RTAI */

IfStepIsToBeRepeated:
		try {   	
			pDM->SetTime(dTime+dCurrTimeStep);
			dTest = pRegularSteps->Advance(dRefTimeStep,
					dCurrTimeStep/dRefTimeStep,
				 	CurrStep, pSM, pNLS,
					qX, qXPrime, iStIter
#ifdef MBDYN_X_CONVSOL				
					, dSolTest
#endif /* MBDYN_X_CONVSOL */				
				);
		}

		catch (NonlinearSolver::NoConvergence) {
			if (dCurrTimeStep > dMinimumTimeStep) {
				/* Riduce il passo */
				CurrStep = StepIntegrator::REPEATSTEP;
				dCurrTimeStep = NewTimeStep(dCurrTimeStep,
						iStIter,
						CurrStep);	       
				DEBUGCOUT("Changing time step"
					" during step " 
					<< iStep << " after "
					<< iStIter << " iterations"
					<< std::endl);
				goto IfStepIsToBeRepeated;
	    		} else {
				std::cerr << std::endl
					<< "Maximum iterations number "
					<< iStIter 
					<< " has been reached during"
					" step " << iStep << ';'
					<< std::endl
					<< "time step dt = "
					<< dCurrTimeStep
					<< " cannot be reduced"
					" further;" << std::endl
					<< "aborting ..." << std::endl;
	       			THROW(ErrMaxIterations());
			}
		}		   
		catch (NonlinearSolver::ErrSimulationDiverged) {
			/*
			 * Mettere qui eventuali azioni speciali 
			 * da intraprendere in caso di errore ...
			 */
			THROW(SimulationDiverged());
		}
		catch (NonlinearSolver::ConvergenceOnSolution) {
#ifdef MBDYN_X_CONVSOL
			bSolConv = true;
#endif /* MBDYN_X_CONVSOL */
		}
		catch (...) {
			THROW(ErrGeneric());
		}

		pDM->AfterConvergence();

	      	dTotErr += dTest;
      		iTotIter += iStIter;

      		pDM->Output();     
	
		if (outputMsg()) {	
      			Out << "Step " << iStep 
				<< " time " << dTime+dCurrTimeStep
				<< " step " << dCurrTimeStep
				<< " iterations " << iStIter
				<< " error " << dTest 
#ifdef MBDYN_X_CONVSOL
				<< " " << dSolTest
				<< " " << bSolConv 
#endif /* MBDYN_X_CONVSOL */
				<< std::endl;
		}
      
     	 	DEBUGCOUT("Step " << iStep
			<< " has been completed successfully in "
			<< iStIter << " iterations" << std::endl);
      
	      	dRefTimeStep = dCurrTimeStep;
      		dTime += dRefTimeStep;

#ifdef MBDYN_X_CONVSOL
		bSolConv = false;
#endif /* MBDYN_X_CONVSOL */

#ifdef __HACK_POD__
		if (fPOD && dTime >= pod.dTime) {
			if (++iPODStep == pod.iSteps) {
				/* output degli stati su di una riga */
#ifdef __HACK_POD_BINARY__
	       			PodOut.write((char *)&qX[0], iNumDofs*sizeof(doublereal));
#else /* !__HACK_POD_BINARY__ */
				PodOut << qX[0]->dGetCoef(1);
				for (integer j = 1; j < iNumDofs; j++) {
					PodOut << "  " << qX[0]->dGetCoef(j+1);
                       		}
                       		PodOut << std::endl;
#endif /* __HACK_POD_BINARY__ */
			}
                     	iPODFrames++;
                      	iPODStep = 0;
		}
		if (iPODFrames >= pod.iFrames){
			fPOD = flag(0);
		}                        
#endif /*__HACK_POD__ */

#ifdef __HACK_EIG__
      		if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
			Eig();
			OneEig.fDone = flag(1);
		}
#endif /* __HACK_EIG__ */
      
      		/* Calcola il nuovo timestep */
      		dCurrTimeStep =
			NewTimeStep(dCurrTimeStep, iStIter, CurrStep);
		DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   	}
}

/* Distruttore */
Solver::~Solver(void)
{
   	DEBUGCOUTFNAME("Solver::~Solver");

   	if (sInputFileName != NULL) {
      		SAFEDELETEARR(sInputFileName);
   	}
   
   	if (sOutputFileName != NULL) {
      		SAFEDELETEARR(sOutputFileName);
   	}
   
	for (int ivec = 0; ivec < iNumPreviousVectors; ivec++) {  
   		if (qX[ivec] != NULL) { 
			SAFEDELETE(qX[ivec]);
		}
	}

   	if (pdWorkSpace != NULL) {	
      		SAFEDELETEARR(pdWorkSpace);
   	}

   	if (pDM != NULL) {	
      		SAFEDELETE(pDM);
	}
}
	
/* Nuovo delta t */
doublereal
Solver::NewTimeStep(doublereal dCurrTimeStep,
				 integer iPerformedIters,
				 StepIntegrator::StepChange Why)
{
   	DEBUGCOUTFNAME("Solver::NewTimeStep");

   	switch (CurrStrategy) {
    	case NOCHANGE:
       		return dCurrTimeStep;
      
	case CHANGE:
		return pStrategyChangeDrive->dGet(dTime);
      
    	case FACTOR:
       		if (Why == StepIntegrator::REPEATSTEP) {
	  		if (dCurrTimeStep*StrategyFactor.dReductionFactor 
	      		    >= dMinimumTimeStep) {
	     			if (fLastChance == flag(1)) {
					fLastChance = flag(0);
	     			}
	     			iStepsAfterReduction = 0;
	     			return dCurrTimeStep*StrategyFactor.dReductionFactor;
	  		} else {
	     			if (fLastChance == flag(0)) {
					fLastChance = flag(1);
					return StrategyFactor.dRaiseFactor*dCurrTimeStep;
	     			} else {
					/*
					 * Fuori viene intercettato
					 * il valore illegale
					 */
					return dCurrTimeStep*StrategyFactor.dReductionFactor;
	     			}
	  		}
       		}
       
       		if (Why == StepIntegrator::NEWSTEP) {
	  		iStepsAfterReduction++;
	  		iStepsAfterRaise++;

#if 0
			/* never raise after iPerformedIters > StrategyFactor.iMinIters */
			if (iPerformedIters > iWeightedPerformedIters) {
				iWeightedPerformedIters = iPerformedIters;
			}
#else
			iWeightedPerformedIters = (10*iPerformedIters + 9*iWeightedPerformedIters)/10;
#endif
	  
	  		if (iPerformedIters <= StrategyFactor.iMinIters
	      		    && iStepsAfterReduction > StrategyFactor.iStepsBeforeReduction
			    && iStepsAfterRaise > StrategyFactor.iStepsBeforeRaise
			    && dCurrTimeStep < dMaxTimeStep) {
	     			iStepsAfterRaise = 0;
				iWeightedPerformedIters = 0;
	     			return dCurrTimeStep*StrategyFactor.dRaiseFactor;
	  		}
	  		return dCurrTimeStep;
       		}
       		break;
      
    	default:
       		std::cerr << "You shouldn't have reached this point!" << std::endl;
       		THROW(Solver::ErrGeneric());
   	}
   
   	return dCurrTimeStep;
}
			
			
					
			
const integer iDefaultMaxIterations = 1;
const doublereal dDefaultFictitiousStepsTolerance = 1.e-6;


/* Dati dell'integratore */
void 
Solver::ReadData(MBDynParser& HP)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::ReadData");

   	/* parole chiave */
   	const char* sKeyWords[] = { 
      		"begin",
		"multistep",
		"end",
	
		"initial" "time",
		"final" "time",
		"time" "step",
		"min" "time" "step",
		"max" "time" "step",
		"tolerance",
		"max" "iterations",
	
		/* DEPRECATED */
		"fictitious" "steps" "number",
		"fictitious" "steps" "ratio",      
		"fictitious" "steps" "tolerance",
		"fictitious" "steps" "max" "iterations",
		/* END OF DEPRECATED */

		"dummy" "steps" "number",
		"dummy" "steps" "ratio",
		"dummy" "steps" "tolerance",
		"dummy" "steps" "max" "iterations",
	
		"abort" "after",
			"input",
			"assembly",
			"derivatives",

			/* DEPRECATED */ "fictitious" "steps" /* END OF DEPRECATED */ ,
			"dummy" "steps",

		"output",
			"none",
			"iterations",
			"residual",
			"solution",
			"jacobian",
			"messages",
	
		"method",
		/* DEPRECATED */ "fictitious" "steps" "method" /* END OF DEPRECATED */ ,
		"dummy" "steps" "method",
	
		"Crank" "Nicholson",
			/* DEPRECATED */ "nostro" /* END OF DEPRECATED */ ,
			"ms",
			"hope",
			"bdf",
	
		"derivatives" "coefficient",
		"derivatives" "tolerance",
		"derivatives" "max" "iterations",
	
		"Newton" "Raphson",
			"true",
			"modified",
		
		"strategy",
			"factor",
			"no" "change",
			"change",
		
		"pod",
		"eigen" "analysis",
		"output" "modes",
		
		"solver",
		"interface" "solver", 
		
			"harwell",
			"meschach",
			"y12",
			"umfpack",
			"umfpack3",
			"empty",
		"matrix" "free",
		"bicgstab",
		"gmres",
		"preconditioner",
		"full" "jacobian",

		/* RTAI stuff */
		"real" "time",

		NULL
   	};
   
   	/* enum delle parole chiave */
   	enum KeyWords {
      		UNKNOWN = -1,
		BEGIN = 0,
		MULTISTEP,
		END,
	
		INITIALTIME,
		FINALTIME,
		TIMESTEP,
		MINTIMESTEP,
		MAXTIMESTEP,
		TOLERANCE,
		MAXITERATIONS,
	
		FICTITIOUSSTEPSNUMBER,
		FICTITIOUSSTEPSRATIO,      
		FICTITIOUSSTEPSTOLERANCE,
		FICTITIOUSSTEPSMAXITERATIONS,

		DUMMYSTEPSNUMBER,
		DUMMYSTEPSRATIO,
		DUMMYSTEPSTOLERANCE,
		DUMMYSTEPSMAXITERATIONS,
	
		ABORTAFTER,
		INPUT,
		ASSEMBLY,
		DERIVATIVES,
		FICTITIOUSSTEPS,
		DUMMYSTEPS,

		OUTPUT,
			NONE,
			ITERATIONS,
			RESIDUAL,
			SOLUTION,
			JACOBIAN,
			MESSAGES,
	
		METHOD,
		FICTITIOUSSTEPSMETHOD,
		DUMMYSTEPSMETHOD,
		CRANKNICHOLSON,
		NOSTRO, 
		MS,
		HOPE,
		BDF,
	
		DERIVATIVESCOEFFICIENT,
		DERIVATIVESTOLERANCE,
		DERIVATIVESMAXITERATIONS,
	
		NEWTONRAPHSON,
		NR_TRUE,
		MODIFIED,

		STRATEGY,
		STRATEGYFACTOR,
		STRATEGYNOCHANGE,
		STRATEGYCHANGE,
	
		POD,
		EIGENANALYSIS,
		OUTPUTMODES,
		
		SOLVER,
		INTERFACESOLVER,
		HARWELL,
		MESCHACH,
		Y12,
		UMFPACK,
		UMFPACK3,
		EMPTY,
		
		MATRIXFREE,
		BICGSTAB,
		GMRES,
		PRECONDITIONER,
		FULLJACOBIAN,

		/* RTAI stuff */
		REALTIME,
	
		LASTKEYWORD
   	};
   
   	/* tabella delle parole chiave */
   	KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   	/* cambia la tabella del parser */
   	HP.PutKeyTable(K);

   	/* legge i dati della simulazione */
   	if (KeyWords(HP.GetDescription()) != BEGIN) {
      		std::cerr << std::endl << "Error: <begin> expected at line " 
			<< HP.GetLineData() << "; aborting ..." << std::endl;
      		THROW(ErrGeneric());
   	}
   
   	if (KeyWords(HP.GetWord()) != MULTISTEP) {
      		std::cerr << std::endl << "Error: <begin: multistep;> expected at line " 
			<< HP.GetLineData() << "; aborting ..." << std::endl;
      		THROW(ErrGeneric());
   	}

   	flag fMethod(0);
   	flag fFictitiousStepsMethod(0);      
	
	/* dati letti qui ma da passare alle classi
	 *	StepIntegration e NonlinearSolver
	 */ 
	
	doublereal dTol = dDefaultTol;
   	doublereal dSolutionTol = 0.;
   	integer iMaxIterations = iDefaultMaxIterations;

        /* Dati dei passi fittizi di trimmaggio iniziale */
   	doublereal dFictitiousStepsTolerance = dDefaultFictitiousStepsTolerance;
   	integer iFictitiousStepsMaxIterations = iDefaultMaxIterations;

   	/* Dati del passo iniziale di calcolo delle derivate */
   	
	doublereal dDerivativesTol = dDefaultTol;
   	doublereal dDerivativesCoef = 1.;   
   	integer iDerivativesMaxIterations = iDefaultMaxIterations;
	

	DriveCaller* pRhoRegular = NULL;
	DriveCaller* pRhoAlgebraicRegular = NULL;	
	DriveCaller* pRhoFictitious = NULL;
	DriveCaller* pRhoAlgebraicFictitious = NULL;
	
	enum StepIntegratorType {
			INT_CRANKNICHOLSON,
			INT_MS2,
			INT_HOPE,
			INT_UNKNOWN
	};
	
	StepIntegratorType RegularType = INT_UNKNOWN, FictitiousType = INT_UNKNOWN; 
	
	
	 	
   	/* Ciclo infinito */
   	while (1) {	
      		KeyWords CurrKeyWord = KeyWords(HP.GetDescription());
      
      		switch (CurrKeyWord) {	 	 
       		case INITIALTIME:
	  		dInitialTime = HP.GetReal();
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Initial time is "
				   << dInitialTime << std::endl);
	  		break;
	 
       		case FINALTIME:
	  		dFinalTime = HP.GetReal();
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Final time is "
				   << dFinalTime << std::endl);
				   
	  		if(dFinalTime <= dInitialTime) {
	     			std::cerr << "warning: final time " << dFinalTime
	       				<< " is less than initial time "
					<< dInitialTime << ';' << std::endl
	       				<< "this will cause the simulation"
					" to abort" << std::endl;
			}
	  		break;
	 
       		case TIMESTEP:
	  		dInitialTimeStep = HP.GetReal();
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Initial time step is "
				   << dInitialTimeStep << std::endl);
	  
	  		if (dInitialTimeStep == 0.) {
	     			std::cerr << "warning, null initial time step"
					" is not allowed" << std::endl;
	  		} else if (dInitialTimeStep < 0.) {
	     			dInitialTimeStep = -dInitialTimeStep;
				std::cerr << "warning, negative initial time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dInitialTimeStep 
					<< " will be considered" << std::endl;
			}
			break;
	 
       		case MINTIMESTEP:
	  		dMinimumTimeStep = HP.GetReal();
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Minimum time step is "
				   << dMinimumTimeStep << std::endl);
	  
	  		if (dMinimumTimeStep == 0.) {
	     			std::cerr << "warning, null minimum time step"
					" is not allowed" << std::endl;
	     			THROW(ErrGeneric());
			} else if (dMinimumTimeStep < 0.) {
				dMinimumTimeStep = -dMinimumTimeStep;
				std::cerr << "warning, negative minimum time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dMinimumTimeStep 
					<< " will be considered" << std::endl;
	  		}
	  		break;
	
       		case MAXTIMESTEP:
	  		dMaxTimeStep = HP.GetReal();
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Max time step is "
				   << dMaxTimeStep << std::endl);
	  
	  		if (dMaxTimeStep == 0.) {
				std::cout << "no max time step limit will be"
					" considered" << std::endl;
			} else if (dMaxTimeStep < 0.) {
				dMaxTimeStep = -dMaxTimeStep;
				std::cerr << "warning, negative max time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dMaxTimeStep 
					<< " will be considered" << std::endl;
	  		}
	  		break;
	 
       		case FICTITIOUSSTEPSNUMBER:
       		case DUMMYSTEPSNUMBER:
	  		iFictitiousStepsNumber = HP.GetInt();
			if (iFictitiousStepsNumber < 0) {
				iFictitiousStepsNumber = 
					iDefaultFictitiousStepsNumber;
				std::cerr << "warning, negative dummy steps number"
					" is illegal;" << std::endl
					<< "resorting to default value "
					<< iDefaultFictitiousStepsNumber
					<< std::endl;		       
			} else if (iFictitiousStepsNumber == 1) {
				std::cerr << "warning, a single dummy step"
					" may be useless" << std::endl;
	  		}
	  
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps number: " 
		     		   << iFictitiousStepsNumber << std::endl);
	  		break;
	 
       		case FICTITIOUSSTEPSRATIO:
       		case DUMMYSTEPSRATIO:
	  		dFictitiousStepsRatio = HP.GetReal();
	  		if (dFictitiousStepsRatio < 0.) {
	     			dFictitiousStepsRatio =
					dDefaultFictitiousStepsRatio;
				std::cerr << "warning, negative dummy steps ratio"
					" is illegal;" << std::endl
					<< "resorting to default value "
					<< dDefaultFictitiousStepsRatio
					<< std::endl;		       
			}
			
	  		if (dFictitiousStepsRatio > 1.) {
				std::cerr << "warning, dummy steps ratio"
					" is larger than one." << std::endl
					<< "Something like 1.e-3 should"
					" be safer ..." << std::endl;
	  		}
	  
	  		DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps ratio: " 
		     		   << dFictitiousStepsRatio << std::endl);
	  		break;
	 
       		case FICTITIOUSSTEPSTOLERANCE:
       		case DUMMYSTEPSTOLERANCE:
	  		dFictitiousStepsTolerance = HP.GetReal();
	  		if (dFictitiousStepsTolerance <= 0.) {
				dFictitiousStepsTolerance =
					dDefaultFictitiousStepsTolerance;
				std::cerr << "warning, negative dummy steps"
					" tolerance is illegal;" << std::endl
					<< "resorting to default value "
					<< dDefaultFictitiousStepsTolerance
					<< std::endl;		       
	  		}
			DEBUGLCOUT(MYDEBUG_INPUT,
				   "Fictitious steps tolerance: "
		     		   << dFictitiousStepsTolerance << std::endl);
	  		break;
	 
       		case ABORTAFTER: {
	  		KeyWords WhenToAbort(KeyWords(HP.GetWord()));
	  		switch (WhenToAbort) {
	   		case INPUT:
	      			fAbortAfterInput = flag(1);
	      			DEBUGLCOUT(MYDEBUG_INPUT, 
			 		"Simulation will abort after"
					" data input" << std::endl);
	      			break;
			
	   		case ASSEMBLY:
	     			fAbortAfterAssembly = flag(1);
	      			DEBUGLCOUT(MYDEBUG_INPUT,
			 		   "Simulation will abort after"
					   " initial assembly" << std::endl);
	      			break;	  
	     
	   		case DERIVATIVES:
	      			fAbortAfterDerivatives = flag(1);
	      			DEBUGLCOUT(MYDEBUG_INPUT, 
			 		   "Simulation will abort after"
					   " derivatives solution" << std::endl);
	      			break;
	     
	   		case FICTITIOUSSTEPS:
	   		case DUMMYSTEPS:
	      			fAbortAfterFictitiousSteps = flag(1);
	      			DEBUGLCOUT(MYDEBUG_INPUT, 
			 		   "Simulation will abort after"
					   " dummy steps solution" << std::endl);
	      			break;
	     
	   		default:
	      			std::cerr << std::endl 
					<< "Don't know when to abort,"
					" so I'm going to abort now" << std::endl;
	      			THROW(ErrGeneric());
	  		}
	  		break;
       		}

		case OUTPUT: {
			unsigned OF = OUTPUT_NONE;

			while (HP.fIsArg()) {
				KeyWords OutputFlag(KeyWords(HP.GetWord()));
				switch (OutputFlag) {
				case NONE:
					OF = OUTPUT_NONE;
					break;

				case ITERATIONS:
					OF |= OUTPUT_ITERS;
					break;

				case RESIDUAL:
					OF |= OUTPUT_RES;
					break;

				case SOLUTION:
					OF |= OUTPUT_SOL;
					break;

				case JACOBIAN:
					OF |= OUTPUT_JAC;
					break;

				case MESSAGES:
					OF |= OUTPUT_MSG;
					break;

				default:
					std::cerr << "Unknown output flag "
						"at line " << HP.GetLineData() 
						<< "; ignored" << std::endl;
					break;
				}
			}

			SetOutputFlags(OF);
			
			break;
		}
	 
       		case METHOD: {
	  		if (fMethod) {
	     			std::cerr << "error: multiple definition"
					" of integration method at line "
					<< HP.GetLineData() << std::endl;
	     			THROW(ErrGeneric());
	  		}
	  		fMethod = flag(1);
	        	  
	  		KeyWords KMethod = KeyWords(HP.GetWord());
	  		switch (KMethod) {
	   		case CRANKNICHOLSON:
				RegularType = INT_CRANKNICHOLSON;				
	      			break;
				
			case BDF: {
				SAFENEWWITHCONSTRUCTOR(pRhoRegular,
						NullDriveCaller, 
						NullDriveCaller(NULL));
				SAFENEWWITHCONSTRUCTOR(pRhoAlgebraicRegular,
						NullDriveCaller, 
						NullDriveCaller(NULL));
				RegularType = INT_MS2;
		  		break;
			}
			
	   		case NOSTRO:
				  silent_cerr("integration method \"nostro\" "
						  "is deprecated; use \"ms\" "
						  "instead at line "
						  << HP.GetLineData()
						  << std::endl);
	   		case MS:
	   		case HOPE: {
	      			pRhoRegular =
					ReadDriveData(NULL, HP, NULL);
				HP.PutKeyTable(K);

	      			pRhoAlgebraicRegular = NULL;
				if (HP.fIsArg()) {
					pRhoAlgebraicRegular = ReadDriveData(NULL, 
							HP, NULL);
					HP.PutKeyTable(K);
				} else {
					pRhoAlgebraicRegular = pRhoRegular->pCopy();
				}
			
	      			switch (KMethod) {
	       			case NOSTRO: 
	       			case MS:
					RegularType = INT_MS2;
		  			break;
		 
	       			case HOPE:
					RegularType = INT_HOPE;	      
		  			break;
					
	       			default:
	          			THROW(ErrGeneric());
	      			}
	      			break;
	   		}
	   		default:
	      			std::cerr << "Unknown integration method at line "
					<< HP.GetLineData() << std::endl;
				THROW(ErrGeneric());
	  		}
	  		break;
       		}

		case FICTITIOUSSTEPSMETHOD:
		case DUMMYSTEPSMETHOD: {
			if (fFictitiousStepsMethod) {
				std::cerr << "error: multiple definition "
					"of dummy steps integration method "
					"at line " << HP.GetLineData()
					<< std::cerr;
				THROW(ErrGeneric());
			}
			fFictitiousStepsMethod = flag(1);	  	
 
			KeyWords KMethod = KeyWords(HP.GetWord());
			switch (KMethod) {
			case CRANKNICHOLSON:
				FictitiousType = INT_CRANKNICHOLSON; 
				break;

			case BDF: 
				SAFENEWWITHCONSTRUCTOR(pRhoFictitious,
					NullDriveCaller,
					NullDriveCaller(NULL));
				SAFENEWWITHCONSTRUCTOR(pRhoAlgebraicFictitious,
					NullDriveCaller,
					NullDriveCaller(NULL));
				FictitiousType = INT_MS2;
				break;
  
			case NOSTRO:
			case MS:
			case HOPE: 	      	     
				pRhoFictitious = ReadDriveData(NULL, HP, NULL);
				HP.PutKeyTable(K);

				if (HP.fIsArg()) {
					pRhoAlgebraicFictitious = ReadDriveData(NULL, HP, NULL);
					HP.PutKeyTable(K);
				} else {
					pRhoAlgebraicFictitious = pRhoFictitious->pCopy();
				}
				HP.PutKeyTable(K);
   
				switch (KMethod) {
				case NOSTRO:
				case MS: 
					FictitiousType = INT_MS2;     
					break;
				
				case HOPE: 
					FictitiousType = INT_HOPE;    	      
					break;

	       			default:
	          			THROW(ErrGeneric());
				}
	      			break;	      
	   		
			default: {
				std::cerr << "Unknown integration method at line " << HP.GetLineData() << std::endl;
				THROW(ErrGeneric());
			}
			}	     
			break;
		} 

		case TOLERANCE: {
			dTol = HP.GetReal();
			if (dTol <= 0.) {
				dTol = dDefaultTol;
				std::cerr << "warning, tolerance <= 0. is illegal; "
					"using default value " << dTol << std::endl;
			}
#ifdef MBDYN_X_CONVSOL
			dSolutionTol = dTol;
			if (HP.fIsArg()) {
				dSolutionTol = HP.GetReal();
			}
			if (dSolutionTol <= 0.) {
				dSolutionTol = 0.;
				std::cerr << "warning, tolerance <= 0. is illegal; "
					"switching to default value " << dSolutionTol 
					<< std::endl;
			}

			DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dTol
					<< ", " << dSolutionTol << std::endl);
#else /* !MBDYN_X_CONVSOL */
			if (HP.fIsArg()) {
				pedantic_cerr("define MBDYN_X_SOLCONV to enable "
						"convergence test on solution" << std::endl);
				(void)HP.GetReal();
			}
#endif /* !MBDYN_X_CONVSOL */
			DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dTol << std::endl);
			break;
		}	

	 
       case DERIVATIVESTOLERANCE: {
	  dDerivativesTol = HP.GetReal();
	  if (dDerivativesTol <= 0.) {
	     dDerivativesTol = 1e-6;
	     std::cerr
	       << "warning, derivatives tolerance <= 0.0 "
	       "is illegal; switching to default value " 
	       << dDerivativesTol
	       << std::endl;
	  }		       		  
	  DEBUGLCOUT(MYDEBUG_INPUT, 
		     "Derivatives tolerance = " << dDerivativesTol
		     << std::endl);
	  break;
       }
	 
       case MAXITERATIONS: {
	  iMaxIterations = HP.GetInt();
	  if (iMaxIterations < 1) {
	     iMaxIterations = iDefaultMaxIterations;
	     std::cerr 
	       << "warning, max iterations < 1 is illegal; "
	       "switching to default value "
	       << iMaxIterations
	       << std::endl;
	  }		       		  
	  DEBUGLCOUT(MYDEBUG_INPUT, 
		     "Max iterations = " << iMaxIterations << std::endl);
	  break;
       }
	 
       case DERIVATIVESMAXITERATIONS: {
	  iDerivativesMaxIterations = HP.GetInt();
	  if (iDerivativesMaxIterations < 1) {
	     iDerivativesMaxIterations = iDefaultMaxIterations;
	     std::cerr 
	       << "warning, derivatives max iterations < 1 is illegal; "
	       "switching to default value "
	       << iDerivativesMaxIterations
	       << std::endl;
	  }		       		  
	  DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives max iterations = " 
		    << iDerivativesMaxIterations << std::endl);
	  break;
       }	     	  	       
	 
       case FICTITIOUSSTEPSMAXITERATIONS:
       case DUMMYSTEPSMAXITERATIONS: {
	  iFictitiousStepsMaxIterations = HP.GetInt();
	  if (iFictitiousStepsMaxIterations < 1) {
	     iFictitiousStepsMaxIterations = iDefaultMaxIterations;
	     std::cerr 
	       << "warning, dummy steps max iterations < 1 is illegal;"
	       " switching to default value "
	       << iFictitiousStepsMaxIterations
	       << std::endl;
	  }
	  DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps max iterations = " 
		    << iFictitiousStepsMaxIterations << std::endl);
	  break;
       }	     	  	       
	 
       case DERIVATIVESCOEFFICIENT: {
	  dDerivativesCoef = HP.GetReal();
	  if (dDerivativesCoef <= 0.) {
	     dDerivativesCoef = 1.;
	     std::cerr 
	       << "warning, derivatives coefficient <= 0. is illegal; "
	       "switching to default value "
	       << dDerivativesCoef
	       << std::endl;
	  }
	  DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives coefficient = "
		    << dDerivativesCoef << std::endl);
	  break;
       }	     	  	       
	 
       case NEWTONRAPHSON: {
	  KeyWords NewRaph(KeyWords(HP.GetWord()));
	  switch(NewRaph) {
	   case MODIFIED: {
	      fTrueNewtonRaphson = 0;
	      if (HP.fIsArg()) {
		 iIterationsBeforeAssembly = HP.GetInt();
	      } else {
		 iIterationsBeforeAssembly = iDefaultIterationsBeforeAssembly;
	      }
	      DEBUGLCOUT(MYDEBUG_INPUT, 
			 "Modified Newton-Raphson will be used;" << std::endl
			 << "matrix will be assembled at most after " 
			 << iIterationsBeforeAssembly
			 << " iterations" << std::endl);
	      break;
	   }
	   default: {
	      std::cerr 
		<< "warning: unknown case; resorting to default" 
		<< std::endl;
	      /* Nota: non c'e' break; 
	       * cosi' esegue anche il caso NR_TRUE */
	   }		       
	   case NR_TRUE: {
	      fTrueNewtonRaphson = 1;
	      iIterationsBeforeAssembly = 0;
	      break;
	   }		 
	  }		  		  
	  break;
       }

       case END: {	     
	  if (KeyWords(HP.GetWord()) != MULTISTEP) {
	     std::cerr << std::endl 
	       << "Error: <end: multistep;> expected at line " 
	       << HP.GetLineData() << "; aborting ..." << std::endl;
	     THROW(ErrGeneric());
	  }
	  goto EndOfCycle;
       }
	 
       case STRATEGY: {
	  switch (KeyWords(HP.GetWord())) {
	     
	   case STRATEGYFACTOR: {
	      CurrStrategy = FACTOR;

	      /*
	       *	strategy: factor ,
	       *		<reduction factor> ,
	       *		<steps before reduction> ,
	       *		<raise factor> ,
	       *		<steps before raise> ,
	       *		<min iterations> ;
	       */
	      
	      StrategyFactor.dReductionFactor = HP.GetReal();
	      if (StrategyFactor.dReductionFactor >= 1.) {
		 std::cerr << "warning, illegal reduction factor at line "
		   << HP.GetLineData() 
		   << "; default value 1. (no reduction) will be used"
		   << std::endl;
		 StrategyFactor.dReductionFactor = 1.;
	      }
	      
	      StrategyFactor.iStepsBeforeReduction = HP.GetInt();
	      if (StrategyFactor.iStepsBeforeReduction <= 0) {
		 std::cerr << "Warning, illegal number of steps "
		   "before reduction at line "
		   << HP.GetLineData() << ';' << std::endl
		   << "default value 1 will be used (it may be dangerous)" 
		   << std::endl;
		 StrategyFactor.iStepsBeforeReduction = 1;
	      }
	      
	      StrategyFactor.dRaiseFactor = HP.GetReal();
	      if (StrategyFactor.dRaiseFactor <= 1.) {
		 std::cerr << "warning, illegal raise factor at line "
		   << HP.GetLineData() 
		   << "; default value 1. (no raise) will be used"
		   << std::endl;
		 StrategyFactor.dRaiseFactor = 1.;
	      }
	      
	      StrategyFactor.iStepsBeforeRaise = HP.GetInt();
	      if (StrategyFactor.iStepsBeforeRaise <= 0) {
		 std::cerr << "Warning, illegal number of steps "
		   "before raise at line "
		   << HP.GetLineData() << ';' << std::endl
		   << "default value 1 will be used (it may be dangerous)" 
		   << std::endl;
		 StrategyFactor.iStepsBeforeRaise = 1;
	      }
	      
	      StrategyFactor.iMinIters = HP.GetInt();
	      if (StrategyFactor.iMinIters <= 0) {
		 std::cerr << "Warning, illegal minimum number "
		   "of iterations at line "
		   << HP.GetLineData() << ';' << std::endl
		   << "default value 0 will be used (never raise)" 
		   << std::endl;
		 StrategyFactor.iMinIters = 1;
	      }
	      
	      DEBUGLCOUT(MYDEBUG_INPUT,
			 "Time step control strategy: Factor" << std::endl
			 << "Reduction factor: " 
			 << StrategyFactor.dReductionFactor
			 << "Steps before reduction: "
			 << StrategyFactor.iStepsBeforeReduction
			 << "Raise factor: "
			 << StrategyFactor.dRaiseFactor
			 << "Steps before raise: "
			 << StrategyFactor.iStepsBeforeRaise 
			 << "Min iterations: "
			 << StrategyFactor.iMinIters << std::endl);
	      break;  
	   }		 
	     
	   case STRATEGYNOCHANGE: {
	      CurrStrategy = NOCHANGE;
	      break;
	   }
	   
	   case STRATEGYCHANGE: {
	      CurrStrategy = CHANGE;
	      pStrategyChangeDrive = ReadDriveData(NULL, HP, NULL);
	      HP.PutKeyTable(K);	      
	      break;
	   }
	     
	   default: {
	      std::cerr << "Unknown time step control strategy at line "
		<< HP.GetLineData() << std::endl;
	      THROW(ErrGeneric());
	   }		 
	  }
	  
	  break;
       }
	 
       case POD:
#ifdef __HACK_POD__
              pod.dTime = HP.GetReal();

	      pod.iSteps = 1;
	      if (HP.fIsArg()) {
              	 pod.iSteps = HP.GetInt();
	      }

	      pod.iFrames = (unsigned int)(-1);
	      if (HP.fIsArg()) {
              	 pod.iFrames = HP.GetInt();
	      }

              fPOD = flag(1);
              DEBUGLCOUT(MYDEBUG_INPUT, "POD analysis will be performed "
			      "since time " << pod.dTime
			      << " for " << pod.iFrames << " frames "
			      << " every " << pod.iSteps << " steps " 
			      << std::endl);
#else/* !__HACK_POD__ */
              std::cerr << "line " << HP.GetLineData()
                   << ": POD analysis not supported (ignored)" << std::endl;
	      for (; HP.fIsArg();) {
		      (void)HP.GetReal();
	      }
#endif /* !__HACK_POD__ */
              break;
 
       case EIGENANALYSIS: {
#ifdef __HACK_EIG__
	  OneEig.dTime = HP.GetReal();
	  if (HP.IsKeyWord("parameter")) {
             dEigParam = HP.GetReal();
	  }
	  OneEig.fDone = flag(0);
	  fEigenAnalysis = flag(1);
	  DEBUGLCOUT(MYDEBUG_INPUT, "Eigenanalysis will be performed at time "
	  	     << OneEig.dTime << " (parameter: " << dEigParam << ")" 
		     << std::endl);
	  if (HP.IsKeyWord("outputmatrices")) {
	     fEigenAnalysis = flag(2);
	  }
#else /* !__HACK_EIG__ */
	  HP.GetReal();
	  if (HP.IsKeyWord("parameter")) {
	     HP.GetReal();
	  }
	  if (HP.IsKeyWord("outputmatrices")) {
	     NO_OP;
	  }
	  std::cerr << HP.GetLineData()
	    << ": eigenanalysis not supported (ignored)" << std::endl;
#endif /* !__HACK_EIG__ */
	  break;
       }

       case OUTPUTMODES:
#ifndef __HACK_EIG__
	  std::cerr << "line " << HP.GetLineData()
	    << ": warning, no eigenvalue support available" << std::endl;
#endif /* !__HACK_EIG__ */
	  if (HP.IsKeyWord("yes") || HP.IsKeyWord("nastran")) {
#ifdef __HACK_EIG__
	     fOutputModes = flag(1);
	     if (HP.fIsArg()) {
		dUpperFreq = HP.GetReal();
		if (dUpperFreq < 0.) {
			std::cerr << "line "<< HP.GetLineData()
				<< ": illegal upper frequency limit " 
				<< dUpperFreq << "; using " 
				<< -dUpperFreq << std::endl;
			dUpperFreq = -dUpperFreq;
		}
		if (HP.fIsArg()) {
			dLowerFreq = HP.GetReal();
			if (dLowerFreq > dUpperFreq) {
				std::cerr << "line "<< HP.GetLineData()
					<< ": illegal lower frequency limit "
					<< dLowerFreq 
					<< " higher than upper frequency limit; using " 
					<< 0. << std::endl;
				dLowerFreq = 0.;
			}
		}
	     }
#endif /* !__HACK_EIG__ */
	  } else if (HP.IsKeyWord("no")) {
#ifdef __HACK_EIG__
	     fOutputModes = flag(0);
#endif /* !__HACK_EIG__ */
	  } else {
	     std::cerr << HP.GetLineData()
	       << ": unknown mode output flag (should be { yes | no })"
	       << std::endl;
	  }
	  break;

       case SOLVER: {
	  switch(KeyWords(HP.GetWord())) {	     
	   case MESCHACH:
#ifdef USE_MESCHACH
	     CurrSolver = LinSol::MESCHACH_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Using meschach sparse LU solver" << std::endl);
	     break;
#endif /* USE_MESCHACH */

	   case Y12:
#ifdef USE_Y12
             CurrSolver = LinSol::Y12_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using y12 sparse LU solver" << std::endl);
	     break;
#endif /* USE_Y12 */
							       
	   case UMFPACK3:
	     pedantic_cerr("\"umfpack3\" is deprecated; "
			     "use \"umfpack\" instead" << std::endl);
	   case UMFPACK:
#ifdef USE_UMFPACK
             CurrSolver = LinSol::UMFPACK_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using umfpack sparse LU solver" << std::endl);
	     break;
#endif /* USE_UMFPACK */

	   case HARWELL: 
#ifdef USE_HARWELL
	     CurrSolver = LinSol::HARWELL_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Using harwell sparse LU solver" << std::endl);	 
	     break;	   
#endif /* USE_HARWELL */

	   case EMPTY: 
	     CurrSolver = LinSol::EMPTY_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"No LU solver" << std::endl);	 
	     break;	   

	   default:
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Unknown solver; switching to default" << std::endl);
	     break;
	  }
	  
	  if (HP.IsKeyWord("workspacesize")) {
	     iWorkSpaceSize = HP.GetInt();
	     if (iWorkSpaceSize < 0) {
		iWorkSpaceSize = 0;
	     }
	  }
	  
	  if (HP.IsKeyWord("pivotfactor")) {
	     dPivotFactor = HP.GetReal();
	     if (dPivotFactor <= 0. || dPivotFactor > 1.) {
		dPivotFactor = 1.;
	     }
	  }
	  
	  DEBUGLCOUT(MYDEBUG_INPUT, "Workspace size: " << iWorkSpaceSize 
		    << ", pivor factor: " << dPivotFactor << std::endl);
	  break;
       }

       case INTERFACESOLVER: {
#ifdef USE_MPI
		switch(KeyWords(HP.GetWord())) {           
		case MESCHACH:
#ifdef USE_MESCHACH
	  		CurrIntSolver = MESCHACH_SOLVER;
	  		DEBUGLCOUT(MYDEBUG_INPUT, 
				"Using meschach sparse LU solver" << std::endl);
	  		break;
#endif /* USE_MESCHACH */

		case Y12:
#ifdef USE_Y12
	  		CurrIntSolver = LinSol::Y12_SOLVER;
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using y12 sparse LU solver" << std::endl);
	  		break;
#endif /* USE_Y12 */

		case UMFPACK3:
	   		pedantic_cerr("\"umfpack3\" is deprecated; "
	   				"use \"umfpack\" instead" 
					<< std::endl);
		case UMFPACK:
#ifdef USE_UMFPACK
	  		CurrIntSolver = LinSol::UMFPACK_SOLVER;
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using umfpack sparse LU solver" << std::endl);
	  		break;
#endif /* USE_UMFPACK */

		case HARWELL: 
#ifdef USE_HARWELL
	  		CurrIntSolver = LinSol::MESCHACH_SOLVER;
	  		std::cerr << "Harwell solver cannot be used "
				"as interface solver; Meschach will be used"
				<< std::endl;
	  		break;
#endif /* USE_HARWELL */                       

	   	case EMPTY: 
	     		CurrSolver = LinSol::EMPTY_SOLVER;
	     		DEBUGLCOUT(MYDEBUG_INPUT, 
				"No LU interface solver" << std::endl);	 
	     		break;	   
	
		default:
	  		DEBUGLCOUT(MYDEBUG_INPUT, 
				"Unknown solver; switching to default" << std::endl);
	  		break;
        	}
	
		if (HP.IsKeyWord("workspacesize")) {
			iIWorkSpaceSize = HP.GetInt();
			if (iIWorkSpaceSize < 0) {
				iIWorkSpaceSize = 0;
			}
		}

		if (HP.IsKeyWord("pivotfactor")) {
			dIPivotFactor = HP.GetReal();
			if (dIPivotFactor <= 0. || dIPivotFactor > 1.) {
				dIPivotFactor = 1.;
			}
		}
	
		DEBUGLCOUT(MYDEBUG_INPUT, "Workspace size: " << iIWorkSpaceSize 
			<< ", pivor factor: " << dIPivotFactor << std::endl);
		break;
#else /* !USE_MPI */
		std::cerr << "Interface solver only allowed when compiled "
			"with MPI support" << std::endl;
		THROW(ErrGeneric());
#endif /* !USE_MPI */
       }
       
       case MATRIXFREE: {
		fMatrixFree = 1;

		switch(KeyWords(HP.GetWord())) {           
		
			case BICGSTAB:
				MFSolverType = MatrixFreeSolver::BICGSTAB;
				break;
			
			case GMRES:
				MFSolverType = MatrixFreeSolver::GMRES;
				break;
		
			default:
				std::cerr << "Unknown iterative solver "
					"at line " << HP.GetLineData()
					<< std::endl;
				THROW(ErrGeneric());
		}	
	       
		if (HP.IsKeyWord("tolerance")) {
			dIterTol = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, "Inner Iterative Solver Tolerance: " 
				<< dIterTol << std::endl);
		}
		if (HP.IsKeyWord("steps")) {
			iIterativeMaxSteps = HP.GetInt();
			DEBUGLCOUT(MYDEBUG_INPUT, "Maximum Number of Inner Steps for Iterative Solver : " 
				<< iIterativeMaxSteps << std::endl);
		}
		if (HP.IsKeyWord("eta")) {
			dIterertiveEtaMax = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, "Maximum Eta Coefficient for Iterative Solver : " 
				<< dIterertiveEtaMax << std::endl);
		}
		break;
	}
	
	case PRECONDITIONER: {
		KeyWords KPrecond = KeyWords(HP.GetWord());
	  	switch (KPrecond) {	
			
		case FULLJACOBIAN: 
			PcType = Preconditioner::FULLJACOBIAN;		
			if (HP.IsKeyWord("steps")) {
				iPrecondSteps = HP.GetInt();
				DEBUGLCOUT(MYDEBUG_INPUT, "Number of Steps before recomputing the preconditioner : " 
					<< iPrecondSteps << std::endl);
			}
			break;
			
		default:
			std::cerr << "Unkown preconditioner at line "
				<< HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}
		
		break;
        }

       case REALTIME: {
#ifdef USE_RTAI
         bRT = true;

	 if (HP.IsKeyWord("allow" "nonroot")) {
	    bRTAllowNonRoot = true;
	 }

	 /* FIXME: use a sae default? */
	 if (HP.IsKeyWord("mode")) {
	    if (HP.IsKeyWord("period")) {
	       RTMode = MBRTAI_WAITPERIOD;
	    } else if (HP.IsKeyWord("semaphore")) {
	       /* FIXME: not implemented yet ... */
	       RTMode = MBRTAI_SEMAPHORE;
	    } else {
		std::cerr << "unknown realtime mode at line "
			<< HP.GetLineData() << std::endl;
		THROW(ErrGeneric());
	    }
	 }

	 if (HP.IsKeyWord("time" "step")) {
	    long long p = HP.GetInt();

	    if (p <= 0) {
               std::cerr << "illegal time step " << p << " at line "
		 << HP.GetLineData() << std::endl;
	       THROW(ErrGeneric());
	    }

	    lRTPeriod = mbdyn_nano2count(p);
	 } else {
	    std::cerr << "need a time step for real time at line "
	      << HP.GetLineData() << std::endl;
	    THROW(ErrGeneric());
	 }

	 if (HP.IsKeyWord("reserve" "stack")) {
	    long size = HP.GetInt();

	    if (size <= 0) {
	       std::cerr << "illegal stack size " << size << " at line "
		       << HP.GetLineData() << std::endl;
	       THROW(ErrGeneric());
	    }

	    RTStackSize = size;
	 }

#else /* !USE_RTAI */
         std::cerr << "need to configure --with-rtai to use realtime" << std::endl;
	 THROW(ErrGeneric());
#endif /* !USE_RTAI */
	 break;
       }

       default: {
	  std::cerr << std::endl << "Unknown description at line " 
	    << HP.GetLineData() << "; aborting ..." << std::endl;
	  THROW(ErrGeneric());
       }	
      }   
   }
   
EndOfCycle: /* esce dal ciclo di lettura */
      
   if (dFinalTime < dInitialTime) {      
      fAbortAfterAssembly = flag(1);
   }   
   
   if (dFinalTime == dInitialTime) {      
      fAbortAfterDerivatives = flag(1);
   }   
      
   /* Metodo di integrazione di default */
   if (!fMethod || pRhoRegular == NULL) {
      ASSERT(RegularType == INT_UNKNOWN);
      
 
      SAFENEWWITHCONSTRUCTOR(pRhoRegular,
		             NullDriveCaller,
			     NullDriveCaller(NULL));
      
      /* DriveCaller per Rho asintotico per variabili algebriche */     
      SAFENEWWITHCONSTRUCTOR(pRhoAlgebraicRegular,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
      RegularType = INT_MS2; 	      
   }

   /* Metodo di integrazione di default */
   if (!fFictitiousStepsMethod || pRhoFictitious == NULL) {
      ASSERT(FictitiousType == INT_UNKNOWN);
      
      SAFENEWWITHCONSTRUCTOR(pRhoFictitious,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
                 
      /* DriveCaller per Rho asintotico per variabili algebriche */     
      SAFENEWWITHCONSTRUCTOR(pRhoAlgebraicFictitious,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
      
      FictitiousType = INT_MS2;
  }
  
  
  /* costruzione dello step solver derivative */
  SAFENEWWITHCONSTRUCTOR(pDerivativeSteps,
			     DerivativeSolver,
			     DerivativeSolver(dDerivativesTol,
		  			      0.,
		  			      dInitialTimeStep*dDerivativesCoef,
		  			      iDerivativesMaxIterations));

     /* First step prediction must always be Crank-Nicholson for accuracy */
  SAFENEWWITHCONSTRUCTOR(pFirstFictitiousStep,
			     CrankNicholsonSolver,
			     CrankNicholsonSolver(dFictitiousStepsTolerance,
		  		0.,
			        iFictitiousStepsMaxIterations));
  
  SAFENEWWITHCONSTRUCTOR(pFirstRegularStep,
			 CrankNicholsonSolver,
			 CrankNicholsonSolver(dTol,
		  			dSolutionTol,
					iMaxIterations));
	
		  
  /* costruzione dello step solver fictitious */
  switch(FictitiousType) {
	  
	  case INT_CRANKNICHOLSON:
  		pFictitiousSteps = pFirstFictitiousStep;
		break;	
		    
	  case INT_MS2:
  		SAFENEWWITHCONSTRUCTOR(pFictitiousSteps,
			     		MultistepSolver,
			     		MultistepSolver(dFictitiousStepsTolerance,
						0.,
						iFictitiousStepsMaxIterations,
					        pRhoFictitious,
						pRhoAlgebraicFictitious));
		break;
		  
	  case INT_HOPE:
  		SAFENEWWITHCONSTRUCTOR(pFictitiousSteps,
			     		HopeSolver,
			     		HopeSolver(dFictitiousStepsTolerance,
		  				dSolutionTol,
						iFictitiousStepsMaxIterations,
					        pRhoFictitious,
						pRhoAlgebraicFictitious));
		break;
		  
	 default:
 	  	std::cerr << "Unknown dummy steps integration method" << std::endl;
	      	THROW(ErrGeneric());
		break;
		
  }	
  
  
  /* costruzione dello step solver per i passi normali */
  switch(RegularType) {
	  
	  case INT_CRANKNICHOLSON:
		pRegularSteps = pFirstRegularStep;
		break;	
		    
	  case INT_MS2:
  		SAFENEWWITHCONSTRUCTOR(pRegularSteps,
			     		MultistepSolver,
			     		MultistepSolver(dTol, 
		  				dSolutionTol,
						iMaxIterations,
					        pRhoRegular,
						pRhoAlgebraicRegular));
		break;
		  
	  case INT_HOPE:
  		SAFENEWWITHCONSTRUCTOR(pRegularSteps,
			     		HopeSolver,
			     		HopeSolver(dTol,
		  				dSolutionTol,
						iMaxIterations,
					        pRhoRegular,
						pRhoAlgebraicRegular));
		break;
		  
	 default:
 	  	std::cerr << "Unknown integration method" << std::endl;
	      	THROW(ErrGeneric());
		break;
		
  }	
}


/* Estrazione autovalori, vincolata alla disponibilita' delle LAPACK */
#ifdef __HACK_EIG__

#include <ac/lapack.h>

static int
do_eig(doublereal b, doublereal re, doublereal im, doublereal h, doublereal& sigma, doublereal& omega, doublereal& csi, doublereal& freq)
{
	int isPi = 0;

	if (fabs(b) > 1.e-16) {
		doublereal d;
		if (im != 0.) {
			d = sqrt(re*re+im*im)/fabs(b);
		} else {
			d = fabs(re/b);
		}
		sigma = log(d)/h;
		omega = atan2(im, re)/h;
		
		isPi = (fabs(im/b) < 1.e-15 && fabs(re/b+1.) < 1.e-15);
		if (isPi) {
			sigma = 0.;
			omega = HUGE_VAL;
		}
		
		d = sqrt(sigma*sigma+omega*omega);
		if (d > 1.e-15 && fabs(sigma)/d > 1.e-15) {
			csi = 100*sigma/d;
		} else {
			csi = 0.;
		}
		
		if (isPi) {
			freq = HUGE_VAL;
		} else {
			freq = omega/(2*M_PI);
		}
	} else {	 
		sigma = 0.;
		omega = 0.;
		csi = 0.;
		freq = 0.;
	}

	return isPi;
}

void 
Solver::Eig(void)
{
   DEBUGCOUTFNAME("Solver::Eig");   

   /*
    * MatA, MatB: MatrixHandlers to eigenanalysis matrices
    * MatL, MatR: MatrixHandlers to eigenvectors, if required
    * AlphaR, AlphaI Beta: eigenvalues
    * WorkVec:    Workspace
    * iWorkSize:  Size of the workspace
    */
   
   DEBUGCOUT("Solver::Eig(): performing eigenanalysis" << std::endl);
   
   char sL[2] = "V";
   char sR[2] = "V";
   
   /* iNumDof is a member, set after dataman constr. */
   integer iSize = iNumDofs;
   /* Minimum workspace size. To be improved */
   integer iWorkSize = 8*iNumDofs;
   integer iInfo = 0;

   /* Workspaces */
   /* 4 matrices iSize x iSize, 3 vectors iSize x 1, 1 vector iWorkSize x 1 */
   doublereal* pd = NULL;
   int iTmpSize = 4*(iSize*iSize)+3*iSize+iWorkSize;
   SAFENEWARR(pd, doublereal, iTmpSize);
   for (int iCnt = iTmpSize; iCnt-- > 0; ) {
      pd[iCnt] = 0.;
   }
      
   /* 4 pointer arrays iSize x 1 for the matrices */
   doublereal** ppd = NULL;
   SAFENEWARR(ppd, doublereal*, 4*iSize);

   /* Data Handlers */
   doublereal* pdTmp = pd;
   doublereal** ppdTmp = ppd;
  
#if 0
   doublereal* pdA = pd;
   doublereal* pdB = pdA+iSize*iSize;
   doublereal* pdVL = pdB+iSize*iSize;
   doublereal* pdVR = pdVL+iSize*iSize;
   doublereal* pdAlphaR = pdVR+iSize*iSize;
   doublereal* pdAlphaI = pdAlphaR+iSize;
   doublereal* pdBeta = pdAlphaI+iSize;
   doublereal* pdWork = pdBeta+iSize;
#endif /* 0 */

   FullMatrixHandler MatA(pdTmp, ppdTmp, iSize*iSize, iSize, iSize);
   MatA.Init();
   pdTmp += iSize*iSize;
   ppdTmp += iSize;
  
   FullMatrixHandler MatB(pdTmp, ppdTmp, iSize*iSize, iSize, iSize);
   MatB.Init();
   pdTmp += iSize*iSize;
   ppdTmp += iSize;
  
   FullMatrixHandler MatL(pdTmp, ppdTmp, iSize*iSize, iSize, iSize);  
   pdTmp += iSize*iSize;
   ppdTmp += iSize;
   
   FullMatrixHandler MatR(pdTmp, ppdTmp, iSize*iSize, iSize, iSize);    
   pdTmp += iSize*iSize;
   
   MyVectorHandler AlphaR(iSize, pdTmp);
   pdTmp += iSize;
   
   MyVectorHandler AlphaI(iSize, pdTmp);
   pdTmp += iSize;
   
   MyVectorHandler Beta(iSize, pdTmp);
   pdTmp += iSize;
   
   MyVectorHandler WorkVec(iWorkSize, pdTmp);
   
   /* Matrices Assembly (vedi eig.ps) */
   doublereal h = dEigParam;
   pDM->AssJac(MatA, -h/2.);
   pDM->AssJac(MatB, h/2.);

#ifdef DEBUG
   DEBUGCOUT(std::endl << "Matrix A:" << std::endl << MatA << std::endl
	     << "Matrix B:" << std::endl << MatB << std::endl);
#endif /* DEBUG */

   if (fEigenAnalysis == 2) {
      char *tmpFileName = NULL;
      const char *srcFileName = NULL;
      if (sOutputFileName == NULL) {
	 srcFileName = sInputFileName;
      } else {
	 srcFileName = sOutputFileName;
      }

      size_t l = strlen(srcFileName);
      SAFENEWARR(tmpFileName, char, l+1+3+1);
      strcpy(tmpFileName, srcFileName);
      strcpy(tmpFileName+l, ".mat");
      
      std::ofstream o(tmpFileName);

      o << MatA << std::endl << MatB << std::endl;

      o.close();
      SAFEDELETEARR(tmpFileName);
   }
   
#ifdef DEBUG_MEMMANAGER
   ASSERT(defaultMemoryManager.fIsValid(MatA.pdGetMat(), 
			   iSize*iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(MatB.pdGetMat(), 
			   iSize*iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(MatL.pdGetMat(), 
			   iSize*iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(MatR.pdGetMat(), 
			   iSize*iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(AlphaR.pdGetVec(), 
			   iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(AlphaI.pdGetVec(), 
			   iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(Beta.pdGetVec(), 
			   iSize*sizeof(doublereal)));
   ASSERT(defaultMemoryManager.fIsValid(WorkVec.pdGetVec(), 
			   iWorkSize*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
     
   
   /* Eigenanalysis */
   __FC_DECL__(dgegv)(sL,
	  sR,
	  &iSize,
	  MatA.pdGetMat(),
	  &iSize,
	  MatB.pdGetMat(),
	  &iSize,	  
	  AlphaR.pdGetVec(),
	  AlphaI.pdGetVec(),
	  Beta.pdGetVec(),
	  MatL.pdGetMat(), 
	  &iSize,
	  MatR.pdGetMat(), 
	  &iSize,
	  WorkVec.pdGetVec(),
	  &iWorkSize,
	  &iInfo);
   
   std::ostream& Out = pDM->GetOutFile();
   Out << "Info: " << iInfo << ", ";
   
   const char* const sErrs[] = {
      "DGGBAL",
	"DGEQRF",
	"DORMQR",
	"DORGQR",
	"DGGHRD",
	"DHGEQZ (other than failed iteration)",
	"DTGEVC",
	"DGGBAK (computing VL)",
	"DGGBAK (computing VR)",
	"DLASCL (various calls)"
   };
   
   if (iInfo == 0) {
      /* = 0:  successful exit */
      Out << "success" << std::endl;

   } else if (iInfo < 0) {
      char *th = "th";

      /* Aaaaah, English! :) */
      if (-iInfo/10 != 10) {
         switch ((-iInfo+20)%10) {
         case 1:
	    th = "st";
	    break;
         case 2:
	    th = "nd";
	    break;
	 case 3:
	    th = "rd";
	    break;
	 }
      }
      /* < 0:  if INFO = -i, the i-th argument had an illegal value. */
      Out << "the " << -iInfo << "-" << th 
	      << " argument had an illegal value" << std::endl;

   } else if (iInfo > 0 && iInfo <= iSize) {
      /* = 1,...,N:   
       * The QZ iteration failed.  No eigenvectors have been   
       * calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)   
       * should be correct for j=INFO+1,...,N. */
      Out << "the QZ iteration failed, but eigenvalues " 
	<< iInfo+1 << "->" << iSize << "should be correct" << std::endl;

   } else if (iInfo > iSize) {
      /* > N:  errors that usually indicate LAPACK problems:   
       * =N+1: error return from DGGBAL
       * =N+2: error return from DGEQRF
       * =N+3: error return from DORMQR
       * =N+4: error return from DORGQR   
       * =N+5: error return from DGGHRD   
       * =N+6: error return from DHGEQZ (other than failed iteration)   
       * =N+7: error return from DTGEVC   
       * =N+8: error return from DGGBAK (computing VL)   
       * =N+9: error return from DGGBAK (computing VR)   
       * =N+10: error return from DLASCL (various calls) */
      Out << "error return from " << sErrs[iInfo-iSize-1] << std::endl;	      
   }
   
   /* Output? */
   Out << "Mode n. " "  " "    Real    " "   " "    Imag    " "  " "    " "   Damp %   " "  Freq Hz" << std::endl;

   for (int iCnt = 1; iCnt <= iSize; iCnt++) {
      Out << std::setw(8) << iCnt << ": ";

      doublereal b = Beta.dGetCoef(iCnt);
      doublereal re = AlphaR.dGetCoef(iCnt);
      doublereal im = AlphaI.dGetCoef(iCnt);
      doublereal sigma;
      doublereal omega;
      doublereal csi;
      doublereal freq;

      int isPi = do_eig(b, re, im, h, sigma, omega, csi, freq);

      if (isPi) {
	      Out << std::setw(12) << 0. << " - " << "          PI j";
      } else {
	      Out << std::setw(12) << sigma << " + " << std::setw(12) << omega << " j";
      }

      if (fabs(csi) > 1.e-15) {
	      Out << "    " << std::setw(12) << csi;
      } else {
	      Out << "    " << std::setw(12) << 0.;
      }

      if (isPi) {
	      Out << "    " << "PI";
      } else {
	      Out << "    " << std::setw(12) << freq;
      }

      Out << std::endl;
   }

#ifdef __HACK_NASTRAN_MODES__
   /* EXPERIMENTAL */
   std::ofstream f06, pch;
   char datebuf[] = "11/14/95"; /* as in the example I used :) */

#if defined(HAVE_STRFTIME) && defined(HAVE_LOCALTIME) && defined(HAVE_TIME)
   time_t currtime = time(NULL);
   struct tm *currtm = localtime(&currtime);
   if (currtm) {
#warning "Your compiler might complain about %y"
#warning "in strftime() yielding only the last"
#warning "two digits of the year; don't worry,"
#warning "it's intended :)"
   	strftime(datebuf, sizeof(datebuf)-1, "%m/%d/%y", currtm);
   }
#endif /* HAVE_STRFTIME && HAVE_LOCALTIME && HAVE_TIME */
   
   if (fOutputModes) {
	   /* crea il file .pch */
	   char *tmp = NULL;

	   SAFENEWARR(tmp, char, strlen(sOutputFileName ? sOutputFileName : sInputFileName) + sizeof(".bdf"));
	   sprintf(tmp, "%s.bdf", sOutputFileName ? sOutputFileName : sInputFileName);
	   pch.open(tmp);
	   
	   pch.setf(std::ios::showpoint);
	   pch 
		   << "BEGIN BULK" << std::endl
		   << "$.......2.......3.......4.......5.......6.......7.......8.......9.......0......." << std::endl;
	   pch << "MAT1           1    1.E0    1.E0            1.E0" << std::endl;
	   pDM->Output_pch(pch);
	   
	   /* crea il file .f06 */
	   sprintf(tmp, "%s.f06", sOutputFileName ? sOutputFileName : sInputFileName);
	   f06.open(tmp);
	   SAFEDELETEARR(tmp);

	   f06.setf(std::ios::showpoint);
	   f06.setf(std::ios::scientific);
   }
#endif /* __HACK_NASTRAN_MODES__ */

   int iPage;
   for (iPage = 1; iPage <= iSize; iPage++) {

#ifdef __HACK_NASTRAN_MODES__
      /* EXPERIMENTAL */
      flag doPlot = 0;
      if (fOutputModes) {

         doublereal b = Beta.dGetCoef(iPage);
	 doublereal re = AlphaR.dGetCoef(iPage);
	 doublereal im = AlphaI.dGetCoef(iPage);
	 doublereal sigma;
	 doublereal omega;
	 doublereal csi;
	 doublereal freq;

	 do_eig(b, re, im, h, sigma, omega, csi, freq);

	 if (freq >= dLowerFreq && freq <= dUpperFreq) {
	    doPlot = 1;
	      
	    f06 
	      << "                                                                                                 CSA/NASTRAN " << datebuf << "    PAGE   " 
	      << std::setw(4) << iPage << std::endl
	      << "MBDyn modal analysis" << std::endl
	      << std::endl
	      << "    LABEL=DISPLACEMENTS, ";
#ifdef HAVE_FMTFLAGS_IN_IOS
	    std::ios::fmtflags iosfl = f06.setf(std::ios::left);
#else /* !HAVE_FMTFLAGS_IN_IOS */
	    long iosfl = f06.setf(std::ios::left);
#endif /* !HAVE_FMTFLAGS_IN_IOS */
	    const char *comment = "(EXPERIMENTAL) MODAL ANALYSIS";
	    int l = strlen(comment);
	    f06 << std::setw(l+1) << comment;
	    f06 << sigma << " " << (omega < 0. ? "-" : "+") << " " << fabs(omega) << " j (" << csi << ", "<< freq << std::setw(80-1-l) << ")";
	    f06.flags(iosfl);
	    f06 << "   SUBCASE " << iPage << std::endl
	      << std::endl
	      << "                                            D I S P L A C E M E N T  V E C T O R" << std::endl
	      << std::endl
	      << "     POINT ID.   TYPE          T1             T2             T3             R1             R2             R3" << std::endl;
	 }
      }
#endif /* __HACK_NASTRAN_MODES__ */
      
      doublereal cmplx = AlphaI.dGetCoef(iPage);
      if (cmplx == 0.) {
         Out << "Mode " << iPage << ":" << std::endl;
         for (int jCnt = 1; jCnt <= iSize; jCnt++) {
            Out << std::setw(12) << jCnt << ": "
	      << std::setw(12) << MatR.dGetCoef(jCnt, iPage) << std::endl;
	 }
	 
	 if (fOutputModes) {
	    /*
	     * per ora sono uguali; in realta' XP e' X * lambda
	     */
	    MyVectorHandler X(iSize, MatR.pdGetMat()+iSize*(iPage-1));
	    MyVectorHandler XP(iSize, MatR.pdGetMat()+iSize*(iPage-1));
	    pDM->Output(X, XP);
	    
#ifdef __HACK_NASTRAN_MODES__
	    /* EXPERIMENTAL */
	    if (doPlot) {
	       pDM->Output_f06(f06, X);
	    }
#endif /* __HACK_NASTRAN_MODES__ */
	 }
      } else {
	 if (cmplx > 0.) {
            Out << "Modes " << iPage << ", " << iPage+1 << ":" << std::endl;
	    for (int jCnt = 1; jCnt <= iSize; jCnt++) {
	       doublereal im = MatR.dGetCoef(jCnt, iPage+1);
	       Out << std::setw(12) << jCnt << ": "
	         << std::setw(12) << MatR.dGetCoef(jCnt, iPage) 
	         << ( im >= 0. ? " + " : " - " ) 
	         << std::setw(12) << fabs(im) << " * j " << std::endl;
            }
	 }

	 if (fOutputModes) {
            /*
	     * uso la parte immaginaria ...
	     */
	    int i = iPage - (cmplx > 0. ? 0 : 1);
	    MyVectorHandler X(iSize, MatR.pdGetMat()+iSize*i);
	    MyVectorHandler XP(iSize, MatR.pdGetMat()+iSize*i);
	    pDM->Output(X, XP);

#ifdef __HACK_NASTRAN_MODES__
	    /* EXPERIMENTAL */
	    if (doPlot) {
	       pDM->Output_f06(f06, X);
	    }
#endif /* __HACK_NASTRAN_MODES__ */
	 }
      }
   }

#ifdef __HACK_NASTRAN_MODES__
   if (fOutputModes) {
      pch 
	      << "ENDDATA" << std::endl;
      f06 
	      << "                                                                                                 CSA/NASTRAN " << datebuf << "    PAGE   " 
	      << std::setw(4) << iPage << std::endl;
      pch.close();
      f06.close();
   }
#endif /* __HACK_NASTRAN_MODES__ */      

   /* Non puo' arrivare qui se le due aree di lavoro non sono definite */
   SAFEDELETEARR(pd);
   SAFEDELETEARR(ppd);
}

#endif /* __HACK_EIG__ */
