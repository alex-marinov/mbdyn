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

#define RTAI_LOG

#include <unistd.h>
#include <ac/float.h>
#include <ac/math.h>

#include <solver.h>
#include "thirdorderstepsol.h"
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
const doublereal defaultIterativeTau = 1.e-7;

/* Costruttore: esegue la simulazione */
Solver::Solver(MBDynParser& HPar,
		const char* sInFName,
		const char* sOutFName,
		bool bPar)
:
CurrStrategy(NOCHANGE),
sInputFileName(NULL),
sOutputFileName(NULL),
HP(HPar),
pStrategyChangeDrive(NULL),
#ifdef __HACK_EIG__
eEigenAnalysis(EIG_NO),
dEigParam(1.),
bOutputModes(false),
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
#ifdef RTAI_LOG
bRTlog(false),
mbxlog(NULL),
LogProcName(NULL),
#endif /*RTAI_LOG*/
#endif /* USE_RTAI */
#ifdef __HACK_POD__
bPOD(0),
iPODStep(0),
iPODFrames(0),
#endif /*__HACK_POD__*/
iNumPreviousVectors(2),
iUnkStates(1),
pdWorkSpace(NULL),
qX(),
qXPrime(),
pX(NULL),
pXPrime(NULL),
dTime(0.),
dInitialTime(0.), 
dFinalTime(0.),
dRefTimeStep(0.),
dInitialTimeStep(1.), 
dMinimumTimeStep(1.),
dMaxTimeStep(1.),
iFictitiousStepsNumber(iDefaultFictitiousStepsNumber),
dFictitiousStepsRatio(dDefaultFictitiousStepsRatio),
eAbortAfter(AFTER_UNKNOWN),
iStepsAfterReduction(0),
iStepsAfterRaise(0),
iWeightedPerformedIters(0),
bLastChance(0),
pDerivativeSteps(NULL),
pFirstRegularStep(NULL),
pRegularSteps(NULL),
pFictitiousSteps(NULL),
ResTest(NonlinearSolverTest::NORM),
SolTest(NonlinearSolverTest::NONE),
bScale(false),
bTrueNewtonRaphson(1),
iIterationsBeforeAssembly(0),
NonlinearSolverType(NonlinearSolver::UNKNOWN),
MFSolverType(MatrixFreeSolver::UNKNOWN),
dIterTol(dDefaultTol),
PcType(Preconditioner::FULLJACOBIAN),
iPrecondSteps(iDefaultPreconditionerSteps),
iIterativeMaxSteps(iDefaultPreconditionerSteps),
dIterertiveEtaMax(defaultIterativeEtaMax),
dIterertiveTau(defaultIterativeTau),
/* FOR PARALLEL SOLVERS*/
bParallel(bPar),
pSDM(NULL),
iNumLocDofs(0),
iNumIntDofs(0),
pLocDofs(NULL),
pIntDofs(NULL),
pDofs(NULL),
pLocalSM(NULL),
/* end of FOR PARALLEL SOLVERS */
pDM(NULL),
iNumDofs(0),
pSM(NULL),
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

#ifdef USE_RTAI
	if (bRT) {
		/* Init RTAI; if init'ed, it will be shut down at exit */
		if (mbdyn_rt_task_init("MBDTSK", 1, 0, 0, &mbdyn_rtai_task)) {
			std::cerr << "unable to init RTAI task" << std::endl;
			THROW(ErrGeneric());
		}
	}
#endif /* USE_RTAI */

#ifdef USE_MPI
	int MyRank = 0;
	if (bParallel) {
		
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

		DEBUGLCOUT(MYDEBUG_MEM, "creating parallel SchurDataManager" 
				<< std::endl);
		
		SAFENEWWITHCONSTRUCTOR(pSDM,
			SchurDataManager,
			SchurDataManager(HP,
				OutputFlags,
				dInitialTime, 
				sInputFileName,
				sNewOutName,
				eAbortAfter == AFTER_INPUT));
		
		pDM = pSDM;

	} else {
#endif /* USE_MPI */

		/* chiama il gestore dei dati generali della simulazione */
		DEBUGLCOUT(MYDEBUG_MEM, "creating DataManager" << std::endl);
		SAFENEWWITHCONSTRUCTOR(pDM,
 				DataManager,
 				DataManager(HP, 
					OutputFlags,
					dInitialTime, 
					sInputFileName,
					sOutputFileName,
					eAbortAfter == AFTER_INPUT));
#ifdef USE_MPI
	}
#endif /* USE_MPI */

	HP.Close();
	
   	/* Si fa dare l'std::ostream al file di output per il log */
   	std::ostream& Out = pDM->GetOutFile();

   	if (eAbortAfter == AFTER_INPUT) {
      		/* Esce */
		pDM->Output(true);
      		Out << "End of Input; no simulation or assembly is required."
			<< std::endl;
      		return;

   	} else if (eAbortAfter == AFTER_ASSEMBLY) {
      		/* Fa l'output dell'assemblaggio iniziale e poi esce */
      		pDM->Output(true);
      		Out << "End of Initial Assembly; no simulation is required."
			<< std::endl;
      		return;
   	}

	/* Qui crea le partizioni: principale fra i processi, se parallelo  */
#ifdef USE_MPI
	if (bParallel) {
		pSDM->CreatePartition();
	}
#endif /* USE_MPI */

   	/* Si fa dare il DriveHandler e linka i drivers di rho ecc. */
   	const DriveHandler* pDH = pDM->pGetDrvHdl();
   	pRegularSteps->SetDriveHandler(pDH);
   	pFictitiousSteps->SetDriveHandler(pDH);
   
   	/* Costruisce i vettori della soluzione ai vari passi */
   	DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << std::endl);

	if (bParallel) {
		iNumDofs = pSDM->HowManyDofs(SchurDataManager::TOTAL);
		pDofs = pSDM->pGetDofsList();
		
		iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
		pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);
		iNumIntDofs = pSDM->HowManyDofs(SchurDataManager::INTERNAL);
		pIntDofs = pSDM->GetDofsList(SchurDataManager::INTERNAL);

	} else {
   		iNumDofs = pDM->iGetNumDofs();
	}
	
   	ASSERT(iNumDofs > 0);        
	
	integer iNSteps = pRegularSteps->GetIntegratorNumPreviousStates();
	integer iFSteps = pFictitiousSteps->GetIntegratorNumPreviousStates();
	iNumPreviousVectors = (iNSteps < iFSteps) ? iFSteps : iNSteps;
	
	iNSteps = pRegularSteps->GetIntegratorNumUnknownStates();
	iFSteps = pFictitiousSteps->GetIntegratorNumUnknownStates();
	iUnkStates = (iNSteps < iFSteps) ? iFSteps : iNSteps;
	
	/* allocate workspace for previous time steps */
	SAFENEWARR(
		pdWorkSpace,
		doublereal, 
		2*(iNumPreviousVectors)*iNumDofs);
	/* allocate MyVectorHandlers for previous time steps: use workspace */
	for (int ivec = 0; ivec < iNumPreviousVectors; ivec++) {  
   		SAFENEWWITHCONSTRUCTOR(pX,
			       	MyVectorHandler,
			       	MyVectorHandler(iNumDofs, pdWorkSpace+ivec*iNumDofs));
		qX.push_back(pX);
   		SAFENEWWITHCONSTRUCTOR(pXPrime,
			       	MyVectorHandler,
			       	MyVectorHandler(iNumDofs,
				pdWorkSpace+((iNumPreviousVectors)+ivec)*iNumDofs));
		qXPrime.push_back(pXPrime);
		pX = NULL;
		pXPrime = NULL;
	}
	/* allocate MyVectorHandlers for unknown time step(s): own memory */
	SAFENEWWITHCONSTRUCTOR(pX,
		       	MyVectorHandler,
		       	MyVectorHandler(iUnkStates*iNumDofs));
	SAFENEWWITHCONSTRUCTOR(pXPrime,
		       	MyVectorHandler,
		       	MyVectorHandler(iUnkStates*iNumDofs));
	

	/* Resetta i vettori */
   	for (int ivec = 0; ivec < iNumPreviousVectors; ivec++) {
		qX[ivec]->Reset(0.);
		qXPrime[ivec]->Reset(0.);
	}
	pX->Reset(0.);
	pXPrime->Reset(0.);

#ifdef __HACK_POD__
	std::ofstream PodOut;
	if (bPOD) {
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


   	/* Immediately link DataManager to current solution              */
	/*                                                               */
	/* this should work as long as the last unknown time step is put */
	/* at the beginning of pX, pXPrime                               */
   	pDM->LinkToSolution(*(pX), *(pXPrime));         

	/* a questo punto si costruisce il nonlinear solver */
	pNLS = AllocateNonlinearSolver();

	MyVectorHandler Scale(iNumDofs);
	VectorHandler *pScale = &Scale;
	if (bScale) {
		/* collects scale factors from data manager */
		pDM->SetScale(*pScale);
	}

	/*
	 * prepare tests for nonlinear solver;
	 *
	 * test on residual may allow pre-scaling;
	 * test on solution (difference between two iterations) does not
	 */
	NonlinearSolverTest *pResTest = NULL;
	if (bScale) {
		NonlinearSolverTestScale *pResTestScale = NULL;

		switch (ResTest) {
		case NonlinearSolverTest::NORM:
			SAFENEW(pResTestScale, NonlinearSolverTestScaleNorm);
			break;

		case NonlinearSolverTest::MINMAX:
			SAFENEW(pResTestScale, NonlinearSolverTestScaleMinMax);
			break;

		default:
			ASSERT(0);
			throw ErrGeneric();
		}

		/* registers scale factors at nonlinear solver */
		pResTestScale->SetScale(pScale);

		pResTest = pResTestScale;
		

	} else {
		switch (ResTest) {
		case NonlinearSolverTest::NONE:
			SAFENEW(pResTest, NonlinearSolverTestNone);
			break;

		case NonlinearSolverTest::NORM:
			SAFENEW(pResTest, NonlinearSolverTestNorm);
			break;

		case NonlinearSolverTest::MINMAX:
			SAFENEW(pResTest, NonlinearSolverTestMinMax);
			break;

		default:
			ASSERT(0);
			throw ErrGeneric();
		}
	}

	NonlinearSolverTest *pSolTest = NULL;
	switch (SolTest) {
	case NonlinearSolverTest::NONE:
		SAFENEW(pSolTest, NonlinearSolverTestNone);
		break;

	case NonlinearSolverTest::NORM:
		SAFENEW(pSolTest, NonlinearSolverTestNorm);
		break;

	case NonlinearSolverTest::MINMAX:
		SAFENEW(pSolTest, NonlinearSolverTestMinMax);
		break;

	default:
		ASSERT(0);
		throw ErrGeneric();
	}

	/* registers tests in nonlinear solver */
	pNLS->SetTest(pResTest, pSolTest);

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
	pDM->SetValue(*(pX), *(pXPrime));
	
#ifdef __HACK_EIG__  
   	if (eEigenAnalysis != EIG_NO && OneEig.dTime <= dTime && !OneEig.bDone) {
	 	Eig();
	 	OneEig.bDone = true;
   	}
#endif /* __HACK_EIG__ */
   
   	integer iTotIter = 0;
	integer iStIter = 0;
   	doublereal dTotErr = 0.;
	doublereal dTest = DBL_MAX;
   	doublereal dSolTest = DBL_MAX;
	bool bSolConv = false;
	/* calcolo delle derivate */
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
   	/* Setup SolutionManager(s) */
	SetupSolmans(pDerivativeSteps->GetIntegratorNumUnknownStates());
	/* Derivative steps */
	try {
		
		dTest = pDerivativeSteps->Advance(this,
				0., 1., StepIntegrator::NEWSTEP,
			 	qX, qXPrime, pX, pXPrime,
				iStIter, dTest, dSolTest);
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
		bSolConv = true;
	}
	catch (...) {
		THROW(ErrGeneric());
	}
	dTotErr  += dTest;
	iTotIter += iStIter;
	  	
	if (outputMsg()) {	
   		Out << "# Derivatives solution step at time " << dInitialTime
     			<< " performed in " << iStIter
     			<< " iterations with " << dTest
     			<< " error" << std::endl;
	}
      
   	DEBUGCOUT("Derivatives solution step has been performed successfully"
		  " in " << iStIter << " iterations" << std::endl);
   
#ifdef USE_EXTERNAL
	/* comunica che gli ultimi dati inviati sono la condizione iniziale */
	if (iFictitiousStepsNumber == 0) {
		External::SendInitial();
	}
#endif /* USE_EXTERNAL */

   	if (eAbortAfter == AFTER_DERIVATIVES) {
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
      		pDM->BeforePredict(*(pX), *(pXPrime),
				   *(qX[0]), *(qXPrime[0]));
      		Flip();

      		dRefTimeStep = dInitialTimeStep*dFictitiousStepsRatio;
      		dCurrTimeStep = dRefTimeStep;
      		pDM->SetTime(dTime+dCurrTimeStep);    
      	
      		DEBUGLCOUT(MYDEBUG_FSTEPS, "Current time step: "
			   << dCurrTimeStep << std::endl);
		
	 	ASSERT(pFirstFictitiousStep != NULL);
		
   		/* Setup SolutionManager(s) */
		SetupSolmans(pFirstFictitiousStep->GetIntegratorNumUnknownStates());
		/* pFirstFictitiousStep */
		try {
	 		dTest = pFirstFictitiousStep->Advance(this,
					dRefTimeStep, 1.,
					StepIntegrator::NEWSTEP,
					qX, qXPrime, pX, pXPrime, 
					iStIter, dTest, dSolTest);
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
			bSolConv = true;
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
		if (iFictitiousStepsNumber > 1) {
   			/* Setup SolutionManager(s) */
			SetupSolmans(pFictitiousSteps->GetIntegratorNumUnknownStates());
		}

      		for (int iSubStep = 2;
		     iSubStep <= iFictitiousStepsNumber;
		     iSubStep++) {
      			pDM->BeforePredict(*(pX), *(pXPrime),
				   	*(qX[0]), *(qXPrime[0]));
	 		Flip();
	 
	 		DEBUGLCOUT(MYDEBUG_FSTEPS, "Fictitious step "
				   << iSubStep 
				   << "; current time step: " << dCurrTimeStep
				   << std::endl);
	 
	 		ASSERT(pFictitiousSteps!= NULL);
			try {
	 			pDM->SetTime(dTime+dCurrTimeStep);
	 			dTest = pFictitiousSteps->Advance(this,
						dRefTimeStep,
						dCurrTimeStep/dRefTimeStep,
					 	StepIntegrator::NEWSTEP, 
						qX, qXPrime, pX, pXPrime,
						iStIter, dTest, dSolTest);
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
				bSolConv = true;
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
      			Out << "# Initial solution after dummy steps "
				"at time " << dTime
				<< " performed in " << iStIter
				<< " iterations with " << dTest 
				<< " error" << std::endl;
		}
			
      		DEBUGLCOUT(MYDEBUG_FSTEPS, 
			   "Fictitious steps have been completed successfully"
			   " in " << iStIter << " iterations" << std::endl);
#ifdef USE_EXTERNAL
	/* comunica che gli ultimi dati inviati sono la condizione iniziale */
		External::SendInitial();
#endif /* USE_EXTERNAL */
   	} /* Fine dei passi fittizi */

   	/* Output delle "condizioni iniziali" */
   	pDM->Output();

	   
        if (outputMsg()) {	
  	 	Out
			<< "# Key for lines starting with \"Step\":" 
				<< std::endl
			<< "# Step Time TStep NIter ResErr SolErr SolConv" 
				<< std::endl
			<< "Step " << 0
     			<< " " << dTime+dCurrTimeStep
     			<< " " << dCurrTimeStep
     			<< " " << iStIter
     			<< " " << dTest
			<< " " << dSolTest
			<< " " << bSolConv
			<< std::endl;
	}
   

   	if (eAbortAfter == AFTER_DUMMY_STEPS) {
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
	/* il prossimo passo e' un regular */
	pNLS->SetExternal(External::REGULAR);
#endif /* USE_EXTERNAL */
	
   	iStep = 1; /* Resetto di nuovo iStep */
      
   	DEBUGCOUT("Step " << iStep << " has been completed successfully in "
		  << iStIter << " iterations" << std::endl);

   			
   	DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   
      	pDM->BeforePredict(*(pX), *(pXPrime),
				*(qX[0]), *(qXPrime[0]));
	
	Flip();
	dRefTimeStep = dInitialTimeStep;   
   	dCurrTimeStep = dRefTimeStep;
	 	
	ASSERT(pFirstRegularStep!= NULL);
	StepIntegrator::StepChange CurrStep 
			= StepIntegrator::NEWSTEP;
   
	/* Setup SolutionManager(s) */
	SetupSolmans(pFirstRegularStep->GetIntegratorNumUnknownStates());
IfFirstStepIsToBeRepeated:
	try {   	
		pDM->SetTime(dTime+dCurrTimeStep);
		dTest = pFirstRegularStep->Advance(this, dRefTimeStep,
				dCurrTimeStep/dRefTimeStep, CurrStep, 
				qX, qXPrime, pX, pXPrime,
				iStIter, dTest, dSolTest);
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
		bSolConv = true;
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
      		Out 
			<< "Step " << iStep
			<< " " << dTime+dCurrTimeStep
			<< " " << dCurrTimeStep
			<< " " << iStIter
			<< " " << dTest
			<< " " << dSolTest
			<< " " << bSolConv
			<< std::endl;
	}

	bSolConv = false;

   	dRefTimeStep = dCurrTimeStep;
   	dTime += dRefTimeStep;
   
   	dTotErr += dTest;
   	iTotIter += iStIter;
   
#ifdef __HACK_EIG__  
   	if (eEigenAnalysis != EIG_NO && OneEig.dTime <= dTime && !OneEig.bDone) {
	 	Eig();
		OneEig.bDone = true;
      	}
#endif /* __HACK_EIG__ */

#ifdef USE_RTAI
	
#ifdef RTAI_LOG		    
	struct {
		int step;
		unsigned long time;
	} msg;
#endif /* RTAI_LOG */

	if (bRT) {
		/* Need timer */
		if (!mbdyn_rt_is_hard_timer_running() ){
			/* FIXME: ??? */
			std::cout << "Hard timer is started by MBDyn" 
				<< std::endl;
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
				<< "; period: " << mbdyn_count2nano(lRTPeriod) << std::endl);
			r = mbdyn_rt_task_make_periodic(mbdyn_rtai_task,
					t,lRTPeriod);

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

		/* FIXME: should check whether RTStackSize is correclty set? */
#ifdef RTAI_LOG	
		if (bRTlog) { 
			char *mbxlogname = "logmb";
			std::cout << "Mbdyn start overruns monitor" << std::endl;

			if (mbdyn_rt_mbx_init(mbxlogname, sizeof(msg)*8, &mbxlog)){
				bRTlog = false;
				std::cerr << "Cannot init mail box log" << std::endl;
			}
			switch (fork()) {
			case 0:
				if (execl(LogProcName, LogProcName, "MBDTSK",
						mbxlogname,NULL) == -1){
				/* error */
				std::cout << "Cannot start log procedure" 
						<< std::endl;
				}
				break;

			case -1:
				std::cerr << "Cannot init log procedure" << std::endl;
				bRTlog = false;

			default:
				mbdyn_rt_task_suspend(mbdyn_rtai_task);
			}
		}
#endif /* RTAI_LOG */

		reserve_stack(RTStackSize);
	}

	int 	RTStpFlag = 0;
	volatile int	RTSteps = 0;

        int t_tot = 0;
	long long t0 = 0, t1;
	int or_counter = 0;
#endif /* USE_RTAI */

    	/* Altri passi regolari */ 
	ASSERT(pRegularSteps != NULL);
	
	/* Setup SolutionManager(s) */
	SetupSolmans(pRegularSteps->GetIntegratorNumUnknownStates());
      	while (true) {
		
		StepIntegrator::StepChange CurrStep 
				= StepIntegrator::NEWSTEP;
	
      		if (dTime >= dFinalTime) {
#ifdef USE_RTAI
			if (bRT && bRTHard) {
				mbdyn_rt_make_soft_real_time();
			}
#endif /* USE_RTAI */
	 		std::cout << "End of simulation at time "
				<< dTime << " after " 
				<< iStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobians: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl;

#ifdef USE_RTAI
			if (bRT){
				std::cout << "total overruns: " << or_counter  << std::endl
					  << "total overrun time: " << t_tot << " micro s" << std::endl;
			}
#endif /* USE_RTAI */

			return;

#ifdef USE_RTAI
		} else if (bRT && RTStpFlag == 1){
			if (bRTHard) {
				mbdyn_rt_make_soft_real_time();
			}
			std::cout << "Simulation is stopped by RTAI task" << std::endl
				<< "Simulation ended at time "
				<< dTime << " after " 
				<< iStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobians: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl;

			std::cout << "Total overruns:" << or_counter  << std::endl
				<< "Total overruns time:" << t_tot << "micro s" << std::endl;

			return;
#endif /* USE_RTAI */

#ifdef HAVE_SIGNAL
      		} else if (!::mbdyn_keep_going) {
#ifdef USE_RTAI
			if (bRT && bRTHard) {
				mbdyn_rt_make_soft_real_time();
			}
#endif /* USE_RTAI */

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
      		pDM->BeforePredict(*(pX), *(pXPrime),
				*(qX[0]), *(qXPrime[0]));
	
		Flip();

#ifdef USE_RTAI
		if (bRT) {
			mbdyn_rt_receive_if(NULL, &RTStpFlag);

			t1 = mbdyn_rt_get_time();
			if ((RTSteps >= 2) && (t1 > (t0 + lRTPeriod))) {
				or_counter++;
				t_tot = t_tot + mbdyn_count2nano(t1 - t0 - lRTPeriod)/1000;

#ifdef RTAI_LOG
				if (bRTlog){
					msg.step = RTSteps;
					msg.time = mbdyn_count2nano(t1 - t0 - lRTPeriod)/1000;

					mbdyn_RT_mbx_send_if(0, 0, mbxlog, &msg, sizeof(msg));
				}
#endif /* RTAI_LOG */
			}


			if (RTWaitPeriod()) {
				mbdyn_rt_task_wait_period();

			} else if (RTSemaphore()) {
				/* FIXME: semaphore must be configurable */
				mbdyn_rt_sem_wait(RTSemPtr);
			}


			t0 = mbdyn_rt_get_time();

			if (bRTHard) {
				if (RTSteps == 2) {
					/* make hard real time */ 
					mbdyn_rt_make_hard_real_time();
				}
			}
			RTSteps++;
		}

#endif /* USE_RTAI */

IfStepIsToBeRepeated:
		try {  
 	
			pDM->SetTime(dTime+dCurrTimeStep);
			dTest = pRegularSteps->Advance(this, dRefTimeStep,
					dCurrTimeStep/dRefTimeStep, CurrStep,
					qX, qXPrime, pX, pXPrime, iStIter,
					dTest, dSolTest);
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
			bSolConv = true;
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
				<< " " << dTime+dCurrTimeStep
				<< " " << dCurrTimeStep
				<< " " << iStIter
				<< " " << dTest 
				<< " " << dSolTest
				<< " " << bSolConv
				<< std::endl;
		}

     	 	DEBUGCOUT("Step " << iStep
			<< " has been completed successfully in "
			<< iStIter << " iterations" << std::endl);
      
	      	dRefTimeStep = dCurrTimeStep;
      		dTime += dRefTimeStep;

		bSolConv = false;

#ifdef __HACK_POD__
		if (bPOD && dTime >= pod.dTime) {
			if (++iPODStep == pod.iSteps) {
				/* output degli stati su di una riga */
#ifdef __HACK_POD_BINARY__
	       			PodOut.write((char *)&pX, iNumDofs*sizeof(doublereal));
	       			PodOut.write((char *)&pXPrime, iNumDofs*sizeof(doublereal));
#else /* !__HACK_POD_BINARY__ */
				PodOut << pX->dGetCoef(1);
				for (integer j = 1; j < iNumDofs; j++) {
					PodOut << "  " << pX->dGetCoef(j+1);
                       		}
                       		PodOut << std::endl;
#endif /* ! __HACK_POD_BINARY__ */
			}
                     	iPODFrames++;
                      	iPODStep = 0;
		}
	
		if (iPODFrames >= pod.iFrames){
			bPOD = false;
		}                        
#endif /*__HACK_POD__ */

#ifdef __HACK_EIG__
      		if (eEigenAnalysis != EIG_NO && OneEig.dTime <= dTime && !OneEig.bDone) {
			Eig();
			OneEig.bDone = true;
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
			SAFEDELETE(qXPrime[ivec]);
		}
	}
	SAFEDELETE(pX);
	SAFEDELETE(pXPrime);

   	if (pdWorkSpace != NULL) {	
      		SAFEDELETEARR(pdWorkSpace);
   	}

   	if (pDM != NULL) {	
      		SAFEDELETE(pDM);
	}
#if defined(USE_RTAI) && defined(RTAI_LOG)
	if (bRTlog&&bRT){
		mbdyn_rt_mbx_delete(&mbxlog);
	}
#endif /* USE_RTAI && RTAI_LOG */
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
	     			if (bLastChance == true) {
					bLastChance = false;
	     			}
	     			iStepsAfterReduction = 0;
	     			return dCurrTimeStep*StrategyFactor.dReductionFactor;
	  		} else {
	     			if (bLastChance == false) {
					bLastChance = true;
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

			iWeightedPerformedIters = (10*iPerformedIters + 9*iWeightedPerformedIters)/10;
	  
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
const doublereal dDefaultFictitiousStepsTolerance = dDefaultTol;


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
		"modify" "residual" "test",
	
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
			"thirdorder",
	
		"derivatives" "coefficient",
		"derivatives" "tolerance",
		"derivatives" "max" "iterations",

		/* DEPRECATED */	
		"true",
		"modified",
		/* END OF DEPRECATED */
		
		"strategy",
			"factor",
			"no" "change",
			"change",
		
		"pod",
		"eigen" "analysis",
		"output" "modes",
		
		"solver",
		"interface" "solver", 
		
		/* DEPRECATED */	
		"preconditioner",
		/* END OF DEPRECATED */

		"nonlinear" "solver",
			"default",
			"newton" "raphson",
			"matrix" "free",
				"bicgstab",
				"gmres",
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
		MODIFY_RES_TEST,
	
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
		THIRDORDER,
	
		DERIVATIVESCOEFFICIENT,
		DERIVATIVESTOLERANCE,
		DERIVATIVESMAXITERATIONS,

		/* DEPRECATED */
		NR_TRUE,
		MODIFIED,
		/* END OF DEPRECATED */

		STRATEGY,
		STRATEGYFACTOR,
		STRATEGYNOCHANGE,
		STRATEGYCHANGE,
	
		POD,
		EIGENANALYSIS,
		OUTPUTMODES,
		
		SOLVER,
		INTERFACESOLVER,
		
		/* DEPRECATED */
		PRECONDITIONER,
		/* END OF DEPRECATED */
	
		NONLINEARSOLVER,
			DEFAULT,
			NEWTONRAPHSON,
			MATRIXFREE,
				BICGSTAB,
				GMRES,
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

   	bool bMethod(false);
   	bool bFictitiousStepsMethod(false);      
	
	/* dati letti qui ma da passare alle classi
	 *	StepIntegration e NonlinearSolver
	 */ 
	
	doublereal dTol = dDefaultTol;
   	doublereal dSolutionTol = 0.;
   	integer iMaxIterations = iDefaultMaxIterations;
	bool bModResTest = false;

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
			INT_THIRDORDER,
			INT_UNKNOWN
	};
	
	StepIntegratorType RegularType = INT_UNKNOWN, FictitiousType = INT_UNKNOWN; 
	
	
	 	
   	/* Ciclo infinito */
   	while (true) {	
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
	      			eAbortAfter = AFTER_INPUT;
	      			DEBUGLCOUT(MYDEBUG_INPUT, 
			 		"Simulation will abort after"
					" data input" << std::endl);
	      			break;
			
	   		case ASSEMBLY:
	     			eAbortAfter = AFTER_ASSEMBLY;
	      			DEBUGLCOUT(MYDEBUG_INPUT,
			 		   "Simulation will abort after"
					   " initial assembly" << std::endl);
	      			break;	  
	     
	   		case DERIVATIVES:
	      			eAbortAfter = AFTER_DERIVATIVES;
	      			DEBUGLCOUT(MYDEBUG_INPUT, 
			 		   "Simulation will abort after"
					   " derivatives solution" << std::endl);
	      			break;
	     
	   		case FICTITIOUSSTEPS:
	   		case DUMMYSTEPS:
	      			eAbortAfter = AFTER_DUMMY_STEPS;
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
			unsigned OF = OUTPUT_DEFAULT;

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
	  		if (bMethod) {
	     			std::cerr << "error: multiple definition"
					" of integration method at line "
					<< HP.GetLineData() << std::endl;
	     			THROW(ErrGeneric());
	  		}
	  		bMethod = true;
	        	  
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
			case THIRDORDER: {
	      			pRhoRegular =
					ReadDriveData(NULL, HP, NULL);
				HP.PutKeyTable(K);
				RegularType = INT_THIRDORDER;
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
			if (bFictitiousStepsMethod) {
				std::cerr << "error: multiple definition "
					"of dummy steps integration method "
					"at line " << HP.GetLineData()
					<< std::cerr;
				THROW(ErrGeneric());
			}
			bFictitiousStepsMethod = true;	  	
 
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
	   		case THIRDORDER:
				pRhoFictitious = ReadDriveData(NULL, HP, NULL);
				HP.PutKeyTable(K);
				FictitiousType = INT_THIRDORDER;
				break;
			default: {
				std::cerr << "Unknown integration method at line " << HP.GetLineData() << std::endl;
				THROW(ErrGeneric());
			}
			}	     
			break;
		} 

		case TOLERANCE: {
			/*
			 * residual tolerance; can be the keyword "null",
			 * which means that the convergence test 
			 * will be computed on the solution, or a number
			 */
			if (HP.IsKeyWord("null")) {
				dTol = 0.;

			} else {
				dTol = HP.GetReal();
				if (dTol < 0.) {
					dTol = dDefaultTol;
					std::cerr << "warning, residual tolerance "
						"< 0. is illegal; "
						"using default value " << dTol
						<< std::endl;
				}
			}

			/* safe default */
			if (dTol == 0.) {
				ResTest = NonlinearSolverTest::NONE;
			}
				
			if (HP.fIsArg()) {
				if (HP.IsKeyWord("test")) {
					if (HP.IsKeyWord("norm")) {
						ResTest = NonlinearSolverTest::NORM;
					} else if (HP.IsKeyWord("minmax")) {
						ResTest = NonlinearSolverTest::MINMAX;
					} else if (HP.IsKeyWord("none")) {
						ResTest = NonlinearSolverTest::NONE;
					} else {
						std::cerr << "unknown test "
							"method at line " 
							<< HP.GetLineData()
							<< std::endl;
						throw ErrGeneric();
					}

					if (HP.IsKeyWord("scale")) {
						if (ResTest == NonlinearSolverTest::NONE) {
							std::cerr << "it's a nonsense "
								"to scale a disabled test; "
								"\"scale\" ignored" 
								<< std::endl;
							bScale = false;
						} else {
							bScale = true;
						}
					}
				}
			}

			if (HP.fIsArg()) {
				if (!HP.IsKeyWord("null")) {
					dSolutionTol = HP.GetReal();
				}

				/* safe default */
				if (dSolutionTol != 0.) {
					SolTest = NonlinearSolverTest::NORM;
				}

				if (HP.fIsArg()) {
					if (HP.IsKeyWord("test")) {
						if (HP.IsKeyWord("norm")) {
							SolTest = NonlinearSolverTest::NORM;
						} else if (HP.IsKeyWord("minmax")) {
							SolTest = NonlinearSolverTest::MINMAX;
						} else if (HP.IsKeyWord("none")) {
							SolTest = NonlinearSolverTest::NONE;
						} else {
							std::cerr << "unknown test "
								"method at line " 
								<< HP.GetLineData() 
								<< std::endl;
							throw ErrGeneric();
						}
					}
				}
				
			} else if (dTol == 0.) {
				std::cerr << "need solution tolerance "
					"with null residual tolerance"
					<< std::endl;
				throw ErrGeneric();
			}

			if (dSolutionTol < 0.) {
				dSolutionTol = 0.;
				std::cerr << "warning, solution tolerance "
					"< 0. is illegal; "
					"solution test is disabled" 
					<< std::endl;
			}

			if (dTol == 0. && dSolutionTol == 0.) {
				std::cerr << "both residual and solution "
					"tolerances are zero" << std::endl;
				throw ErrGeneric();
			}

			DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dTol
					<< ", " << dSolutionTol << std::endl);
			break;
		}	

	 
		case DERIVATIVESTOLERANCE: {
			dDerivativesTol = HP.GetReal();
			if (dDerivativesTol <= 0.) {
				dDerivativesTol = dDefaultTol;
				std::cerr << "warning, derivatives "
					"tolerance <= 0.0 is illegal; "
					"using default value " 
					<< dDerivativesTol
					<< std::endl;
			}
			DEBUGLCOUT(MYDEBUG_INPUT, 
					"Derivatives tolerance = " 
					<< dDerivativesTol
					<< std::endl);
			break;
		}
	 
		case MAXITERATIONS: {
			iMaxIterations = HP.GetInt();
			if (iMaxIterations < 1) {
				iMaxIterations = iDefaultMaxIterations;
				std::cerr << "warning, max iterations "
					"< 1 is illegal; using default value "
					<< iMaxIterations
					<< std::endl;
			}
			DEBUGLCOUT(MYDEBUG_INPUT, 
					"Max iterations = "
					<< iMaxIterations << std::endl);
			break;
		}

		case MODIFY_RES_TEST: {
			bModResTest = true;
			DEBUGLCOUT(MYDEBUG_INPUT, 
					"Modify residual test" << std::endl);
			break;
		}

		case DERIVATIVESMAXITERATIONS: {
			iDerivativesMaxIterations = HP.GetInt();
			if (iDerivativesMaxIterations < 1) {
				iDerivativesMaxIterations = iDefaultMaxIterations;
				std::cerr << "warning, derivatives "
					"max iterations < 1 is illegal; "
					"using default value "
					<< iDerivativesMaxIterations
					<< std::endl;
			}
			DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives "
					"max iterations = "
					<< iDerivativesMaxIterations
					<< std::endl);
			break;
		}

		case FICTITIOUSSTEPSMAXITERATIONS:
		case DUMMYSTEPSMAXITERATIONS: {
			iFictitiousStepsMaxIterations = HP.GetInt();
			if (iFictitiousStepsMaxIterations < 1) {
				iFictitiousStepsMaxIterations = iDefaultMaxIterations;
				std::cerr << "warning, dummy steps "
					"max iterations < 1 is illegal; "
					"using default value "
					<< iFictitiousStepsMaxIterations
					<< std::endl;
			}
			DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps "
					"max iterations = " 
					<< iFictitiousStepsMaxIterations
					<< std::endl);
			break;
		}

		case DERIVATIVESCOEFFICIENT: {
			dDerivativesCoef = HP.GetReal();
			if (dDerivativesCoef <= 0.) {
				dDerivativesCoef = 1.;
				std::cerr << "warning, derivatives "
					"coefficient <= 0. is illegal; "
					"using default value "
					<< dDerivativesCoef
					<< std::endl;
			}
			DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives coefficient = "
					<< dDerivativesCoef << std::endl);
			break;
		}

		case NEWTONRAPHSON: {
			pedantic_cout("Newton Raphson is deprecated; use "
					"\"nonlinear solver: newton raphson "
					"[ , modified, <steps> ]\" instead"
					<< std::endl);
			KeyWords NewRaph(KeyWords(HP.GetWord()));
			switch(NewRaph) {
			case MODIFIED:
				bTrueNewtonRaphson = 0;
				if (HP.fIsArg()) {
					iIterationsBeforeAssembly = HP.GetInt();
		  		} else {
		       			iIterationsBeforeAssembly = iDefaultIterationsBeforeAssembly;
				}
				DEBUGLCOUT(MYDEBUG_INPUT, "Modified "
						"Newton-Raphson will be used; "
						"matrix will be assembled "
						"at most after "
						<< iIterationsBeforeAssembly
						<< " iterations" << std::endl);
				break;

			default:
				std::cerr << "warning: unknown case; "
					"using default" << std::endl;
			
			/* no break: fall-thru to next case */
			case NR_TRUE:
				bTrueNewtonRaphson = 1;
				iIterationsBeforeAssembly = 0;
				break;
			}
			break;
		}

		case END:
			if (KeyWords(HP.GetWord()) != MULTISTEP) {
				std::cerr << "<end: multistep;> expected "
					"at line " << HP.GetLineData() 
					<< "; aborting ..." << std::endl;
				THROW(ErrGeneric());
			}
			goto EndOfCycle;

		case STRATEGY: {
			switch (KeyWords(HP.GetWord())) {
			case STRATEGYFACTOR: {
				CurrStrategy = FACTOR;

				/*
				 * strategy: factor ,
				 *     <reduction factor> ,
				 *     <steps before reduction> ,
				 *     <raise factor> ,
				 *     <steps before raise> ,
				 *     <min iterations> ;
				 */

				StrategyFactor.dReductionFactor = HP.GetReal();
				if (StrategyFactor.dReductionFactor >= 1.) {
					std::cerr << "warning, "
						"illegal reduction factor "
						"at line " << HP.GetLineData() 
						<< "; default value 1. "
						"(no reduction) will be used"
						<< std::endl;
					StrategyFactor.dReductionFactor = 1.;
				}

				StrategyFactor.iStepsBeforeReduction = HP.GetInt();
				if (StrategyFactor.iStepsBeforeReduction <= 0) {
					std::cerr << "warning, "
						"illegal number of steps "
						"before reduction at line "
						<< HP.GetLineData()
						<< "; default value 1 will be "
						"used (it may be dangerous)" 
						<< std::endl;
					StrategyFactor.iStepsBeforeReduction = 1;
				}

				StrategyFactor.dRaiseFactor = HP.GetReal();
				if (StrategyFactor.dRaiseFactor <= 1.) {
					std::cerr << "warning, "
						"illegal raise factor at line "
						<< HP.GetLineData() 
						<< "; default value 1. "
						"(no raise) will be used"
						<< std::endl;
					StrategyFactor.dRaiseFactor = 1.;
				}

				StrategyFactor.iStepsBeforeRaise = HP.GetInt();
				if (StrategyFactor.iStepsBeforeRaise <= 0) {
					std::cerr << "warning, "
						"illegal number of steps "
						"before raise at line "
						<< HP.GetLineData()
						<< "; default value 1 will be "
						"used (it may be dangerous)" 
						<< std::endl;
					StrategyFactor.iStepsBeforeRaise = 1;
				}

				StrategyFactor.iMinIters = HP.GetInt();
				if (StrategyFactor.iMinIters <= 0) {
					std::cerr << "warning, "
						"illegal minimum number "
						"of iterations at line "
						<< HP.GetLineData()
						<< "; default value 0 will be "
						"used (never raise)" 
						<< std::endl;
					StrategyFactor.iMinIters = 1;
				}

				DEBUGLCOUT(MYDEBUG_INPUT,
						"Time step control strategy: "
						"Factor" << std::endl
						<< "Reduction factor: " 
						<< StrategyFactor.dReductionFactor
						<< "Steps before reduction: "
						<< StrategyFactor.iStepsBeforeReduction
						<< "Raise factor: "
						<< StrategyFactor.dRaiseFactor
						<< "Steps before raise: "
						<< StrategyFactor.iStepsBeforeRaise 
						<< "Min iterations: "
						<< StrategyFactor.iMinIters
						<< std::endl);
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

			default:
				std::cerr << "unknown time step control "
					"strategy at line "
					<< HP.GetLineData() << std::endl;
				THROW(ErrGeneric());
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

			bPOD = true;
			DEBUGLCOUT(MYDEBUG_INPUT, "POD analysis will be "
					"performed since time " << pod.dTime
					<< " for " << pod.iFrames 
					<< " frames  every " << pod.iSteps 
					<< " steps" << std::endl);
#else/* !__HACK_POD__ */
			std::cerr << "line " << HP.GetLineData()
				<< ": POD analysis not supported (ignored)" 
				<< std::endl;
			for (; HP.fIsArg();) {
				(void)HP.GetReal();
			}
#endif /* !__HACK_POD__ */
			break;
		
		case EIGENANALYSIS:
#ifdef __HACK_EIG__
			OneEig.dTime = HP.GetReal();
			if (HP.IsKeyWord("parameter")) {
				dEigParam = HP.GetReal();
			}
			OneEig.bDone = false;
			eEigenAnalysis = EIG_YES;
			DEBUGLCOUT(MYDEBUG_INPUT, "Eigenanalysis will be "
					"performed at time " << OneEig.dTime
					<< " (parameter: " << dEigParam << ")" 
					<< std::endl);
			if (HP.IsKeyWord("output" "matrices")) {
				eEigenAnalysis = EIG_OUTPUTMATRICES;
			}
#else /* !__HACK_EIG__ */
			HP.GetReal();
			if (HP.IsKeyWord("parameter")) {
				HP.GetReal();
			}
			
			if (HP.IsKeyWord("output" "matrices")) {
				NO_OP;
			}
			
			std::cerr << HP.GetLineData() 
				<< ": eigenanalysis not supported (ignored)"
				<< std::endl;
#endif /* !__HACK_EIG__ */
			break;

		case OUTPUTMODES:
#ifndef __HACK_EIG__
			std::cerr << "line " << HP.GetLineData()
				<< ": warning, no eigenvalue support available"
				<< std::endl;
#endif /* !__HACK_EIG__ */
			if (HP.IsKeyWord("yes") || HP.IsKeyWord("nastran")) {
#ifdef __HACK_EIG__
				bOutputModes = true;
				if (HP.fIsArg()) {
					dUpperFreq = HP.GetReal();
					if (dUpperFreq < 0.) {
						std::cerr << "line "
							<< HP.GetLineData()
							<< ": illegal upper "
							"frequency limit " 
							<< dUpperFreq 
							<< "; using " 
							<< -dUpperFreq
							<< std::endl;
						dUpperFreq = -dUpperFreq;
					}
					
					if (HP.fIsArg()) {
						dLowerFreq = HP.GetReal();
						if (dLowerFreq > dUpperFreq) {
							std::cerr << "line "
								<< HP.GetLineData()
								<< ": illegal lower "
								"frequency limit "
								<< dLowerFreq 
								<< " higher than upper "
								"frequency limit; using " 
								<< 0. << std::endl;
							dLowerFreq = 0.;
						}
					}
				}
#endif /* !__HACK_EIG__ */
			} else if (HP.IsKeyWord("no")) {
#ifdef __HACK_EIG__
				bOutputModes = false;
#endif /* !__HACK_EIG__ */
			} else {
				std::cerr << "line " << HP.GetLineData()
					<< ": unknown mode output flag "
					"(should be { yes | no })"
					<< std::endl;
			}
			break;

		case SOLVER:
			CurrSolver.Read(HP);
			
			/* ripristina la tabella del parser */
			HP.PutKeyTable(K);

			break;

		case INTERFACESOLVER:
			CurrIntSolver.Read(HP, true);
			
			/* ripristina la tabella del parser */
			HP.PutKeyTable(K);

#ifndef USE_MPI
			std::cerr << "Interface solver only allowed "
				"when compiled with MPI support" << std::endl;
#endif /* ! USE_MPI */
			break;

		case NONLINEARSOLVER:
			switch (KeyWords(HP.GetWord())) {
			case DEFAULT:
				NonlinearSolverType = NonlinearSolver::DEFAULT;
				break;

			case NEWTONRAPHSON:
				NonlinearSolverType = NonlinearSolver::NEWTONRAPHSON;
				break;

			case MATRIXFREE:
				NonlinearSolverType = NonlinearSolver::MATRIXFREE;
				break;

			default:
				std::cerr << "unknown nonlinear solver "
					"at line " << HP.GetLineData()
					<< std::endl;
				THROW(ErrGeneric());
				break;
			}

			switch (NonlinearSolverType) {
			case NonlinearSolver::NEWTONRAPHSON:
    				bTrueNewtonRaphson = true;
				bKeepJac = false;
    				iIterationsBeforeAssembly = 0;

				if (HP.IsKeyWord("modified")) {
					bTrueNewtonRaphson = false;
					iIterationsBeforeAssembly = HP.GetInt();

					if (HP.IsKeyWord("keep" "jacobian")) {
						bKeepJac = true;
					}

					DEBUGLCOUT(MYDEBUG_INPUT, "modified "
							"Newton-Raphson "
							"will be used; "
							"matrix will be "
							"assembled at most "
							"after "
							<< iIterationsBeforeAssembly
							<< " iterations"
							<< std::endl);
				}
				break;

			case NonlinearSolver::MATRIXFREE:
				switch (KeyWords(HP.GetWord())) {
				case DEFAULT:
					MFSolverType = MatrixFreeSolver::DEFAULT;
					break;

 
				case BICGSTAB:
					MFSolverType = MatrixFreeSolver::BICGSTAB;
					break;

				case GMRES:
					MFSolverType = MatrixFreeSolver::GMRES;
					break;

				default:
					std::cerr << "unknown iterative "
						"solver at line " 
						<< HP.GetLineData()
						<< std::endl;
					THROW(ErrGeneric());
				}

				if (HP.IsKeyWord("tolerance")) {
					dIterTol = HP.GetReal();
					DEBUGLCOUT(MYDEBUG_INPUT,"inner "
							"iterative solver "
							"tolerance: "
							<< dIterTol
							<< std::endl);
				}
				
				if (HP.IsKeyWord("steps")) {
					iIterativeMaxSteps = HP.GetInt();
					DEBUGLCOUT(MYDEBUG_INPUT, "maximum "
							"number of inner "
							"steps for iterative "
							"solver: " 
							<< iIterativeMaxSteps
							<< std::endl);
				}
				
				if (HP.IsKeyWord("tau")) {
					dIterertiveTau = HP.GetReal();
					DEBUGLCOUT(MYDEBUG_INPUT,
							"tau scaling "
							"coefficient "
							"for iterative "
							"solver: "
							<< dIterertiveTau
							<< std::endl);
				}

				if (HP.IsKeyWord("eta")) {
					dIterertiveEtaMax = HP.GetReal();
					DEBUGLCOUT(MYDEBUG_INPUT, "maximum "
							"eta coefficient "
							"for iterative "
							"solver: " 
	  						<< dIterertiveEtaMax
							<< std::endl);
				}

				if (HP.IsKeyWord("preconditioner")) {
					KeyWords KPrecond = KeyWords(HP.GetWord());
					switch (KPrecond) {
					case FULLJACOBIAN:
						PcType = Preconditioner::FULLJACOBIAN;		
						if (HP.IsKeyWord("steps")) {
							iPrecondSteps = HP.GetInt();
							DEBUGLCOUT(MYDEBUG_INPUT,
									"number of steps "
									"before recomputing "
									"the preconditioner: "
									<< iPrecondSteps
									<< std::endl);
						}
						break;

						/* add other preconditioners
						 * here */

					default:
						std::cerr << "unknown "
							"preconditioner "
							"at line "
							<< HP.GetLineData()
							<< std::endl;
						THROW(ErrGeneric());
					}
					break;
				}
				break;

			default:
				ASSERT(0);
				THROW(ErrGeneric());
			}
			break;

		case REALTIME:
#ifdef USE_RTAI
			bRT = true;

			if (HP.IsKeyWord("allow" "nonroot")) {
				bRTAllowNonRoot = true;
			}

			/* FIXME: use a safe default? */
			if (HP.IsKeyWord("mode")) {
				if (HP.IsKeyWord("period")) {
					RTMode = MBRTAI_WAITPERIOD;
					
				} else if (HP.IsKeyWord("semaphore")) {
					/* FIXME: not implemented yet ... */
					RTMode = MBRTAI_SEMAPHORE;
				} else {
					std::cerr << "unknown realtime mode "
						"at line " << HP.GetLineData()
						<< std::endl;
					THROW(ErrGeneric());
				}
			}
			
			if (HP.IsKeyWord("time" "step")) {
				long long p = HP.GetInt();

				if (p <= 0) {
					std::cerr << "illegal time step " 
						<< p << " at line " 
						<< HP.GetLineData()
						<< std::endl;
					THROW(ErrGeneric());
				}

				lRTPeriod = mbdyn_nano2count(p);
				
			} else {
				std::cerr << "need a time step for real time "
					"at line " << HP.GetLineData()
					<< std::endl;
				THROW(ErrGeneric());
			}

			if (HP.IsKeyWord("reserve" "stack")) {
				long size = HP.GetInt();

				if (size <= 0) {
					std::cerr << "illegal stack size "
						<< size << " at line "
						<< HP.GetLineData()
						<< std::endl;
					THROW(ErrGeneric());
				}

				RTStackSize = size;
			}
			
			if (HP.IsKeyWord("hard" "real" "time")) {
				bRTHard = true;
			}
#ifdef RTAI_LOG
			if (HP.IsKeyWord("real" "time" "log")) {
				if (HP.IsKeyWord("file" "name")){
					const char *m = HP.GetFileName();
					SAFESTRDUP(LogProcName, m);
				} else {
					/* FIXME */
					SAFESTRDUP(LogProcName, "./logproc");
				}
				bRTlog = true;
			}
#endif /*RTAI_LOG*/

#else /* !USE_RTAI */
			std::cerr << "need to configure --with-rtai "
				"to use realtime" << std::endl;
			THROW(ErrGeneric());
#endif /* !USE_RTAI */
			break;

		default:
			std::cerr << "unknown description at line " 
				<< HP.GetLineData() << "; aborting ..."
				<< std::endl;
			THROW(ErrGeneric());
		}
	}
   
EndOfCycle: /* esce dal ciclo di lettura */

	if (dFinalTime < dInitialTime) {      
		eAbortAfter = AFTER_ASSEMBLY;
	}

	if (dFinalTime == dInitialTime) {      
		eAbortAfter = AFTER_DERIVATIVES;
	}

	/* Metodo di integrazione di default */
	if (!bMethod || pRhoRegular == NULL) {
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
	if (!bFictitiousStepsMethod || pRhoFictitious == NULL) {
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
				iDerivativesMaxIterations,
				bModResTest));

	/* First step prediction must always be Crank-Nicholson for accuracy */
	SAFENEWWITHCONSTRUCTOR(pFirstFictitiousStep,
			CrankNicholsonSolver,
			CrankNicholsonSolver(dFictitiousStepsTolerance,
				0.,
				iFictitiousStepsMaxIterations,
				bModResTest));

	SAFENEWWITHCONSTRUCTOR(pFirstRegularStep,
			CrankNicholsonSolver,
			CrankNicholsonSolver(dTol,
				dSolutionTol,
				iMaxIterations,
				bModResTest));

	/* costruzione dello step solver fictitious */
	switch (FictitiousType) {
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
					pRhoAlgebraicFictitious,
					bModResTest));
		break;
		  
	case INT_HOPE:
		SAFENEWWITHCONSTRUCTOR(pFictitiousSteps,
				HopeSolver,
				HopeSolver(dFictitiousStepsTolerance,
					dSolutionTol,
					iFictitiousStepsMaxIterations,
					pRhoFictitious,
					pRhoAlgebraicFictitious,
					bModResTest));
		break;
		  
	case INT_THIRDORDER:
		SAFENEWWITHCONSTRUCTOR(pFictitiousSteps,
				ThirdOrderIntegrator,
				ThirdOrderIntegrator(dFictitiousStepsTolerance,
					dSolutionTol,
					iFictitiousStepsMaxIterations,
					pRhoFictitious,
					bModResTest));
		break;
		
	default:
		std::cerr << "unknown dummy steps integration method"
			<< std::endl;
		THROW(ErrGeneric());
		break;
	}

	/* costruzione dello step solver per i passi normali */
	switch (RegularType) {
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
					pRhoAlgebraicRegular,
					bModResTest));
		break;
		  
	case INT_HOPE:
		SAFENEWWITHCONSTRUCTOR(pRegularSteps,
				HopeSolver,
				HopeSolver(dTol,
					dSolutionTol,
					iMaxIterations,
					pRhoRegular,
					pRhoAlgebraicRegular,
					bModResTest));
		break;
		  
	case INT_THIRDORDER:
		SAFENEWWITHCONSTRUCTOR(pRegularSteps,
				ThirdOrderIntegrator,
				ThirdOrderIntegrator(dTol,
					dSolutionTol,
					iMaxIterations,
					pRhoRegular,
					bModResTest));
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
do_eig(const doublereal& b, const doublereal& re,
		const doublereal& im, const doublereal& h,
		doublereal& sigma, doublereal& omega,
		doublereal& csi, doublereal& freq)
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

   if (eEigenAnalysis == EIG_OUTPUTMATRICES) {
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
   
   if (bOutputModes) {
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
      bool doPlot = false;
      if (bOutputModes) {

         doublereal b = Beta.dGetCoef(iPage);
	 doublereal re = AlphaR.dGetCoef(iPage);
	 doublereal im = AlphaI.dGetCoef(iPage);
	 doublereal sigma;
	 doublereal omega;
	 doublereal csi;
	 doublereal freq;

	 do_eig(b, re, im, h, sigma, omega, csi, freq);

	 if (freq >= dLowerFreq && freq <= dUpperFreq) {
	    doPlot = true;
	      
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
	 
	 if (bOutputModes) {
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

	 if (bOutputModes) {
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
   if (bOutputModes) {
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


SolutionManager *const
Solver::AllocateSolman(integer iNLD, integer iLWS)
{
	SolutionManager *pCurrSM = CurrSolver.GetSolutionManager(iNLD, iLWS);
	
	/* special extra parameters if required */
	switch (CurrSolver.GetSolver()) {
	case LinSol::UMFPACK_SOLVER:
#if defined(USE_RTAI) && defined(HAVE_UMFPACK_TIC_DISABLE)
		if (bRT) {
			/* disable profiling, to avoid times() system call */
			umfpack_tic_disable();
		}
#endif /* USE_RTAI && HAVE_UMFPACK_TIC_DISABLE */
		break;

	default:
		break;
	}

	return pCurrSM;
};


SolutionManager *const
Solver::AllocateSchurSolman(integer iStates)
{
	SolutionManager *pSSM(NULL);

#ifdef USE_MPI
	doublereal dIPivotFactor = CurrIntSolver.dGetPivotFactor();
	integer iIWorkSpaceSize = CurrIntSolver.iGetWorkSpaceSize();

	switch (CurrIntSolver.GetSolver()) {
	case LinSol::Y12_SOLVER:
#ifdef USE_Y12
		SAFENEWWITHCONSTRUCTOR(pSSM, 
			SchurSolutionManager,
			SchurSolutionManager(iNumDofs, iStates, pLocDofs,
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
				SchurSolutionManager(iNumDofs, iStates,
					pLocDofs,
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
				SchurSolutionManager(iNumDofs, iStates,
					pLocDofs,
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

	default: 
		std::cerr << "Unknown interface solver. Aborting ..."
			<< std::endl;
		THROW(ErrGeneric());

		break;
	}
#else /* !USE_MPI */
		std::cerr << "Configure --with-mpi "
			"to enable schur solver" << std::endl;
		THROW(ErrGeneric());
#endif /* !USE_MPI */
	return pSSM;
};

NonlinearSolver *const
Solver::AllocateNonlinearSolver()
{
	NonlinearSolver *pNLS;

	switch (NonlinearSolverType) {
	case NonlinearSolver::MATRIXFREE:
		switch (MFSolverType) {
		case MatrixFreeSolver::BICGSTAB:
			SAFENEWWITHCONSTRUCTOR(pNLS,
					BiCGStab,
					BiCGStab(PcType, 
						iPrecondSteps,
						dIterTol, 
						iIterativeMaxSteps,
						dIterertiveEtaMax,
						dIterertiveTau));
			break;

		default:
			pedantic_cout("unknown matrix free solver type; "
					"using default" << std::endl);
			/* warning: should be unreachable */
			
		case MatrixFreeSolver::GMRES:
			SAFENEWWITHCONSTRUCTOR(pNLS,
					Gmres,
					Gmres(PcType, 
						iPrecondSteps,
						dIterTol, 
						iIterativeMaxSteps,
						dIterertiveEtaMax,
						dIterertiveTau));
			break;
		}
		break;

	default:
		pedantic_cout("unknown nonlinear solver type; using default"
				<< std::endl);

	case NonlinearSolver::NEWTONRAPHSON:
		SAFENEWWITHCONSTRUCTOR(pNLS,
				NewtonRaphsonSolver,
				NewtonRaphsonSolver(bTrueNewtonRaphson,
					bKeepJac,
					iIterationsBeforeAssembly));  
		break;
	}
	return pNLS;
}

void
Solver::SetupSolmans(integer iStates)
{
   	DEBUGLCOUT(MYDEBUG_MEM, "creating SolutionManager\n\tsize = "
		   << iNumDofs*iUnkStates << 
		   "\n\tnumdofs = " << iNumDofs <<
		   "\n\tnumstates = " iStates << std::endl);

	/*delete previous solmans*/
	if (pSM != 0) {
		SAFEDELETE(pSM);
		pSM = 0;
	}
	if (pLocalSM != 0) {
		SAFEDELETE(pLocalSM);
		pLocalSM = 0;
	}
	
	integer iWorkSpaceSize = CurrSolver.iGetWorkSpaceSize();
	integer iLWS = iWorkSpaceSize;
	integer iNLD = iNumDofs*iStates;
	if (bParallel) {
		/*FIXME BEPPE!*/
		iLWS = iWorkSpaceSize*iNumLocDofs/(iNumDofs*iNumDofs);
		/*FIXME: GIUSTO QUESTO?*/
		iNLD = iNumLocDofs*iStates;
	}

	SolutionManager *pCurrSM = AllocateSolman(iNLD, iLWS);

	/*
	 * This is the LOCAL solver if instantiating a parallel
	 * integrator; otherwise it is the MAIN solver
	 */
	if (bParallel) {
		pLocalSM = pCurrSM;

		/* Crea il solutore di Schur globale */
		pSM = AllocateSchurSolman(iStates);

	} else {
		pSM = pCurrSM;
	}
	/*
	 * FIXME: at present there MUST be a pSM
	 * (even for matrix-free nonlinear solvers)
	 */
	if (pSM == NULL) {
		std::cerr << "No linear solver defined" << std::endl;
		THROW(ErrGeneric());
	}
};
