/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

/* metodo multistep */

/* 
 * Copyright 1999-2000 Giuseppe Quaranta <giuquaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include <ac/float.h>
#include <ac/math.h>

#ifdef HAVE_SIGNAL
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif /* HAVE_SIGNAL_H */
#endif /* HAVE_SIGNAL */

#ifdef USE_MPI
#include <schur.h>
#include <mynewmem.h>

#ifdef MPI_PROFILING
extern "C" {
#include <mpe.h>
#include <stdio.h>
}
#endif /* MPI_PROFILING */

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

#ifdef HAVE_SIGNAL
extern volatile sig_atomic_t keep_going;
extern __sighandler_t sh_term;
extern __sighandler_t sh_int;
extern __sighandler_t sh_hup;

static void
modify_final_time_handler(int signum)
{
        ::keep_going = 0;
        switch (signum) {
        case SIGTERM:
                signal(signum, ::sh_term);
                break;

	case SIGINT:
		signal(signum, ::sh_int);
		break;

	case SIGHUP:
		signal(signum, ::sh_hup);
		break;
	}
}
#endif /* HAVE_SIGNAL */

/* SchurMultiStepIntegrator - begin */

/* Parametri locali */
const integer iDefaultMaxIterations = 1;

const integer iDefaultFictitiousStepsNumber = 0;
const doublereal dDefaultFictitiousStepsRatio = 1.e-3;
const doublereal dDefaultFictitiousStepsRho = 0.;
const doublereal dDefaultFictitiousStepsTolerance = 1.e-6;
const doublereal dDefaultToll = 1e-6;
const integer iDefaultIterationsBeforeAssembly = 2;
/*
 * Default solver
 */
SchurMultiStepIntegrator::SolverType
#if defined(USE_Y12)
SchurMultiStepIntegrator::defaultSolver = SchurMultiStepIntegrator::Y12_SOLVER;
#elif /* !USE_Y12 */ defined(USE_UMFPACK3)
SchurMultiStepIntegrator::defaultSolver = SchurMultiStepIntegrator::UMFPACK3_SOLVER;
#elif /* !USE_UMFPACK3 */ defined(USE_HARWELL)
SchurMultiStepIntegrator::defaultSolver = SchurMultiStepIntegrator::HARWELL_SOLVER;
#elif /* !USE_HARWELL */ defined(USE_MESCHACH)
SchurMultiStepIntegrator::defaultSolver = SchurMultiStepIntegrator::MESCHACH_SOLVER;
#else /* !USE_MESCHACH */
#error "need a solver!"
#endif /* !USE_MESCHACH */

/* Costruttore: esegue la simulazione */
SchurMultiStepIntegrator::SchurMultiStepIntegrator(MBDynParser& HPar,
                     const char* sInFName,
                     const char* sOutFName,
                     flag fPar)
: 
CurrStrategy(NOCHANGE),
CurrSolver(defaultSolver),
CurrIntSolver(defaultSolver),
sInputFileName(NULL),
sOutputFileName(NULL),
HP(HPar),
#ifdef __HACK_EIG__
fEigenAnalysis(0),
dEigParam(1.),
fOutputModes(0),
dUpperFreq(FLT_MAX),
dLowerFreq(0.),
#endif /* __HACK_EIG__ */
pdWorkSpace(NULL),
pXCurr(NULL),
pXPrimeCurr(NULL),
pXPrev(NULL),
pXPrimePrev(NULL),
pXPrev2(NULL),
pXPrimePrev2(NULL),
pSM(NULL),
pIntSM(NULL),
pDM(NULL),
pSDM(NULL),
DofIterator(), 
iNumDofs(0),
pDofs(NULL),
iNumLocDofs(0),
pLocDofs(NULL),
iNumIntDofs(0),
pIntDofs(NULL),
dTime(0.),
dInitialTime(0.),
dFinalTime(0.),
dRefTimeStep(0.),
dInitialTimeStep(1.),
dMinimumTimeStep(1.),
dToll(dDefaultToll),
iMaxIterations(iDefaultMaxIterations),
iFictitiousStepsNumber(iDefaultFictitiousStepsNumber),
dFictitiousStepsRatio(dDefaultFictitiousStepsRatio),
dFictitiousStepsRho(dDefaultFictitiousStepsRho),
dFictitiousStepsTolerance(dDefaultFictitiousStepsTolerance),
iFictitiousStepsMaxIterations(iDefaultMaxIterations),
dDerivativesToll(1e-6),
dDerivativesCoef(1.),
iDerivativesMaxIterations(iDefaultMaxIterations),
fAbortAfterInput(0),
fAbortAfterAssembly(0),
fAbortAfterDerivatives(0),
fAbortAfterFictitiousSteps(0),
fTrueNewtonRaphson(1),
iIterationsBeforeAssembly(0),
iPerformedIterations(0),
iStepsAfterReduction(0),
iStepsAfterRaise(0),
fLastChance(0),
pMethod(NULL),
pFictitiousStepsMethod(NULL),
db0Differential(0.),
db0Algebraic(0.),
iLWorkSpaceSize(0),
dLPivotFactor(1.),
iIWorkSpaceSize(0),
dIPivotFactor(1.),
fParallel(fPar)
{
	DEBUGCOUTFNAME(sClassName() << "::SchurMultiStepIntegrator");

	if (sInFName != NULL) {
		SAFESTRDUP(sInputFileName, sInFName);
	}
	if (sOutFName != NULL) {
		SAFESTRDUP(sOutputFileName, sOutFName);
	}

   	/* Legge i dati relativi al metodo di integrazione */
   	ReadData(HP);
}

void
SchurMultiStepIntegrator::Run(void)
{
	DEBUGCOUTFNAME("SchurMultiStepIntegrator::Run");
#ifdef MPI_PROFILING
   	int  err = MPE_Init_log();
   	if (err) {
     		std::cerr << "Impossible to write jumpshot log file" << std::endl;
   	}
   	if (MyRank == 0) {
     		MPE_Describe_state(1, 2, "Create Partition", "yellow");
     		MPE_Describe_state(3, 4, "Test Broadcast", "steel blue");
     		MPE_Describe_state(5, 6, "Initialize Communications", "blue");
     		MPE_Describe_state(7, 8, "Prediction", "midnight blue");
     		MPE_Describe_state(13, 14, "ISend", "light salmon");
     		MPE_Describe_state(15, 16, "Local Update", "green");
     		MPE_Describe_state(19, 20, "IRecv", "cyan");
     		MPE_Describe_state(23, 24, "Compute Jac", "DarkGreen");
     		MPE_Describe_state(27, 28, "Compute Res", "orchid");
     		MPE_Describe_state(29, 30, "Output", "khaki1");
     		MPE_Describe_state(31, 32, "Solve Local", "firebrick");
     		MPE_Describe_state(35, 36, "Solve Schur", "orange");
     		MPE_Describe_state(33, 34, "Rotor Trust Exchange", "grey");
     		MPE_Describe_state(41,42, "Ass. Schur", "blue");
   	}
   
   	MPE_Start_log();
#endif /* MPI_PROFILING */
   	/*
    	 * E' necessario poter determinare in questa routine
	 * quale e' il master in modo da far calcolare la soluzione
	 * solo su di esso
	 */
   int MyRank = MPI::COMM_WORLD.Get_rank();
   /* chiama il gestore dei dati generali della simulazione */


   /* I file di output vengono stampati localmente da ogni processo aggiungendo
      al termine del OutputFileName il rank del processo */
   int iOutLen;
   char* sNewOutName = NULL;
   if (sOutputFileName == NULL) {
     iOutLen = strlen(sInputFileName);
     SAFENEWARR(sNewOutName, char, iOutLen+4+1);
     strcpy(sNewOutName, sInputFileName);
   } else {
     iOutLen = strlen(sOutputFileName);
     SAFENEWARR(sNewOutName, char, iOutLen+4+1);
     strcpy(sNewOutName, sOutputFileName);
   }
   sprintf(sNewOutName+iOutLen,".%.3d", MyRank);
   
   DEBUGLCOUT(MYDEBUG_MEM, "creating ParallelDataManager" << std::endl);
   SAFENEWWITHCONSTRUCTOR(pSDM,
			  SchurDataManager,
			  SchurDataManager(HP,
					     dInitialTime,
					     sInputFileName,
					     sNewOutName,
					     fAbortAfterInput));
   pDM = pSDM;

   /* Si fa dare l'ostream al file di output per il log */
   std::ostream& Out = pDM->GetOutFile();
   
   if (fAbortAfterInput) {
     /* Esce */
     Out << "End of Input; no simulation or assembly is required." 
     		<< std::endl;
     return;
   } else if (fAbortAfterAssembly) {
     /* Fa l'output dell'assemblaggio iniziale e poi esce */
     pDM->Output();
     Out << "End of Initial Assembly; no simulation is required." 
     		<< std::endl;
     return;
   }

    /* Qui crea le partizioni: principale fra i processi, se parallelo  */
#ifdef MPI_PROFILING
   MPE_Log_event(1, 0, "start");
#endif /* MPI_PROFILING */ 
   pSDM->CreatePartition();
#ifdef MPI_PROFILING
   MPE_Log_event(2, 0, "end");
#endif /* MPI_PROFILING */ 


   /* Si fa dare il DriveHandler e linka i drivers di rho ecc. */
   const DriveHandler* pDH = pDM->pGetDrvHdl();
   pMethod->SetDriveHandler(pDH);
   pFictitiousStepsMethod->SetDriveHandler(pDH);

   /* Costruisce i vettori della soluzione ai vari passi */
   DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << std::endl);
   
   iNumDofs = pSDM->HowManyDofs(SchurDataManager::TOTAL);
   pDofs = pSDM->pGetDofsList();

   iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
   pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);
   iNumIntDofs = pSDM->HowManyDofs(SchurDataManager::INTERNAL);
   pIntDofs = pSDM->GetDofsList(SchurDataManager::INTERNAL);

   ASSERT(iNumDofs > 0);
   
   SAFENEWARR(pdWorkSpace, doublereal, 6*iNumDofs);
   SAFENEWWITHCONSTRUCTOR(pXCurr,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace));
   SAFENEWWITHCONSTRUCTOR(pXPrimeCurr,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace+iNumDofs));
   SAFENEWWITHCONSTRUCTOR(pXPrev,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace+2*iNumDofs));
   SAFENEWWITHCONSTRUCTOR(pXPrimePrev,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace+3*iNumDofs));
   SAFENEWWITHCONSTRUCTOR(pXPrev2,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace+4*iNumDofs));
   SAFENEWWITHCONSTRUCTOR(pXPrimePrev2,
			  MyVectorHandler,
			  MyVectorHandler(iNumDofs, pdWorkSpace+5*iNumDofs));

   /* Resetta i vettori */
   pXCurr->Reset(0.);
   pXPrimeCurr->Reset(0.);
   pXPrev->Reset(0.);
   pXPrimePrev->Reset(0.);
   pXPrev2->Reset(0.);
   pXPrimePrev2->Reset(0.);

   /* Subito collega il DataManager alla soluzione corrente */
   pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);

   /* costruisce il SolutionManager */
   DEBUGLCOUT(MYDEBUG_MEM, "creating SolutionManager, size = " << iNumDofs << std::endl);
   

    /* Crea il solutore locale */
    integer iLocWorkSpaceSize = iLWorkSpaceSize/(iNumDofs*iNumDofs)* iNumLocDofs;
    switch (CurrSolver) {
     	case Y12_SOLVER: 
#ifdef USE_Y12
      		SAFENEWWITHCONSTRUCTOR(pSM,
			Y12SparseLUSolutionManager,
			Y12SparseLUSolutionManager(iNumLocDofs,
				iLocWorkSpaceSize,
				dLPivotFactor == -1. ? 1. : dLPivotFactor));
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

    	case MESCHACH_SOLVER:
		std::cerr << "Sorry Meschach cannot be used as a local parallel "
			<< "solver. Switching to Harwell...." << std::endl;
   	case HARWELL_SOLVER:
#ifdef USE_HARWELL
      		SAFENEWWITHCONSTRUCTOR(pSM,
			HarwellSparseLUSolutionManager,
			HarwellSparseLUSolutionManager(iNumLocDofs,
				iLocWorkSpaceSize,
				dLPivotFactor == -1. ? 1. : dLPivotFactor));
      		break;
#else /* !USE_HARWELL */
      		std::cerr << "Configure with --with-harwell "
			"to enable Harwell solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_HARWELL */

   	case UMFPACK3_SOLVER:
#ifdef USE_UMFPACK3
      		SAFENEWWITHCONSTRUCTOR(pSM,
			Umfpack3SparseLUSolutionManager,
			Umfpack3SparseLUSolutionManager(iNumLocDofs, 
				0, dLPivotFactor));
      		break;
#else /* !USE_UMFPACK3 */
      		std::cerr << "Configure with --with-umfpack3 "
			"to enable Umfpack3 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_UMFPACK3 */
   	}

   
    /* Crea il solutore di Schu globale */
    switch (CurrIntSolver) {
     	case Y12_SOLVER: 
#ifdef USE_Y12
		{ 
		  Y12SparseLUSolutionManager* pSM;
   		  SAFENEWWITHCONSTRUCTOR(pIntSM,
			  SchurSolutionManager,
			  SchurSolutionManager(iNumDofs, pLocDofs, 
						iNumLocDofs, 
						pIntDofs, iNumIntDofs,
						pSM,
						pSM,
						iIWorkSpaceSize, 
						dIPivotFactor== -1.? 1. : dIPivotFactor));
		}
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

   	case HARWELL_SOLVER:
      		std::cerr << "Harwell solver cannot be used as interface"
			"solver. Switching to Meschach..." << std::endl;

    	case MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		{ 
		MeschachSparseLUSolutionManager* pSM;
   		SAFENEWWITHCONSTRUCTOR(pIntSM,
			  SchurSolutionManager,
			  SchurSolutionManager(iNumDofs, pLocDofs, 
						iNumLocDofs, 
						pIntDofs, iNumIntDofs,
						pSM,
						pSM,
						iIWorkSpaceSize, 
						dIPivotFactor== -1.? 1. : dIPivotFactor));
		}
      		break;
#else /* !USE_MESCHACH */
      		std::cerr << "Configure with --with-meschach "
			"to enable Meschach solver"
			<< std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

   	case UMFPACK3_SOLVER:
#ifdef USE_UMFPACK3
		{
		Umfpack3SparseLUSolutionManager* pSM;
   		SAFENEWWITHCONSTRUCTOR(pIntSM,
			  SchurSolutionManager,
			  SchurSolutionManager(iNumDofs, pLocDofs, 
						iNumLocDofs, 
						pIntDofs, iNumIntDofs,
						pSM,
						pSM,
						0, dIPivotFactor));
      		}
		break;
#else /* !USE_UMFPACK3 */
      		std::cerr << "Configure with --with-umfpack3 "
			"to enable Umfpack3 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_UMFPACK3 */
   	}


   /* Puntatori agli handlers del solution manager */
   VectorHandler* pRes = pIntSM->pResHdl();
   VectorHandler* pSol = pIntSM->pSolHdl();
   MatrixHandler* pJac = pIntSM->pMatHdl();

   /* Legenda:
    *   MS - SchurMultiStepIntegrator
    *   DM - DataManager
    *   OM - DofManager
    *   NM - NodeManager
    *   EM - ElemManager
    *   SM - SolutionManager */


   /* cicli vari */


   /* Dell'assemblaggio iniziale dei vincoli se ne occupa il DataManager
    * in quanto e' lui il responsabile dei dati della simulazione,
    * e quindi anche della loro coerenza. Inoltre e' lui a sapere
    * quali equazioni sono di vincolo o meno. */

   /* Dialoga con il DataManager per dargli il tempo iniziale
    * e per farsi inizializzare i vettori di soluzione e derivata */
   dTime = dInitialTime;
   pDM->SetTime(dTime);
   pDM->SetValue(*pXCurr, *pXPrimeCurr);
   DofIterator = pDM->GetDofIterator();

#ifdef __HACK_EIG__
   if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
	 Eig();
	 OneEig.fDone = flag(1);
   }
#endif /* __HACK_EIG__ */

   integer iTotIter = 0;
   doublereal dTotErr = 0.;

    /* calcolo delle derivate */
   DEBUGLCOUT(MYDEBUG_DERIVATIVES, "derivatives solution step" << std::endl);
   doublereal dTest = 1.;

#ifdef HAVE_SIGNAL
   	::sh_term = signal(SIGTERM, modify_final_time_handler);
   	::sh_int = signal(SIGINT, modify_final_time_handler);
   	::sh_hup = signal(SIGHUP, modify_final_time_handler);
   
   	if (::sh_term == SIG_IGN) {
      		signal (SIGTERM, SIG_IGN);
   	}
   	if (::sh_int == SIG_IGN) {
      		signal (SIGINT, SIG_IGN);
   	}
   	if (::sh_hup == SIG_IGN) {
      		signal (SIGHUP, SIG_IGN);
   	}
#endif /* HAVE_SIGNAL */

   int iIterCnt = 0;
   while (1) {
   
#ifdef MPI_PROFILING
	MPE_Log_event(27, 0, "start");
#endif /* MPI_PROFILING */ 

     	pRes->Reset(0.);
     	pDM->AssRes(*pRes, dDerivativesCoef);

#ifdef MPI_PROFILING
     	MPE_Log_event(28, 0, "end");
     	MPE_Log_event(3, 0, "start");
#endif /* MPI_PROFILING */

     	dTest = this->MakeTest(*pRes, *pXPrimeCurr);

#ifdef MPI_PROFILING
     	MPE_Log_event(4, 0, "end");
#endif /* MPI_PROFILING */

#ifdef DEBUG
     	if (DEBUG_LEVEL_MATCH(MYDEBUG_DERIVATIVES|MYDEBUG_RESIDUAL)) {
       		std::cout << "Residual:" << std::endl;
       		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
	   			<< pRes->dGetCoef(iTmpCnt) << std::endl;
       		}
     	}
#endif /* DEBUG */
     	if (dTest < dDerivativesToll) {
       		goto EndOfDerivatives;
     	}
     
     	iIterCnt++;
     
     	DEBUGLCOUT(MYDEBUG_DERIVATIVES, 
		"calculating derivatives, iteration "
		<< std::setw(4) << iIterCnt << ", test = " << dTest << std::endl);

     	if (iIterCnt > iDerivativesMaxIterations || !isfinite(dTest)) {
       		std::cerr << std::endl 
			<< "Maximum iterations number " << iIterCnt
	    		<< " has been reached during initial derivatives calculation;"
	    		<< std::endl << "aborting ..." << std::endl;
       		pDM->Output();

       		THROW(SchurMultiStepIntegrator::ErrMaxIterations());
     	}
#ifdef MPI_PROFILING
     	MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */
 
     	pIntSM->MatrInit(0.);
     	pDM->AssJac(*pJac, dDerivativesCoef);
    
#ifdef MPI_PROFILING
     	MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */
    
     	pIntSM->Solve();
  

#ifdef DEBUG
     /* Output della soluzione */
     	if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_DERIVATIVES)) {
       		std::cout << "Solution:" << std::endl;
       		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
        		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
	     			<< pSol->dGetCoef(iTmpCnt) << std::endl;
       		}
     	}
#endif /* debug */
     
#ifdef MPI_PROFILING
     	MPE_Log_event(15, 0, "start local Update");
#endif /* MPI-PROFILING */
     
     	this->Update(*pSol);
     	pDM->DerivativesUpdate();
     
#ifdef MPI_PROFILING
     	MPE_Log_event(16, 0, "end local Update");
#endif /* MPI-PROFILING */
   }
   
EndOfDerivatives:
   
   dTotErr += dTest;
   iTotIter += iIterCnt;
   
   Out << "Derivatives solution step at time " << dInitialTime
       << " performed in " << iIterCnt
       << " iterations with " << dTest
       << " error" << std::endl;
   
   DEBUGCOUT("Derivatives solution step has been performed successfully in "
	     << iIterCnt << " iterations" << std::endl);
   
   if (fAbortAfterDerivatives) {
     	/* Fa l'output della soluzione delle derivate iniziali ed esce */
     	pDM->Output();
     	Out << "End of derivatives; no simulation is required." << std::endl;
     	return;
#ifdef HAVE_SIGNAL
   } else if (!::keep_going) {
   	/*
	 * Fa l'output della soluzione delle derivate iniziali ed esce
	 */
      	pDM->Output();
      	Out << "Interrupted during derivatives computation." << std::endl;
      	return;
#endif /* HAVE_SIGNAL */
   }
   
   /* Dati comuni a passi fittizi e normali */
   integer iStep = 1;
   doublereal dCurrTimeStep = 0.;
   /* First step prediction must always be Crank-Nicholson for accuracy */
   CrankNicholson cn;
   
   if (iFictitiousStepsNumber > 0) {
      	/* passi fittizi */
    
     	/* inizio integrazione: primo passo a predizione lineare con sottopassi
      	 * di correzione delle accelerazioni e delle reazioni vincolari */
     	pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *pXPrev, *pXPrimePrev);
     	this->Flip();
     
     	/***********************************************************************
      	* primo passo fittizio
      	***********************************************************************/
     
     	/* Passo ridotto per step fittizi di messa a punto */
     	dRefTimeStep = dInitialTimeStep*dFictitiousStepsRatio;
     	dCurrTimeStep = dRefTimeStep;
     	pDM->SetTime(dTime+dCurrTimeStep);
     
     	DEBUGLCOUT(MYDEBUG_FSTEPS, "Current time step: " << dCurrTimeStep << std::endl);

     	/* First step prediction must always be Crank-Nicholson for accuracy */
     	cn.SetCoef(dRefTimeStep, 1., MultiStepIntegrationMethod::NEWSTEP, db0Differential, db0Algebraic);
#ifdef MPI_PROFILING
     	MPE_Log_event(7, 0, "start Predict");
#endif

     	FirstStepPredict(&cn);

#ifdef MPI_PROFILING
     	MPE_Log_event(8, 0, "end Predict");
#endif

     	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
     	pDM->AfterPredict();
     
#ifdef DEBUG
     	if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_PRED)) {
       		std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	    		<< std::endl;
       		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 		std::cout << std::setw(4) << iTmpCnt << ": "
	      			<< std::setw(12) << pXCurr->dGetCoef(iTmpCnt)
	      			<< std::setw(12) << pXPrev->dGetCoef(iTmpCnt)
	      			<< std::setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	      			<< std::setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	      			<< std::setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	      			<< std::setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	      			<< std::endl;
       		}
     }
#endif /* DEBUG */

     iIterCnt = 0;
     while (1) {
       
     	/* l02: EM calcolo del residuo */
       
#ifdef MPI_PROFILING
       	MPE_Log_event(27, 0, "start");
#endif /* MPI_PROFILING */ 
       
       	pRes->Reset(0.);
       	pDM->AssRes(*pRes, db0Differential);

#ifdef MPI_PROFILING
     	MPE_Log_event(28, 0, "end");
        MPE_Log_event(3, 0, "start");
#endif /* MPI_PROFILING */ 

       dTest = this->MakeTest(*pRes, *pXPrimeCurr);

#ifdef MPI_PROFILING
        MPE_Log_event(4, 0, "end");
#endif /* MPI_PROFILING */  
       
       
#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_RESIDUAL)) {
	 	std::cout << "Residual:" << std::endl;
	 	for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
				<< pRes->dGetCoef(iTmpCnt) << std::endl;
	 	}
       }
#endif
       
       if (dTest < dFictitiousStepsTolerance) {
	 goto EndOfFirstFictitiousStep;
       }
       
       iIterCnt++;
       if (iIterCnt > iFictitiousStepsMaxIterations || !isfinite(dTest)) {
	 	std::cerr << std::endl << "Maximum iterations number " << iIterCnt
	      		<< " has been reached during first dummy step;"
	      		<< std::endl << "time step dt = " << dCurrTimeStep
	      		<< " cannot be reduced further;" << std::endl
	      		<< "aborting ..." << std::endl;
	 	pDM->Output(); 
	 	THROW(SchurMultiStepIntegrator::ErrMaxIterations());
       }

#ifdef MPI_PROFILING
       MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */ 
       
       pIntSM->MatrInit(0.);
       pDM->AssJac(*pJac, db0Differential);

#ifdef MPI_PROFILING
        MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */ 
    
       pIntSM->Solve();
     
#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	 	std::cout << "Solution:" << std::endl;
	 	for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
				<< pSol->dGetCoef(iTmpCnt) << std::endl;
	 	}
       }
#endif
       
       
#ifdef MPI_PROFILING
       MPE_Log_event(15, 0, "start local Update");
#endif /* MPI-PROFILING */
       
       this->Update(*pSol);
       pDM->Update();

#ifdef MPI_PROFILING
       MPE_Log_event(16, 0, "end local Update");
#endif /* MPI-PROFILING */
    }
   
EndOfFirstFictitiousStep:
     
    dRefTimeStep = dCurrTimeStep;
    dTime += dRefTimeStep;
   
    dTotErr += dTest;
    iTotIter += iIterCnt;

#ifdef HAVE_SIGNAL
     if (!::keep_going) {
	 /*
	  * Fa l'output della soluzione delle derivate iniziali
	  * ed esce
	  */
#ifdef DEBUG_FICTITIOUS
	 pDM->Output();
#endif /* DEBUG_FICTITIOUS */
	 Out << "Interrupted during first dummy step." << std::endl;
	 return;
      }
#endif /* HAVE_SIGNAL */   

#ifdef DEBUG_FICTITIOUS
     pDM->Output();
#endif
     
     /**********************************************************************
      * Passi fittizi successivi
      **********************************************************************/
     for (int iSubStep = 2; iSubStep <= iFictitiousStepsNumber; iSubStep++) {
       	pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *pXPrev, *pXPrimePrev);
       	this->Flip();
     
       	DEBUGLCOUT(MYDEBUG_FSTEPS, "Fictitious step " << iSubStep
		  << "; current time step: " << dCurrTimeStep << std::endl);
     
       	pDM->SetTime(dTime+dCurrTimeStep);
       	ASSERT(pFictitiousStepsMethod != NULL);
       	pFictitiousStepsMethod->SetCoef(dRefTimeStep,
				       dCurrTimeStep/dRefTimeStep,
				       MultiStepIntegrationMethod::NEWSTEP,
				       db0Differential,
				       db0Algebraic);
#ifdef MPI_PROFILING
       	MPE_Log_event(7, 0, "start Predict");
#endif

       	Predict(pFictitiousStepsMethod);

#ifdef MPI_PROFILING
       	MPE_Log_event(8, 0, "end Predict");
#endif

	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
       	pDM->AfterPredict();

#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_PRED)) {
	 	std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	      		<< std::endl;
	 	for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   		std::cout << std::setw(4) << iTmpCnt << ": "
				<< std::setw(12) << pXCurr->dGetCoef(iTmpCnt)
				<< std::setw(12) << pXPrev->dGetCoef(iTmpCnt)
				<< std::setw(12) << pXPrev2->dGetCoef(iTmpCnt)
				<< std::setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
				<< std::setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
				<< std::setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
				<< std::endl;
	 	}
      	}
#endif /* DEBUG */
       
       iIterCnt = 0;
       while (1) {
#ifdef MPI_PROFILING
         MPE_Log_event(27, 0, "start");
#endif /* MPI_PROFILING */ 
	 
	 pRes->Reset(0.);
	 pDM->AssRes(*pRes, db0Differential);
	 
#ifdef MPI_PROFILING
         MPE_Log_event(28, 0, "end");
#endif /* MPI_PROFILING */

#ifdef MPI_PROFILING
         MPE_Log_event(3, 0, "start");
#endif /* MPI_PROFILING */
 	 	 
	 dTest = this->MakeTest(*pRes, *pXPrimeCurr);

#ifdef MPI_PROFILING
         MPE_Log_event(4, 0, "end");
#endif /* MPI_PROFILING */	 
	 
#ifdef DEBUG
	 if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_RESIDUAL)) {
	   	std::cout << "Residual:" << std::endl;
	   	for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	     		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
		  		<< pRes->dGetCoef(iTmpCnt) << std::endl;
	   	}
	 }
#endif /* DEBUG */
	 
	 if (dTest < dFictitiousStepsTolerance) {
	   	goto EndOfFictitiousStep;
	 }
	 
	 iIterCnt++;
	 if (iIterCnt > iFictitiousStepsMaxIterations || !isfinite(dTest)) {
	   	std::cerr << std::endl << "Maximum iterations number " << iIterCnt
			<< " has been reached during dummy step "
			<< iSubStep << ';' << std::endl
			<< "aborting ..." << std::endl;
	   	pDM->Output();
	   	THROW(SchurMultiStepIntegrator::ErrMaxIterations());
	 }

#ifdef MPI_PROFILING
         MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */  
	 
	 pIntSM->MatrInit(0.);
	 pDM->AssJac(*pJac, db0Differential);
	 
#ifdef MPI_PROFILING
	 MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */
	 
	 pIntSM->Solve();
	 
#ifdef DEBUG
	 if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	   	std::cout << "Solution:" << std::endl;
	   	for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	     		std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
		  		<< pSol->dGetCoef(iTmpCnt) << std::endl;
	   	}
	 }
#endif /* DEBUG */
	 
#ifdef MPI_PROFILING
	 MPE_Log_event(15, 0, "start local Update");
#endif /* MPI-PROFILING */
	 
	 this->Update(*pSol);
	 pDM->Update();

#ifdef MPI_PROFILING
	 MPE_Log_event(16, 0, "end local Update");
#endif /* MPI-PROFILING */
       }
     
EndOfFictitiousStep:
     
       dTotErr += dTest;
       iTotIter += iIterCnt;
     
#ifdef DEBUG
       if (DEBUG_LEVEL(MYDEBUG_FSTEPS)) {
		 pDM->Output();
	 
	 	Out << "Step " << iStep
	     		<< " at time " << dTime+dCurrTimeStep
	     		<< " time step " << dCurrTimeStep
	     		<< " performed in " << iIterCnt
	     		<< " iterations with " << dTest
	     		<< " error" << std::endl;
       }
#endif /* DEBUG */
       
       DEBUGLCOUT(MYDEBUG_FSTEPS, "Substep " << iSubStep
		  << " of step " << iStep
		  << " has been completed successfully in "
		  << iIterCnt << " iterations" << std::endl);

#ifdef HAVE_SIGNAL
	if (!::keep_going) {
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
     
     Out << "Initial solution after dummy steps at time "
	 << dTime+dCurrTimeStep
	 << " performed in " << iIterCnt
	 << " iterations with " << dTest
	 << " error" << std::endl;
     
     DEBUGLCOUT(MYDEBUG_FSTEPS,
		"Fictitious steps have been completed successfully in "
		<< iIterCnt << " iterations" << std::endl);
   } /* Fine dei passi fittizi */


   /* Output delle "condizioni iniziali" */
   pDM->Output();

   Out << "Step " << 0
       << " at time " << dTime+dCurrTimeStep
       << " with time step " << dCurrTimeStep
       << " performed in " << iIterCnt
       << " iterations with " << dTest
       << " error" << std::endl;
   
   if (fAbortAfterFictitiousSteps) {
     Out << "End of dummy steps; no simulation is required." << std::endl;
     return;
#ifdef HAVE_SIGNAL
   } else if (!::keep_going) {
      	/* Fa l'output della soluzione ed esce */
      	Out << "Interrupted during dummy steps." << std::endl;
      	return;
#endif /* HAVE_SIGNAL */
   }

   iStep = 1; /* Resetto di nuovo iStep */
   
   DEBUGCOUT("Step " << iStep << " has been completed successfully in "
	     << iIterCnt << " iterations" << std::endl);
   
   
   pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *pXPrev, *pXPrimePrev);
   this->Flip();

   dRefTimeStep = dInitialTimeStep;
   dCurrTimeStep = dRefTimeStep;
   
   DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   
   /**************************************************************************
    * Primo passo regolare
    **************************************************************************/

 IfFirstStepIsToBeRepeated:
   pDM->SetTime(dTime+dCurrTimeStep);
   
   /* doublereal dRho = Rho.pGetDriveCaller()->dGet();
    * doublereal dRhoAlgebraic = RhoAlgebraic.pGetDriveCaller()->dGet();
    * metto rho = 1 perche' cosi' rispetto certi teoremi sulla precisione di
    * questo passo (Petzold, 89) (ricado nella regola dei trapezi) */
   cn.SetCoef(dRefTimeStep, 1., MultiStepIntegrationMethod::NEWSTEP, db0Differential, db0Algebraic);
#ifdef MPI_PROFILING
   MPE_Log_event(7, 0, "start Predict");
#endif

   FirstStepPredict(&cn);

#ifdef MPI_PROFILING
   MPE_Log_event(8, 0, "end Predict");
#endif

   pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
   pDM->AfterPredict();
   
#ifdef DEBUG
   if (DEBUG_LEVEL(MYDEBUG_PRED)) {
     std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	  << std::endl;
     for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
       std::cout << std::setw(4) << iTmpCnt << ": "
	    << std::setw(12) << pXCurr->dGetCoef(iTmpCnt)
	    << std::setw(12) << pXPrev->dGetCoef(iTmpCnt)
	    << std::setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	    << std::setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	    << std::setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	    << std::setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	    << std::endl;
     }
   }
#endif /* DEBUG */
   
   iIterCnt = 0;
   iPerformedIterations = iIterationsBeforeAssembly;
   while (1) {
#ifdef MPI_PROFILING
     MPE_Log_event(27,0,"start Residual");
#endif
     
     pRes->Reset(0.);
     pDM->AssRes(*pRes, db0Differential);

#ifdef MPI_PROFILING
     MPE_Log_event(28,0,"end Residual");
#endif

#ifdef MPI_PROFILING
     MPE_Log_event(3,0,"start");
#endif

     dTest = this->MakeTest(*pRes, *pXPrimeCurr);

#ifdef MPI_PROFILING
     MPE_Log_event(4,0,"end");
#endif


#ifdef DEBUG
     if (DEBUG_LEVEL(MYDEBUG_RESIDUAL)) {
       std::cout << "Residual:" << std::endl;
       std::cout << iStep  << "   " << iIterCnt <<std::endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
	      << pRes->dGetCoef(iTmpCnt) << std::endl;
       }
     }
#endif /* DEBUG */

     if (dTest < dToll) {
       goto EndOfFirstStep;
     }
     

     iIterCnt++;
     if (iIterCnt > iMaxIterations || !isfinite(dTest)) {
     if (dCurrTimeStep > dMinimumTimeStep) {
       /* Riduce il passo */
       dCurrTimeStep = this->NewTimeStep(dCurrTimeStep,
					 iIterCnt,
					 MultiStepIntegrationMethod::REPEATSTEP);
       dRefTimeStep = dCurrTimeStep;
       DEBUGCOUT("Changing time step during first step after "
		 << iIterCnt << " iterations" << std::endl);
       goto IfFirstStepIsToBeRepeated;
       
     } else {
       std::cerr << std::endl << "Maximum iterations number " << iIterCnt
	    << " has been reached during first step;"
	    << std::endl << "time step dt = " << dCurrTimeStep
	    << " cannot be reduced further;" << std::endl
	    << "aborting ..." << std::endl;
       
       pDM->Output();
       
       THROW(SchurMultiStepIntegrator::ErrMaxIterations());
     }
     }

     if (iPerformedIterations < iIterationsBeforeAssembly) {
       iPerformedIterations++;
     } else {
       iPerformedIterations = 0;
#ifdef MPI_PROFILING
       MPE_Log_event(23,0,"start Jacobian");
#endif

       pIntSM->MatrInit(0.);
       pDM->AssJac(*pJac, db0Differential);

#ifdef MPI_PROFILING
       MPE_Log_event(24,0,"end Jacobian");
#endif
      }

     pIntSM->Solve();
     
#ifdef DEBUG
     if (DEBUG_LEVEL(MYDEBUG_SOL)) {
       std::cout << "Solution:" << std::endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 std::cout << "Dof" << setw(4) << iTmpCnt << ": "
	      << pSol->dGetCoef(iTmpCnt) << std::endl;
       }
     }
#endif /* DEBUG */


#ifdef MPI_PROFILING
     MPE_Log_event(15, 0, "start local Update");
#endif /* MPI-PROFILING */

     this->Update(*pSol);
     pDM->Update();

#ifdef MPI_PROFILING
     MPE_Log_event(16, 0, "end local Update");
#endif /* MPI-PROFILING */
   }
   
EndOfFirstStep:
   
   pDM->Output();

#ifdef HAVE_SIGNAL
   	if (!::keep_going) {
      		/* Fa l'output della soluzione al primo passo ed esce */
      		Out << "Interrupted during first dummy step." << std::endl;
      		return;
   	} else {
#endif /* HAVE_SIGNAL */
   
   Out << "Step " << iStep
       << " at time " << dTime+dCurrTimeStep
       << " with time step " << dCurrTimeStep
       << " performed in " << iIterCnt
       << " iterations with " << dTest
       << " error" << std::endl;
       
#ifdef HAVE_SIGNAL
   	}
#endif /* HAVE_SIGNAL */   
   
   dRefTimeStep = dCurrTimeStep;
   dTime += dRefTimeStep;
   
   dTotErr += dTest;
   iTotIter += iIterCnt;

#ifdef __HACK_EIG__
   if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
	 Eig();
	OneEig.fDone = flag(1);
   }
#endif /* __HACK_EIG__ */
   
   /*************************************************************************
    * Altri passi regolari
    *************************************************************************/
   while (1) {
     MultiStepIntegrationMethod::StepChange CurrStep
       = MultiStepIntegrationMethod::NEWSTEP;
     
     if (dTime >= dFinalTime) {
       std::cout << "End of simulation at time " << dTime << " after "
	    << iStep << " steps;" << std::endl
	    << "total iterations: " << iTotIter << std::endl
	    << "total error: " << dTotErr << std::endl;
       pDM->MakeRestart();
       return;
#ifdef HAVE_SIGNAL
      } else if (!::keep_going) {
	 std::cout << "Interrupted!" << std::endl
	   	<< "Simulation ended at time "
		<< dTime << " after " 
		<< iStep << " steps;" << std::endl
		<< "total iterations: " << iTotIter << std::endl
		<< "total error: " << dTotErr << std::endl;
	 return;
#endif /* HAVE_SIGNAL */
     }
     
     iStep++;

     pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *pXPrev, *pXPrimePrev);
     this->Flip();
     
IfStepIsToBeRepeated:
     pDM->SetTime(dTime+dCurrTimeStep);
     pMethod->SetCoef(dRefTimeStep,
		      dCurrTimeStep/dRefTimeStep,
		      CurrStep,
		      db0Differential,
		      db0Algebraic);
#ifdef MPI_PROFILING
     MPE_Log_event(7, 0, "start Predict");
#endif

     Predict(pMethod);

#ifdef MPI_PROFILING
     MPE_Log_event(8, 0, "end Predict");
#endif

     pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
     pDM->AfterPredict();
     
#ifdef DEBUG
     if (DEBUG_LEVEL(MYDEBUG_PRED)) {
       std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	    << std::endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 std::cout << std::setw(4) << iTmpCnt << ": "
	      << std::setw(12) << pXCurr->dGetCoef(iTmpCnt)
	      << std::setw(12) << pXPrev->dGetCoef(iTmpCnt)
	      << std::setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	      << std::setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	      << std::setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	      << std::setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	      << std::endl;
       }
     }
#endif /* DEBUG */
     
     iIterCnt = 0;
     iPerformedIterations = iIterationsBeforeAssembly;
     while (1) {
#ifdef MPI_PROFILING
      MPE_Log_event(27,0,"start Residual");
#endif

     pRes->Reset(0.);
     pDM->AssRes(*pRes, db0Differential);
#ifdef MPI_PROFILING
     		MPE_Log_event(28,0,"end Residual");
#endif

#ifdef USE_EXCEPTIONS
	try {
#endif /* USE_EXCEPTIONS */

#ifdef MPI_PROFILING
     		MPE_Log_event(3,0,"start");
#endif

     		dTest = this->MakeTest(*pRes, *pXPrimeCurr);
              
#ifdef MPI_PROFILING
     		MPE_Log_event(4,0,"end");
#endif
#ifdef USE_EXCEPTIONS
	}
	catch (SchurMultiStepIntegrator::ErrSimulationDiverged) {
		/*
		 * Mettere qui eventuali azioni speciali 
		 * da intraprendere in caso di errore ...
		 */
#if 0
		std::cerr << *pJac << std::endl;
#endif /* 0 */
		throw;
	}
#endif /* USE_EXCEPTIONS */       

#ifdef DEBUG
       if (DEBUG_LEVEL(MYDEBUG_RESIDUAL)) {
	 std::cout << "Residual:" <<std::endl;
	 std::cout << iStep  << "   " << iIterCnt <<std::endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
           std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
		<< pRes->dGetCoef(iTmpCnt) << std::endl;
	 }
       }
#endif /* DEBUG */
       
       if (dTest < dToll) {
	 CurrStep = MultiStepIntegrationMethod::NEWSTEP;
	 goto EndOfStep;
       }
       
       iIterCnt++;
       if (iIterCnt > iMaxIterations || !isfinite(dTest)) {
	 if (dCurrTimeStep > dMinimumTimeStep) {
	   /* Riduce il passo */
	   CurrStep = MultiStepIntegrationMethod::REPEATSTEP;
	   dCurrTimeStep = this->NewTimeStep(dCurrTimeStep,
					     iIterCnt, CurrStep);
	   DEBUGCOUT("Changing time step during step "
		     << iStep << " after "
		     << iIterCnt << " iterations" << std::endl);
	   goto IfStepIsToBeRepeated;
	   
	 } else {
	   std::cerr << std::endl << "Maximum iterations number " << iIterCnt
		<< " has been reached during step " << iStep << ';'
		<< std::endl << "time step dt = " << dCurrTimeStep
		<< " cannot be reduced further;" << std::endl
		<< "aborting ..." << std::endl;
	   
	   /* pDM->Output();  temporarily commented */
	   THROW(SchurMultiStepIntegrator::ErrMaxIterations());
	 }
       }
       
       if (iPerformedIterations < iIterationsBeforeAssembly) {
	 iPerformedIterations++;
       } else {
	 iPerformedIterations = 0;
	 
#ifdef MPI_PROFILING
	 MPE_Log_event(23,0,"start Jacobian");
#endif
	 pIntSM->MatrInit(0.);
	 pDM->AssJac(*pJac, db0Differential);

#ifdef MPI_PROFILING
	 MPE_Log_event(24,0,"end Jacobian");
#endif

       }
       pIntSM->Solve();
       
#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_MPI)) {
	 std::cout << "Solution:" << std::endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
		<< pSol->dGetCoef(iTmpCnt) << std::endl;
	 }
       }
#endif /* DEBUG */
       
       
#ifdef MPI_PROFILING
       MPE_Log_event(15, 0, "start local Update");
#endif /* MPI-PROFILING */

       this->Update(*pSol);
       pDM->Update();

#ifdef MPI_PROFILING
       MPE_Log_event(16, 0, "end local Update");
#endif /* MPI-PROFILING */
     }
     
EndOfStep:
     
     dTotErr += dTest;
     iTotIter += iIterCnt;
     
     pDM->Output();
     
     Out << "Step " << iStep
	 << " at time " << dTime+dCurrTimeStep
	 << " with time step " << dCurrTimeStep
	 << " performed in " << iIterCnt
	 << " iterations with " << dTest
	 << " error" << std::endl;
     
     DEBUGCOUT("Step " << iStep << " has been completed successfully in "
	       << iIterCnt << " iterations" << std::endl);
     
     dRefTimeStep = dCurrTimeStep;
     dTime += dRefTimeStep;
     
     
#ifdef __HACK_EIG__
      if (fEigenAnalysis && OneEig.dTime <= dTime && !OneEig.fDone) {
	Eig();
	OneEig.fDone = flag(1);
      }
#endif /* __HACK_EIG__ */
     
     
     /* Calcola il nuovo timestep */
     dCurrTimeStep = this->NewTimeStep(dCurrTimeStep, iIterCnt, CurrStep);
     DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   }
}


/* Distruttore */
SchurMultiStepIntegrator::~SchurMultiStepIntegrator(void)
{
   DEBUGCOUTFNAME(sClassName() << "::~SchurMultiStepIntegrator");
#ifdef MPI_PROFILING
   MPE_Finish_log("mbdyn.mpi");
#endif /* MPI_PROFILING */

   if (sInputFileName != NULL) {
      SAFEDELETEARR(sInputFileName);
   }

   if (sOutputFileName != NULL) {
      SAFEDELETEARR(sOutputFileName);
   }

   if (pIntSM != NULL)  {
      SAFEDELETE(pIntSM);
   }

   if (pXPrimePrev2 != NULL) {
      SAFEDELETE(pXPrimePrev2);
   }

   if (pXPrev2 != NULL) {
      SAFEDELETE(pXPrev2);
   }

   if (pXPrimePrev != NULL) {
      SAFEDELETE(pXPrimePrev);
   }

   if (pXPrev != NULL) {
      SAFEDELETE(pXPrev);
   }

   if (pXPrimeCurr != NULL) {
      SAFEDELETE(pXPrimeCurr);
   }

   if (pXCurr != NULL) {
      SAFEDELETE(pXCurr);
   }

   if (pdWorkSpace != NULL) {
      SAFEDELETEARR(pdWorkSpace);
   }
   
   if (pDM != NULL) {
     SAFEDELETE(pDM);
   }
}


/* Test sul residuo */
doublereal
SchurMultiStepIntegrator::MakeTest(const VectorHandler& Res,
				     const VectorHandler& XP)
{
  DEBUGCOUTFNAME(sClassName() << "::MakeTest");
   
  Dof CurrDof;
   
  doublereal dRes = 0.;
  doublereal dXPr = 0.;
  
  /*
   * chiama la routine di comunicazione per la trasmissione del residuo
   * delle interfacce
   */
  pIntSM->StartExchInt();
 
  /* calcola il test per i dofs locali */
  int DCount = 0;
  for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
    DCount = pLocDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    doublereal d = Res.dGetCoef(DCount);
    dRes += d*d;
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      d = XP.dGetCoef(DCount);
      dXPr += d*d;
    }
    /* else if ALGEBRAIC: non aggiunge nulla */
  }

  for (int iCntp1 = 0; 
		  iCntp1 < pSDM->HowManyDofs(SchurDataManager::MYINTERNAL); 
		  iCntp1++) {
    DCount = (pSDM->GetDofsList(SchurDataManager::MYINTERNAL))[iCntp1];
    CurrDof = pDofs[DCount-1];
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      doublereal d = XP.dGetCoef(DCount);
      dXPr += d*d;
    }
    /* else if ALGEBRAIC: non aggiunge nulla */
  }
  
  /* verifica completamento trasmissioni */
  pIntSM->ComplExchInt(dRes, dXPr);

  
  dRes /= (1.+dXPr);

  if (!isfinite(dRes)) {
    std::cerr << "The simulation diverged; aborting ..." << std::endl;
    THROW(SchurMultiStepIntegrator::ErrSimulationDiverged());
  }
  return sqrt(dRes);
}


/* Predizione al primo passo */
void
SchurMultiStepIntegrator::FirstStepPredict(MultiStepIntegrationMethod* pM)
{
  DEBUGCOUTFNAME(sClassName() << "::FirstStepPredict");
  
  Dof CurrDof;
  

  /* Combinazione lineare di stato e derivata al passo precedente ecc. */
  /* dofs locali */
  int DCount =  0;
  for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
    DCount = pLocDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      doublereal dXnm1 = pXPrev->dGetCoef(DCount);
      doublereal dXPnm1 = pXPrimePrev->dGetCoef(DCount);

      doublereal dXPn = pM->dPredDer(dXnm1, 0., dXPnm1, 0.);
      pXPrimeCurr->fPutCoef(DCount, dXPn);
      pXCurr->fPutCoef(DCount, pM->dPredState(dXnm1, 0., dXPn, dXPnm1, 0.));

    } else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
      doublereal dXnm1 = pXPrev->dGetCoef(DCount);
      
      doublereal dXn = pM->dPredDerAlg(0., dXnm1, 0.);
      pXCurr->fPutCoef(DCount, dXn);
      pXPrimeCurr->fPutCoef(DCount, pM->dPredStateAlg(0., dXn, dXnm1, 0.));
      
    } else {
      std::cerr << sClassName() << "::FirstStepPredict(): unknown dof order" << std::endl;
      THROW(ErrGeneric());
    }
  }

  /* Combinazione lineare di stato e derivata al passo precedente ecc. */
  /* dofs interfaccia */
 
  DCount =  0;
  for (int iCntp1 = 0; iCntp1 < iNumIntDofs; iCntp1++) {
    DCount = pIntDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      doublereal dXnm1 = pXPrev->dGetCoef(DCount);
      doublereal dXPnm1 = pXPrimePrev->dGetCoef(DCount);

      doublereal dXPn = pM->dPredDer(dXnm1, 0., dXPnm1, 0.);
      pXPrimeCurr->fPutCoef(DCount, dXPn);
      pXCurr->fPutCoef(DCount, pM->dPredState(dXnm1, 0., dXPn, dXPnm1, 0.));

    } else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
      doublereal dXnm1 = pXPrev->dGetCoef(DCount);
      
      doublereal dXn = pM->dPredDerAlg(0., dXnm1, 0.);
      pXCurr->fPutCoef(DCount, dXn);
      pXPrimeCurr->fPutCoef(DCount, pM->dPredStateAlg(0., dXn, dXnm1, 0.));
      
    } else {
      std::cerr << sClassName() << "::FirstStepPredict(): unknown dof order" << std::endl;
      THROW(ErrGeneric());
    }
  }
}


/* Predizione al passo generico */
void
SchurMultiStepIntegrator::Predict(MultiStepIntegrationMethod* pM)
{
   // Note: pM must be initialised prior to calling Predict()

   DEBUGCOUTFNAME(sClassName() << "::Predict");

   Dof CurrDof;
  

   /* Combinazione lineare di stato e derivata al passo precedente ecc. */
   /* Dofs locali */
   int DCount = 0;
   for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
     DCount = pLocDofs[iCntp1];
     CurrDof = pDofs[DCount-1];
     if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
       doublereal dXnm1 = pXPrev->dGetCoef(DCount);
       doublereal dXnm2 = pXPrev2->dGetCoef(DCount);
       doublereal dXPnm1 = pXPrimePrev->dGetCoef(DCount);
       doublereal dXPnm2 = pXPrimePrev2->dGetCoef(DCount);
       
     doublereal dXPn = pM->dPredDer(dXnm1, dXnm2, dXPnm1, dXPnm2);
     pXPrimeCurr->fPutCoef(DCount, dXPn);
     pXCurr->fPutCoef(DCount, pM->dPredState(dXnm1, dXnm2, dXPn, dXPnm1, dXPnm2));

     } else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
       doublereal dXnm1 = pXPrev->dGetCoef(DCount);
       doublereal dXnm2 = pXPrev2->dGetCoef(DCount);
       doublereal dXInm1 = pXPrimePrev->dGetCoef(DCount);
       
       doublereal dXn = pM->dPredDerAlg(dXInm1, dXnm1, dXnm2);
       pXCurr->fPutCoef(DCount, dXn);
       pXPrimeCurr->fPutCoef(DCount, pM->dPredStateAlg(dXInm1, dXn, dXnm1, dXnm2));
       
      } else {
#ifdef __GNUC__
     std::cerr << __FUNCTION__ << ": ";
#else
     std::cerr << sClassName() << "::Predict(): ";
#endif /* __GNUC__ */
     std::cerr << "unknown dof order" << std::endl;
     THROW(ErrGeneric());
      }
   }


   /* Dofs interfaccia */
   DCount = 0;
   for (int iCntp1 = 0; iCntp1 < iNumIntDofs; iCntp1++) {
     DCount = pIntDofs[iCntp1];
     CurrDof = pDofs[DCount-1];
     if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
       doublereal dXnm1 = pXPrev->dGetCoef(DCount);
       doublereal dXnm2 = pXPrev2->dGetCoef(DCount);
       doublereal dXPnm1 = pXPrimePrev->dGetCoef(DCount);
       doublereal dXPnm2 = pXPrimePrev2->dGetCoef(DCount);
       
     doublereal dXPn = pM->dPredDer(dXnm1, dXnm2, dXPnm1, dXPnm2);
     pXPrimeCurr->fPutCoef(DCount, dXPn);
     pXCurr->fPutCoef(DCount, pM->dPredState(dXnm1, dXnm2, dXPn, dXPnm1, dXPnm2));

      } else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
     doublereal dXnm1 = pXPrev->dGetCoef(DCount);
     doublereal dXnm2 = pXPrev2->dGetCoef(DCount);
     doublereal dXInm1 = pXPrimePrev->dGetCoef(DCount);

     doublereal dXn = pM->dPredDerAlg(dXInm1, dXnm1, dXnm2);
     pXCurr->fPutCoef(DCount, dXn);
     pXPrimeCurr->fPutCoef(DCount, pM->dPredStateAlg(dXInm1, dXn, dXnm1, dXnm2));

      } else {
#ifdef __GNUC__
     std::cerr << __FUNCTION__ << ": ";
#else
     std::cerr << sClassName() << "::Predict(): ";
#endif /* __GNUC__ */
     std::cerr << "unknown dof order" << std::endl;
     THROW(ErrGeneric());
      }
   }
}


/* Nuovo delta t */
doublereal SchurMultiStepIntegrator::NewTimeStep(doublereal dCurrTimeStep,
						   integer iPerformedIters,
						   MultiStepIntegrationMethod::StepChange Why )
{
   DEBUGCOUTFNAME(sClassName() << "::NewTimeStep");

   switch (CurrStrategy) {
    case FACTOR: {
       if (Why == MultiStepIntegrationMethod::REPEATSTEP) {
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
        /* Fuori viene intercettato il valore illegale */
        return dCurrTimeStep*StrategyFactor.dReductionFactor;
         }
      }
       }

       if (Why == MultiStepIntegrationMethod::NEWSTEP) {
      iStepsAfterReduction++;
      iStepsAfterRaise++;

      if (iPerformedIters <= StrategyFactor.iMinIters
          && iStepsAfterReduction > StrategyFactor.iStepsBeforeReduction
          && iStepsAfterRaise > StrategyFactor.iStepsBeforeRaise
          && dCurrTimeStep < dMaxTimeStep) {
         iStepsAfterRaise = 0;
         return dCurrTimeStep*StrategyFactor.dRaiseFactor;
      }
      return dCurrTimeStep;
       }

       break;
    }

    case NOCHANGE: {
       return dCurrTimeStep;
    }

    default: {
       std::cerr << "You shouldn't have reached this point!" << std::endl;
       THROW(SchurMultiStepIntegrator::ErrGeneric());
    }
   }

   return dCurrTimeStep;
}


/* Aggiornamento della soluzione nel passo fittizio */
void
SchurMultiStepIntegrator::DerivativesUpdate(const VectorHandler& Sol)
{
  DEBUGCOUTFNAME(sClassName() << "::DerivativesUpdate");

  Dof CurrDof;
  
  /* dofs locali */
  int DCount = 0;
  for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
    DCount = pLocDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    doublereal d = Sol.dGetCoef(DCount);
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
       pXPrimeCurr->fIncCoef(DCount, d);
    } else {
      pXCurr->fIncCoef(DCount, d);
    }
  }
  
  /* dofs interfaccia */
  DCount = 0;
  for (int iCntp1 = 0; iCntp1 < iNumIntDofs; iCntp1++) {
    DCount = pIntDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    doublereal d = Sol.dGetCoef(DCount);
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      pXPrimeCurr->fIncCoef(DCount, d);
    } else {
      pXCurr->fIncCoef(DCount, d);
    }
  }
}


/* Aggiornamento normale */
void
SchurMultiStepIntegrator::Update(const VectorHandler& Sol)
{
  DEBUGCOUTFNAME(sClassName() << "::Update");
  
  Dof CurrDof;
  
  /* dofs locali */
  int DCount = 0;
  for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
    DCount = pLocDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    doublereal d = Sol.dGetCoef(DCount);
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      pXPrimeCurr->fIncCoef(DCount, d);
      /* Nota: b0Differential e b0Algebraic possono essere distinti;
       * in ogni caso sono calcolati dalle funzioni di predizione
       * e sono dati globali */
      pXCurr->fIncCoef(DCount, db0Differential*d);
    } else {
      pXCurr->fIncCoef(DCount, d);
      pXPrimeCurr->fIncCoef(DCount, db0Algebraic*d);
    }
  }
  
  /* dofs interfaccia locale */
  DCount = 0;
  for (int iCntp1 = 0; iCntp1 < iNumIntDofs; iCntp1++) {
    DCount = pIntDofs[iCntp1];
    CurrDof = pDofs[DCount-1];
    doublereal d = Sol.dGetCoef(DCount);
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      pXPrimeCurr->fIncCoef(DCount, d);
      /* Nota: b0Differential e b0Algebraic possono essere distinti;
       * in ogni caso sono calcolati dalle funzioni di predizione
       * e sono dati globali */
      pXCurr->fIncCoef(DCount, db0Differential*d);
    } else {
      pXCurr->fIncCoef(DCount, d);
      pXPrimeCurr->fIncCoef(DCount, db0Algebraic*d);
    }
  }
}
/* Dati dell'integratore */
void 
SchurMultiStepIntegrator::ReadData(MBDynParser& HP)
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

		"eigen" "analysis",
		"output" "modes",
		
		"solver",
		"interface" "solver", 
		"harwell",
		"meschach",
		"y12",
		"umfpack3"
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
	
		EIGENANALYSIS,
		OUTPUTMODES,
		
		SOLVER,
		INTERFACESOLVER,
		HARWELL,
		MESCHACH,
		Y12,
		UMFPACK3,
	
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
      		THROW(SchurMultiStepIntegrator::ErrGeneric());
   	}
   
   	if (KeyWords(HP.GetWord()) != MULTISTEP) {
      		std::cerr << std::endl << "Error: <begin: multistep;> expected at line " 
			<< HP.GetLineData() << "; aborting ..." << std::endl;
      		THROW(SchurMultiStepIntegrator::ErrGeneric());
   	}

   	flag fMethod(0);
   	flag fFictitiousStepsMethod(0);      
     
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
	     			THROW(SchurMultiStepIntegrator::ErrGeneric());
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
	      			THROW(SchurMultiStepIntegrator::ErrGeneric());
	  		}
	  		break;
       		}
	 
       		case METHOD: {
	  		if (fMethod) {
	     			std::cerr << "error: multiple definition"
					" of integration method at line "
					<< HP.GetLineData() << std::endl;
	     			THROW(SchurMultiStepIntegrator::ErrGeneric());
	  		}
	  		fMethod = flag(1);
	        	  
	  		KeyWords KMethod = KeyWords(HP.GetWord());
	  		switch (KMethod) {
	   		case CRANKNICHOLSON:
	      			SAFENEW(pMethod,
		      			CrankNicholson); /* no constructor */
	      			break;
				
			case BDF: {
				DriveCaller* pRho = NULL;
				SAFENEWWITHCONSTRUCTOR(pRho,
						NullDriveCaller, 
						NullDriveCaller(NULL));
				DriveCaller* pRhoAlgebraic = NULL;
				SAFENEWWITHCONSTRUCTOR(pRhoAlgebraic,
						NullDriveCaller, 
						NullDriveCaller(NULL));

		  		SAFENEWWITHCONSTRUCTOR(pMethod,
				 	NostroMetodo,
				 	NostroMetodo(pRho, pRhoAlgebraic));
		  		break;
			}
			
	   		case NOSTRO:
	   		case MS:
	   		case HOPE: {
	      			DriveCaller* pRho =
					ReadDriveData(NULL, HP, NULL);
				HP.PutKeyTable(K);

	      			DriveCaller* pRhoAlgebraic = NULL;
				if (HP.fIsArg()) {
					pRhoAlgebraic = ReadDriveData(NULL, 
							HP, NULL);
					HP.PutKeyTable(K);
				} else {
					pRhoAlgebraic = pRho->pCopy();
				}
			
	      			switch (KMethod) {
	       			case NOSTRO: 
	       			case MS:
		  			SAFENEWWITHCONSTRUCTOR(pMethod,
					 	NostroMetodo,
					 	NostroMetodo(pRho,
							     pRhoAlgebraic));
		  			break;
		 
	       			case HOPE:	      
		  			SAFENEWWITHCONSTRUCTOR(pMethod,
					 	Hope,
						Hope(pRho, pRhoAlgebraic));
		  			break;
					
	       			default:
	          			THROW(ErrGeneric());
	      			}
	      			break;
	   		}
	   		default:
	      			std::cerr << "Unknown integration method at line "
					<< HP.GetLineData() << std::endl;
				THROW(SchurMultiStepIntegrator::ErrGeneric());
	  		}
	  		break;
       		}

       case FICTITIOUSSTEPSMETHOD:
       case DUMMYSTEPSMETHOD: {
	  if (fFictitiousStepsMethod) {
	     std::cerr << "error: multiple definition of dummy steps integration"
	       " method at line "
	       << HP.GetLineData();
	     THROW(SchurMultiStepIntegrator::ErrGeneric());
	  }
	  fFictitiousStepsMethod = flag(1);	  	
	  
	  KeyWords KMethod = KeyWords(HP.GetWord());
	  switch (KMethod) {
	   case CRANKNICHOLSON:
	      SAFENEW(pFictitiousStepsMethod,
		      CrankNicholson); /* no constructor */
	      break;

           case BDF: {
	      DriveCaller* pRho = NULL;
	      SAFENEWWITHCONSTRUCTOR(pRho,
			NullDriveCaller,
			NullDriveCaller(NULL));
  	      DriveCaller* pRhoAlgebraic = NULL;
	      SAFENEWWITHCONSTRUCTOR(pRhoAlgebraic,
			NullDriveCaller,
			NullDriveCaller(NULL));

	      SAFENEWWITHCONSTRUCTOR(pFictitiousStepsMethod,
			NostroMetodo,
			NostroMetodo(pRho, pRhoAlgebraic));
	      break;
	   }
	   
	   case NOSTRO:
	   case MS:
	   case HOPE: {	      	     
	      DriveCaller* pRho = ReadDriveData(NULL, HP, NULL);
	      HP.PutKeyTable(K);

	      DriveCaller* pRhoAlgebraic;
	      if (HP.fIsArg()) {
	      	 pRhoAlgebraic = ReadDriveData(NULL, HP, NULL);
		 HP.PutKeyTable(K);
	      } else {
		 pRhoAlgebraic = pRho->pCopy();
	      }
	      HP.PutKeyTable(K);
	      
	      switch (KMethod) {
	       case NOSTRO:
	       case MS: {
		  SAFENEWWITHCONSTRUCTOR(pFictitiousStepsMethod,
					 NostroMetodo,
					 NostroMetodo(pRho, pRhoAlgebraic));
		  break;
	       }
		 
	       case HOPE: {	      
		  SAFENEWWITHCONSTRUCTOR(pFictitiousStepsMethod,
					 Hope,
					 Hope(pRho, pRhoAlgebraic));
		  break;
	       }
	       default:
	          THROW(ErrGeneric());
	      }
	      break;	      
	   }
	   default: {
	      std::cerr << "Unknown integration method at line " << HP.GetLineData() << std::endl;
	      THROW(SchurMultiStepIntegrator::ErrGeneric());
	   }	     
	  }	  
	  break;
       }	 

       case TOLERANCE: {
	  dToll = HP.GetReal();
	  if (dToll <= 0.) {
	     dToll = dDefaultToll;
	     std::cerr 
	       << "warning, tolerance <= 0. is illegal; switching to default value "
	       << dToll << std::endl;
	  }		       		  
	  DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dToll << std::endl);
	  break;
       }
	 
       case DERIVATIVESTOLERANCE: {
	  dDerivativesToll = HP.GetReal();
	  if (dDerivativesToll <= 0.) {
	     dDerivativesToll = 1e-6;
	     std::cerr 
	       << "warning, derivatives tolerance <= 0. is illegal; switching to default value " 
	       << dDerivativesToll
	       << std::endl;
	  }		       		  
	  DEBUGLCOUT(MYDEBUG_INPUT, 
		     "Derivatives toll = " << dDerivativesToll << std::endl);
	  break;
       }
	 
       case MAXITERATIONS: {
	  iMaxIterations = HP.GetInt();
	  if (iMaxIterations < 1) {
	     iMaxIterations = iDefaultMaxIterations;
	     std::cerr 
	       << "warning, max iterations < 1 is illegal; switching to default value "
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
	       << "warning, derivatives max iterations < 1 is illegal; switching to default value "
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
	       << "warning, derivatives coefficient <= 0. is illegal; switching to default value "
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
	     THROW(SchurMultiStepIntegrator::ErrGeneric());
	  }
	  goto EndOfCycle;
       }
	 
       case STRATEGY: {
	  switch (KeyWords(HP.GetWord())) {
	     
	   case STRATEGYFACTOR: {
	      CurrStrategy = FACTOR;
	      
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
		 std::cerr << "Warning, illegal number of steps before reduction at line "
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
		 std::cerr << "Warning, illegal number of steps before raise at line "
		   << HP.GetLineData() << ';' << std::endl
		   << "default value 1 will be used (it may be dangerous)" 
		   << std::endl;
		 StrategyFactor.iStepsBeforeRaise = 1;
	      }
	      
	      StrategyFactor.iMinIters = HP.GetInt();
	      if (StrategyFactor.iMinIters <= 0) {
		 std::cerr << "Warning, illegal minimum number of iterations at line "
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
	     
	   default: {
	      std::cerr << "Unknown time step control strategy at line "
		<< HP.GetLineData() << std::endl;
	      THROW(SchurMultiStepIntegrator::ErrGeneric());
	   }		 
	  }
	  
	  break;
       }
	 
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
#else /* !__HACK_EIG__ */
	  HP.GetReal();
	  if (HP.IsKeyWord("parameter")) {
	     HP.GetReal();
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
	  if (HP.IsKeyWord("yes")) {
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
	     CurrSolver = MESCHACH_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Using meschach sparse LU solver" << std::endl);
	     break;
#endif /* USE_MESCHACH */

	   case Y12:
#ifdef USE_Y12
             CurrSolver = Y12_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using y12 sparse LU solver" << std::endl);
	     break;
#endif /* USE_Y12 */
							       
	   case UMFPACK3:
#ifdef USE_UMFPACK3
             CurrSolver = UMFPACK3_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using umfpack3 sparse LU solver" << std::endl);
	     break;
#endif /* USE_UMFPACK3 */

#ifdef USE_HARWELL
	   case HARWELL: 
	     CurrSolver = HARWELL_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Using harwell sparse LU solver" << std::endl);	 
	     break;	   
#endif /* USE_HARWELL */

	   default:
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Unknown solver; switching to default" << std::endl);
	     break;
	  }
  
	  if (HP.IsKeyWord("workspacesize")) {
	     iLWorkSpaceSize = HP.GetInt();
	     if (iLWorkSpaceSize < 0) {
		iLWorkSpaceSize = 0;
	     }
	  }
	  
	  if (HP.IsKeyWord("pivotfactor")) {
	     dLPivotFactor = HP.GetReal();
	     if (dLPivotFactor <= 0. || dLPivotFactor > 1.) {
		dLPivotFactor = 1.;
	     }
	  }
	  
	  DEBUGLCOUT(MYDEBUG_INPUT, "Workspace size: " << iLWorkSpaceSize 
		    << ", pivor factor: " << dLPivotFactor << std::endl);
	  break;
       }
       case INTERFACESOLVER: {
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
             CurrIntSolver = Y12_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using y12 sparse LU solver" << std::endl);
	     break;
#endif /* USE_Y12 */
							       
	   case UMFPACK3:
#ifdef USE_UMFPACK3
             CurrIntSolver = UMFPACK3_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT,
			"Using umfpack3 sparse LU solver" << std::endl);
	     break;
#endif /* USE_UMFPACK3 */

#ifdef USE_HARWELL
	   case HARWELL: 
	        CurrIntSolver = MESCHACH_SOLVER;
		std::cerr << "Harwell solver cannot be used as interface "
		"solver. Meschach will be used ..." << std::endl;
		break;
#endif /* USE_HARWELL */		 	
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
       }
	 	 
       default: {
	  std::cerr << std::endl << "Unknown description at line " 
	    << HP.GetLineData() << "; aborting ..." << std::endl;
	  THROW(SchurMultiStepIntegrator::ErrGeneric());
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
   if (!fMethod) {
      ASSERT(pMethod == NULL);
      
      DriveCaller* pRho = NULL;
      SAFENEWWITHCONSTRUCTOR(pRho,
		             NullDriveCaller,
			     NullDriveCaller(NULL));
      
      /* DriveCaller per Rho asintotico per variabili algebriche */     
      DriveCaller* pRhoAlgebraic = NULL;
      SAFENEWWITHCONSTRUCTOR(pRhoAlgebraic,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
      
      SAFENEWWITHCONSTRUCTOR(pMethod,
			     NostroMetodo,
			     NostroMetodo(pRho, pRhoAlgebraic));
   }

   /* Metodo di integrazione di default */
   if (!fFictitiousStepsMethod) {
      ASSERT(pFictitiousStepsMethod == NULL);
      
      DriveCaller* pRho = NULL;
      SAFENEWWITHCONSTRUCTOR(pRho,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
                 
      /* DriveCaller per Rho asintotico per variabili algebriche */     
      DriveCaller* pRhoAlgebraic = NULL;
      SAFENEWWITHCONSTRUCTOR(pRhoAlgebraic,
			     NullDriveCaller,
			     NullDriveCaller(NULL));
      
      SAFENEWWITHCONSTRUCTOR(pFictitiousStepsMethod,
			     NostroMetodo,
			     NostroMetodo(pRho, pRhoAlgebraic));
   }
   
   return;
}
   

#ifdef __HACK_EIG__
/* Estrazione autovalori, vincolata alla disponibilita' delle LAPACK */
void
SchurMultiStepIntegrator::Eig(void)
{
	std::cerr << "SchurMultiStepIntegrator::Eig() not available" << std::endl;
}
#endif /* __HACK_EIG__ */

/* SchurMultiStepIntegrator - end */

#endif /* !USE_MPI */

