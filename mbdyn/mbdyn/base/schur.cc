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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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
 * Copyright 1999-2000 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_MPI

#include <schur.h>
#include <mynewmem.h>
#include <memmans.h>

#ifdef USE_MESCHACH
#include <mschwrap.h>
#endif /* USE_MESCHACH */

#ifdef MPI_PROFILING
extern "C" {
#include <mpe.h>
#include <stdio.h>
}
#endif /* MPI_PROFILING */

/* SchurMultiStepIntegrator - begin */

#ifdef DEBUG_MEMMANAGER
clMemMan MSmm(SchurMultiStepIntegrator::sClassName());
#endif /* DEBUG_MEMMANAGER */

/* Parametri locali */
const integer iDefaultMaxIterations = 1;

const integer iDefaultFictitiousStepsNumber = 0;
const doublereal dDefaultFictitiousStepsRatio = 1.e-3;
const doublereal dDefaultFictitiousStepsRho = 0.;
const doublereal dDefaultFictitiousStepsTolerance = 1.e-6;
const doublereal dDefaultToll = 1e-6;
const integer iDefaultIterationsBeforeAssembly = 2;

/* Costruttore: esegue la simulazione */
SchurMultiStepIntegrator::SchurMultiStepIntegrator(MBDynParser& HPar,
                     const char* sInFName,
                     const char* sOutFName,
                     flag fPar)
: sInputFileName(NULL),
sOutputFileName(NULL),
HP(HPar),
CurrStrategy(NOCHANGE),
CurrSolver(HARWELL_SOLVER),
pdWorkSpace(NULL),
pXCurr(NULL),
pXPrimeCurr(NULL),
pXPrev(NULL),
pXPrimePrev(NULL),
pXPrev2(NULL),
pXPrimePrev2(NULL),
pSM(NULL),
pDM(NULL),
DofIterator(), 
iNumDofs(0),
pDofs(NULL),
iNumLocDofs(0),
pLocDofs(NULL),
iNumIntDofs(0),
pIntDofs(NULL),
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
#ifdef __HACK_EIG__
fEigenAnalysis(0),		/********* TEMPORARY *******/
#endif /* __HACK_EIG__ */
iWorkSpaceSize(0),
dPivotFactor(1.),
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
   
   DEBUGLCOUT(MYDEBUG_MEM, "creating ParallelDataManager" << endl);
   SAFENEWWITHCONSTRUCTOR(pDM,
			  SchurDataManager,
			  SchurDataManager(HP,
					     dInitialTime,
					     sInputFileName,
					     sNewOutName,
					     fAbortAfterInput));
   
#ifdef MPI_PROFILING
   int  err = MPE_Init_log();
   if (err) {
     cerr << "Impossible to write jumpshot log file" << endl;
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

   /* Si fa dare l'ostream al file di output per il log */
   ostream& Out = pDM->GetOutFile();
   
   if (fAbortAfterInput) {
     /* Esce */
     Out << "End of Input; no simulation or assembly is required." << endl;
     return;
   } else if (fAbortAfterAssembly) {
     /* Fa l'output dell'assemblaggio iniziale e poi esce */
     pDM->Output();
     Out << "End of Initial Assembly; no simulation is required." << endl;
     return;
   }

    /* Qui crea le partizioni: principale fra i processi, se parallelo  */
#ifdef MPI_PROFILING
   MPE_Log_event(1, 0, "start");
#endif /* MPI_PROFILING */ 
   pDM->CreatePartition();
#ifdef MPI_PROFILING
   MPE_Log_event(2, 0, "end");
#endif /* MPI_PROFILING */ 


   /* Si fa dare il DriveHandler e linka i drivers di rho ecc. */
   const DriveHandler* pDH = pDM->pGetDrvHdl();
   pMethod->SetDriveHandler(pDH);
   pFictitiousStepsMethod->SetDriveHandler(pDH);

   /* Costruisce i vettori della soluzione ai vari passi */
   DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << endl);
   
   iNumDofs = pDM->HowManyDofs(1);
   pDofs = pDM->pGetDofsList();

   iNumLocDofs = pDM->HowManyDofs(2);
   pLocDofs = pDM->GetDofsList(2);
   iNumIntDofs = pDM->HowManyDofs(3);
   pIntDofs = pDM->GetDofsList(3);

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
   DEBUGLCOUT(MYDEBUG_MEM, "creating SolutionManager, size = " << iNumDofs << endl);
   
#ifdef USE_MESCHACH
   if (CurrSolver == MESCHACH_SOLVER) {
     cerr << "Sorry ! Schur method currently implemented with Harwell Solver only." << endl;
     THROW(ErrGeneric());
   } else if (CurrSolver == HARWELL_SOLVER) {
#endif
     SAFENEWWITHCONSTRUCTOR(pSM,
			    SchurSolutionManager,
			    SchurSolutionManager(iNumDofs, pLocDofs, 
						 iNumLocDofs, 
						 pIntDofs, iNumIntDofs,
						 iWorkSpaceSize, dPivotFactor));
#ifdef USE_MESCHACH
   } else {
     cerr << "don't know which solver to use!" << endl;
     THROW(ErrGeneric());
   }
#endif

   /* Puntatori agli handlers del solution manager */
   VectorHandler* pRes = pSM->pResHdl();
   VectorHandler* pSol = pSM->pSolHdl();
   MatrixHandler* pJac = pSM->pMatHdl();

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
   doublereal dTime = dInitialTime;
   pDM->SetTime(dTime);
   pDM->SetValue(*pXCurr, *pXPrimeCurr);
   DofIterator = pDM->GetDofIterator();

#ifdef __HACK_EIG__
   /************************
    * TEMPORANEO: FA UN'EIGENANALYSIS */
   if (fEigenAnalysis) {
     if (OneEig.dTime <= dTime && !OneEig.fDone) {
       Eig();
       OneEig.fDone = flag(1);
     }
   }
#endif /* __HACK_EIG__ */

   integer iTotIter = 0;
   doublereal dTotErr = 0.;

   /**************************************************************************
    * calcolo delle derivate
    **************************************************************************/
   DEBUGLCOUT(MYDEBUG_DERIVATIVES, "derivatives solution step" << endl);
   doublereal dTest = 1.;

   int iIterCnt = 0;
   while (1) {
   
#ifdef MPI_PROFILING
     MPE_Log_event(27, 0, "start");
#endif /* MPI_PROFILING */ 
     pRes->Reset(0.);
     pDM->AssRes(*pRes, dDerivativesCoef);
#ifdef MPI_PROFILING
     MPE_Log_event(28, 0, "end");
#endif /* MPI_PROFILING */

#ifdef DEBUG
     if (DEBUG_LEVEL_MATCH(MYDEBUG_DERIVATIVES|MYDEBUG_RESIDUAL)) {
       cout << "Residual:" << endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 cout << "Dof" << setw(4) << iTmpCnt << ": " 
	   << pRes->dGetCoef(iTmpCnt) << endl;
       }
     }
#endif /* DEBUG */
     
#ifdef MPI_PROFILING
     MPE_Log_event(3, 0, "start");
#endif /* MPI_PROFILING */
     dTest = this->MakeTest(*pRes, *pXPrimeCurr);
#ifdef MPI_PROFILING
     MPE_Log_event(4, 0, "end");
#endif /* MPI_PROFILING */

     if (dTest < dDerivativesToll) {
       goto EndOfDerivatives;
     }
     
     iIterCnt++;
     
     DEBUGLCOUT(MYDEBUG_DERIVATIVES, "calculating derivatives, iteration "
		<< setw(4) << iIterCnt << ", test = " << dTest << endl);

     if (iIterCnt > iDerivativesMaxIterations || !isfinite(dTest)) {
       cerr << endl << "Maximum iterations number " << iIterCnt
	    << " has been reached during initial derivatives calculation;"
	    << endl << "aborting ..." << endl;
       
       pDM->Output();

       THROW(SchurMultiStepIntegrator::ErrMaxIterations());
     }

#ifdef MPI_PROFILING
     MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */ 
     pSM->MatrInit(0.);
     pDM->AssJac(*pJac, dDerivativesCoef);
    
#ifdef MPI_PROFILING
     MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */
    
     pSM->Solve();
  

#ifdef DEBUG
     /* Output della soluzione */
     if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_DERIVATIVES)) {
       cout << "Solution:" << endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
        cout << "Dof" << setw(4) << iTmpCnt << ": "
	     << pSol->dGetCoef(iTmpCnt) << endl;
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
       << " error" << endl;
   
   DEBUGCOUT("Derivatives solution step has been performed successfully in "
	     << iIterCnt << " iterations" << endl);
   
   if (fAbortAfterDerivatives) {
     /* Fa l'output della soluzione delle derivate iniziali ed esce */
     pDM->Output();
     Out << "End of derivatives; no simulation is required." << endl;
     return;
   }
   
   /* Dati comuni a passi fittizi e normali */
   integer iStep = 1;
   doublereal dCurrTimeStep = 0.;
   /* First step prediction must always be Crank-Nicholson for accuracy */
   CrankNicholson cn;
   
   if (iFictitiousStepsNumber > 0) {
      /***********************************************************************
       * passi fittizi
       ***********************************************************************/
     
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
     
     DEBUGLCOUT(MYDEBUG_FSTEPS, "Current time step: " << dCurrTimeStep << endl);

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
       cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	    << endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 cout << setw(4) << iTmpCnt << ": "
	      << setw(12) << pXCurr->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrev->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	      << endl;
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
	 cout << "Residual:" << endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   cout << "Dof" << setw(4) << iTmpCnt << ": "
		<< pRes->dGetCoef(iTmpCnt) << endl;
	 }
       }
#endif
       
       if (dTest < dFictitiousStepsTolerance) {
	 goto EndOfFirstFictitiousStep;
       }
       
       iIterCnt++;
       if (iIterCnt > iFictitiousStepsMaxIterations || !isfinite(dTest)) {
	 cerr << endl << "Maximum iterations number " << iIterCnt
	      << " has been reached during first dummy step;"
	      << endl << "time step dt = " << dCurrTimeStep
	      << " cannot be reduced further;" << endl
	      << "aborting ..." << endl;
	 
	 pDM->Output();
	 
	 THROW(SchurMultiStepIntegrator::ErrMaxIterations());
       }

#ifdef MPI_PROFILING
       MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */ 
       pSM->MatrInit(0.);
       pDM->AssJac(*pJac, db0Differential);
#ifdef MPI_PROFILING
        MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */ 
    
       pSM->Solve();
     
#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	 cout << "Solution:" << endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   cout << "Dof" << setw(4) << iTmpCnt << ": "
		<< pSol->dGetCoef(iTmpCnt) << endl;
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
		  << "; current time step: " << dCurrTimeStep << endl);
     
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
	 cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	      << endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   cout << setw(4) << iTmpCnt << ": "
		<< setw(12) << pXCurr->dGetCoef(iTmpCnt)
		<< setw(12) << pXPrev->dGetCoef(iTmpCnt)
		<< setw(12) << pXPrev2->dGetCoef(iTmpCnt)
		<< setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
		<< setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
		<< setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
		<< endl;
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
	   cout << "Residual:" << endl;
	   for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	     cout << "Dof" << setw(4) << iTmpCnt << ": "
		  << pRes->dGetCoef(iTmpCnt) << endl;
	   }
	 }
#endif /* DEBUG */
	 
	 if (dTest < dFictitiousStepsTolerance) {
	   goto EndOfFictitiousStep;
	 }
	 
	 iIterCnt++;
	 if (iIterCnt > iFictitiousStepsMaxIterations || !isfinite(dTest)) {
	   cerr << endl << "Maximum iterations number " << iIterCnt
		<< " has been reached during dummy step "
		<< iSubStep << ';' << endl
		<< "aborting ..." << endl;
	 
	   pDM->Output();
	   
	   THROW(SchurMultiStepIntegrator::ErrMaxIterations());
	 }

#ifdef MPI_PROFILING
         MPE_Log_event(23, 0, "start");
#endif /* MPI_PROFILING */  
	 pSM->MatrInit(0.);
	 pDM->AssJac(*pJac, db0Differential);
#ifdef MPI_PROFILING
	 MPE_Log_event(24, 0, "end");
#endif /* MPI_PROFILING */
	 
	 pSM->Solve();
	 
#ifdef DEBUG
	 if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	   cout << "Solution:" << endl;
	   for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	     cout << "Dof" << setw(4) << iTmpCnt << ": "
		  << pSol->dGetCoef(iTmpCnt) << endl;
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
	     << " error" << endl;
       }
#endif /* DEBUG */
       
       DEBUGLCOUT(MYDEBUG_FSTEPS, "Substep " << iSubStep
		  << " of step " << iStep
		  << " has been completed successfully in "
		  << iIterCnt << " iterations" << endl);
       
       dTime += dRefTimeStep;
     }
     
     Out << "Initial solution after dummy steps at time "
	 << dTime+dCurrTimeStep
	 << " performed in " << iIterCnt
	 << " iterations with " << dTest
	 << " error" << endl;
     
     DEBUGLCOUT(MYDEBUG_FSTEPS,
		"Fictitious steps have been completed successfully in "
		<< iIterCnt << " iterations" << endl);
   } /* Fine dei passi fittizi */


   /* Output delle "condizioni iniziali" */
   pDM->Output();

   Out << "Step " << 0
       << " at time " << dTime+dCurrTimeStep
       << " with time step " << dCurrTimeStep
       << " performed in " << iIterCnt
       << " iterations with " << dTest
       << " error" << endl;
   
   if (fAbortAfterFictitiousSteps) {
     Out << "End of dummy steps; no simulation is required." << endl;
     return;
   }

   iStep = 1; /* Resetto di nuovo iStep */
   
   DEBUGCOUT("Step " << iStep << " has been completed successfully in "
	     << iIterCnt << " iterations" << endl);
   
   
   pDM->BeforePredict(*pXCurr, *pXPrimeCurr, *pXPrev, *pXPrimePrev);
   this->Flip();

   dRefTimeStep = dInitialTimeStep;
   dCurrTimeStep = dRefTimeStep;
   
   DEBUGCOUT("Current time step: " << dCurrTimeStep << endl);
   
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
     cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	  << endl;
     for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
       cout << setw(4) << iTmpCnt << ": "
	    << setw(12) << pXCurr->dGetCoef(iTmpCnt)
	    << setw(12) << pXPrev->dGetCoef(iTmpCnt)
	    << setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	    << setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	    << setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	    << setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	    << endl;
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
       cout << "Residual:" << endl;
       cout << iStep  << "   " << iIterCnt <<endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 cout << "Dof" << setw(4) << iTmpCnt << ": "
	      << pRes->dGetCoef(iTmpCnt) << endl;
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
		 << iIterCnt << " iterations" << endl);
       goto IfFirstStepIsToBeRepeated;
       
     } else {
       cerr << endl << "Maximum iterations number " << iIterCnt
	    << " has been reached during first step;"
	    << endl << "time step dt = " << dCurrTimeStep
	    << " cannot be reduced further;" << endl
	    << "aborting ..." << endl;
       
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
       pSM->MatrInit(0.);
       pDM->AssJac(*pJac, db0Differential);
#ifdef MPI_PROFILING
       MPE_Log_event(24,0,"end Jacobian");
#endif
      }

     pSM->Solve();
     
#ifdef DEBUG
     if (DEBUG_LEVEL(MYDEBUG_SOL)) {
       cout << "Solution:" << endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 cout << "Dof" << setw(4) << iTmpCnt << ": "
	      << pSol->dGetCoef(iTmpCnt) << endl;
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
   
   Out << "Step " << iStep
       << " at time " << dTime+dCurrTimeStep
       << " with time step " << dCurrTimeStep
       << " performed in " << iIterCnt
       << " iterations with " << dTest
       << " error" << endl;
   
   dRefTimeStep = dCurrTimeStep;
   dTime += dRefTimeStep;
   
   dTotErr += dTest;
   iTotIter += iIterCnt;

#ifdef __HACK_EIG__
   /************************
    * TEMPORANEO: FA UN'EIGENANALYSIS */
   if (fEigenAnalysis) {
     if (OneEig.dTime <= dTime && !OneEig.fDone) {
       Eig();
       OneEig.fDone = flag(1);
     }
   }
#endif /* __HACK_EIG__ */
   
   /*************************************************************************
    * Altri passi regolari
    *************************************************************************/
   while (1) {
     MultiStepIntegrationMethod::StepChange CurrStep
       = MultiStepIntegrationMethod::NEWSTEP;
     
     if (dTime >= dFinalTime) {
       cout << "End of simulation at time " << dTime << " after "
	    << iStep << " steps;" << endl
	    << "total iterations: " << iTotIter << endl
	    << "total error: " << dTotErr << endl;
       pDM->MakeRestart();
       return;
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
       cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  ,   XPrime  ,   XPPrev  ,   XPPrev2"
	    << endl;
       for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	 cout << setw(4) << iTmpCnt << ": "
	      << setw(12) << pXCurr->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrev->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrev2->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimeCurr->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimePrev->dGetCoef(iTmpCnt)
	      << setw(12) << pXPrimePrev2->dGetCoef(iTmpCnt)
	      << endl;
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
	 cout << "Residual:" <<endl;
	 cout << iStep  << "   " << iIterCnt <<endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
           cout << "Dof" << setw(4) << iTmpCnt << ": "
		<< pRes->dGetCoef(iTmpCnt) << endl;
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
		     << iIterCnt << " iterations" << endl);
	   goto IfStepIsToBeRepeated;
	   
	 } else {
	   cerr << endl << "Maximum iterations number " << iIterCnt
		<< " has been reached during step " << iStep << ';'
		<< endl << "time step dt = " << dCurrTimeStep
		<< " cannot be reduced further;" << endl
		<< "aborting ..." << endl;
	   
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
	 pSM->MatrInit(0.);
	 pDM->AssJac(*pJac, db0Differential);
#ifdef MPI_PROFILING
	 MPE_Log_event(24,0,"end Jacobian");
#endif
       }
       
       pSM->Solve();
       
#ifdef DEBUG
       if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_MPI)) {
	 cout << "Solution:" << endl;
	 for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	   cout << "Dof" << setw(4) << iTmpCnt << ": "
		<< pSol->dGetCoef(iTmpCnt) << endl;
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
	 << " error" << endl;
     
     DEBUGCOUT("Step " << iStep << " has been completed successfully in "
	       << iIterCnt << " iterations" << endl);
     
     dRefTimeStep = dCurrTimeStep;
     dTime += dRefTimeStep;
     
     
#ifdef __HACK_EIG__
     /************************
      * TEMPORANEO: FA UN'EIGENANALYSIS */
     if (fEigenAnalysis) {
       if (OneEig.dTime <= dTime && !OneEig.fDone) {
	 Eig();
	 OneEig.fDone = flag(1);
       }
     }
#endif /* __HACK_EIG__ */
     
     
     /* Calcola il nuovo timestep */
     dCurrTimeStep = this->NewTimeStep(dCurrTimeStep, iIterCnt, CurrStep);
     DEBUGCOUT("Current time step: " << dCurrTimeStep << endl);
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

   if (pSM != NULL)  {
      SAFEDELETE(pSM);
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
  //DofIterator.fGetFirst(CurrDof);
   
  doublereal dRes = 0.;
  doublereal dXPr = 0.;
  
  /*
   * chiama la routine di comunicazione per la trasmissione del residuo
   * delle interfacce
   */
  pSM->StartExchInt();
 
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

  for (int iCntp1 = 0; iCntp1 < pDM->HowManyDofs(4); iCntp1++) {
    DCount = (pDM->GetDofsList(4))[iCntp1];
    CurrDof = pDofs[DCount-1];
    if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
      doublereal d = XP.dGetCoef(DCount);
      dXPr += d*d;
    }
    /* else if ALGEBRAIC: non aggiunge nulla */
  }
  
  /* verifica completamento trasmissioni */
  pSM->ComplExchInt(dRes, dXPr);

  
  dRes /= (1.+dXPr);

  if (!isfinite(dRes)) {
    cerr << "The simulation diverged; aborting ..." << endl;
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
      cerr << sClassName() << "::FirstStepPredict(): unknown dof order" << endl;
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
      cerr << sClassName() << "::FirstStepPredict(): unknown dof order" << endl;
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
     cerr << __FUNCTION__ << ": ";
#else
     cerr << sClassName() << "::Predict(): ";
#endif /* __GNUC__ */
     cerr << "unknown dof order" << endl;
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
     cerr << __FUNCTION__ << ": ";
#else
     cerr << sClassName() << "::Predict(): ";
#endif /* __GNUC__ */
     cerr << "unknown dof order" << endl;
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
       cerr << "You shouldn't have reached this point!" << endl;
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
   DEBUGCOUTFNAME(sClassName() << "::ReadData");

   /* parole chiave */
   const char* sKeyWords[] = {
      "begin",
    "multistep",
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
    "nostro",
    "ms",
    "hope",
    "rho",
    "algebraic" "rho",
    "derivatives" "coefficient",
    "derivatives" "tolerance",
    "derivatives" "max" "iterations",
    "Newton" "Raphson",
    "true",
    "modified",
    "end",
    "strategy",
    "factor",
    "no" "change",
    "eigen" "analysis",
    "solver",
    "harwell",
    "meschach"
   };

   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
    BEGIN = 0,
    MULTISTEP,
    INITIALTIME,
    FINALTIME,
    TIMESTEP,
    MINTIMESTEP,
    MAXTIMESTEP,
    TOLERANCE,
    MAXITERATIONS,

    /* DEPRECATED */
    FICTITIOUSSTEPSNUMBER,
    FICTITIOUSSTEPSRATIO,
    FICTITIOUSSTEPSTOLERANCE,
    FICTITIOUSSTEPSMAXITERATIONS,
    /* END OF DEPRECATED */

    DUMMYSTEPSNUMBER,
    DUMMYSTEPSRATIO,
    DUMMYSTEPSTOLERANCE,
    DUMMYSTEPSMAXITERATIONS,
    
    ABORTAFTER,
    INPUT,
    ASSEMBLY,
    DERIVATIVES,
    /* DEPRECATED */ FICTITIOUSSTEPS /* END OF DEPRECATED */ ,
    DUMMYSTEPS,
    METHOD,
    /* DEPRECATED */ FICTITIOUSSTEPSMETHOD /* END OF DEPRECATED */ ,
    DUMMYSTEPSMETHOD,
    CRANKNICHOLSON,
    NOSTRO,
    MS,
    HOPE,
    RHO,
    ALGEBRAICRHO,
    DERIVATIVESCOEFFICIENT,
    DERIVATIVESTOLERANCE,
    DERIVATIVESMAXITERATIONS,
    NEWTONRAPHSON,
    NR_TRUE,
    MODIFIED,
    END,
    STRATEGY,
    STRATEGYFACTOR,
    STRATEGYNOCHANGE,
    EIGENANALYSIS,
    SOLVER,
    HARWELL,
    MESCHACH,
    LASTKEYWORD
   };

   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);

   /* cambia la tabella del parser */
   HP.PutKeyTable(K);

   /* legge i dati della simulazione */
   if (KeyWords(HP.GetDescription()) != BEGIN) {
      cerr << endl << "Error: <begin> expected at line "
    << HP.GetLineData() << "; aborting ..." << endl;
      THROW(SchurMultiStepIntegrator::ErrGeneric());
   }

   if (KeyWords(HP.GetWord()) != MULTISTEP) {
      cerr << endl << "Error: <begin: multistep;> expected at line "
    << HP.GetLineData() << "; aborting ..." << endl;
      THROW(SchurMultiStepIntegrator::ErrGeneric());
   }


   flag fMethod(0);
   flag fFictitiousStepsMethod(0);

   /* Ciclo infinito */
   while (1) {
      KeyWords CurrKeyWord = KeyWords(HP.GetDescription());

      switch (CurrKeyWord) {

       case INITIALTIME: {
      dInitialTime = HP.GetReal();
      DEBUGLCOUT(MYDEBUG_INPUT, "Initial time is " << dInitialTime << endl);
      break;
       }

       case FINALTIME: {
      dFinalTime = HP.GetReal();
      DEBUGLCOUT(MYDEBUG_INPUT, "Final time is " << dFinalTime << endl);

      if(dFinalTime <= dInitialTime) {
         cerr << "warning: final time " << dFinalTime
           << " is less than initial time " << dInitialTime
           << ';' << endl
           << "this will cause the simulation to abort" << endl;
      }
      break;
       }

       case TIMESTEP: {
      dInitialTimeStep = HP.GetReal();
      DEBUGLCOUT(MYDEBUG_INPUT, "Initial time step is " << dInitialTimeStep << endl);

      if (dInitialTimeStep == 0.) {
         cerr
           << "warning, null initial time step is not allowed"
           << endl;
      } else if (dInitialTimeStep < 0.) {
         dInitialTimeStep = -dInitialTimeStep;
         cerr
           << "warning, negative initial time step is not allowed;"
           << endl << "its modulus " << dInitialTimeStep
           << " will be considered" << endl;
      }
      break;
       }

       case MINTIMESTEP: {
      dMinimumTimeStep = HP.GetReal();
      DEBUGLCOUT(MYDEBUG_INPUT, "Minimum time step is " << dMinimumTimeStep << endl);

      if (dMinimumTimeStep == 0.) {
         cerr
           << "warning, null minimum time step is not allowed"
           << endl;
         THROW(SchurMultiStepIntegrator::ErrGeneric());
      } else if (dMinimumTimeStep < 0.) {
         dMinimumTimeStep = -dMinimumTimeStep;
         cerr
           << "warning, negative minimum time step is not allowed;"
           << endl << "its modulus " << dMinimumTimeStep
           << " will be considered" << endl;
      }
      break;
       }

       case MAXTIMESTEP: {
      dMaxTimeStep = HP.GetReal();
      DEBUGLCOUT(MYDEBUG_INPUT, "Max time step is " << dMaxTimeStep << endl);

      if (dMaxTimeStep == 0.) {
         cout << "no max time step limit will be considered" << endl;
      } else if (dMaxTimeStep < 0.) {
         dMaxTimeStep = -dMaxTimeStep;
         cerr
           << "warning, negative max time step is not allowed;"
           << endl << "its modulus " << dMaxTimeStep
           << " will be considered" << endl;
      }
      break;
       }

       case FICTITIOUSSTEPSNUMBER:
         cerr << "warning: deprecated keyword \"fictitious steps number\""
	   " at line " << HP.GetLineData() << ";" << endl
	   << "use \"dummy steps number\" instead" << endl;
       case DUMMYSTEPSNUMBER: {
      iFictitiousStepsNumber = HP.GetInt();
      if (iFictitiousStepsNumber < 0) {
         iFictitiousStepsNumber = iDefaultFictitiousStepsNumber;
         cerr << "warning, negative dummy steps number is illegal;" << endl
           << "resorting to default value " << iDefaultFictitiousStepsNumber << endl;
      } else if (iFictitiousStepsNumber == 1) {
         cerr << "warning, a single dummy step may be useless" << endl;
      }

      DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps number: "
             << iFictitiousStepsNumber << endl);
      break;
       }

       case FICTITIOUSSTEPSRATIO:
         cerr << "warning: deprecated keyword \"fictitious steps ratio\""
           " at line " << HP.GetLineData() << ";" << endl
	   << "use \"dummy steps ratio\" instead" << endl;
       case DUMMYSTEPSRATIO: {
      dFictitiousStepsRatio = HP.GetReal();
      if (dFictitiousStepsRatio < 0.) {
         dFictitiousStepsRatio = dDefaultFictitiousStepsRatio;
         cerr << "warning, negative dummy steps ratio is illegal;" << endl
           << "resorting to default value " << dDefaultFictitiousStepsRatio << endl;
      }

      if (dFictitiousStepsRatio > 1.) {
         cerr << "warning, dummy steps ratio is larger than one." << endl
           << "Something like 1.e-3 should be safer ..." << endl;
      }

      DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps ratio: "
             << dFictitiousStepsRatio << endl);
      break;
       }

       case FICTITIOUSSTEPSTOLERANCE:
         cerr << "warning: deprecated keyword \"fictitious steps tolerance\""
	   " at line " << HP.GetLineData() << ";" << endl
	   << "use \"dummy steps tolerance\" instead" << endl;
       case DUMMYSTEPSTOLERANCE: {
      dFictitiousStepsTolerance = HP.GetReal();
      if (dFictitiousStepsTolerance <= 0.) {
         dFictitiousStepsTolerance = dDefaultFictitiousStepsTolerance;
         cerr << "warning, negative dummy steps tolerance is illegal;" << endl
           << "resorting to default value " << dDefaultFictitiousStepsTolerance << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps tolerance: "
             << dFictitiousStepsTolerance << endl);
      break;
       }

       case ABORTAFTER: {
      KeyWords WhenToAbort(KeyWords(HP.GetWord()));
      switch(WhenToAbort) {
       case INPUT: {
          fAbortAfterInput = flag(1);
          DEBUGLCOUT(MYDEBUG_INPUT,
             "Simulation will abort after data input" << endl);
          break;
       }
       case ASSEMBLY: {
          fAbortAfterAssembly = flag(1);
          DEBUGLCOUT(MYDEBUG_INPUT,
             "Simulation will abort after initial assembly" << endl);
          break;
       }

       case DERIVATIVES: {
          fAbortAfterDerivatives = flag(1);
          DEBUGLCOUT(MYDEBUG_INPUT,
             "Simulation will abort after derivatives solution" << endl);
          break;
       }

       case FICTITIOUSSTEPS:
         cerr << "warning: deprecated keyword \"fictitious steps\""
	   " at line " << HP.GetLineData() << ";" << endl
	   << "use \"dummy steps\" instead" << endl;
       case DUMMYSTEPS: {
          fAbortAfterFictitiousSteps = flag(1);
          DEBUGLCOUT(MYDEBUG_INPUT,
             "Simulation will abort after dummy steps solution" << endl);
          break;
       }

       default: {
          cerr << endl
        << "Don't know when to abort, so I'm going to abort now" << endl;
          THROW(SchurMultiStepIntegrator::ErrGeneric());
       }
      }
      break;
       }

       case METHOD: {
      if (fMethod) {
         cerr << "error: multiple definition of integration method at line "
           << HP.GetLineData();
         THROW(SchurMultiStepIntegrator::ErrGeneric());
      }
      fMethod = flag(1);

      KeyWords KMethod = KeyWords(HP.GetWord());
      switch (KMethod) {
       case CRANKNICHOLSON: {
          SAFENEW(pMethod, CrankNicholson); // no constructor
          break;
       }
       case NOSTRO:
       case MS:
       case HOPE: {
          DriveCaller* pRho = ReadDriveData(NULL, HP, NULL);
          HP.PutKeyTable(K);

          DriveCaller* pRhoAlgebraic = ReadDriveData(NULL, HP, NULL);
          HP.PutKeyTable(K);

          switch (KMethod) {
           case NOSTRO:
	   case MS: {
          SAFENEWWITHCONSTRUCTOR(pMethod,
                     NostroMetodo,
                     NostroMetodo(pRho, pRhoAlgebraic));
          break;
           }

           case HOPE: {
          SAFENEWWITHCONSTRUCTOR(pMethod,
                     Hope,
                     Hope(pRho, pRhoAlgebraic));
          break;
           }
	   default:
          THROW(SchurMultiStepIntegrator::ErrGeneric());
          }
          break;
       }
       default: {
          cerr << "Unknown integration method at line " << HP.GetLineData() << endl;
          THROW(SchurMultiStepIntegrator::ErrGeneric());
       }
      }
      break;
       }

       case FICTITIOUSSTEPSMETHOD:
         cerr << "warning: deprecated keyword \"fictitious steps method\""
	   " at line " << HP.GetLineData() << ";" << endl
           << "use \"dummy steps method\" instead" << endl;
       case DUMMYSTEPSMETHOD: {
      if (fFictitiousStepsMethod) {
         cerr << "error: multiple definition of dummy steps"
	   " integration method at line " << HP.GetLineData() << endl;
         THROW(SchurMultiStepIntegrator::ErrGeneric());
      }
      fFictitiousStepsMethod = flag(1);

      KeyWords KMethod = KeyWords(HP.GetWord());
      switch (KMethod) {
       case CRANKNICHOLSON: {
          SAFENEW(pFictitiousStepsMethod, CrankNicholson); // no constructor
          break;
       }
       case NOSTRO:
       case MS:
       case HOPE: {
          DriveCaller* pRho = ReadDriveData(NULL, HP, NULL);
          HP.PutKeyTable(K);

          DriveCaller* pRhoAlgebraic = ReadDriveData(NULL, HP, NULL);
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
          cerr << "Unknown integration method at line " << HP.GetLineData() << endl;
          THROW(SchurMultiStepIntegrator::ErrGeneric());
       }
      }
      break;
       }

       case TOLERANCE: {
      dToll = HP.GetReal();
      if (dToll <= 0.) {
         dToll = dDefaultToll;
         cerr
           << "warning, tolerance <= 0. is illegal; switching to default value "
           << dToll << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dToll << endl);
      break;
       }

       case DERIVATIVESTOLERANCE: {
      dDerivativesToll = HP.GetReal();
      if (dDerivativesToll <= 0.) {
         dDerivativesToll = 1e-6;
         cerr
           << "warning, derivatives tolerance <= 0. is illegal; switching to default value "
           << dDerivativesToll
           << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT,
             "Derivatives toll = " << dDerivativesToll << endl);
      break;
       }

       case MAXITERATIONS: {
      iMaxIterations = HP.GetInt();
      if (iMaxIterations < 1) {
         iMaxIterations = iDefaultMaxIterations;
         cerr
           << "warning, max iterations < 1 is illegal; switching to default value "
           << iMaxIterations
           << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT,
             "Max iterations = " << iMaxIterations << endl);
      break;
       }

       case DERIVATIVESMAXITERATIONS: {
      iDerivativesMaxIterations = HP.GetInt();
      if (iDerivativesMaxIterations < 1) {
         iDerivativesMaxIterations = iDefaultMaxIterations;
         cerr
           << "warning, derivatives max iterations < 1 is illegal; switching to default value "
           << iDerivativesMaxIterations
           << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives max iterations = "
            << iDerivativesMaxIterations << endl);
      break;
       }

       case FICTITIOUSSTEPSMAXITERATIONS:
         cerr << "warning: deprecated keyword"
	   " \"fictitious steps max iterations\""
	   " at line " << HP.GetLineData() << ";" << endl
	   << "use \"dummy steps max iterations\" instead" << endl;
       case DUMMYSTEPSMAXITERATIONS: {
      iFictitiousStepsMaxIterations = HP.GetInt();
      if (iFictitiousStepsMaxIterations < 1) {
         iFictitiousStepsMaxIterations = iDefaultMaxIterations;
         cerr
           << "warning, dummy steps max iterations < 1 is illegal;"
	   " switching to default value " << iFictitiousStepsMaxIterations
           << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "Fictitious steps max iterations = "
            << iFictitiousStepsMaxIterations << endl);
      break;
       }

       case DERIVATIVESCOEFFICIENT: {
      dDerivativesCoef = HP.GetReal();
      if (dDerivativesCoef <= 0.) {
         dDerivativesCoef = 1.;
         cerr
           << "warning, derivatives coefficient <= 0. is illegal; switching to default value "
           << dDerivativesCoef
           << endl;
      }
      DEBUGLCOUT(MYDEBUG_INPUT, "Derivatives coefficient = "
            << dDerivativesCoef << endl);
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
             "Modified Newton-Raphson will be used;" << endl
             << "matrix will be assembled at most after "
             << iIterationsBeforeAssembly
             << " iterations" << endl);
          break;
       }
       default: {
          cerr
        << "warning: unknown case; resorting to default"
        << endl;
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
         cerr << endl
           << "Error: <end: multistep;> expected at line "
           << HP.GetLineData() << "; aborting ..." << endl;
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
         cerr << "warning, illegal reduction factor at line "
           << HP.GetLineData()
           << "; default value 1. (no reduction) will be used"
           << endl;
         StrategyFactor.dReductionFactor = 1.;
          }

          StrategyFactor.iStepsBeforeReduction = HP.GetInt();
          if (StrategyFactor.iStepsBeforeReduction <= 0) {
         cerr << "Warning, illegal number of steps before reduction at line "
           << HP.GetLineData() << ';' << endl
           << "default value 1 will be used (it may be dangerous)"
           << endl;
         StrategyFactor.iStepsBeforeReduction = 1;
          }

          StrategyFactor.dRaiseFactor = HP.GetReal();
          if (StrategyFactor.dRaiseFactor <= 1.) {
         cerr << "warning, illegal raise factor at line "
           << HP.GetLineData()
           << "; default value 1. (no raise) will be used"
           << endl;
         StrategyFactor.dRaiseFactor = 1.;
          }

          StrategyFactor.iStepsBeforeRaise = HP.GetInt();
          if (StrategyFactor.iStepsBeforeRaise <= 0) {
         cerr << "Warning, illegal number of steps before raise at line "
           << HP.GetLineData() << ';' << endl
           << "default value 1 will be used (it may be dangerous)"
           << endl;
         StrategyFactor.iStepsBeforeRaise = 1;
          }

          StrategyFactor.iMinIters = HP.GetInt();
          if (StrategyFactor.iMinIters <= 0) {
         cerr << "Warning, illegal minimum number of iterations at line "
           << HP.GetLineData() << ';' << endl
           << "default value 0 will be used (never raise)"
           << endl;
         StrategyFactor.iMinIters = 1;
          }



          DEBUGLCOUT(MYDEBUG_INPUT,
             "Time step control strategy: Factor" << endl
             << "Reduction factor: "
             << StrategyFactor.dReductionFactor
             << "Steps before reduction: "
             << StrategyFactor.iStepsBeforeReduction
             << "Raise factor: "
             << StrategyFactor.dRaiseFactor
             << "Steps before raise: "
             << StrategyFactor.iStepsBeforeRaise
             << "Min iterations: "
             << StrategyFactor.iMinIters << endl);

          break;
       }

       case STRATEGYNOCHANGE: {
          CurrStrategy = NOCHANGE;
          break;
       }

       default: {
          cerr << "Unknown time step control strategy at line "
        << HP.GetLineData() << endl;
          THROW(SchurMultiStepIntegrator::ErrGeneric());
       }
      }

      break;
       }

#ifdef __HACK_EIG__
       case EIGENANALYSIS: {
      OneEig.dTime = HP.GetReal();
      OneEig.fDone = flag(0);
      fEigenAnalysis = flag(1);

      DEBUGLCOUT(MYDEBUG_INPUT, "Eigenanalysis will be performed at time "
             << OneEig.dTime << endl);
      break;
       }
#endif /* __HACK_EIG__ */

       case SOLVER: {
      switch(KeyWords(HP.GetWord())) {
       case MESCHACH:
#if defined(USE_MESCHACH)
         CurrSolver = MESCHACH_SOLVER;
         DEBUGLCOUT(MYDEBUG_INPUT,
            "Using meschach sparse LU solver" << endl);
         break;
#endif // USE_MESCHACH
       default:
         DEBUGLCOUT(MYDEBUG_INPUT,
            "Unknown solver; switching to default" << endl);
       case HARWELL:
         CurrSolver = HARWELL_SOLVER;
         DEBUGLCOUT(MYDEBUG_INPUT,
            "Using harwell sparse LU solver" << endl);
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
            << ", pivor factor: " << dPivotFactor << endl);
      break;
       }

       default: {
      cerr << endl << "Unknown description at line "
        << HP.GetLineData() << "; aborting ..." << endl;
      THROW(SchurMultiStepIntegrator::ErrGeneric());
       }
      }
   }


   if (dFinalTime < dInitialTime) {
      fAbortAfterAssembly = flag(1);
   }

   if (dFinalTime == dInitialTime) {
      fAbortAfterDerivatives = flag(1);
   }


EndOfCycle:

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
	cerr << "SchurMultiStepIntegrator::Eig() not available" << endl;
}
#endif /* __HACK_EIG__ */

/* SchurMultiStepIntegrator - end */

#endif /* USE_MPI */

