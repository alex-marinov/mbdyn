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

#include <multistp.h>
#include <mynewmem.h>

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

#ifdef HAVE_SIGNAL
volatile sig_atomic_t keep_going = 1;
__sighandler_t sh_term = SIG_DFL;
__sighandler_t sh_int = SIG_DFL;
__sighandler_t sh_hup = SIG_DFL;

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

/* MultiStepIntegrator - begin */

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
MultiStepIntegrator::SolverType
#if defined(USE_Y12)
MultiStepIntegrator::defaultSolver = MultiStepIntegrator::Y12_SOLVER;
#elif /* !USE_Y12 */ defined(USE_UMFPACK3)
MultiStepIntegrator::defaultSolver = MultiStepIntegrator::UMFPACK3_SOLVER;
#elif /* !USE_UMFPACK3 */ defined(USE_HARWELL)
MultiStepIntegrator::defaultSolver = MultiStepIntegrator::HARWELL_SOLVER;
#elif /* !USE_HARWELL */ defined(USE_MESCHACH)
MultiStepIntegrator::defaultSolver = MultiStepIntegrator::MESCHACH_SOLVER;
#else /* !USE_MESCHACH */
#error "need a solver!"
#endif /* !USE_MESCHACH */

/* Costruttore: esegue la simulazione */
MultiStepIntegrator::MultiStepIntegrator(MBDynParser& HPar,
					 const char* sInFName,
					 const char* sOutFName)
:
CurrStrategy(NOCHANGE),
CurrSolver(defaultSolver),
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
pDM(NULL),
DofIterator(), iNumDofs(0),
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
iOutputFlags(0),
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
iWorkSpaceSize(0),
dPivotFactor(-1.)
{
	DEBUGCOUTFNAME("MultiStepIntegrator::MultiStepIntegrator");

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
MultiStepIntegrator::Run(void)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::Run");
   
   	/* chiama il gestore dei dati generali della simulazione */
   	DEBUGLCOUT(MYDEBUG_MEM, "creating DataManager" << std::endl);
   	SAFENEWWITHCONSTRUCTOR(pDM,
			       DataManager,
			       DataManager(HP, 
			       		   dInitialTime, 
					   sInputFileName,
					   sOutputFileName,
					   fAbortAfterInput));

   	/* Si fa dare l'std::ostream al file di output per il log */
   	std::ostream& Out = pDM->GetOutFile();

   	if (fAbortAfterInput) {
      		/* Esce */
		pDM->Output();
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

   	/* Si fa dare il DriveHandler e linka i drivers di rho ecc. */
   	const DriveHandler* pDH = pDM->pGetDrvHdl();
   	pMethod->SetDriveHandler(pDH);
   	pFictitiousStepsMethod->SetDriveHandler(pDH);
   
   	/* Costruisce i vettori della soluzione ai vari passi */
   	DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << std::endl);
   
   	iNumDofs = pDM->iGetNumDofs();
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
			       MyVectorHandler(iNumDofs,
			       		       pdWorkSpace+2*iNumDofs));
   	SAFENEWWITHCONSTRUCTOR(pXPrimePrev,
			       MyVectorHandler,
			       MyVectorHandler(iNumDofs,
			       		       pdWorkSpace+3*iNumDofs));
   	SAFENEWWITHCONSTRUCTOR(pXPrev2, 
			       MyVectorHandler,
			       MyVectorHandler(iNumDofs,
			       		       pdWorkSpace+4*iNumDofs));
   	SAFENEWWITHCONSTRUCTOR(pXPrimePrev2, 
			       MyVectorHandler,
			       MyVectorHandler(iNumDofs,
			       		       pdWorkSpace+5*iNumDofs));
   
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
   	DEBUGLCOUT(MYDEBUG_MEM, "creating SolutionManager, size = "
		   << iNumDofs << std::endl);
   
   	switch (CurrSolver) {
     	case Y12_SOLVER: 
#ifdef USE_Y12
      		SAFENEWWITHCONSTRUCTOR(pSM,
			Y12SparseLUSolutionManager,
			Y12SparseLUSolutionManager(iNumDofs,
				iWorkSpaceSize,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

    	case MESCHACH_SOLVER:
#ifdef USE_MESCHACH
      		SAFENEWWITHCONSTRUCTOR(pSM,
			MeschachSparseLUSolutionManager,
			MeschachSparseLUSolutionManager(iNumDofs,
				iWorkSpaceSize,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_MESCHACH */
      		std::cerr << "Configure with --with-meschach "
			"to enable Meschach solver"
			<< std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

   	case HARWELL_SOLVER:
#ifdef USE_HARWELL
      		SAFENEWWITHCONSTRUCTOR(pSM,
			HarwellSparseLUSolutionManager,
			HarwellSparseLUSolutionManager(iNumDofs,
				iWorkSpaceSize,
				dPivotFactor == -1. ? 1. : dPivotFactor));
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
			Umfpack3SparseLUSolutionManager(iNumDofs, 
				0, dPivotFactor));
      		break;
#else /* !USE_UMFPACK3 */
      		std::cerr << "Configure with --with-umfpack3 "
			"to enable Umfpack3 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_UMFPACK3 */
   	}
   
   	/* Puntatori agli handlers del solution manager */
   	VectorHandler* pRes = pSM->pResHdl();
   	VectorHandler* pSol = pSM->pSolHdl(); 
   	MatrixHandler* pJac = pSM->pMatHdl();

   	/*
	 * Legenda:
	 *   MS - MultiStepIntegrator 
	 *   DM - DataManager
	 *   OM - DofManager
	 *   NM - NodeManager
	 *   EM - ElemManager
	 *   SM - SolutionManager
	 */

   
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
      		pRes->Reset(0.);
      		pDM->AssRes(*pRes, dDerivativesCoef);
      		dTest = MakeTest(*pRes, *pXPrimeCurr);
      
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
			   << std::setw(4) << iIterCnt
			   << ", test = " << dTest << std::endl);
	
      		if (iIterCnt > iDerivativesMaxIterations || !isfinite(dTest)) {
	 		std::cerr << std::endl
				<< "Maximum iterations number " << iIterCnt 
				<< " has been reached during initial"
				" derivatives calculation;" << std::endl
				<< "aborting ..." << std::endl;	 
	 		pDM->Output();	 
	 		THROW(MultiStepIntegrator::ErrMaxIterations());
      		}
      
      		pSM->MatrInit(0.);
      		pDM->AssJac(*pJac, dDerivativesCoef);
      		pSM->Solve();
      
#ifdef DEBUG
      		/* Output della soluzione */
      		if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_DERIVATIVES)) {
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
					<< pSol->dGetCoef(iTmpCnt) << std::endl;
	 		}
      		}
#endif /* DEBUG */
      	
      		Update(*pSol);
      		pDM->DerivativesUpdate();	
   	}
   
EndOfDerivatives:
   
   	dTotErr += dTest;	
   	iTotIter += iIterCnt;
   
   	Out << "Derivatives solution step at time " << dInitialTime
     		<< " performed in " << iIterCnt
     		<< " iterations with " << dTest
     		<< " error" << std::endl;
      
   	DEBUGCOUT("Derivatives solution step has been performed successfully"
		  " in " << iIterCnt << " iterations" << std::endl);
   
   	if (fAbortAfterDerivatives) {
      		/*
		 * Fa l'output della soluzione delle derivate iniziali ed esce
		 */
      		pDM->Output();
      		Out << "End of derivatives; no simulation is required."
			<< std::endl;
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
      
      		/*
		 * inizio integrazione: primo passo a predizione lineare
		 * con sottopassi di correzione delle accelerazioni
		 * e delle reazioni vincolari
		 */
      		pDM->BeforePredict(*pXCurr, *pXPrimeCurr,
				   *pXPrev, *pXPrimePrev);
      		Flip();
         
       		/* primo passo fittizio */
      		/* Passo ridotto per step fittizi di messa a punto */
      		dRefTimeStep = dInitialTimeStep*dFictitiousStepsRatio;
      		dCurrTimeStep = dRefTimeStep;
      		pDM->SetTime(dTime+dCurrTimeStep);    
      
      		DEBUGLCOUT(MYDEBUG_FSTEPS, "Current time step: "
			   << dCurrTimeStep << std::endl);
      
      		/*
		 * First step prediction must always be Crank-Nicholson
		 * for accuracy
		 */
      		cn.SetCoef(dRefTimeStep, 1.,
			   MultiStepIntegrationMethod::NEWSTEP,
			   db0Differential, db0Algebraic);
      		FirstStepPredict(&cn);
      		pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      		pDM->AfterPredict();
      
#ifdef DEBUG
      		if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_PRED)) {
	 		std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  "
				",   XPrime  ,   XPPrev  ,   XPPrev2" 
		   		<< std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	    			std::cout << std::setw(4) << iTmpCnt << ": "
					<< std::setw(12)
					<< pXCurr->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrev->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrev2->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimeCurr->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimePrev->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimePrev2->dGetCoef(iTmpCnt) 
					<< std::endl;
	 		}
      		}
#endif /* DEBUG */
      
      		iIterCnt = 0;   
      		while (1) {
	 		/* l02: EM calcolo del residuo */
	 		pRes->Reset(0.);
			pDM->AssRes(*pRes, db0Differential);
			dTest = MakeTest(*pRes, *pXPrimeCurr);
	 
#ifdef DEBUG
			if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_RESIDUAL)) {
				std::cout << "Residual:" << std::endl;
				for (int iTmpCnt = 1;
				     iTmpCnt <= iNumDofs;
				     iTmpCnt++) {	 
       					std::cout << "Dof"
						<< std::setw(4) << iTmpCnt << ": " 
						<< pRes->dGetCoef(iTmpCnt) << std::endl;
	    			}
	 		}
#endif /* DEBUG */
	 
	 		if (dTest < dFictitiousStepsTolerance) {
	   			goto EndOfFirstFictitiousStep;
	 		}
	 
	 		iIterCnt++;
	 		if (iIterCnt > iFictitiousStepsMaxIterations
			    || !isfinite(dTest)) {
			    	std::cerr << std::endl
					<< "Maximum iterations number "
					<< iIterCnt 
					<< " has been reached during"
					" first dummy step;" << std::endl
					<< "time step dt = " << dCurrTimeStep 
					<< " cannot be reduced further;"
					<< std::endl
					<< "aborting ..." << std::endl;
	    			pDM->Output();
	    			THROW(MultiStepIntegrator::ErrMaxIterations());
	 		}
	 
	 		pSM->MatrInit(0.);
	 		pDM->AssJac(*pJac, db0Differential);
	 		pSM->Solve();
	 
#ifdef DEBUG
	 		if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	    			std::cout << "Solution:" << std::endl;
	    			for (int iTmpCnt = 1;
				     iTmpCnt <= iNumDofs;
				     iTmpCnt++) {
				     	std::cout << "Dof"
						<< std::setw(4) << iTmpCnt << ": "
						<< pSol->dGetCoef(iTmpCnt)
						<< std::endl;
	    			}
	 		}
#endif /* DEBUG */
	 
	 		Update(*pSol); 
	 		pDM->Update();
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
#endif /* DEBUG_FICTITIOUS */
            
       		/* Passi fittizi successivi */
      		for (int iSubStep = 2;
		     iSubStep <= iFictitiousStepsNumber;
		     iSubStep++) {
	 		pDM->BeforePredict(*pXCurr, *pXPrimeCurr,
					   *pXPrev, *pXPrimePrev);
	 		Flip();
	 
	 		DEBUGLCOUT(MYDEBUG_FSTEPS, "Fictitious step "
				   << iSubStep 
				   << "; current time step: " << dCurrTimeStep
				   << std::endl);
	 
	 		pDM->SetTime(dTime+dCurrTimeStep);
	 		ASSERT(pFictitiousStepsMethod != NULL);
	 		pFictitiousStepsMethod->SetCoef(dRefTimeStep, 
				dCurrTimeStep/dRefTimeStep, 
				MultiStepIntegrationMethod::NEWSTEP,
				db0Differential, db0Algebraic);	
	 		Predict(pFictitiousStepsMethod);
	 		pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
	 		pDM->AfterPredict();          
	 
#ifdef DEBUG
	 		if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_PRED)) {
	    			std::cout << "Dof:      XCurr  ,    XPrev  "
					",   XPrev2  ,   XPrime  "
					",   XPPrev  ,   XPPrev2" << std::endl;
	    			for (int iTmpCnt = 1;
				     iTmpCnt <= iNumDofs;
				     iTmpCnt++) {	 
				     	std::cout << std::setw(4) << iTmpCnt << ": " 
						<< std::setw(12)
						<< pXCurr->dGetCoef(iTmpCnt) 
						<< std::setw(12)
						<< pXPrev->dGetCoef(iTmpCnt) 
						<< std::setw(12)
						<< pXPrev2->dGetCoef(iTmpCnt) 
						<< std::setw(12)
						<< pXPrimeCurr->dGetCoef(iTmpCnt) 
						<< std::setw(12)
						<< pXPrimePrev->dGetCoef(iTmpCnt) 
						<< std::setw(12)
						<< pXPrimePrev2->dGetCoef(iTmpCnt) 
						<< std::endl;
	    			}
	 		}
#endif /* DEBUG */
	 
	 		iIterCnt = 0;
	 		while (1) { 
	    			pRes->Reset(0.);
	    			pDM->AssRes(*pRes, db0Differential);
	    			dTest = MakeTest(*pRes, *pXPrimeCurr);
	    
#ifdef DEBUG
	    			if (DEBUG_LEVEL_MATCH(MYDEBUG_FSTEPS|MYDEBUG_RESIDUAL)) {
	       				std::cout << "Residual:" << std::endl;
	       				for (int iTmpCnt = 1;
					     iTmpCnt <= iNumDofs;
					     iTmpCnt++) {	    
		  				std::cout << "Dof"
							<< std::setw(4) << iTmpCnt
							<< ": " 
							<< pRes->dGetCoef(iTmpCnt)
							<< std::endl;
	       				}
	    			}
#endif /* DEBUG */

	    
	    			if (dTest < dFictitiousStepsTolerance) {
	       				goto EndOfFictitiousStep;
	    			}
	    
	    			iIterCnt++;
	    			if (iIterCnt > iFictitiousStepsMaxIterations
				    || !isfinite(dTest) ) {
				    	std::cerr << std::endl
						<< "Maximum iterations number "
						<< iIterCnt 
						<< " has been reached"
						" during dummy step " 
						<< iSubStep << ';' << std::endl
						<< "aborting ..." << std::endl;
						
	       				pDM->Output();
					THROW(MultiStepIntegrator::ErrMaxIterations());
	    			}
	    
				if (iOutputFlags & OUTPUT_ITERS) {
					Out << "\tIteration " << iIterCnt
						<< " " << dTest << " J"
						<< std::endl;
				}

	    			pSM->MatrInit(0.);
	    			pDM->AssJac(*pJac, db0Differential);
	    			pSM->Solve();
	    
#ifdef DEBUG
	    			if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_FSTEPS)) {
	       				std::cout << "Solution:" << std::endl;
	       				for (int iTmpCnt = 1;
					     iTmpCnt <= iNumDofs;
					     iTmpCnt++) {
					     	std::cout << "Dof"
							<< std::setw(4) << iTmpCnt
							<< ": " 
							<< pSol->dGetCoef(iTmpCnt)
							<< std::endl;
	       				}
	    			}
#endif /* DEBUG */
	    
	    			Update(*pSol); 
	    			pDM->Update();
	 		}
	 
EndOfFictitiousStep:
	 
	 		dTotErr += dTest;	
	 		iTotIter += iIterCnt;
	 
#ifdef DEBUG
	 		if (DEBUG_LEVEL(MYDEBUG_FSTEPS)) {
	    			pDM->Output();
				
	    			Out << "Step " << iStep 
					<< " time " << dTime+dCurrTimeStep
					<< " step " << dCurrTimeStep
					<< " iterations " << iIterCnt
					<< " error " << dTest << std::endl;
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
      
      		Out << "Initial solution after dummy steps at time " << dTime
			<< " performed in " << iIterCnt
			<< " iterations with " << dTest 
			<< " error" << std::endl;
			
      		DEBUGLCOUT(MYDEBUG_FSTEPS, 
			   "Fictitious steps have been completed successfully"
			   " in " << iIterCnt << " iterations" << std::endl);
   	} /* Fine dei passi fittizi */
   
   
   	/* Output delle "condizioni iniziali" */
   	pDM->Output();
   
   	Out << "Step " << 0
     		<< " time " << dTime+dCurrTimeStep
     		<< " step " << dCurrTimeStep
     		<< " iterations " << iIterCnt
     		<< " error " << dTest << std::endl;
   
   	if (fAbortAfterFictitiousSteps) {
      		Out << "End of dummy steps; no simulation is required."
			<< std::endl;
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
   	Flip();
   
   	dRefTimeStep = dInitialTimeStep;   
   	dCurrTimeStep = dRefTimeStep;
   
   	DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   
    	/* Primo passo regolare */
   
IfFirstStepIsToBeRepeated:
   	pDM->SetTime(dTime+dCurrTimeStep);

    	/*
	 * metto rho = 1 perche' cosi' rispetto certi teoremi
	 * sulla precisione di questo passo (Petzold, 89)
	 * (ricado nella regola dei trapezi)
	 */
   	cn.SetCoef(dRefTimeStep, 1.,
		   MultiStepIntegrationMethod::NEWSTEP,
		   db0Differential, db0Algebraic);
   	FirstStepPredict(&cn);
   	pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
   	pDM->AfterPredict();    

#ifdef DEBUG
   	if (DEBUG_LEVEL(MYDEBUG_PRED)) {
      		std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  "
			",   XPrime  ,   XPPrev  ,   XPPrev2" << std::endl;
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
      		pRes->Reset(0.);
      		pDM->AssRes(*pRes, db0Differential);      
      		dTest = MakeTest(*pRes, *pXPrimeCurr);
      
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
				dCurrTimeStep =
					NewTimeStep(dCurrTimeStep, 
							  iIterCnt, 
							  MultiStepIntegrationMethod::REPEATSTEP);
				dRefTimeStep = dCurrTimeStep;
				DEBUGCOUT("Changing time step during"
					  " first step after "
					  << iIterCnt << " iterations"
					  << std::endl);
	    			goto IfFirstStepIsToBeRepeated;
	 		} else {
	    			std::cerr << std::endl
					<< "Maximum iterations number "
					<< iIterCnt
					<< " has been reached during"
					" first step;"
					<< std::endl
					<< "time step dt = " << dCurrTimeStep 
					<< " cannot be reduced further;"
					<< std::endl
					<< "aborting ..." << std::endl;
	    			pDM->Output();
				THROW(MultiStepIntegrator::ErrMaxIterations());
	 		}
      		}
      
      		/* Modified Newton-Raphson ... */
      		if (iPerformedIterations < iIterationsBeforeAssembly) {
	 		iPerformedIterations++;
			
			if (iOutputFlags & OUTPUT_ITERS) {
				Out << "\tIteration " << iIterCnt
					<< " " << dTest << " M"
					<< std::endl;
			}

      		} else {
	 		iPerformedIterations = 0;
	 		pSM->MatrInit(0.);
	 		pDM->AssJac(*pJac, db0Differential);

			if (iOutputFlags & OUTPUT_ITERS) {
				Out << "\tIteration " << iIterCnt
					<< " " << dTest << " J"
					<< std::endl;
			}

      		}
      		pSM->Solve();   
      
#ifdef DEBUG
      		if (DEBUG_LEVEL(MYDEBUG_SOL)) {      
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
					<< pSol->dGetCoef(iTmpCnt) << std::endl;
			}
		}
#endif /* DEBUG */


      		Update(*pSol); 
      		pDM->Update();
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
			<< " time " << dTime+dCurrTimeStep
			<< " step " << dCurrTimeStep
			<< " iterations " << iIterCnt
			<< " error " << dTest << std::endl;
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
   
    	/* Altri passi regolari */
   	while (1) {
      		MultiStepIntegrationMethod::StepChange CurrStep 
			= MultiStepIntegrationMethod::NEWSTEP;
      
      		if (dTime >= dFinalTime) {
	 		std::cout << "End of simulation at time "
				<< dTime << " after " 
				<< iStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total error: " << dTotErr << std::endl;
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
	   
      		pDM->BeforePredict(*pXCurr, *pXPrimeCurr,
				   *pXPrev, *pXPrimePrev);
      		Flip();
      		
IfStepIsToBeRepeated:
      		pDM->SetTime(dTime+dCurrTimeStep);
      		pMethod->SetCoef(dRefTimeStep, 
				 dCurrTimeStep/dRefTimeStep, 
				 CurrStep,
				 db0Differential, 
				 db0Algebraic);	
      		Predict(pMethod);
      		pDM->LinkToSolution(*pXCurr, *pXPrimeCurr);
      		pDM->AfterPredict();        
   
#ifdef DEBUG
      		if (DEBUG_LEVEL(MYDEBUG_PRED)) {
			std::cout << "Dof:      XCurr  ,    XPrev  ,   XPrev2  "
				",   XPrime  ,   XPPrev  ,   XPPrev2" << std::endl;
			for (int iTmpCnt = 1; iTmpCnt <= iNumDofs; iTmpCnt++) {
	    			std::cout << std::setw(4) << iTmpCnt << ": " 
					<< std::setw(12)
					<< pXCurr->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrev->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrev2->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimeCurr->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimePrev->dGetCoef(iTmpCnt) 
					<< std::setw(12)
					<< pXPrimePrev2->dGetCoef(iTmpCnt) 
					<< std::endl;
	 		}
      		}
#endif /* DEBUG */

      		iIterCnt = 0;
      		iPerformedIterations = iIterationsBeforeAssembly;
      		while (1) {
	 		pRes->Reset(0.);
	 		pDM->AssRes(*pRes, db0Differential);
#ifdef USE_EXCEPTIONS
			try {
#endif /* USE_EXCEPTIONS */
	 			dTest = MakeTest(*pRes, *pXPrimeCurr);
#ifdef USE_EXCEPTIONS
			} catch (MultiStepIntegrator::ErrSimulationDiverged) {
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
	    			std::cout << "Residual:" << std::endl;
	    			std::cout << iStep  << "   " << iIterCnt << std::endl;
	    			for (int iTmpCnt = 1;
				     iTmpCnt <= iNumDofs;
				     iTmpCnt++) {
	       				std::cout << "Dof"
						<< std::setw(4) << iTmpCnt << ": " 
						<< pRes->dGetCoef(iTmpCnt)
						<< std::endl;
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
					dCurrTimeStep =
						NewTimeStep(dCurrTimeStep,
								  iIterCnt,
								  CurrStep);	       
					DEBUGCOUT("Changing time step"
						  " during step " 
						  << iStep << " after "
						  << iIterCnt << " iterations"
						  << std::endl);
					goto IfStepIsToBeRepeated;
	    			} else {
					std::cerr << std::endl
						<< "Maximum iterations number "
						<< iIterCnt 
						<< " has been reached during"
						" step " << iStep << ';'
						<< std::endl
						<< "time step dt = "
						<< dCurrTimeStep
						<< " cannot be reduced"
						" further;" << std::endl
						<< "aborting ..." << std::endl;
#if 0
					std::cerr << *pJac << std::endl;
#endif /* 0 */
					
	       				THROW(MultiStepIntegrator::ErrMaxIterations());
		    		}
	 		}

			/* Modified Newton-Raphson ... */
	 		if (iPerformedIterations < iIterationsBeforeAssembly) {
	    			iPerformedIterations++;

				if (iOutputFlags & OUTPUT_ITERS) {
					Out << "\tIteration " << iIterCnt
						<< " " << dTest << " M"
						<< std::endl;
				}

	 		} else {
	    			iPerformedIterations = 0;
	    			pSM->MatrInit(0.);
	    			pDM->AssJac(*pJac, db0Differential);

				if (iOutputFlags & OUTPUT_ITERS) {
					Out << "\tIteration " << iIterCnt
						<< " " << dTest << " J"
						<< std::endl;
				}

	 		}
	 		pSM->Solve();

#ifdef DEBUG
	 		if (DEBUG_LEVEL_MATCH(MYDEBUG_SOL|MYDEBUG_MPI)) {
	    			std::cout << "Solution:" << std::endl;
	    			for (int iTmpCnt = 1;
				     iTmpCnt <= iNumDofs;
				     iTmpCnt++) {	    
	       				std::cout << "Dof"
						<< std::setw(4) << iTmpCnt << ": " 
						<< pSol->dGetCoef(iTmpCnt)
						<< std::endl;
	    			}
	 		}
#endif /* DEBUG */

	 		Update(*pSol);
	 		pDM->Update();
      		}
	 
EndOfStep:

      		dTotErr += dTest;	
      		iTotIter += iIterCnt;
      
      		pDM->Output();     
	 
      		Out << "Step " << iStep 
			<< " time " << dTime+dCurrTimeStep
			<< " step " << dCurrTimeStep
			<< " iterations " << iIterCnt
			<< " error " << dTest << std::endl;
      
      		DEBUGCOUT("Step " << iStep
			  << " has been completed successfully in "
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
      		dCurrTimeStep =
			NewTimeStep(dCurrTimeStep, iIterCnt, CurrStep);
		DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
   	}
}

/* Distruttore */
MultiStepIntegrator::~MultiStepIntegrator(void)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::~MultiStepIntegrator");

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
MultiStepIntegrator::MakeTest(const VectorHandler& Res, 
			      const VectorHandler& XP)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::MakeTest");
   
   	Dof CurrDof;
   	DofIterator.fGetFirst(CurrDof);
      
   	doublereal dRes = 0.;
   	doublereal dXPr = 0.;
   
   	for (int iCntp1 = 1;
	     iCntp1 <= iNumDofs; 
	     iCntp1++, DofIterator.fGetNext(CurrDof)) {
      		doublereal d = Res.dGetCoef(iCntp1);
      		dRes += d*d;
      		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
	 		d = XP.dGetCoef(iCntp1);
	 		dXPr += d*d;
      		}
      		/* else if ALGEBRAIC: non aggiunge nulla */
   	}
   	dRes /= (1.+dXPr);
   	if (!isfinite(dRes)) {      
      		std::cerr << "The simulation diverged; aborting ..." << std::endl;       
      		THROW(MultiStepIntegrator::ErrSimulationDiverged());
   	}
   	return sqrt(dRes);
}

/* Predizione al primo passo */
void 
MultiStepIntegrator::FirstStepPredict(MultiStepIntegrationMethod* pM)
{      
   	DEBUGCOUTFNAME("MultiStepIntegrator::FirstStepPredict");
   
   	Dof CurrDof;
   	DofIterator.fGetFirst(CurrDof);
   
   	/* Combinazione lineare di stato e derivata al passo precedente ecc. */
   	for (int iCntp1 = 1;
	     iCntp1 <= iNumDofs;
	     iCntp1++, DofIterator.fGetNext(CurrDof)) {
      		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXPnm1 = pXPrimePrev->dGetCoef(iCntp1);

	 		doublereal dXPn = pM->dPredDer(dXnm1, 0., dXPnm1, 0.);
			doublereal dXn = pM->dPredState(dXnm1, 0.,
							dXPn, dXPnm1, 0.);
			
	 		pXPrimeCurr->fPutCoef(iCntp1, dXPn);	 
	 		pXCurr->fPutCoef(iCntp1, dXn);
      		} else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	
	 		doublereal dXn = pM->dPredDerAlg(0., dXnm1, 0.);
			doublereal dXIn = pM->dPredStateAlg(0., dXn,
							    dXnm1, 0.);
			
	 		pXCurr->fPutCoef(iCntp1, dXn);
	 		pXPrimeCurr->fPutCoef(iCntp1, dXIn);

      		} else {
	 		std::cerr << "MultiStepIntegrator::FirstStepPredict():"
				" unknown dof order" << std::endl;
	 		THROW(ErrGeneric());
      		}
   	}
}

/* Predizione al passo generico */
void 
MultiStepIntegrator::Predict(MultiStepIntegrationMethod* pM)
{   
   	/* Note: pM must be initialised prior to calling Predict() */
   
   	DEBUGCOUTFNAME("MultiStepIntegrator::Predict");
   
   	Dof CurrDof;
   	DofIterator.fGetFirst(CurrDof);
   
   	/* Combinazione lineare di stato e derivata al passo precedente ecc. */
   	for (int iCntp1 = 1;
	     iCntp1 <= iNumDofs;
	     iCntp1++, DofIterator.fGetNext(CurrDof)) {
      		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXnm2 = pXPrev2->dGetCoef(iCntp1);
	 		doublereal dXPnm1 = pXPrimePrev->dGetCoef(iCntp1);
	 		doublereal dXPnm2 = pXPrimePrev2->dGetCoef(iCntp1);

	 		doublereal dXPn = pM->dPredDer(dXnm1, dXnm2,
						       dXPnm1, dXPnm2);
			doublereal dXn = pM->dPredState(dXnm1, dXnm2, 
							dXPn, dXPnm1, dXPnm2);
			
	 		pXPrimeCurr->fPutCoef(iCntp1, dXPn);
	 		pXCurr->fPutCoef(iCntp1, dXn);
			
      		} else if (CurrDof.Order == DofOrder::ALGEBRAIC) {
	 		doublereal dXnm1 = pXPrev->dGetCoef(iCntp1);
	 		doublereal dXnm2 = pXPrev2->dGetCoef(iCntp1);
	 		doublereal dXInm1 = pXPrimePrev->dGetCoef(iCntp1);

	 		doublereal dXn = pM->dPredDerAlg(dXInm1, dXnm1, dXnm2);
			doublereal dXIn = pM->dPredStateAlg(dXInm1, dXn,
							    dXnm1, dXnm2);
			
	 		pXCurr->fPutCoef(iCntp1, dXn);
	 		pXPrimeCurr->fPutCoef(iCntp1, dXIn);

      		} else {
	 		std::cerr << "unknown dof order" << std::endl;
	 		THROW(ErrGeneric());
      		}
   	}
}

/* Nuovo delta t */
doublereal
MultiStepIntegrator::NewTimeStep(doublereal dCurrTimeStep,
				 integer iPerformedIters,
				 MultiStepIntegrationMethod::StepChange Why)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::NewTimeStep");

   	switch (CurrStrategy) {
    	case FACTOR:
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
					/*
					 * Fuori viene intercettato
					 * il valore illegale
					 */
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
      
    	case NOCHANGE:
       		return dCurrTimeStep;
      
    	default:
       		std::cerr << "You shouldn't have reached this point!" << std::endl;
       		THROW(MultiStepIntegrator::ErrGeneric());
   	}
   
   	return dCurrTimeStep;
}

/* Aggiornamento della soluzione nel passo fittizio */
void 
MultiStepIntegrator::DerivativesUpdate(const VectorHandler& Sol)
{
   	DEBUGCOUTFNAME("MultiStepIntegrator::DerivativesUpdate");
   
   	Dof CurrDof;
   	DofIterator.fGetFirst(CurrDof);

   	for (int iCntp1 = 1;
	     iCntp1 <= iNumDofs;
	     iCntp1++, DofIterator.fGetNext(CurrDof)) {		
      		doublereal d = Sol.dGetCoef(iCntp1);
      		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
	 		pXPrimeCurr->fIncCoef(iCntp1, d);
      		} else {
	 		pXCurr->fIncCoef(iCntp1, d);
      		}
   	}
}

/* Aggiornamento normale */
void 
MultiStepIntegrator::Update(const VectorHandler& Sol)
{   
   	DEBUGCOUTFNAME("MultiStepIntegrator::Update");

   	Dof CurrDof;
   	DofIterator.fGetFirst(CurrDof);
 
   	for (int iCntp1 = 1;
	     iCntp1 <= iNumDofs;
	     iCntp1++, DofIterator.fGetNext(CurrDof)) {		
	     	doublereal d = Sol.dGetCoef(iCntp1);
		if (CurrDof.Order == DofOrder::DIFFERENTIAL) {
			pXPrimeCurr->fIncCoef(iCntp1, d);
			/*
			 * Nota: b0Differential e b0Algebraic
			 * possono essere distinti;
			 * in ogni caso sono calcolati dalle funzioni
			 * di predizione e sono dati globali
			 */
	 		pXCurr->fIncCoef(iCntp1, db0Differential*d);
      		} else {
	 		pXCurr->fIncCoef(iCntp1, d);
	 		pXPrimeCurr->fIncCoef(iCntp1, db0Algebraic*d);
      		}
   	}
}

/* Dati dell'integratore */
void 
MultiStepIntegrator::ReadData(MBDynParser& HP)
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
			"iterations",
	
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

		OUTPUT,
		ITERATIONS,
	
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
      		THROW(MultiStepIntegrator::ErrGeneric());
   	}
   
   	if (KeyWords(HP.GetWord()) != MULTISTEP) {
      		std::cerr << std::endl << "Error: <begin: multistep;> expected at line " 
			<< HP.GetLineData() << "; aborting ..." << std::endl;
      		THROW(MultiStepIntegrator::ErrGeneric());
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
	     			THROW(MultiStepIntegrator::ErrGeneric());
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
	      			THROW(MultiStepIntegrator::ErrGeneric());
	  		}
	  		break;
       		}

		case OUTPUT: {
			while (HP.fIsArg()) {
				KeyWords OutputFlag(KeyWords(HP.GetWord()));
				switch (OutputFlag) {
				case ITERATIONS:
					iOutputFlags |= OUTPUT_ITERS;
					break;

				default:
					std::cerr << "Unknown output flag "
						"at line " << HP.GetLineData() 
						<< "; ignored" << std::endl;
					break;
				}
			}
			
			break;
		}
	 
       		case METHOD: {
	  		if (fMethod) {
	     			std::cerr << "error: multiple definition"
					" of integration method at line "
					<< HP.GetLineData() << std::endl;
	     			THROW(MultiStepIntegrator::ErrGeneric());
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
				THROW(MultiStepIntegrator::ErrGeneric());
	  		}
	  		break;
       		}

       case FICTITIOUSSTEPSMETHOD:
       case DUMMYSTEPSMETHOD: {
	  if (fFictitiousStepsMethod) {
	     std::cerr << "error: multiple definition of dummy steps integration"
	       " method at line "
	       << HP.GetLineData();
	     THROW(MultiStepIntegrator::ErrGeneric());
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
	      THROW(MultiStepIntegrator::ErrGeneric());
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
	     THROW(MultiStepIntegrator::ErrGeneric());
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
	      THROW(MultiStepIntegrator::ErrGeneric());
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

	   case HARWELL: 
	     CurrSolver = HARWELL_SOLVER;
	     DEBUGLCOUT(MYDEBUG_INPUT, 
			"Using harwell sparse LU solver" << std::endl);	 
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
	 	 
       default: {
	  std::cerr << std::endl << "Unknown description at line " 
	    << HP.GetLineData() << "; aborting ..." << std::endl;
	  THROW(MultiStepIntegrator::ErrGeneric());
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
   
   
/* Estrazione autovalori, vincolata alla disponibilita' delle LAPACK */
#ifdef __HACK_EIG__

#include <../libraries/libobjs/dgegv.h>

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
MultiStepIntegrator::Eig(void)
{
   DEBUGCOUTFNAME("MultiStepIntegrator::Eig");   

   /*
    * MatA, MatB: MatrixHandlers to eigenanalysis matrices
    * MatL, MatR: MatrixHandlers to eigenvectors, if required
    * AlphaR, AlphaI Beta: eigenvalues
    * WorkVec:    Workspace
    * iWorkSize:  Size of the workspace
    */
   
   DEBUGCOUT("MultiStepIntegrator::Eig(): performing eigenanalysis" << std::endl);
   
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
      /* < 0:  if INFO = -i, the i-th argument had an illegal value. */
      Out << "the " << -iInfo << "-th argument had illegal value" << std::endl;
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
   if (fOutputModes) {
	   /* crea il file .pch */
	   pch.open("mbdyn.bdf");
	   pch.setf(std::ios::showpoint);
	   pch << "$.......2.......3.......4.......5.......6.......7.......8.......9.......0......." << std::endl;
	   pch << "MAT1           1    1.E0    1.E0            1.E0" << std::endl;
	   pDM->Output_pch(pch);
	   
	   /* crea il file .f06 */
	   f06.open("mbdyn.f06");
	   f06.setf(std::ios::showpoint);
	   f06.setf(std::ios::scientific);
   }
#endif /* __HACK_NASTRAN_MODES__ */

   for (int iCnt = 1; iCnt <= iSize; iCnt++) {

#ifdef __HACK_NASTRAN_MODES__
      /* EXPERIMENTAL */
      flag doPlot = 0;
      if (fOutputModes) {

         doublereal b = Beta.dGetCoef(iCnt);
	 doublereal re = AlphaR.dGetCoef(iCnt);
	 doublereal im = AlphaI.dGetCoef(iCnt);
	 doublereal sigma;
	 doublereal omega;
	 doublereal csi;
	 doublereal freq;

	 do_eig(b, re, im, h, sigma, omega, csi, freq);

	 if (freq >= dLowerFreq && freq <= dUpperFreq) {
	    doPlot = 1;
	      
	    f06 
	      << "                                                                                                 CSA/NASTRAN 11/14/95    PAGE   " 
	      << std::setw(4) << iCnt << std::endl
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
	    f06 << "   SUBCASE " << iCnt << std::endl
	      << std::endl
	      << "                                            D I S P L A C E M E N T  V E C T O R" << std::endl
	      << std::endl
	      << "     POINT ID.   TYPE          T1             T2             T3             R1             R2             R3" << std::endl;
	 }
      }
#endif /* __HACK_NASTRAN_MODES__ */
      
      doublereal cmplx = AlphaI.dGetCoef(iCnt);
      if (cmplx == 0.) {
         Out << "Mode " << iCnt << ":" << std::endl;
         for (int jCnt = 1; jCnt <= iSize; jCnt++) {
            Out << std::setw(12) << jCnt << ": "
	      << std::setw(12) << MatR.dGetCoef(jCnt, iCnt) << std::endl;
	 }
	 
	 if (fOutputModes) {
	    /*
	     * per ora sono uguali; in realta' XP e' X * lambda
	     */
	    MyVectorHandler X(iSize, MatR.pdGetMat()+iSize*(iCnt-1));
	    MyVectorHandler XP(iSize, MatR.pdGetMat()+iSize*(iCnt-1));
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
            Out << "Modes " << iCnt << ", " << iCnt+1 << ":" << std::endl;
	    for (int jCnt = 1; jCnt <= iSize; jCnt++) {
	       doublereal im = MatR.dGetCoef(jCnt, iCnt+1);
	       Out << std::setw(12) << jCnt << ": "
	         << std::setw(12) << MatR.dGetCoef(jCnt, iCnt) 
	         << ( im >= 0. ? " + " : " - " ) 
	         << std::setw(12) << fabs(im) << " * j " << std::endl;
            }
	 }

	 if (fOutputModes) {
            /*
	     * uso la parte immaginaria ...
	     */
	    int i = iCnt - (cmplx > 0. ? 0 : 1);
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
   pch.close();
   f06.close();
#endif /* __HACK_NASTRAN_MODES__ */      

   /* Non puo' arrivare qui se le due aree di lavoro non sono definite */
   SAFEDELETEARR(pd);
   SAFEDELETEARR(ppd);
}

#endif /* __HACK_EIG__ */

/* MultiStepIntegrator - end */

