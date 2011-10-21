/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
 * Copyright 1999-2011 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * This copyright statement applies to the MPI related code, which was
 * merged from files schur.h/schur.cc
 */

/*
 *
 * Copyright (C) 2008
 * Alessandro Fumagalli <alessandro.fumagalli@polimi.it>
 *
 */

/* metodo per la soluzione del modello */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

/* required for configure time macros with paths */
#include "mbdefs.h"

#define RTAI_LOG

#include <cstring>
#include <limits>
#include <unistd.h>
#include <cerrno>

#include <cfloat>
#include <cmath>
#include "ac/sys_sysinfo.h"

#include "invsolver.h"
#include "dataman.h"
#include "mtdataman.h"
#include "thirdorderstepsol.h"
#include "nr.h"
#include "bicg.h"
#include "gmres.h"
#include "solman.h"
#include "stepsol.h"
#include <vector>
#include "readlinsol.h"

#include "solver_impl.h"

InverseSolver::InverseSolver(MBDynParser& HPar,
		const std::string& sInFName,
		const std::string& sOutFName,
		bool bPar)
: Solver(HPar, sInFName, sOutFName, bPar),
ProblemType(InverseDynamics::FULLY_ACTUATED_COLLOCATED),
pXPrimePrime(NULL), pLambda(NULL)
{
	DEBUGCOUTFNAME("InverseSolver::InverseSolver");
}

void
InverseSolver::Run(void)
{
	DEBUGCOUTFNAME("InverseSolver::Run");

	mbdyn_signal_init(1);

	/* Legge i dati relativi al metodo di integrazione */
	ReadData(HP);

/*FIXME:*/
//	bParallel = false;

#ifdef USE_MULTITHREAD
	/* check for thread potential */
	if (nThreads == 0) {
		int n = get_nprocs();

		if (n > 1) {
			silent_cout("no multithread requested "
					"with a potential of " << n
					<< " CPUs" << std::endl);
			nThreads = n;

		} else {
			nThreads = 1;
		}
	}
#endif /* USE_MULTITHREAD */

	if (pRTSolver) {
		pRTSolver->Setup();
	}

#ifdef USE_SCHUR
	int mpi_finalize = 0;

	if (bParallel) {
		DEBUGLCOUT(MYDEBUG_MEM, "creating parallel SchurDataManager"
				<< std::endl);

		SAFENEWWITHCONSTRUCTOR(pSDM,
			SchurDataManager,
			SchurDataManager(HP,
				OutputFlags,
				this,
				dInitialTime,
				sOutputFileName.c_str(),
				sInputFileName.c_str(),
				eAbortAfter == AFTER_INPUT));

		/* FIXME: who frees sNewOutname? */

		pDM = pSDM;

	} else
#endif // USE_SCHUR
	{
		/* chiama il gestore dei dati generali della simulazione */
#ifdef USE_MULTITHREAD
		if (nThreads > 1) {
			if (!(CurrLinearSolver.GetSolverFlags() & LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS)) {
				/* conservative: dir may use too much memory */
				if (!CurrLinearSolver.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS)) {
					bool b;

#if defined(USE_UMFPACK)
					b = CurrLinearSolver.SetSolver(LinSol::UMFPACK_SOLVER,
							LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS);
#elif defined(USE_Y12)
					b = CurrLinearSolver.SetSolver(LinSol::Y12_SOLVER,
							LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS);
#else
					b = false;
#endif
					if (!b) {
						silent_cerr("unable to select a CC-capable solver"
									<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}

			silent_cout("Creating multithread solver "
					"with " << nThreads << " threads "
					"and "
					<< CurrLinearSolver.GetSolverName()
					<< " linear solver"
					<< std::endl);

			SAFENEWWITHCONSTRUCTOR(pDM,
					MultiThreadDataManager,
					MultiThreadDataManager(HP,
						OutputFlags,
						this,
						dInitialTime,
						sOutputFileName.c_str(),
						sInputFileName.c_str(),
						eAbortAfter == AFTER_INPUT,
						nThreads));

		} else
#endif /* USE_MULTITHREAD */
		{
			DEBUGLCOUT(MYDEBUG_MEM, "creating DataManager"
					<< std::endl);

			silent_cout("Creating scalar solver "
					"with "
					<< CurrLinearSolver.GetSolverName()
					<< " linear solver"
					<< std::endl);

			SAFENEWWITHCONSTRUCTOR(pDM,
					DataManager,
					DataManager(HP,
						OutputFlags,
						this,
						dInitialTime,
						sOutputFileName.c_str(),
						sInputFileName.c_str(),
						eAbortAfter == AFTER_INPUT));
		}
	}

	{ // log of symbol table
		std::ostream& out = pDM->GetLogFile();
		out << HP.GetMathParser().GetSymbolTable();
	}
	HP.Close();

	/* Si fa dare l'std::ostream al file di output per il log */
	std::ostream& Out = pDM->GetOutFile();

	/* Qui crea le partizioni: principale fra i processi, se parallelo  */
#ifdef USE_SCHUR
	if (bParallel) {
		pSDM->CreatePartition();
	}
#endif // USE_SCHUR

	const DriveHandler* pDH = pDM->pGetDrvHdl();
	pRegularSteps->SetDriveHandler(pDH);
	/* Costruisce i vettori della soluzione ai vari passi */
	DEBUGLCOUT(MYDEBUG_MEM, "creating solution vectors" << std::endl);

#ifdef USE_SCHUR
	if (bParallel) {
		iNumDofs = pSDM->HowManyDofs(SchurDataManager::TOTAL);
		pDofs = pSDM->pGetDofsList();

		iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
		pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);

		iNumIntDofs = pSDM->HowManyDofs(SchurDataManager::INTERNAL);
		pIntDofs = pSDM->GetDofsList(SchurDataManager::INTERNAL);

	} else
#endif // USE_SCHUR
	{
		iNumDofs = pDM->iGetNumDofs();
	}

	/* relink those known drive callers that might need
	 * the data manager, but were verated ahead of it */
	if (pStrategyChangeDrive) {
		pStrategyChangeDrive->SetDrvHdl(pDM->pGetDrvHdl());
	}

	ASSERT(iNumDofs > 0);

#if 0
	/* sono i passi precedenti usati dall'integratore */
	integer iRSteps = pRegularSteps->GetIntegratorNumPreviousStates();
	integer iRUnkStates = pRegularSteps->GetIntegratorNumUnknownStates();
#endif // 0

	/* FIXME: pdWorkspace?*/

#if 1
	/* allocate workspace for previous time steps */
	SAFENEWARR(pdWorkSpace, doublereal,
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
#endif

	SAFENEWWITHCONSTRUCTOR(pX,
		MyVectorHandler,
		MyVectorHandler(iUnkStates*iNumDofs));
	SAFENEWWITHCONSTRUCTOR(pXPrime,
		MyVectorHandler,
		MyVectorHandler(iUnkStates*iNumDofs));
	SAFENEWWITHCONSTRUCTOR(pXPrimePrime,
		MyVectorHandler,
		MyVectorHandler(iUnkStates*iNumDofs));
	SAFENEWWITHCONSTRUCTOR(pLambda,
		MyVectorHandler,
		MyVectorHandler(iUnkStates*iNumDofs));

	pX->Reset();
	pXPrime->Reset();
	pXPrimePrime->Reset();
	pLambda->Reset();

	/*
	 * Immediately link DataManager to current solution
	 *
	 * this should work as long as the last unknown time step is put
	 * at the beginning of pX, pXPrime
	 */

	pDM->LinkToSolution(*pX, *pXPrime, *pXPrimePrime, *pLambda);

	/* a questo punto si costruisce il nonlinear solver */
	pNLS = AllocateNonlinearSolver();

	/* FIXME: Serve?*/
	MyVectorHandler Scale(iNumDofs);
	if (bScale) {
		/* collects scale factors from data manager */
		pDM->SetScale(Scale);
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
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* registers scale factors at nonlinear solver */
		pResTestScale->SetScale(&Scale);

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
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// pippero
	switch (GetProblemType()) {
	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED:
		{
			NonlinearSolverTestRange *pRT = new NonlinearSolverTestRange(pResTest);
			NonlinearSolverTestRange *pST = new NonlinearSolverTestRange(pSolTest);
			pDM->IDSetTest(pRT, pST);
			pResTest = pRT;
			pSolTest = pST;
		} break;

	default:
		break;
	}

	/* registers tests in nonlinear solver */
	pNLS->SetTest(pResTest, pSolTest);

	/*
	 * Dell'assemblaggio iniziale dei vincoli se ne occupa il DataManager
	 * in quanto e' lui il responsabile dei dati della simulazione,
	 * e quindi anche della loro coerenza. Inoltre e' lui a sapere
	 * quali equazioni sono di vincolo o meno.
	 */

	pDM->SetValue(*pX, *pXPrime);


	/*
	 * Prepare output
	 */
	pDM->OutputPrepare();

	/*
	 * Dialoga con il DataManager per dargli il tempo iniziale
	 * e per farsi inizializzare i vettori di soluzione e derivata
	 */
	/* FIXME: the time is already set by DataManager, but FileDrivers
	 * have not been ServePending'd
	 */
	dTime = dInitialTime - dInitialTimeStep;
	pDM->SetTime(dTime + dInitialTimeStep, dInitialTimeStep, 0);

	integer iTotIter = 0;
	integer iStIter = 0;
	doublereal dTotErr = 0.;
	doublereal dTest = std::numeric_limits<double>::max();
	doublereal dSolTest = std::numeric_limits<double>::max();
	bool bSolConv = false;
	bool bOut = false;

	dTest = 0.;
	dSolTest = 0.;

	mbdyn_signal_init(0);

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
	if (pOutputMeter) {
		pOutputMeter->SetDrvHdl(pDM->pGetDrvHdl());
		pNLS->SetOutputMeter(pOutputMeter->pCopy());
	}

	pRegularSteps->SetDataManager(pDM);
	pRegularSteps->OutputTypes(DEBUG_LEVEL_MATCH(MYDEBUG_PRED));

#ifdef USE_EXTERNAL
	pNLS->SetExternal(External::EMPTY);
#endif /* USE_EXTERNAL */

	doublereal dCurrTimeStep = 0.;

#ifdef USE_EXTERNAL
	/* comunica che gli ultimi dati inviati sono la condizione iniziale */
	External::SendInitial();
#endif /* USE_EXTERNAL */

#ifdef USE_EXTERNAL
	/* il prossimo passo e' un regular */
	pNLS->SetExternal(External::REGULAR);
#endif /* USE_EXTERNAL */

	long lStep = -1;
	DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);

	if (mbdyn_stop_at_end_of_time_step()) {
		/* Fa l'output della soluzione al primo passo ed esce */
		Out << "Interrupted during first step." << std::endl;
		throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
	}

	bSolConv = false;
	dRefTimeStep = dInitialTimeStep;
	dCurrTimeStep = dRefTimeStep;

	if (pRTSolver) {
		pRTSolver->Init();
	}

	bool bOutputCounter = outputCounter() && isatty(fileno(stderr));
	const char *outputCounterPrefix = bOutputCounter ? "\n" : "";
	const char *outputCounterPostfix = outputStep() ? "\n" : "\r";

	ASSERT(pRegularSteps != NULL);

	/* Setup SolutionManager(s) */
	SetupSolmans(pRegularSteps->GetIntegratorNumUnknownStates(), true);

	/* this is here to force AfterConvergence()
	 * since we need to explicitly SetTime() to the next time step */
	pDM->SetValue(*pX, *pXPrime);

	pCurrStepIntegrator = pRegularSteps;
	while (true) {
		StepIntegrator::StepChange CurrStep
				= StepIntegrator::NEWSTEP;

		if (dTime >= dFinalTime) {
			if (pRTSolver) {
				pRTSolver->StopCommanded();
			}

			silent_cout(outputCounterPrefix
				<< "End of simulation at time "
				<< dTime << " after "
				<< lStep << " steps;" << std::endl
				<< "output in file \"" << sOutputFileName << "\"" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobian matrices: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl);

			if (pRTSolver) {
				pRTSolver->Log();
			}

			return;

		} else if (pRTSolver && pRTSolver->IsStopCommanded()) {
			silent_cout(outputCounterPrefix
				<< "Simulation is stopped by RTAI task" << std::endl
				<< "Simulation ended at time "
				<< dTime << " after "
				<< lStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobian matrices: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl);
			pRTSolver->Log();
			return;

		} else if (mbdyn_stop_at_end_of_time_step()
#ifdef USE_MPI
				|| (MPI_Finalized(&mpi_finalize), mpi_finalize)
#endif /* USE_MPI */
				)
		{
			if (pRTSolver) {
				pRTSolver->StopCommanded();
			}

			silent_cout(outputCounterPrefix
				<< "Interrupted!" << std::endl
				<< "Simulation ended at time "
				<< dTime << " after "
				<< lStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobian matrices: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl);

			if (pRTSolver) {
				pRTSolver->Log();
			}

			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		}

		lStep++;

		if (pRTSolver) {
			pRTSolver->Wait();
		}

		int retries = -1;
IfStepIsToBeRepeated:
		try {
			retries++;
			pDM->SetTime(dTime + dCurrTimeStep, dCurrTimeStep, lStep);
			if (outputStep()) {
				silent_cout("Step(" << lStep << ':' << retries << ") t=" << dTime + dCurrTimeStep << " dt=" << dCurrTimeStep << std::endl);
			}
			dTest = dynamic_cast<InverseDynamicsStepSolver *>(pRegularSteps)->Advance(this, dRefTimeStep,
					CurrStep, pX, pXPrime, pXPrimePrime, pLambda,
					iStIter, dTest, dSolTest);
		}

		catch (NonlinearSolver::NoConvergence) {
			if (dCurrTimeStep > dMinTimeStep) {
				/* Riduce il passo */
				CurrStep = StepIntegrator::REPEATSTEP;
				doublereal dOldCurrTimeStep = dCurrTimeStep;
				dCurrTimeStep = NewTimeStep(dCurrTimeStep,
						iStIter,
						CurrStep);
				if (dCurrTimeStep < dOldCurrTimeStep) {
					DEBUGCOUT("Changing time step"
						" during step "
						<< lStep << " after "
						<< iStIter << " iterations"
						<< std::endl);
					goto IfStepIsToBeRepeated;
				}
			}

			silent_cerr(outputCounterPrefix
				<< "Max iterations number "
				<< pRegularSteps->GetIntegratorMaxIters()
				<< " has been reached during"
				" step " << lStep << ';'
				<< std::endl
				<< "time step dt="
				<< dCurrTimeStep
				<< " cannot be reduced"
				" further;" << std::endl
				<< "aborting..." << std::endl);
			throw ErrMaxIterations(MBDYN_EXCEPT_ARGS);
		}

		catch (NonlinearSolver::ErrSimulationDiverged) {
			/*
			 * Mettere qui eventuali azioni speciali
			 * da intraprendere in caso di errore ...
			 */
			silent_cerr(outputCounterPrefix
				<< "Simulation diverged after "
				<< iStIter << " iterations, before "
				"reaching max iteration number "
				<< pRegularSteps->GetIntegratorMaxIters()
				<< " during step " << lStep << ';'
				<< std::endl
				<< "time step dt="
				<< dCurrTimeStep
				<< " cannot be reduced"
				" further;" << std::endl
				<< "aborting..." << std::endl);
			throw SimulationDiverged(MBDYN_EXCEPT_ARGS);
		}
		catch (NonlinearSolver::ConvergenceOnSolution) {
			bSolConv = true;
		}

		catch (EndOfSimulation& eos) {
			silent_cerr("Simulation ended during a regular step:\n"
				<< eos.what() << "\n");
			mbdyn_set_stop_at_end_of_time_step();
#ifdef USE_MPI
			MBDynComm.Abort(0);
#endif /* USE_MPI */

			if (pRTSolver) {
				pRTSolver->StopCommanded();
			}

			silent_cout(outputCounterPrefix
				<< "Simulation ended at time "
				<< dTime << " after "
				<< lStep << " steps;" << std::endl
				<< "total iterations: " << iTotIter << std::endl
				<< "total Jacobian matrices: " << pNLS->TotalAssembledJacobian() << std::endl
				<< "total error: " << dTotErr << std::endl);

			if (pRTSolver) {
				pRTSolver->Log();
			}

			return;
		}
		catch (...) {
			throw;
		}

		dTotErr += dTest;
		iTotIter += iStIter;

		bOut = pDM->Output(lStep, dTime + dCurrTimeStep, dCurrTimeStep);

		if (outputMsg()) {
			Out << "Step " << lStep
				<< " " << dTime+dCurrTimeStep
				<< " " << dCurrTimeStep
				<< " " << iStIter
				<< " " << dTest
				<< " " << dSolTest
				<< " " << bSolConv
				<< " " << bOut
				<< std::endl;
		}

		if (bOutputCounter) {
			silent_cout("Step " << std::setw(5) << lStep
				<< " " << std::setw(13) << dTime + dCurrTimeStep
				<< " " << std::setw(13) << dCurrTimeStep
				<< " " << std::setw(4) << iStIter
				<< " " << std::setw(13) << dTest
				<< " " << std::setw(13) << dSolTest
				<< " " << bSolConv
				<< " " << bOut
				<< outputCounterPostfix);
		}

		DEBUGCOUT("Step " << lStep
			<< " has been successfully completed "
			"in " << iStIter << " iterations" << std::endl);

		dRefTimeStep = dCurrTimeStep;
		dTime += dRefTimeStep;

		bSolConv = false;

		/* Calcola il nuovo timestep */
		dCurrTimeStep =
			NewTimeStep(dCurrTimeStep, iStIter, CurrStep);
		DEBUGCOUT("Current time step: " << dCurrTimeStep << std::endl);
	} // while (true)  END OF ENDLESS-LOOP
}  // InverseSolver::Run()

/* Distruttore */
InverseSolver::~InverseSolver(void)
{
	DEBUGCOUTFNAME("Solver::~Solver");

	SAFEDELETE(pX);
	SAFEDELETE(pXPrime);
	SAFEDELETE(pXPrimePrime);
	SAFEDELETE(pLambda);

	if (pdWorkSpace != NULL) {
		SAFEDELETEARR(pdWorkSpace);
	}

	if (pDM != NULL) {
		SAFEDELETE(pDM);
	}

	if (pRegularSteps) {
		SAFEDELETE(pRegularSteps);
	}

	if (pSM) {
		SAFEDELETE(pSM);
	}

	if (pNLS) {
		SAFEDELETE(pNLS);
	}
}

/* scrive il contributo al file di restart */
std::ostream &
InverseSolver::Restart(std::ostream& out,DataManager::eRestart type) const
{
#if 0
	out << "begin: inverse dynamics;" << std::endl;
	switch(type) {
	case DataManager::ATEND:
		out << "  #  initial time: " << pDM->dGetTime() << ";"
			<< std::endl
			<< "  #  final time: " << dFinalTime << ";"
			<< std::endl
			<< "  #  time step: " << dInitialTimeStep << ";"
			<< std::endl;
		break;
	case DataManager::ITERATIONS:
	case DataManager::TIME:
	case DataManager::TIMES:
		out << "  initial time: " << pDM->dGetTime()<< ";" << std::endl
			<< "  final time: " << dFinalTime << ";" << std::endl
			<< "  time step: " << dInitialTimeStep << ";"
			<< std::endl;
		break;
	default:
		ASSERT(0);
	}

	out << "  max iterations: " << pRegularSteps->GetIntegratorMaxIters()
		<< ";" << std::endl
		<< "  tolerance: " << pRegularSteps->GetIntegratorDTol();
	switch(ResTest) {
	case NonlinearSolverTest::NORM:
		out << ", test, norm" ;
		break;
	case NonlinearSolverTest::MINMAX:
		out << ", test, minmax" ;
		break;
	case NonlinearSolverTest::NONE:
		NO_OP;
	default:
		silent_cerr("unhandled nonlinear solver test type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (bScale) {
		out << ", scale"
			<< ", " << pRegularSteps->GetIntegratorDSolTol();
	}

	switch (SolTest) {
	case NonlinearSolverTest::NORM:
		out << ", test, norm" ;
		break;
	case NonlinearSolverTest::MINMAX:
		out << ", test, minmax" ;
		break;
	case NonlinearSolverTest::NONE:
		NO_OP;
	default:
		ASSERT(0);
	}
	out
		<< ";" << std::endl
		<< "  derivatives max iterations: " << pDerivativeSteps->GetIntegratorMaxIters() << ";" << std::endl
		<< "  derivatives tolerance: " << pDerivativeSteps->GetIntegratorDTol() << ";" << std::endl
		<< "  derivatives coefficient: " << dDerivativesCoef << ";" << std::endl;
	if (iDummyStepsNumber) {
		out
			<< "  dummy steps max iterations: " << pDummySteps->GetIntegratorMaxIters() << ";" << std::endl
			<< "  dummy steps tolerance: " << pDummySteps->GetIntegratorDTol() << ";" << std::endl;
	}
	out
		<< "  dummy steps number: " << iDummyStepsNumber << ";" << std::endl
		<< "  dummy steps ratio: " << dDummyStepsRatio << ";" << std::endl;
	switch (NonlinearSolverType) {
	case NonlinearSolver::MATRIXFREE:
		out << "  #  nonlinear solver: matrix free;" << std::endl;
		break;
	case NonlinearSolver::NEWTONRAPHSON:
	default :
		out << "  nonlinear solver: newton raphson";
		if (!bTrueNewtonRaphson) {
			out << ", modified, " << iIterationsBeforeAssembly;
			if (bKeepJac) {
				out << ", keep jacobian matrix";
			}
			if (bHonorJacRequest) {
				out << ", honor element requests";
			}
		}
		out << ";" << std::endl;
	}
	out << "  solver: ";
	RestartLinSol(out, CurrLinearSolver);
	out << "end: multistep;" << std::endl << std::endl;
#endif
	return out;
}

InverseDynamics::Type
InverseSolver::GetProblemType(void) const
{
	return ProblemType;
}

void
InverseSolver::GetWeight(InverseDynamics::Order iOrder, doublereal& dw1, doublereal& dw2) const
{
	switch (iOrder) {
	case InverseDynamics::POSITION:
		dw1 = this->dw1[0];
		dw2 = this->dw2[0];
		break;

	case InverseDynamics::VELOCITY:
		dw1 = this->dw1[1];
		dw2 = this->dw2[1];
		break;

	case InverseDynamics::ACCELERATION:
		dw1 = 0.;
		dw2 = this->dw1[2] + this->dw2[2];
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Dati dell'integratore */
void
InverseSolver::ReadData(MBDynParser& HP)
{
	DEBUGCOUTFNAME("InverseDynamics::ReadData");

	/* parole chiave */
	const char* sKeyWords[] = {
		"begin",
		"inverse" "dynamics",
		"end",

		"initial" "time",
		"final" "time",
		"time" "step",
		"min" "time" "step",
		"max" "time" "step",
		"tolerance",
		"max" "iterations",
		"modify" "residual" "test",

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
			/* DEPRECATED */ "jacobian" /* END OF DEPRECATED */ ,
			"jacobian" "matrix",
			"bailout",
			"messages",
		"output" "meter",

		"method",
			"inverse" "default",


		/* DEPRECATED */
		"true",
		"modified",
		/* END OF DEPRECATED */

		"strategy",
			"factor",
			"no" "change",
			"change",


		/* DEPRECATED */
		"solver",
		"interface" "solver",
		/* END OF DEPRECATED */
		"linear" "solver",
		"interface" "linear" "solver",

		/* DEPRECATED */
		"preconditioner",
		/* END OF DEPRECATED */

		"nonlinear" "solver",
			"default",
			"newton" "raphson",
			"matrix" "free",
				"bicgstab",
				"gmres",
					/* DEPRECATED */ "full" "jacobian" /* END OF DEPRECATED */ ,
					"full" "jacobian" "matrix",

		/* RTAI stuff */
		"real" "time",

		/* multithread stuff */
		"threads",

		"problem" "type",
			"fully" "actuated" "collocated",
			"fully" "actuated" "non" "collocated",
			"underdetermined" "underactuated" "collocated",
			"underdetermined" "fully" "actuated",

		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		BEGIN = 0,
		INVERSEDYNAMICS,
		END,

		INITIALTIME,
		FINALTIME,
		TIMESTEP,
		MINTIMESTEP,
		MAXTIMESTEP,
		TOLERANCE,
		MAXITERATIONS,
		MODIFY_RES_TEST,

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
			JACOBIANMATRIX,
			BAILOUT,
			MESSAGES,
		OUTPUTMETER,

		METHOD,
			INVERSEDEFAULT,

		/* DEPRECATED */
		NR_TRUE,
		MODIFIED,
		/* END OF DEPRECATED */

		STRATEGY,
		STRATEGYFACTOR,
		STRATEGYNOCHANGE,
		STRATEGYCHANGE,

		/* DEPRECATED */
		SOLVER,
		INTERFACESOLVER,
		/* END OF DEPRECATED */
		LINEARSOLVER,
		INTERFACELINEARSOLVER,

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
					FULLJACOBIANMATRIX,

		/* RTAI stuff */
		REALTIME,

		THREADS,

		PROBLEMTYPE,
			FAC,
			FANC,
			UDUAC,
			UDFA,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* legge i dati della simulazione */
	if (KeyWords(HP.GetDescription()) != BEGIN) {
		silent_cerr("Error: <begin> expected at line "
			<< HP.GetLineData() << "; aborting..." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (KeyWords(HP.GetWord()) != INVERSEDYNAMICS) {
		silent_cerr("Error: <begin: inverse dynamics;> expected at line "
			<< HP.GetLineData() << "; aborting..." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* dati letti qui ma da passare alle classi
	 *	StepIntegration e NonlinearSolver
	 */

	doublereal dTol = ::dDefaultTol;
	doublereal dSolutionTol = 0.;
	integer iMaxIterations = ::iDefaultMaxIterations;
	bool bModResTest = false;

#ifdef USE_MULTITHREAD
	bool bSolverThreads(false);
	unsigned nSolverThreads(0);
#endif /* USE_MULTITHREAD */


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
				silent_cerr("warning: final time " << dFinalTime
					<< " is less than initial time "
					<< dInitialTime << ';' << std::endl
					<< "this will cause the simulation"
					" to abort" << std::endl);
			}
			break;

		case TIMESTEP:
			dInitialTimeStep = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, "Initial time step is "
				<< dInitialTimeStep << std::endl);

			if (dInitialTimeStep == 0.) {
				silent_cerr("warning, null initial time step"
					" is not allowed" << std::endl);
			} else if (dInitialTimeStep < 0.) {
				dInitialTimeStep = -dInitialTimeStep;
				silent_cerr("warning, negative initial time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dInitialTimeStep
					<< " will be considered" << std::endl);
			}
			break;

		case MINTIMESTEP:
			dMinTimeStep = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, "Min time step is "
				<< dMinTimeStep << std::endl);

			if (dMinTimeStep == 0.) {
				silent_cerr("warning, null minimum time step"
					" is not allowed" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			} else if (dMinTimeStep < 0.) {
				dMinTimeStep = -dMinTimeStep;
				silent_cerr("warning, negative minimum time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dMinTimeStep
					<< " will be considered" << std::endl);
			}
			break;

		case MAXTIMESTEP:
			dMaxTimeStep = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT, "Max time step is "
				<< dMaxTimeStep << std::endl);

			if (dMaxTimeStep == 0.) {
				silent_cout("no max time step limit will be"
					" considered" << std::endl);
			} else if (dMaxTimeStep < 0.) {
				dMaxTimeStep = -dMaxTimeStep;
				silent_cerr("warning, negative max time step"
					" is not allowed;" << std::endl
					<< "its modulus " << dMaxTimeStep
					<< " will be considered" << std::endl);
			}
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
				silent_cerr("Don't know when to abort,"
					" so I'm going to abort now" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;
		}

		case OUTPUT: {
			unsigned OF = OUTPUT_DEFAULT;
			bool setOutput = false;

			while (HP.IsArg()) {
				KeyWords OutputFlag(KeyWords(HP.GetWord()));
				switch (OutputFlag) {
				case NONE:
					OF = OUTPUT_NONE;
					setOutput = true;
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
				case JACOBIANMATRIX:
					OF |= OUTPUT_JAC;
					break;

				case BAILOUT:
					OF |= OUTPUT_BAILOUT;
					break;

				case MESSAGES:
					OF |= OUTPUT_MSG;
					break;

				default:
					silent_cerr("Unknown output flag "
						"at line " << HP.GetLineData()
						<< "; ignored" << std::endl);
					break;
				}
			}

			if (setOutput) {
				SetOutputFlags(OF);
			} else {
				AddOutputFlags(OF);
			}

			break;
		}

		case OUTPUTMETER:
			SetOutputMeter(HP.GetDriveCaller(true));
			break;

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
					dTol = ::dDefaultTol;
					silent_cerr("warning, residual tolerance "
						"< 0. is illegal; "
						"using default value " << dTol
						<< std::endl);
				}
			}

			/* safe default */
			if (dTol == 0.) {
				ResTest = NonlinearSolverTest::NONE;
			}

			if (HP.IsArg()) {
				if (HP.IsKeyWord("test")) {
					if (HP.IsKeyWord("norm")) {
						ResTest = NonlinearSolverTest::NORM;
					} else if (HP.IsKeyWord("minmax")) {
						ResTest = NonlinearSolverTest::MINMAX;
					} else if (HP.IsKeyWord("none")) {
						ResTest = NonlinearSolverTest::NONE;
					} else {
						silent_cerr("unknown test "
							"method at line "
							<< HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					if (HP.IsKeyWord("scale")) {
						if (ResTest == NonlinearSolverTest::NONE) {
							silent_cerr("it's a nonsense "
								"to scale a disabled test; "
								"\"scale\" ignored"
								<< std::endl);
							bScale = false;
						} else {
							bScale = true;
						}
					}
				}
			}

			if (HP.IsArg()) {
				if (!HP.IsKeyWord("null")) {
					dSolutionTol = HP.GetReal();
				}

				/* safe default */
				if (dSolutionTol != 0.) {
					SolTest = NonlinearSolverTest::NORM;
				}

				if (HP.IsArg()) {
					if (HP.IsKeyWord("test")) {
						if (HP.IsKeyWord("norm")) {
							SolTest = NonlinearSolverTest::NORM;
						} else if (HP.IsKeyWord("minmax")) {
							SolTest = NonlinearSolverTest::MINMAX;
						} else if (HP.IsKeyWord("none")) {
							SolTest = NonlinearSolverTest::NONE;
						} else {
							silent_cerr("unknown test "
								"method at line "
								<< HP.GetLineData()
								<< std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
					}
				}

			} else if (dTol == 0.) {
				silent_cerr("need solution tolerance "
					"with null residual tolerance"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (dSolutionTol < 0.) {
				dSolutionTol = 0.;
				silent_cerr("warning, solution tolerance "
					"< 0. is illegal; "
					"solution test is disabled"
					<< std::endl);
			}

			if (dTol == 0. && dSolutionTol == 0.) {
				silent_cerr("both residual and solution "
					"tolerances are zero" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			DEBUGLCOUT(MYDEBUG_INPUT, "tolerance = " << dTol
					<< ", " << dSolutionTol << std::endl);
			break;
		}


		case MAXITERATIONS: {
			iMaxIterations = HP.GetInt();
			if (iMaxIterations < 1) {
				iMaxIterations = ::iDefaultMaxIterations;
				silent_cerr("warning, max iterations "
					"< 1 is illegal; using default value "
					<< iMaxIterations
					<< std::endl);
			}
			DEBUGLCOUT(MYDEBUG_INPUT,
					"Max iterations = "
					<< iMaxIterations << std::endl);
			break;
		}

		case MODIFY_RES_TEST:
			if (bParallel) {
				silent_cerr("\"modify residual test\" "
					"not supported by schur data manager "
					"at line " << HP.GetLineData()
					<< "; ignored" << std::endl);
			} else {
				bModResTest = true;
				DEBUGLCOUT(MYDEBUG_INPUT,
					"Modify residual test" << std::endl);
			}
			break;




		case NEWTONRAPHSON: {
			pedantic_cout("Newton Raphson is deprecated; use "
					"\"nonlinear solver: newton raphson "
					"[ , modified, <steps> ]\" instead"
					<< std::endl);
			KeyWords NewRaph(KeyWords(HP.GetWord()));
			switch(NewRaph) {
			case MODIFIED:
				bTrueNewtonRaphson = 0;
				if (HP.IsArg()) {
					iIterationsBeforeAssembly = HP.GetInt();
				} else {
					iIterationsBeforeAssembly = ::iDefaultIterationsBeforeAssembly;
				}
				DEBUGLCOUT(MYDEBUG_INPUT, "Modified "
						"Newton-Raphson will be used; "
						"matrix will be assembled "
						"at most after "
						<< iIterationsBeforeAssembly
						<< " iterations" << std::endl);
				break;

			default:
				silent_cerr("warning: unknown case; "
					"using default" << std::endl);

			/* no break: fall-thru to next case */
			case NR_TRUE:
				bTrueNewtonRaphson = 1;
				iIterationsBeforeAssembly = 0;
				break;
			}
			break;
		}

		case END:
			if (KeyWords(HP.GetWord()) != INVERSEDYNAMICS) {
				silent_cerr("\"end: inverse dynamics;\" expected "
					"at line " << HP.GetLineData()
					<< "; aborting..." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
				 *     <min iterations> ,
				 *     <max iterations> ;
				 */

				StrategyFactor.dReductionFactor = HP.GetReal();
				if (StrategyFactor.dReductionFactor >= 1.) {
					silent_cerr("warning, "
						"illegal reduction factor "
						"at line " << HP.GetLineData()
						<< "; default value 1. "
						"(no reduction) will be used"
						<< std::endl);
					StrategyFactor.dReductionFactor = 1.;
				}

				StrategyFactor.iStepsBeforeReduction = HP.GetInt();
				if (StrategyFactor.iStepsBeforeReduction <= 0) {
					silent_cerr("warning, "
						"illegal number of steps "
						"before reduction at line "
						<< HP.GetLineData()
						<< "; default value 1 will be "
						"used (it may be dangerous)"
						<< std::endl);
					StrategyFactor.iStepsBeforeReduction = 1;
				}

				StrategyFactor.dRaiseFactor = HP.GetReal();
				if (StrategyFactor.dRaiseFactor <= 1.) {
					silent_cerr("warning, "
						"illegal raise factor at line "
						<< HP.GetLineData()
						<< "; default value 1. "
						"(no raise) will be used"
						<< std::endl);
					StrategyFactor.dRaiseFactor = 1.;
				}

				StrategyFactor.iStepsBeforeRaise = HP.GetInt();
				if (StrategyFactor.iStepsBeforeRaise <= 0) {
					silent_cerr("warning, "
						"illegal number of steps "
						"before raise at line "
						<< HP.GetLineData()
						<< "; default value 1 will be "
						"used (it may be dangerous)"
						<< std::endl);
					StrategyFactor.iStepsBeforeRaise = 1;
				}

				StrategyFactor.iMinIters = HP.GetInt();
				if (StrategyFactor.iMinIters <= 0) {
					silent_cerr("warning, "
						"illegal minimum number "
						"of iterations at line "
						<< HP.GetLineData()
						<< "; default value 0 will be "
						"used (never raise)"
						<< std::endl);
					StrategyFactor.iMinIters = 1;
				}

				StrategyFactor.iMaxIters = 0;
				if (HP.IsArg()) {
					StrategyFactor.iMaxIters = HP.GetInt();
					if (StrategyFactor.iMaxIters <= 0) {
						silent_cerr("warning, "
							"illegal mmaximim number "
							"of iterations at line "
							<< HP.GetLineData()
							<< "; default value will be "
							"used"
							<< std::endl);
						StrategyFactor.iMaxIters = 0;
					}
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
						<< "Max iterations: "
						<< StrategyFactor.iMaxIters
						<< std::endl);
				break;
			}

			case STRATEGYNOCHANGE: {
				CurrStrategy = NOCHANGE;
				break;
			}

			case STRATEGYCHANGE: {
				CurrStrategy = CHANGE;
				pStrategyChangeDrive = HP.GetDriveCaller(true);
				break;
			}

			default:
				silent_cerr("unknown time step control "
					"strategy at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;
		}

		case SOLVER:
			silent_cerr("\"solver\" keyword at line "
					<< HP.GetLineData()
					<< " is deprecated; "
					"use \"linear solver\" instead"
					<< std::endl);
		case LINEARSOLVER:
			ReadLinSol(CurrLinearSolver, HP);
			break;

		case INTERFACESOLVER:
			silent_cerr("\"interface solver\" keyword at line "
					<< HP.GetLineData()
					<< " is deprecated; "
					"use \"interface linear solver\" "
					"instead" << std::endl);
		case INTERFACELINEARSOLVER:
			ReadLinSol(CurrIntSolver, HP, true);

#ifndef USE_MPI
			silent_cerr("Interface solver only allowed "
				"when compiled with MPI support" << std::endl);
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
				silent_cerr("unknown nonlinear solver "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
						pedantic_cout("Use of deprecated \"keep jacobian\" at line " << HP.GetLineData() << std::endl);
						bKeepJac = true;

					} else if (HP.IsKeyWord("keep" "jacobian" "matrix")) {
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
					if (HP.IsKeyWord("honor" "element" "requests")) {
						bHonorJacRequest = true;
						DEBUGLCOUT(MYDEBUG_INPUT,
								"honor elements' "
								"request to update "
								"the preconditioner"
								<< std::endl);
					}
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
					silent_cerr("unknown iterative "
						"solver at line "
						<< HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
					case FULLJACOBIANMATRIX:
						PcType = Preconditioner::FULLJACOBIANMATRIX;
						if (HP.IsKeyWord("steps")) {
							iPrecondSteps = HP.GetInt();
							DEBUGLCOUT(MYDEBUG_INPUT,
									"number of steps "
									"before recomputing "
									"the preconditioner: "
									<< iPrecondSteps
									<< std::endl);
						}
						if (HP.IsKeyWord("honor" "element" "requests")) {
							bHonorJacRequest = true;
							DEBUGLCOUT(MYDEBUG_INPUT,
									"honor elements' "
									"request to update "
									"the preconditioner"
									<< std::endl);
						}
						break;

						/* add other preconditioners
						 * here */

					default:
						silent_cerr("unknown "
							"preconditioner "
							"at line "
							<< HP.GetLineData()
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					break;
				}
				break;

			default:
				ASSERT(0);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case REALTIME:
			pRTSolver = ReadRTSolver(this, HP);
			break;

		case THREADS:
			if (HP.IsKeyWord("auto")) {
#ifdef USE_MULTITHREAD
				int n = get_nprocs();
				/* sanity checks ... */
				if (n <= 0) {
					silent_cerr("got " << n << " CPUs "
							"at line "
							<< HP.GetLineData()
							<< std::endl);
					nThreads = 1;
				} else {
					nThreads = n;
				}
#else /* ! USE_MULTITHREAD */
				silent_cerr("configure with "
						"--enable-multithread "
						"for multithreaded assembly"
						<< std::endl);
#endif /* ! USE_MULTITHREAD */

			} else if (HP.IsKeyWord("disable")) {
#ifdef USE_MULTITHREAD
				nThreads = 1;
#endif /* USE_MULTITHREAD */

			} else {
#ifdef USE_MULTITHREAD
				bool bAssembly = false;
				bool bSolver = false;
				bool bAll = true;
				unsigned nt;
#endif // USE_MULTITHREAD

				if (HP.IsKeyWord("assembly")) {
#ifdef USE_MULTITHREAD
					bAll = false;
					bAssembly = true;
#endif // USE_MULTITHREAD

				} else if (HP.IsKeyWord("solver")) {
#ifdef USE_MULTITHREAD
					bAll = false;
					bSolver = true;
#endif // USE_MULTITHREAD
				}

#ifdef USE_MULTITHREAD
				nt =
#endif // USE_MULTITHREAD
				HP.GetInt();

#ifdef USE_MULTITHREAD
				if (bAll || bAssembly) {
					nThreads = nt;
				}

				if (bAll || bSolver) {
					bSolverThreads = true;
					nSolverThreads = nt;
				}
#else /* ! USE_MULTITHREAD */
				silent_cerr("configure with "
						"--enable-multithread "
						"for multithreaded assembly"
						<< std::endl);
#endif /* ! USE_MULTITHREAD */
			}
			break;

		case PROBLEMTYPE:
			switch (KeyWords(HP.GetWord())) {
			case FAC:
				// default
				ProblemType = InverseDynamics::FULLY_ACTUATED_COLLOCATED;
				break;

			case FANC:
				ProblemType = InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED;
				break;

			case UDUAC:
				ProblemType = InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED;
				if (HP.IsKeyWord("weights")) {
					for (int iCnt = 0; iCnt <= InverseDynamics::VELOCITY; iCnt++) {
						dw1[iCnt] = HP.GetReal();
						dw2[iCnt] = HP.GetReal();
					}
					dw1[InverseDynamics::ACCELERATION] = 0.;
					dw2[InverseDynamics::ACCELERATION] = 0.;

				} else {
					for (int iCnt = 0; iCnt <= InverseDynamics::ACCELERATION; iCnt++) {
						dw1[iCnt] = 1.;
						dw2[iCnt] = 0.;
					}
				}

				break;

			case UDFA:
				ProblemType = InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED;
				if (HP.IsKeyWord("weights")) {
					for (int iCnt = 0; iCnt <= InverseDynamics::ACCELERATION; iCnt++) {
						dw1[iCnt] = HP.GetReal();
						dw2[iCnt] = HP.GetReal();
					}

				} else {
					for (int iCnt = 0; iCnt <= InverseDynamics::ACCELERATION; iCnt++) {
						dw1[iCnt] = 1.;
						dw2[iCnt] = 0.;
					}
				}
				break;

			default:
				silent_cerr("inverse dynamics: unrecognized problem type at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (HP.IsArg()) {
				silent_cerr("inverse dynamics: semicolon expected at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		default:
			silent_cerr("unknown description at line "
				<< HP.GetLineData() << "; aborting..."
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

EndOfCycle: /* esce dal ciclo di lettura */

	switch (CurrStrategy) {
	case FACTOR:
		if (StrategyFactor.iMaxIters <= StrategyFactor.iMinIters) {
			silent_cerr("warning, "
				<< "strategy maximum number "
				<< "of iterations "
				<< "is <= minimum: "
				<< StrategyFactor.iMaxIters << " <= "
				<< StrategyFactor.iMinIters << "; "
				<< "the maximum global iteration value "
				<< iMaxIterations << " "
				<< "will be used"
				<< std::endl);
			StrategyFactor.iMaxIters = iMaxIterations;
		}
		break;

	default:
		break;
	}

	if (dFinalTime < dInitialTime) {
		eAbortAfter = AFTER_ASSEMBLY;
	}

	if (dFinalTime == dInitialTime) {
		eAbortAfter = AFTER_DERIVATIVES;
	}

	SAFENEWWITHCONSTRUCTOR(pRegularSteps,
		InverseDynamicsStepSolver,
		InverseDynamicsStepSolver(iMaxIterations,
			dTol,
			dSolutionTol,
			0, /* Previous Steps to be used */
			1, /* Unknown States: FIXME: 2? 3? */
			bModResTest));

#ifdef USE_MULTITHREAD
	if (bSolverThreads) {
		if (CurrLinearSolver.SetNumThreads(nSolverThreads)) {
			silent_cerr("linear solver "
					<< CurrLinearSolver.GetSolverName()
					<< " does not support "
					"threaded solution" << std::endl);
		}
	}
#endif /* USE_MULTITHREAD */
}

