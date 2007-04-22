/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
  * Copyright (C) 2007
  * Alessandro Fumagalli	<alessandro.fumagalli@polimi.it>
  *
  * Handler of inverse dynamics problem:
  * 
  *   - Initialize data (nodes and elements)
  *   - Allocate memory for solution vectors 
  *   - Handle outputs 
  *
  */
  
#ifndef INVERSESOLVER_H
#define INVERSESOLVER_H  

#define RTAI_LOG

#include <unistd.h>
#include "ac/float.h"
#include "ac/math.h"

class InverseSolver;
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "dataman.h"
#include "schurdataman.h"
#include "schsolman.h"
#include "solver.h"
#include <deque>
#include "linsol.h"
#include "stepsol.h"
#include "nonlin.h"
#include "mfree.h"
#include "precond.h"

class InverseSolver : public Solver {
public:
	class ErrGeneric {};
	class ErrMaxIterations{};
	class SimulationDiverged{};

private:
#ifdef USE_MULTITHREAD
	unsigned nThreads;
#endif /* USE_MULTITHREAD */

   	enum Strategy {
		NOCHANGE,
		CHANGE,
		FACTOR
	} CurrStrategy;

   	const char *sInputFileName;
   	const char *sOutputFileName;
   	MBDynParser& HP;
	 
   	/* Dati per strategia FACTOR */
   	struct {
      		doublereal dReductionFactor;
      		doublereal dRaiseFactor;
      		integer iStepsBeforeReduction;
      		integer iStepsBeforeRaise;
      		integer iMinIters;
      		integer iMaxIters;
   	} StrategyFactor;

   	/* Dati per strategia DRIVER_CHANGE */
	DriveCaller* pStrategyChangeDrive;

#ifdef USE_RTAI
	bool bRT;
	bool bRTAllowNonRoot;
	enum {
		MBRTAI_UNKNOWN,
		MBRTAI_WAITPERIOD,
		MBRTAI_SEMAPHORE,
		MBRTAI_LASTMODE
	} RTMode;
	bool bRTHard;
	long long lRTPeriod;		/* if RTMode == MBRTAI_WAITPERIOD */
	void *RTSemPtr_in;			/* if RTMode == MBRTAI_SEMAPHORE */
	void *RTSemPtr_out;
	unsigned long RTStackSize;
	int RTCpuMap;
#ifdef RTAI_LOG
	bool bRTlog;
	void *mbxlog;
	char *LogProcName;
#endif /*RTAI_LOG*/

	bool RTWaitPeriod(void) const {
		return (RTMode == MBRTAI_WAITPERIOD);
	};

	bool RTSemaphore(void) const {
		return (RTMode == MBRTAI_SEMAPHORE);
	};
#endif /* USE_RTAI */

   	/* Strutture di gestione dei dati */
	integer iNumPreviousVectors;
	integer iUnkStates;
	doublereal* pdWorkSpace;
	
	MyVectorHandler* pX;                  /* queque vett. stati inc. */
  	MyVectorHandler* pXPrime;             /* queque vett. stati d. inc. */ 
	MyVectorHandler* pXPrimePrime;

   	/* Dati della simulazione */
   	doublereal dTime;
   	doublereal dInitialTime;
   	doublereal dFinalTime;
   	doublereal dRefTimeStep;
   	doublereal dInitialTimeStep;
   	doublereal dMinimumTimeStep;
   	doublereal dMaxTimeStep;

   	/* Dati dei passi fittizi di trimmaggio iniziale */
   	integer iFictitiousStepsNumber;
   	doublereal dFictitiousStepsRatio;
   
   	/* Flags vari */
	enum AbortAfter {
		AFTER_UNKNOWN,
		AFTER_INPUT,
		AFTER_ASSEMBLY,
		AFTER_DERIVATIVES,
		AFTER_DUMMY_STEPS
	};
	AbortAfter eAbortAfter;

	/* Parametri per la variazione passo */
   	integer iStepsAfterReduction;
   	integer iStepsAfterRaise;
	integer iWeightedPerformedIters;
   	flag bLastChance;

   	/* Parametri per il metodo */
	enum StepIntegratorType {
			INT_CRANKNICHOLSON,
			INT_MS2,
			INT_HOPE,
			INT_THIRDORDER,
			INT_IMPLICITEULER,
			INV_DEFAULT,
			INT_UNKNOWN
	};
	StepIntegratorType RegularType,  FictitiousType;
	
   	StepIntegrator* pRegularSteps;
	
	doublereal dDerivativesCoef;
	
	/* Type of linear solver */
	LinSol CurrLinearSolver;

	/* Parameters for convergence tests */
	NonlinearSolverTest::Type ResTest;
	NonlinearSolverTest::Type SolTest;
	bool bScale;

   	/* Parametri per solutore nonlineare */
   	bool bTrueNewtonRaphson;
   	bool bKeepJac;
   	integer iIterationsBeforeAssembly;
	NonlinearSolver::Type NonlinearSolverType;
	MatrixFreeSolver::SolverType MFSolverType;
	doublereal dIterTol;
	Preconditioner::PrecondType PcType;
	integer iPrecondSteps;
	integer iIterativeMaxSteps;
	doublereal dIterertiveEtaMax;
	doublereal dIterertiveTau;
	bool bHonorJacRequest;

/* FOR PARALLEL SOLVERS */
	bool bParallel;
	SchurDataManager *pSDM;
	LinSol CurrIntSolver;
	integer iNumLocDofs;		/* Dimensioni problema locale */
	integer iNumIntDofs;		/* Dimensioni interfaccia locale */
	integer* pLocDofs;		/* Lista dof locali (stile FORTRAN) */
	integer* pIntDofs;		/* Lista dei dofs di interfaccia */
	Dof* pDofs;
	SolutionManager *pLocalSM;
/* end of FOR PARALLEL SOLVERS */

	/* gestore dei dati */
	DataManager* pDM;
     	/* Dimensioni del problema; FIXME: serve ancora? */
   	integer iNumDofs;

	/* il solution manager v*/
	SolutionManager *pSM;
	NonlinearSolver* pNLS;
	
	/* corregge i puntatori per un nuovo passo */
   	inline void Flip(void);

   	/* Lettura dati */
   	void ReadData(MBDynParser& HP);

   	/* Nuovo delta t */
   	doublereal NewTimeStep(doublereal dCurrTimeStep, 
			       integer iPerformedIters,
			       StepIntegrator::StepChange Dmy 
			       = StepIntegrator::NEWSTEP);
	/* Alloca Solman */
	SolutionManager *const AllocateSolman(integer iNLD, integer iLWS = 0);
	/* Alloca SchurSolman */
	SolutionManager *const AllocateSchurSolman(integer iStates);
	/* Alloca Nonlinear Solver */
	NonlinearSolver *const AllocateNonlinearSolver();
	/* Alloca tutti i solman*/
	void SetupSolmans(integer iStates, bool bCanBeParallel = false);
   

public:   
   	/* costruttore */
   	InverseSolver(MBDynParser& HP, 
		       	    const char* sInputFileName, 
			    const char* sOutputFileName,
			    bool bParallel = false);

   	/* esegue la simulazione */
   	void Run(void);

   	/* distruttore: esegue tutti i distruttori e libera la memoria */
   	~InverseSolver(void);
	std::ostream & Restart(std::ostream& out, DataManager::eRestart type) const;

	/* EXPERIMENTAL */
	/* FIXME: better const'ify? */
	DataManager *pGetDataManager(void) const {
		return pDM;
	};
	SolutionManager *pGetSolutionManager(void) const {
		return pSM;
	};
	const LinSol& GetLinearSolver(void) const {
		return CurrLinearSolver;
	};
	NonlinearSolver *pGetNonlinearSolver(void) const {
		return pNLS;
	};
	doublereal GetDInitialTimeStep(void) const {
		return dInitialTimeStep;
	};
};

/* InverseSolver - end */

#endif /* INVERSESOLVER_H */

