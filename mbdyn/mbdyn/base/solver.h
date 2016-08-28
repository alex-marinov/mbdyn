/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
  * Copyright (C) 2003-2015
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * Classe che gestisce la soluzione del problema:
  * 
  *   - Inizializza i dati (nodi ed elementi)
  *   - Alloca i vettori degli stati necessari
  *   - Stabililisce il tipo di integratore e quanti passi fare con esso
  *   - determina il passo temporale
  *   - gestisce le operazioni di output
  *
  */
  
#ifndef SOLVER_H
#define SOLVER_H  

#include <unistd.h>
#include <cfloat>
#include <cmath>
#include <limits>

class Solver;
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "dataman.h"
#include "schurdataman.h"
#include "schsolman.h"
#include <deque>
#include "linsol.h"
#include "stepsol.h"
#include "nonlin.h"
#include "linesearch.h"
#include "mfree.h"
#include "cstring"
#include "precond.h"
#include "rtsolver.h"
#include "TimeStepControl.h"
extern "C" int mbdyn_stop_at_end_of_iteration(void);
extern "C" int mbdyn_stop_at_end_of_time_step(void);
extern "C" void mbdyn_set_stop_at_end_of_iteration(void);
extern "C" void mbdyn_set_stop_at_end_of_time_step(void);

class Solver : public SolverDiagnostics, protected NonlinearSolverOptions {
public:
 	class ErrGeneric : public MBDynErrBase {
  	public:
 		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
  	};
 	class ErrMaxIterations : public MBDynErrBase {
  	public:
 		ErrMaxIterations(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
  	};
 	class SimulationDiverged : public MBDynErrBase {
  	public:
 		SimulationDiverged(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
  	};
 	class EndOfSimulation : public MBDynErrBase {
  	private:
  		int EndCode;
  	public:
 		EndOfSimulation(const int e, MBDYN_EXCEPT_ARGS_DECL_NODEF) : 
 		MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU), EndCode(e) {};
  	};
  
protected:
#ifdef USE_MULTITHREAD
	unsigned nThreads;

	void ThreadPrepare(void);
#endif /* USE_MULTITHREAD */
	TimeStepControl * timeStepPtr;
	doublereal dCurrTimeStep;
	integer iStIter;
	doublereal dTime;
	DriveOwner MaxTimeStep;
	doublereal dMinTimeStep;
	StepIntegrator::StepChange CurrStep = StepIntegrator::NEWSTEP;
	std::string sInputFileName;
	std::string sOutputFileName;
   	MBDynParser& HP;
   	integer iMaxIterations;
//***** Method for StepTime ****/
public:
	
	doublereal dGetdMinTimeStep() {
		return dMinTimeStep;
	}

	void dSetMinTimeStep(doublereal val) {
		dMinTimeStep = val;
	}
	
public:
   	/* Dati per esecuzione di eigenanalysis */
	
	struct EigenAnalysis {
		bool bAnalysis;
		enum {
			EIG_NONE			= 0x0U,

			EIG_OUTPUT_MATRICES		= 0x1U,
			EIG_OUTPUT_FULL_MATRICES	= 0x2U,
			EIG_OUTPUT_SPARSE_MATRICES	= 0x4U,

			EIG_OUTPUT_EIGENVECTORS		= 0x8U,

			EIG_OUTPUT_MATRICES_MASK	= (EIG_OUTPUT_MATRICES|EIG_OUTPUT_FULL_MATRICES|EIG_OUTPUT_SPARSE_MATRICES),

			EIG_OUTPUT_GEOMETRY		= 0x10U,

			EIG_OUTPUT			= (EIG_OUTPUT_MATRICES_MASK|EIG_OUTPUT_EIGENVECTORS|EIG_OUTPUT_GEOMETRY),

			EIG_SOLVE			= 0x100U,

			EIG_PERMUTE			= 0x200U,
			EIG_SCALE			= 0x400U,
			EIG_BALANCE			= (EIG_PERMUTE|EIG_SCALE),

			EIG_USE_LAPACK			= 0x1000U,
			EIG_USE_ARPACK			= 0x2000U,
			EIG_USE_JDQZ			= 0x4000U,
			EIG_USE_MASK			= (EIG_USE_LAPACK|EIG_USE_ARPACK|EIG_USE_JDQZ),

			EIG_LAST
		};
		unsigned uFlags;

		// TODO: allow to specify eigenanalysis parameters
		// for each analysis in the list
		std::vector<doublereal> Analyses;
		std::vector<doublereal>::iterator currAnalysis;

		doublereal dParam;
		bool bOutputModes;

		doublereal dUpperFreq;
		doublereal dLowerFreq;

		// text output precision control
		int iMatrixPrecision;
		int iResultsPrecision;

		// ARPACK specific
		struct ARPACK {
			integer iNEV;
			integer iNCV;
			doublereal dTOL;
			ARPACK(void) : iNEV(0), iNCV(0), dTOL(0.) { NO_OP; };
		} arpack;

		// JDQZ specific
		struct JDQZ {
			doublereal eps;
			integer kmax;
			integer jmax;
			integer jmin;
			integer method;
			integer m;
			integer l;
			integer mxmv;
			integer maxstep;
			doublereal lock;
			integer order;
			integer testspace;

			enum Method {
				GMRES = 1,
				BICGSTAB = 2
			};

			JDQZ(void)
			: eps(std::numeric_limits<doublereal>::epsilon()),
			method(BICGSTAB), m(30), l(2), mxmv(100),
			maxstep(1000), lock(1e-9), order(0), testspace(3)
			{ NO_OP; };
		} jdqz;

		EigenAnalysis(void)
		: bAnalysis(false),
		uFlags(EIG_NONE),
		dParam(1.),
		bOutputModes(false),
		dUpperFreq(std::numeric_limits<doublereal>::max()),
		dLowerFreq(-1.)
		{
			currAnalysis = Analyses.end();
		};
	};

protected:
	struct EigenAnalysis EigAn;

	void Eig(bool bNewLine = false);

	RTSolverBase *pRTSolver;

   	/* Strutture di gestione dei dati */
	integer iNumPreviousVectors;
	integer iUnkStates;
	doublereal* pdWorkSpace;
	std::deque<MyVectorHandler*> qX;      /* queque vett. stati */
  	std::deque<MyVectorHandler*> qXPrime; /* queque vett. stati der. */ 
	MyVectorHandler* pX;                  /* queque vett. stati inc. */
  	MyVectorHandler* pXPrime;             /* queque vett. stati d. inc. */ 

   	/* Dati della simulazione */
   	
   	doublereal dInitialTime;
   	doublereal dFinalTime;
   	doublereal dRefTimeStep;
   	doublereal dInitialTimeStep;
   	
   	
   	doublereal dMaxResidual;
   	doublereal dMaxResidualDiff;

   	enum TimeStepFlags {
   		TS_SOFT_LIMIT = 0,
   		TS_HARD_LIMIT = 1
   	} eTimeStepLimit;

   	/* Dati dei passi fittizi di trimmaggio iniziale */
   	integer iDummyStepsNumber;
   	doublereal dDummyStepsRatio;
   
   	/* Flags vari */
	enum AbortAfter {
		AFTER_UNKNOWN,
		AFTER_INPUT,
		AFTER_ASSEMBLY,
		AFTER_DERIVATIVES,
		AFTER_DUMMY_STEPS
	};
	AbortAfter eAbortAfter;

   	/* Parametri per il metodo */
	enum StepIntegratorType {
			INT_CRANKNICOLSON,
			INT_MODCRANKNICOLSON,
			INT_MS2,
			INT_HOPE,
			INT_THIRDORDER,
			INT_IMPLICITEULER,
			INT_UNKNOWN
	};
	StepIntegratorType RegularType, DummyType;
	
   	StepIntegrator* pDerivativeSteps;
   	StepIntegrator* pFirstDummyStep;
	StepIntegrator* pDummySteps;
   	StepIntegrator* pFirstRegularStep;
   	StepIntegrator* pRegularSteps;
	StepIntegrator* pCurrStepIntegrator;
	
	DriveCaller* pRhoRegular;
	DriveCaller* pRhoAlgebraicRegular;
	DriveCaller* pRhoDummy;
	DriveCaller* pRhoAlgebraicDummy;
	
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
	struct LineSearchParameters LineSearch;

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

	/* Workaround: call this function instead of MaxTimeStep.dGet()
	 * until all postponed drive callers have been instantiated */
	doublereal dGetInitialMaxTimeStep() const;

protected:
	enum {
		SOLVER_STATUS_UNINITIALIZED,
		SOLVER_STATUS_PREPARED,
		SOLVER_STATUS_STARTED
	} eStatus;
	
	bool bOutputCounter;
	std::string outputCounterPrefix;
	std::string outputCounterPostfix;
	integer iTotIter;
	
	doublereal dTotErr;
	doublereal dTest;
	doublereal dSolTest;
	bool bSolConv;
	bool bOut;
	long lStep;

public:   
   	/* costruttore */
   	Solver(MBDynParser& HP, 
		const std::string& sInputFileName, 
		const std::string& sOutputFileName,
		unsigned int nThreads,
		bool bParallel = false);

   	/* distruttore: esegue tutti i distruttori e libera la memoria */
   	virtual ~Solver(void);

	virtual bool Prepare(void);
	virtual bool Start(void);
	virtual bool Advance(void);

   	/* esegue la simulazione */
   	void Run(void);

	std::ostream & Restart(std::ostream& out, DataManager::eRestart type) const;

	/* EXPERIMENTAL */
	/* FIXME: better const'ify? */
	virtual DataManager *pGetDataManager(void) const {
		return pDM;
	};
	virtual SolutionManager *pGetSolutionManager(void) const {
		return pSM;
	};
	virtual const LinSol& GetLinearSolver(void) const {
		return CurrLinearSolver;
	};
	virtual NonlinearSolver *pGetNonlinearSolver(void) const {
		return pNLS;
	};
	virtual StepIntegrator *pGetStepIntegrator(void) const {
		// if (pCurrStepIntegrator == NULL)
		// 	return pRegularSteps;
		ASSERT(pCurrStepIntegrator != 0);

		return pCurrStepIntegrator;
	};
	virtual doublereal dGetFinalTime(void) const {
		return dFinalTime;
	};
	virtual doublereal dGetInitialTimeStep(void) const {
		return dInitialTimeStep;
	};

	virtual clock_t GetCPUTime(void) const;

	virtual void PrintResidual(const VectorHandler& Res, integer iIterCnt) const;
	virtual void PrintSolution(const VectorHandler& Sol, integer iIterCnt) const;
	virtual void CheckTimeStepLimit(doublereal dErr, doublereal dErrDiff) const throw(NonlinearSolver::TimeStepLimitExceeded, NonlinearSolver::MaxResidualExceeded);
};

inline void
Solver::Flip(void)
{
	/*
	 * switcha i puntatori; in questo modo non e' necessario
	 * copiare i vettori per cambiare passo
	 */
	qX.push_front(qX.back()); 
	qX.pop_back();
	qXPrime.push_front(qXPrime.back());
	qXPrime.pop_back();

	/* copy from pX, pXPrime to qx[0], qxPrime[0] */
	MyVectorHandler* x = qX[0];
	MyVectorHandler* xp = qXPrime[0];
	for (integer i = 1; i <= iNumDofs; i++) {
		x->PutCoef(i, pX->operator()(i));
		xp->PutCoef(i, pXPrime->operator()(i));
	}			
}

/* Solver - end */

#endif /* SOLVER_H */

