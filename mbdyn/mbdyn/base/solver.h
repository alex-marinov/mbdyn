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
#include <ac/float.h>
#include <ac/math.h>


#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <dataman.h>
#ifdef USE_MPI
#include <schurdataman.h>
#include <schsolman.h>
#endif /* USE_MPI */
#include <deque>
#include <integr.h>
#include <stepsol.h>
#include <nonlin.h>
#include <mfree.h>
#include <precond.h>

class Solver : public SolverDiagnostics {
public:
	class ErrGeneric {};
	class ErrMaxIterations{};
	class SimulationDiverged{};

private:
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
   	} StrategyFactor;

   	/* Dati per strategia DRIVER_CHANGE */
	DriveCaller* pStrategyChangeDrive;
 
#ifdef __HACK_EIG__
   	/* Dati per esecuzione di eigenanalysis */
   	struct WhenEigen {
      		doublereal dTime;
      		flag fDone;
   	} OneEig;
	flag fEigenAnalysis;
	doublereal dEigParam;
	flag fOutputModes;
	doublereal dUpperFreq;
	doublereal dLowerFreq;
	void Eig(void);
#endif /* __HACK_EIG__ */

#ifdef USE_RTAI
	bool bRT;
	bool bRTWaitPeriod;
	bool bRTHard;
	long long lRTPeriod;
	unsigned long RTStackSize;
#endif /* USE_RTAI */

#ifdef __HACK_POD__
       struct PODData {
               doublereal dTime;
               unsigned int iSteps;
               unsigned int iFrames;
	       PODData(void) : dTime(0.), iSteps(0), iFrames(0) {};
       } pod;

       flag fPOD;
       unsigned int iPODStep;
       unsigned int iPODFrames;
#endif /*__HACK_POD__*/


   	/* Strutture di gestione dei dati */
	integer iNumPreviousVectors;
	doublereal* pdWorkSpace;
	std::deque<MyVectorHandler*> qX;      /* queque di vettori degli stati */
  	std::deque<MyVectorHandler*> qXPrime; /* queque di vettori degli stati derivati */ 

	DataManager* pDM;		/* gestore dei dati */
   	
   	integer iNumDofs;     		/* Dimensioni del problema */



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
   	flag fAbortAfterInput;
   	flag fAbortAfterAssembly;
   	flag fAbortAfterDerivatives;
   	flag fAbortAfterFictitiousSteps;

	/* Parametri per la variazione passo */
   	integer iStepsAfterReduction;
   	integer iStepsAfterRaise;
	integer iWeightedPerformedIters;
   	flag fLastChance;

   	/* Parametri per il metodo */
   	StepIntegrator* pDerivativeSteps;
   	StepIntegrator* pFirstRegularStep;
   	StepIntegrator* pRegularSteps;
   	StepIntegrator* pFirstFictitiousStep;
	StepIntegrator* pFictitiousSteps;
	
	LinSol::SolverType CurrSolver;
	integer iWorkSpaceSize;
   	doublereal dPivotFactor;
   	/* Parametri per Newton-Raphson modificato */
   	flag fTrueNewtonRaphson;
   	integer iIterationsBeforeAssembly;
	flag fMatrixFree;
	MatrixFreeSolver::SolverType MFSolverType;
	doublereal dIterTol;
	Preconditioner::PrecondType PcType;
	integer iPrecondSteps;
	integer iIterativeMaxSteps;
	doublereal dIterertiveEtaMax;

#ifdef USE_MPI
	flag fParallel;
	SchurDataManager *pSDM;
	LinSol::SolverType CurrIntSolver;
	integer iNumLocDofs;		/* Dimensioni problema locale */
	integer iNumIntDofs;		/* Dimensioni interfaccia locale */
	integer* pLocDofs;		/* Lista dof locali (stile FORTRAN) */
	integer* pIntDofs;		/* Lista dei dofs di interfaccia */
	Dof* pDofs;
	integer iIWorkSpaceSize;
	doublereal dIPivotFactor;
	SolutionManager *pLocalSM;
        SchurSolutionManager *pSSM;
#endif /* USE_MPI */  

	/* il solution manager v*/
	SolutionManager *pSM;
	SolutionManager *pCurrSM;

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
   

public:   
   	/* costruttore */
   	Solver(MBDynParser& HP, 
		       	    const char* sInputFileName, 
			    const char* sOutputFileName,
			    flag fParallel = 0);

   	/* esegue la simulazione */
   	void Run(void);

   	/* distruttore: esegue tutti i distruttori e libera la memoria */
   	~Solver(void);
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
}

/* Solver - end */

#endif /* SOLVER_H */




