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
 * il metodo e' strutturato nl modo seguente:
 * il driver costruisce un oggetto della classe MultiStepIntegrator;
 * questo comporta: la lettura dei dati, la creazione del gestore dei dati
 * e del solution manager, che sono virtualmente indipendenti dal metodo.
 * Quindi vengono preparate tutte le strutture piu' strettamente dipendenti
 * dal metodo. Infine viene avviata la soluzione vera e propria (sempre da
 * parte del costruttore!). Al termine, quando il programma si conclude, 
 * viene invocato il distruttore della classe, che si occupa di invocare 
 * i distruttori di tutti gli oggetti dinamici contenuti. In modo del tutto
 * simile e' possibile scrivere oggetti per metodi di integrazione diversi
 * sfruttando gli stessi DataManager e SolutionManager usati per questo,
 * oppure sostituire uno di essi (tipicamente il SolutionManager, che ora 
 * e' un HarwellSparse) senza alterare la struttura del programma.
 */


#ifndef MULTISTP_H
#define MULTISTP_H


#include <myassert.h>
#include <mynewmem.h>
#include <solman.h>
// #include <fullmh.h>
#include <dataman.h>
#include <integr.h>
#include <intmeth.h>


/* MultiStepIntegrator - begin */

class MultiStepIntegrator : public Integrator {
public:
   	class ErrGeneric {};
   	class ErrMaxIterations {};
   	class ErrSimulationDiverged {};
   
private:   
   	enum Strategy {
		NOCHANGE,
		FACTOR
	} CurrStrategy;
	
   	enum SolverType {
		HARWELL_SOLVER,
		MESCHACH_SOLVER,
		Y12_SOLVER,
                UMFPACK3_SOLVER   
	} CurrSolver;
   
private:
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
#endif /* __HACK_EIG__ */
   
   	/* 
	 * puntatori alle strutture di gestione delle soluzioni ai vari passi
	 */
   	doublereal* pdWorkSpace;
  	VectorHandler* pXCurr;		/* stato al passo corrente */
   	VectorHandler* pXPrimeCurr;  	/* derivata al passo corrente */
   	VectorHandler* pXPrev;       	/* stato al passo precedente */
   	VectorHandler* pXPrimePrev;  	/* derivata al passo precedente */
   	VectorHandler* pXPrev2;      	/* stato due passi prima */
   	VectorHandler* pXPrimePrev2; 	/* derivata due passi prima */
   
   	/* Strutture di gestione dei dati e della soluzione */
   	SolutionManager* pSM; 		/* SolutionManager generico */
   	DataManager* pDM;     		/* gestore dei dati */
   	VecIter<Dof> DofIterator; 	/* Iteratore per la struttura dei Dof,
					 * passato da DM */
   	integer iNumDofs;     		/* Dimensioni del problema */
   
   	/* Dati della simulazione */
   	doublereal dTime;
   	doublereal dInitialTime;
   	doublereal dFinalTime;
   	doublereal dRefTimeStep;
   	doublereal dInitialTimeStep;
   	doublereal dMinimumTimeStep;
   	doublereal dMaxTimeStep;
   	doublereal dToll;
   	integer iMaxIterations;

   	/* Dati dei passi fittizi di trimmaggio iniziale */
   	integer iFictitiousStepsNumber;
   	doublereal dFictitiousStepsRatio;
   	doublereal dFictitiousStepsRho;
   	doublereal dFictitiousStepsTolerance;
   	integer iFictitiousStepsMaxIterations;

   	/* Dati del passo iniziale di calcolo delle derivate */
   	doublereal dDerivativesToll;
   	doublereal dDerivativesCoef;   
   	integer iDerivativesMaxIterations;
   
   	/* Flags vari */
   	flag fAbortAfterInput;
   	flag fAbortAfterAssembly;
   	flag fAbortAfterDerivatives;
   	flag fAbortAfterFictitiousSteps;
   
   	/* Parametri per Newton-Raphson modificato */
   	flag fTrueNewtonRaphson;
   	integer iIterationsBeforeAssembly;
   	integer iPerformedIterations;
   
   	/* Parametri per la variazione passo */
   	integer iStepsAfterReduction;
   	integer iStepsAfterRaise;
   	flag fLastChance;
   
   	/* Parametri per il metodo */
   	MultiStepIntegrationMethod* pMethod;
   	MultiStepIntegrationMethod* pFictitiousStepsMethod;  
   
   	/* Parametri di correzione (globali) */
   	doublereal db0Differential;
   	doublereal db0Algebraic;
   
   	/* Dimensioni del workspace (se 0, su misura per la matrice) */
   	integer iWorkSpaceSize;
   	doublereal dPivotFactor;

   	/* Test sul residuo */
   	doublereal MakeTest(const VectorHandler& Res, const VectorHandler& XP);

   	/* corregge i puntatori per un nuovo passo */
   	inline void Flip(void);
    
   	/* Lettura dati */
   	void ReadData(MBDynParser& HP);
   
   	/* Predizione al primo passo */
   	void FirstStepPredict(MultiStepIntegrationMethod* pM);
   
   	/* Predizione al passo generico */
   	void Predict(MultiStepIntegrationMethod* pM);
   
   	/* Nuovo delta t */
   	doublereal NewTimeStep(doublereal dCurrTimeStep, 
			       integer iPerformedIters,
			       MultiStepIntegrationMethod::StepChange Dmy 
			       = MultiStepIntegrationMethod::NEWSTEP);
   
   	/* Aggiornamento della soluzione nel passo fittizio */
   	void DerivativesUpdate(const VectorHandler& Sol);
   
   	/* Aggiornamento normale */
   	void Update(const VectorHandler& Sol);

#ifdef __HACK_EIG__
   	/* Estrazione Autovalori (sperimentale) */
   	void Eig(void);
#endif /* __HACK_EIG__ */
   
public:   
   	/* costruttore */
   	MultiStepIntegrator(MBDynParser& HP, 
		       	    const char* sInputFileName, 
			    const char* sOutputFileName);

   	/* esegue la simulazione */
   	virtual void Run(void);
   
   	/* distruttore: esegue tutti i distruttori e libera la memoria */
   	~MultiStepIntegrator(void);
};

inline void
MultiStepIntegrator::Flip(void)
{
	/*
	 * switcha i puntatori; in questo modo non e' necessario
	 * copiare i vettori per cambiare passo
	 */
	VectorHandler* p = pXPrev2;
	pXPrev2 = pXPrev;
	pXPrev = pXCurr;
	pXCurr = p;
	p = pXPrimePrev2;
	pXPrimePrev2 = pXPrimePrev;
	pXPrimePrev = pXPrimeCurr;
	pXPrimeCurr = p;
}

/* MultiStepIntegrator - end */

#endif /* MULTISTP_H */

