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
 * il metodo e' strutturato nel modo seguente:
 * il driver costruisce un oggetto della classe MultiStepIntegrator;
 * questo comporta: la lettura dei dati, la creazione del gestore dei dati
 * e del solution manager, che sono virtualmente indipendenti dal metodo.
 * Quindi vengono preparate tutte le strutture piu' strettamente dipendenti
 * dal metodo. Infine viene avviata la soluzione vera e propria.
 * Al termine, quando il programma si conclude, 
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
#include <dataman.h>
#ifdef USE_MPI
#include <schsolman.h>
#include <schurdataman.h>
#endif /* USE_MPI */
#include <integr.h>
#include <intmeth.h>

#ifdef __HACK_POD__
#include <ac/fstream>
 
class PODMat {
       doublereal* A;
       unsigned long Rows, Cols, Cnt;
       mutable bool bOutput;
       char *sfname;
 
 public:
 
       PODMat(unsigned long rows, unsigned long cols, const char *sfn = NULL)
       : Rows(rows), Cols(cols), Cnt(0), bOutput(false), sfname(NULL) {
               SAFENEWARR(A, doublereal, rows*cols);

	       if (sfn == NULL) {
		       SAFESTRDUP(sfname, "MBDyn.POD");
	       } else {
		       size_t l = strlen(sfn);

		       SAFENEWARR(sfname, char, l+sizeof(".POD"));

		       memcpy(sfname, sfn, l);
		       memcpy(sfname+l, ".POD", sizeof(".POD"));
	       }
       };
 
       ~PODMat(void) {
	       if (!bOutput) {
		       Output();
	       }
               SAFEDELETEARR(A);
	       SAFEDELETEARR(sfname);
       }
 
       void AddTVec(VectorHandler* Vec, unsigned long pos) {
	       doublereal *d = A + pos*Rows;
               for (unsigned long i = 0; i < Rows; i++) {
                       d[i] = Vec->dGetCoef(i+1);
               }
	       Cnt++;
       };
 
       void Output(void) const {
	       /*
	        * Not required
	        */
	       if (Cols == 0) {
		       return;
	       }

	       /*
	        * Simulation ended earlier than expected
	        * -- do not waste CPU, get what's available!
	        */
	       if (Cnt < Cols) {
		       pedantic_cout("warning, only " << Cnt 
				       << " out of " << Cols
				       << " frames computed" << std::endl);
	       }

	       std::ofstream out(sfname);
	       if (!out) {
		       std::cerr << "unable to open file \"" << sfname
			       << "\"" << std::endl;
		       THROW(ErrGeneric());
	       }

#ifdef __HACK_POD_BINARY__
	       /* matrix size is coded at the beginning */
	       out.write((char *)&Rows, sizeof(unsigned long));
	       out.write((char *)&Cnt, sizeof(unsigned long));
	       out.write((char *)A, Cnt*Rows*sizeof(doublereal));
#else /* !__HACK_POD_BINARY__ */
	       doublereal *d = A;
               for (unsigned long i = 0; i < Cnt; i++) {         
                       out << d[0];
                       for (unsigned long j = 1; j < Rows; j++) {
                               out << "  " << d[j];
                       }
                       out << std::endl;
		       d += Rows;
               }
#endif /* __HACK_POD_BINARY__ */
	       bOutput = true;

	       out.close();
       };
};
 
#endif /* __HACK_POD__ */

/* MultiStepIntegrator - begin */

class MultiStepIntegrator : public Integrator {

public:
   	class ErrGeneric {};
   	class ErrMaxIterations {};
   	class ErrSimulationDiverged {};
 
private:
   	enum Strategy {
		NOCHANGE,
		CHANGE,
		FACTOR
	} CurrStrategy;

	Integrator::SolverType CurrSolver;
	
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
#endif /* __HACK_EIG__ */

#ifdef __HACK_POD__
        /* Dati per il cacole delle matrici delle covarianze */
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
   	DataManager* pDM;		/* gestore dei dati */
   	VecIter<Dof> DofIterator; 	/* Iteratore per la struttura dei Dof,
					 * passato da DM */
   	integer iNumDofs;     		/* Dimensioni del problema */
        
        flag fIterative;
	doublereal  dIterTol;          /* tolleranza per la convergenza della
	                                   soluzione di Ax = b con metodi
					   iterativi
					 */
					  
	integer iIterativeMaxSteps;     /* Iterazioni massime per la soluzione 
	                                   di un sistema linear con i metodi
					   iterativi 
					 */
					    
#ifdef USE_MPI
	flag fParallel;

	SolutionManager* pLocalSM;

	SchurDataManager *pSDM;
	SchurSolutionManager *pSSM;

	/* Strutture gestione parallelo */
	integer iNumLocDofs;		/* Dimensioni problema locale */
	integer iNumIntDofs;		/* Dimensioni interfaccia locale */
	integer* pLocDofs;		/* Lista dof locali (stile FORTRAN) */
	integer* pIntDofs;		/* Lista dei dofs di interfaccia */
	Dof* pDofs;

	integer iIWorkSpaceSize;
	doublereal dIPivotFactor;

	Integrator::SolverType CurrIntSolver;
#endif /* USE_MPI */
   
   	/* Dati della simulazione */
   	doublereal dTime;
   	doublereal dInitialTime;
   	doublereal dFinalTime;
   	doublereal dRefTimeStep;
   	doublereal dInitialTimeStep;
   	doublereal dMinimumTimeStep;
   	doublereal dMaxTimeStep;
	doublereal dTol;
#ifdef MBDYN_X_CONVSOL
   	doublereal dSolutionTol;
#endif /* MBDYN_X_CONVSOL */
   	integer iMaxIterations;

   	/* Dati dei passi fittizi di trimmaggio iniziale */
   	integer iFictitiousStepsNumber;
   	doublereal dFictitiousStepsRatio;
   	doublereal dFictitiousStepsRho;
   	doublereal dFictitiousStepsTolerance;
   	integer iFictitiousStepsMaxIterations;

   	/* Dati del passo iniziale di calcolo delle derivate */
   	doublereal dDerivativesTol;
   	doublereal dDerivativesCoef;   
   	integer iDerivativesMaxIterations;
   
   	/* Flags vari */
   	flag fAbortAfterInput;
   	flag fAbortAfterAssembly;
   	flag fAbortAfterDerivatives;
   	flag fAbortAfterFictitiousSteps;

	integer iOutputFlags;
	enum {
		OUTPUT_NONE		= 0x0000,
		OUTPUT_ITERS		= 0x0001,
		OUTPUT_RES		= 0x0002,
		OUTPUT_JAC		= 0x0004,

		OUTPUT_MASK		= 0x00FF
	};

	inline bool outputIters(void) const {
		return (iOutputFlags & OUTPUT_ITERS);
	};
 
	inline bool outputRes(void) const {
		return (iOutputFlags & OUTPUT_RES);
	};
 
	inline bool outputJac(void) const {
		return (iOutputFlags & OUTPUT_JAC);
	};
 
   	/* Parametri per Newton-Raphson modificato */
   	flag fTrueNewtonRaphson;
   	integer iIterationsBeforeAssembly;
   	integer iPerformedIterations;

   	/* Parametri per la variazione passo */
   	integer iStepsAfterReduction;
   	integer iStepsAfterRaise;
	integer iWeightedPerformedIters;
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
#ifdef __HACK_SCALE_RES__
   	doublereal MakeTest(const VectorHandler& Res, const VectorHandler& XP, const VectorHandler& Scale);
#else /* !__HACK_SCALE_RES__ */
   	doublereal MakeTest(const VectorHandler& Res, const VectorHandler& XP);
#endif /* !__HACK_SCALE_RES__ */

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
			    const char* sOutputFileName,
			    flag fParallel = 0,
			    flag fIter = 0);

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

