/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/*****************************************************************************
 *                                                                           *
 *                            HARWLIB C++ WRAPPER                            *
 *                                                                           *
 *****************************************************************************/

/*
 * Wrapper for Harwell sparse LU solution 
 * http://www.netlib.org/harwell/
 */

/*
 * Uso:
 * l'HarwellSparseSolutionManager rende disponibile lo spazio per la matrice, 
 * la soluzione ed il termine noto. Inoltre rende disponibili
 * metodi per l'handling della matrice e dei vettori 
 * (inserimento e lettura dei coefficienti). Infine consente la soluzione del 
 * problema lineare mediante fattorizzazione LU e successiva sostituzone.
 * 
 * Uso consigliato:
 * // creare un oggetto HarwellSparseSolutionManager:
 *         HarwellSparseSolutionManager SM(size, workspacesize, pivotfactor);
 * 
 * // ottenere l'handling alla matrice ed al termine noto
 *         SparseMatrixHandler* pMH = SM.pMatHdl();
 *         VectorHandler* pRH = SM.pResHdl();
 * 
 * // Ogni volta che si desidera riassemblare la matrice:
 * // inizializzare lo SparseMatrixHandler
 *         SM.MatrReset();
 * 
 * // Operare sulla matrice e sui vettori
 *         pMH->PutCoef(row, col, coef);
 *         coef = pMH->operator()(row, col);
 *         pRH->PutCoef(row, coef);
 *         coef = pRH->operator()(row);
 * 
 * // Risolvere il problema; in questa fase si assume che: 
 * //  - la matrice sia stata completata;
 * //  - il termine noto sia stato completato. 
 * // In tale caso viene eseguita la fattorizzazione e quindi la soluzione.
 * // Se la fattorizzazione e' gia' stata compiuta, la matrice non e' piu'
 * // stata modificata e si e' modificato il termine noto, la chiamata a 
 * // HarwellSparseSolutionManager::Solve()
 * // semplicemente esegue una nuova sostituzione all'indietro.
 *         SM.Solve();
 * 
 * // Se si desidera modificare la matrice per una nuova soluzione, occorre
 * // inizializzare di nuovo lo SparseMatrixHandler, con:
 *         SM.MatrReset();
 * 
 * // Per i parametri di inizializzazione e per eventuali modifiche
 * // fare riferimento al codice sorgente ed alla libreria originaria
 * // in FORTRAN od al suo porting in C, files <harwlib.f>, <harwlib.c>
 * // Il pacchetto e' costituito dai files:
 * 
 *         <harwlib.f>  libreria originaria in FORTRAN
 *         <harwlib.c>  conversione in C mediante f2c
 *         <harwlib.h>  header per <harwlib.c> con dichiarazione delle funzioni 
 *                      e dei common
 *         <harwrap.cc> wrapper in C++
 *         <harwrap.h>  header per <harwrap.cc>. 
 *                      E' sufficiente includere questo per usare il wrapper.
 */

#ifndef HARWRAP_H
#define HARWRAP_H

#ifdef USE_HARWELL

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <harwlib.h>
#include <solman.h>
#include <ls.h>
#include <submat.h>
#include <spmapmh.h>

/* classi dichiarate in questo file */
class HarwellSolver;      	/* solutore */
class HarwellSparseSolutionManager;  /* gestore della soluzione */


/* HarwellSolver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

static char sLUClassName[] = "HarwellSolver";

class HarwellSolver {
   	friend class HarwellSparseSolutionManager;


public:
   	class ErrFactorization : public MBDynErrBase {
    	private: 
      		integer iErrCode;
    	public:
      		ErrFactorization(integer i, MBDYN_EXCEPT_ARGS_DECL) : 
			MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU), iErrCode(i) {
			NO_OP;
		};
      		integer iGetErrCode(void) const {
			return iErrCode;
		};
   	};
   
private:
   	integer iMatSize;

   	std::vector<integer>*const piRow;
   	std::vector<integer>*const piCol;
   	std::vector<doublereal>*const pdMat;
   
   	integer iN;         	/* ordine della matrice */
   	integer iNonZeroes; 	/* coeff. non nulli */
   	doublereal* pdRhs;  	/* Soluzione e termine noto */
   
   	doublereal dU;      	/* parametro di pivoting */
   	integer* piKeep;    	/* vettore di lavoro */
   	integer* piW;       	/* vettore di lavoro */
   	doublereal* pdW;    	/* vettore di lavoro */
   
protected:
   	/* Costruttore: si limita ad allocare la memoria */
   	HarwellSolver(integer iMatOrd, integer iSize,
		    std::vector<integer>*const piTmpRow, 
		    std::vector<integer>*const piTmpCol, 
		    std::vector<doublereal>*const  pdTmpMat,
		   	doublereal* pdTmpRhs, doublereal dPivotFact)
	: iMatSize(iSize),
	piRow(piTmpRow),
	piCol(piTmpCol),
	pdMat(pdTmpMat),
     	iN(iMatOrd),
	iNonZeroes(0),
	pdRhs(pdTmpRhs), 
     	dU(dPivotFact),
	piKeep(NULL),
	piW(NULL),
	pdW(NULL) {		
		ASSERT(iMatSize > 0);
		ASSERT(piRow != NULL);
		ASSERT(piCol != NULL);
		ASSERT(pdMat != NULL);   
		ASSERT(*ppiRow != NULL);
		ASSERT(*ppiCol != NULL);
		ASSERT(*ppdMat != NULL);
		ASSERT(pdRhs != NULL);
		ASSERT(iN > 0);
	
#ifdef DEBUG_MEMMANAGER	
		ASSERT(piRow->Size()==iMatSize);
		ASSERT(piCol->Size()==iMatSize);
		ASSERT(pdMat->Size()==iMatSize);
		ASSERT(defaultMemoryManager.fIsValid(pdRhs, 
					iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

		SAFENEWARR(piKeep, integer, 5*iN);
		SAFENEWARR(piW, integer, 8*iN);
		SAFENEWARR(pdW, doublereal, iN);
	
#ifdef DEBUG	
		for (integer iCnt = 0; iCnt < 5*iN; iCnt++) {
	   		piKeep[iCnt] = 0;
		}
		for (integer iCnt = 0; iCnt < 8*iN; iCnt++) {
	   		piW[iCnt] = 0;
		}
		for (integer iCnt = 0; iCnt < 1*iN; iCnt++) {
	   		pdW[iCnt] = 0.;
		}
#endif /* DEBUG */
     	};
   
   	/* Distruttore */
   	~HarwellSolver(void) {
      		if (pdW != NULL) {	     
	 		SAFEDELETEARR(pdW);
      		}
      		if (piW != NULL) {	     
	 		SAFEDELETEARR(piW);
      		}
      		if (piKeep != NULL) {	     
	 		SAFEDELETEARR(piKeep);
      		}
   	};

#ifdef DEBUG	
   	void IsValid(void) const {
      		ASSERT(iMatSize > 0);
      		ASSERT(piRow != NULL);
      		ASSERT(piCol != NULL);
      		ASSERT(pdMat != NULL);   
      		ASSERT(pdRhs != NULL);
      		ASSERT(iN > 0);
      
#ifdef DEBUG_MEMMANAGER	
		ASSERT(piRow->Size()==iMatSize);
		ASSERT(piCol->Size()==iMatSize);
		ASSERT(pdMat->Size()==iMatSize);
      		ASSERT(defaultMemoryManager.fIsValid(pdRhs, 
					iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
      
      		ASSERT(piKeep != NULL);
      		ASSERT(piW != NULL);
      		ASSERT(pdW != NULL);
      
#ifdef DEBUG_MEMMANAGER
      		ASSERT(defaultMemoryManager.fIsBlock(piKeep, 
					5*iN*sizeof(integer)));
      		ASSERT(defaultMemoryManager.fIsBlock(piW, 
					8*iN*sizeof(integer)));
      		ASSERT(defaultMemoryManager.fIsBlock(pdW, 
					1*iN*sizeof(doublereal)));	
#endif /* DEBUG_MEMMANAGER */
   	};
#endif /* DEBUG */
   
   	/* Fattorizza la matrice */
   	bool bLUFactor(void) {
#ifdef DEBUG
      		IsValid();
#endif /* DEBUG */
      
      		ASSERT(iNonZeroes > 0);
      
      		integer iLicn = iMatSize;
      		integer iLirn = iMatSize;
      		integer iFlag = 0;
      
      		DEBUGCOUT("Calling ma28ad_()," << std::endl
			<< "iN         = " << iN << std::endl
			<< "iNonZeroes = " << iNonZeroes << std::endl
			<< "pdMat      = " << *ppdMat << std::endl
			<< "iLicn      = " << iLicn << std::endl
			<< "piRow      = " << *ppiRow << std::endl
			<< "iLirn      = " << iLirn << std::endl
			<< "piCol      = " << *ppiCol << std::endl
			<< "dU         = " << dU << std::endl
			<< "piKeep     = " << piKeep << std::endl
			<< "piW        = " << piW << std::endl
			<< "iFlag      = " << iFlag << std::endl);
      
      		__FC_DECL__(ma28ad)(&iN, &iNonZeroes, &((*pdMat)[0]),
					&iLicn, &((*piRow)[0]),
					&iLirn, &((*piCol)[0]), 
					&dU, piKeep, piW, pdW,
					&iFlag);
      
      		if (iFlag < 0) { 
	 		silent_cerr(sLUClassName 
	   			<< ": error during factorization, code "
				<< iFlag << std::endl);
	 		throw ErrFactorization(iFlag);
      		}
      
		/* FIXME: handle iFlag > 0 ??? */
      		return true;
   	};
   
   	/* Risolve */
   	void Solve(void) {
#ifdef DEBUG
      		IsValid();
#endif /* DEBUG */
      
      		integer iLicn = iMatSize;
      		integer iMtype = 1;
      		__FC_DECL__(ma28cd)(&iN, &((*pdMat)[0]), &iLicn, &((*piCol)[0]),
				    piKeep, pdRhs, pdW, &iMtype);
   	};
};

/* HarwellSolver - end */


/* HarwellSparseSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class HarwellSparseSolutionManager : public SolutionManager {
protected:
   	integer iMatMaxSize;  /* Dimensione max della matrice (per resize) */
   	integer iMatSize;     /* ordine della matrice */
   	std::vector<integer> iRow;       /* array di interi con:
			       * tabella di SparseData/indici di riga
			       * di HarwellSolver */
   	std::vector<integer> iCol;       /* array di interi con:
	                       * keys di SparseData/indici di colonna
			       * di HarwellSolver */
   	std::vector<integer> iColStart;
   	std::vector<doublereal> dMat;    /* array di reali con la matrice */
   	std::vector<doublereal> dVec;    /* array di reali con residuo/soluzione */
   
	mutable SpMapMatrixHandler MH; /* SparseMatrixHandler */
/*   	SparseMatrixHandler* pMH; puntatore a SparseMatrixHandler */
   	VectorHandler* pVH;   /* puntatore a VectorHandler */
   	HarwellSolver* pLU; /* puntatore a HarwellSolver */
   
   	flag fHasBeenReset;   /* flag di matrice resettata */
   
   	/* Prepara i vettori e la matrice per il solutore */
   	void PacVec(void); 
   
   	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
    	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
   	void Factor(void);
   
public:   
   	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
   	HarwellSparseSolutionManager(integer iSize,
				       integer iWorkSpaceSize = 0,
				       const doublereal& dPivotFactor = 1.0);
   
   	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
   	~HarwellSparseSolutionManager(void);
   
   	/* Usata per il debug */
#ifdef DEBUG
   	void IsValid(void) const;
#endif /* DEBUG */
   
   	/* Inizializza il gestore delle matrici */
   	void MatrReset();
   
   	/* Risolve il sistema */
   	void Solve(void);
   
   	/* sposta il puntatore al vettore del residuo */
   	void pdSetResVec(doublereal* pRes){
		pLU->pdRhs = pRes;
	};
   
   	/* sposta il puntatore al vettore del residuo */
   	void pdSetSolVec(doublereal* pSol){
		pLU->pdRhs = pSol;
	};

   	/* Rende disponibile l'handler per la matrice */
   	MatrixHandler* pMatHdl(void) const {
      		return &MH;
   	};
   
   	/* Rende disponibile l'handler per il termine noto */
   	VectorHandler* pResHdl(void) const {
      		ASSERT(pVH != NULL);	
      		return pVH;
   	};

   	/* Rende disponibile l'handler per la soluzione (e' lo stesso 
    	 * del termine noto, ma concettualmente sono separati) */
   	VectorHandler* pSolHdl(void) const {
      		ASSERT(pVH != NULL);	
      		return pVH;
   	};
};

/* HarwellSparseSolutionManager - end */

#endif /* USE_HARWELL */

#endif /* HARWRAP_H */

