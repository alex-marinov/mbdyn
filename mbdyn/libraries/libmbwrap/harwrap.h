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
 * l'HarwellSparseLUSolutionManager rende disponibile lo spazio per la matrice, 
 * la soluzione ed il termine noto. Inoltre rende disponibili
 * metodi per l'handling della matrice e dei vettori 
 * (inserimento e lettura dei coefficienti). Infine consente la soluzione del 
 * problema lineare mediante fattorizzazione LU e successiva sostituzone.
 * 
 * Uso consigliato:
 * // creare un oggetto HarwellSparseLUSolutionManager:
 *         HarwellSparseLUSolutionManager SM(size, workspacesize, pivotfactor);
 * 
 * // ottenere l'handling alla matrice ed al termine noto
 *         SparseMatrixHandler* pMH = SM.pMatHdl();
 *         VectorHandler* pRH = SM.pResHdl();
 * 
 * // Ogni volta che si desidera riassemblare la matrice:
 * // inizializzare lo SparseMatrixHandler
 *         SM.MatrInit();
 * 
 * // Operare sulla matrice e sui vettori
 *         pMH->fPutCoef(row, col, coef);
 *         coef = pMH->dGetCoef(row, col);
 *         pRH->fPutCoef(row, coef);
 *         coef = pRH->dGetCoef(row);
 * 
 * // Risolvere il problema; in questa fase si assume che: 
 * //  - la matrice sia stata completata;
 * //  - il termine noto sia stato completato. 
 * // In tale caso viene eseguita la fattorizzazione e quindi la soluzione.
 * // Se la fattorizzazione e' gia' stata compiuta, la matrice non e' piu'
 * // stata modificata e si e' modificato il termine noto, la chiamata a 
 * // HarwellSparseLUSolutionManager::Solve()
 * // semplicemente esegue una nuova sostituzione all'indietro.
 *         SM.Solve();
 * 
 * // Se si desidera modificare la matrice per una nuova soluzione, occorre
 * // inizializzare di nuovo lo SparseMatrixHandler, con:
 *         SM.MatrInit();
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
#include <spdata.h>
#include <harwlib.h>
#include <solman.h>
#include <submat.h>
#include <sparsemh.h>

/* classi dichiarate in questo file */
class HarwellLUSolver;      	/* solutore */
class HarwellSparseLUSolutionManager;  /* gestore della soluzione */


/* HarwellLUSolver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

static char sLUClassName[] = "HarwellLUSolver";

class HarwellLUSolver {
   	friend class HarwellSparseLUSolutionManager;


public:
   	class ErrFactorisation {
    	private: 
      		int iErrCode;
    	public:
      		ErrFactorisation(int i) : iErrCode(i) {
			NO_OP;
		};
      		int iGetErrCode(void) const {
			return iErrCode;
		};
   	};
   
private:
   	integer iMatSize;
   	integer** ppiRow;
   	integer** ppiCol;
   	doublereal** ppdMat;
   
   	integer iN;         	/* ordine della matrice */
   	integer iNonZeroes; 	/* coeff. non nulli */
   	doublereal* pdRhs;  	/* Soluzione e termine noto */
   
   	doublereal dU;      	/* parametro di pivoting */
   	integer* piKeep;    	/* vettore di lavoro */
   	integer* piW;       	/* vettore di lavoro */
   	doublereal* pdW;    	/* vettore di lavoro */
   
protected:
   	/* Costruttore: si limita ad allocare la memoria */
   	HarwellLUSolver(integer iMatOrd, integer iSize,
		   	integer** ppiTmpRow, integer** ppiTmpCol, 
		   	doublereal** ppdTmpMat,
		   	doublereal* pdTmpRhs, doublereal dPivotFact)
	: iMatSize(iSize),
	ppiRow(ppiTmpRow),
	ppiCol(ppiTmpCol),
	ppdMat(ppdTmpMat),
     	iN(iMatOrd),
	iNonZeroes(0),
	pdRhs(pdTmpRhs), 
     	dU(dPivotFact),
	piKeep(NULL),
	piW(NULL),
	pdW(NULL) {		
		ASSERT(iMatSize > 0);
		ASSERT(ppiRow != NULL);
		ASSERT(ppiCol != NULL);
		ASSERT(ppdMat != NULL);   
		ASSERT(*ppiRow != NULL);
		ASSERT(*ppiCol != NULL);
		ASSERT(*ppdMat != NULL);
		ASSERT(pdRhs != NULL);
		ASSERT(iN > 0);
	
#ifdef DEBUG_MEMMANAGER	
		ASSERT(defaultMemoryManager.fIsValid(*ppiRow, 
					iMatSize*sizeof(integer)));
		ASSERT(defaultMemoryManager.fIsValid(*ppiCol, 
					iMatSize*sizeof(integer)));
		ASSERT(defaultMemoryManager.fIsValid(*ppdMat, 
					iMatSize*sizeof(doublereal)));
		ASSERT(defaultMemoryManager.fIsValid(pdRhs, 
					iN*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */

		SAFENEWARR(piKeep, integer, 5*iN);
		SAFENEWARR(piW, integer, 8*iN);
		SAFENEWARR(pdW, doublereal, iN);
	
#ifdef DEBUG	
		for (int iCnt = 0; iCnt < 5*iN; iCnt++) {
	   		piKeep[iCnt] = 0;
		}
		for (int iCnt = 0; iCnt < 8*iN; iCnt++) {
	   		piW[iCnt] = 0;
		}
		for (int iCnt = 0; iCnt < 1*iN; iCnt++) {
	   		pdW[iCnt] = 0.;
		}
#endif /* DEBUG */
     	};
   
   	/* Distruttore */
   	~HarwellLUSolver(void) {
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
   
   	void IsValid(void) const {
      		ASSERT(iMatSize > 0);
      		ASSERT(ppiRow != NULL);
      		ASSERT(ppiCol != NULL);
      		ASSERT(ppdMat != NULL);   
      		ASSERT(*ppiRow != NULL);
     		ASSERT(*ppiCol != NULL);
      		ASSERT(*ppdMat != NULL);
      		ASSERT(pdRhs != NULL);
      		ASSERT(iN > 0);
      
#ifdef DEBUG_MEMMANAGER	
      		ASSERT(defaultMemoryManager.fIsValid(*ppiRow, 
					iMatSize*sizeof(integer)));
      		ASSERT(defaultMemoryManager.fIsValid(*ppiCol, 
					iMatSize*sizeof(integer)));
      		ASSERT(defaultMemoryManager.fIsValid(*ppdMat, 
					iMatSize*sizeof(doublereal)));
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
   
   	/* Fattorizza la matrice */
   	flag fLUFactor(void) {
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
      
      		__FC_DECL__(ma28ad)(&iN, &iNonZeroes, *ppdMat, &iLicn, *ppiRow,
				    &iLirn, *ppiCol, &dU, piKeep, piW, pdW,
				    &iFlag);
      
      		if (iFlag < 0) { 
	 		std::cerr << sLUClassName 
	   			<< ": error during factorization, code "
				<< iFlag << std::endl;	 
	 		THROW(HarwellLUSolver::ErrFactorisation(iFlag));
      		}
      
      		return iFlag;		
   	};
   
   	/* Risolve */
   	void Solve(void) {
#ifdef DEBUG
      		IsValid();
#endif
      
      		integer iLicn = iMatSize;
      		integer iMtype = 1;
      		__FC_DECL__(ma28cd)(&iN, *ppdMat, &iLicn, *ppiCol,
				    piKeep, pdRhs, pdW, &iMtype);
   	};
};

/* HarwellLUSolver - end */


/* HarwellSparseLUSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class HarwellSparseLUSolutionManager : public SolutionManager {
public: 
   	class ErrGeneric {};
   
private:
   
protected:
   	integer iMatMaxSize;  /* Dimensione max della matrice (per resize) */
   	integer iMatSize;     /* ordine della matrice */
   	integer* piRow;       /* puntatore ad array di interi con:
			       * tabella di SparseData/indici di riga
			       * di HarwellLUSolver */
   	integer* piCol;       /* puntatore ad array di interi con:
	                       * keys di SparseData/indici di colonna
			       * di HarwellLUSolver */
   	doublereal* pdMat;    /* puntatore ad array di reali con la matrice */
   	doublereal* pdVec;    /* p. ad array di reali con residuo/soluzione */
   
   	SparseMatrixHandler* pMH; /* puntatore a SparseMatrixHandler */
   	VectorHandler* pVH;   /* puntatore a VectorHandler */
   	HarwellLUSolver* pLU; /* puntatore a HarwellLUSolver */
   
   	flag fHasBeenReset;   /* flag di matrice resettata */
   
   	/* Prepara i vettori e la matrice per il solutore */
   	void PacVec(void); 
   
   	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
    	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
   	void Factor(void);
   
public:   
   	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
   	HarwellSparseLUSolutionManager(integer iSize,
				       integer iWorkSpaceSize = 0,
				       const doublereal& dPivotFactor = 1.0);
   
   	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
   	~HarwellSparseLUSolutionManager(void);
   
   	/* Usata per il debug */
   	void IsValid(void) const;
   
   	/* Ridimensiona le matrici */
   	void MatrResize(integer iNewSize) {
   		THROW(ErrNotImplementedYet());
   	};
   
   	/* Inizializza il gestore delle matrici */
   	void MatrInit(const doublereal& dResetVal = 0.);
   
   	/* Risolve il sistema */
   	void Solve(void);
   
   	/* sposta il puntatore al vettore del residuo */
   	void ChangeResPoint(doublereal* pRes){
		pLU->pdRhs = pRes;
	};
   
   	/* sposta il puntatore al vettore del residuo */
   	void ChangeSolPoint(doublereal* pSol){
		pLU->pdRhs = pSol;
	};

   	/* Rende disponibile l'handler per la matrice */
   	MatrixHandler* pMatHdl(void) const {
      		ASSERT(pMH != NULL);	
      		return pMH;
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

/* HarwellSparseLUSolutionManager - end */

#endif /* USE_HARWELL */

#endif /* HARWRAP_H */

