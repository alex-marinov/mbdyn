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

/*****************************************************************************
 *                                                                           *
 *                            Y12 C++ WRAPPER                                *
 *                                                                           *
 *****************************************************************************/

/*
 * Wrapper for Y12 sparse LU solution 
 * http://www.netlib.org/y12/
 */

/*
 * Uso:
 * il Y12SparseLUSolutionManager rende disponibile lo spazio per la matrice, 
 * la soluzione ed il termine noto. Inoltre rende disponibili
 * metodi per l'handling della matrice e dei vettori 
 * (inserimento e lettura dei coefficienti). Infine consente la soluzione del 
 * problema lineare mediante fattorizzazione LU e successiva sostituzone.
 * 
 * Uso consigliato:
 * // creare un oggetto Y12SparseLUSolutionManager:
 *         Y12SparseLUSolutionManager SM(size, workspacesize, pivotfactor);
 * 
 * // ottenere l'handling alla matrice ed al termine noto
 *         MatrixHandler* pMH = SM.pMatHdl();
 *         VectorHandler* pRH = SM.pResHdl();
 * 
 * // Ogni volta che si desidera riassemblare la matrice:
 * // inizializzare il MatrixHandler
 *         SM.MatrInit();
 * 
 * // Operare sulla matrice e sui vettori
 *         pMH->PutCoef(row, col, coef);
 *         coef = pMH->dGetCoef(row, col);
 *         pRH->PutCoef(row, coef);
 *         coef = pRH->dGetCoef(row);
 * 
 * // Risolvere il problema; in questa fase si assume che: 
 * //  - la matrice sia stata completata;
 * //  - il termine noto sia stato completato. 
 * // In tale caso viene eseguita la fattorizzazione e quindi la soluzione.
 * // Se la fattorizzazione e' gia' stata compiuta, la matrice non e' piu'
 * // stata modificata e si e' modificato il termine noto, la chiamata a 
 * // Y12SparseLUSolutionManager::Solve()
 * // semplicemente esegue una nuova sostituzione all'indietro.
 *         SM.Solve();
 * 
 * // Se si desidera modificare la matrice per una nuova soluzione, occorre
 * // inizializzare di nuovo il MatrixHandler, con:
 *         SM.MatrInit();
 * 
 * // Per i parametri di inizializzazione e per eventuali modifiche
 * // fare riferimento al codice sorgente ed alla libreria originaria
 * // in FORTRAN od al suo porting in C, files <harwlib.f>, <harwlib.c>
 * // Il pacchetto e' costituito dai files:
 * 
 *         <y12lib.f>  libreria originaria in FORTRAN
 *         <y12lib.h>  header per <y12lib> con dichiarazione delle funzioni 
 *                      e dei common
 *         <y12wrap.cc> wrapper in C++
 *         <y12wrap.h>  header per <y12wrap.cc>. 
 *                      E' sufficiente includere questo per usare il wrapper.
 */

#ifndef Y12WRAP_H
#define Y12WRAP_H

#ifdef USE_Y12

#include <ac/iostream>
#include <vector>

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>
#include <solman.h>
#include <submat.h>
#include <spmapmh.h>

/* classi dichiarate in questo file */
class Y12LUSolver;      	/* solutore */
class Y12SparseLUSolutionManager;  /* gestore della soluzione */


/* Y12LUSolver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

class Y12LUSolver {
   	friend class Y12SparseLUSolutionManager;


public:
   	class ErrFactorization {
    	private: 
      		int iErrCode;
    	public:
      		ErrFactorization(int i) : iErrCode(i) {
			NO_OP;
		};
      		int iGetErrCode(void) const {
			return iErrCode;
		};
   	};
   
private:

 	/*
	 * indices of array iIFLAG
	 */
	enum {
		I_1 = 0,
		I_2 = 1,
		I_3 = 2,
		I_4 = 3,
		I_5 = 4,
		I_6 = 5,
		I_7 = 6,
		I_8 = 7,
		I_9 = 8,
		I10 = 9
	};

   	integer iMatSize;
   	integer iCurSize;

   	std::vector<integer>*const piRow;
   	std::vector<integer>*const piCol;
   	std::vector<doublereal>*const pdMat;
   
   	integer iN;         	/* ordine della matrice */
   	integer iNonZeroes; 	/* coeff. non nulli */
   	doublereal* pdRhs;  /* Soluzione e termine noto */

   	integer* piHA;    	/* vettore di lavoro */
   	integer  iIFLAG[10];    /* vettore di lavoro */
	doublereal dAFLAG[8];	/* vettore di lavoro */
   	doublereal* pdPIVOT;    /* vettore di lavoro */

	integer iFirstSol;      /* 1 se prima backsubst, else 0 */

	void PutError(std::ostream& out, int rc) const; /* scrive l'errore */

	bool SetCurSize(integer i);
	integer iGetCurSize(void) const;
   
protected:
   	/* Costruttore: si limita ad allocare la memoria */
   	Y12LUSolver(integer iMatOrd, integer iSize,
		    std::vector<integer>*const piTmpRow, 
		    std::vector<integer>*const piTmpCol, 
		    std::vector<doublereal>*const  pdTmpMat,
		    doublereal*  pdTmpRhs, 
		    integer iPivotParam);
   	/* Distruttore */
   	~Y12LUSolver(void);

#ifdef DEBUG	
   	void IsValid(void) const;
#endif /* DEBUG */
   	
	/* Fattorizza la matrice */
   	bool bLUFactor(void);
   
   	/* Risolve */
   	void Solve(void);
};

/* Y12LUSolver - end */


/* Y12SparseLUSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class Y12SparseLUSolutionManager : public SolutionManager {
public: 
   	class ErrGeneric {};
   
private:
   
protected:
   	integer iMatMaxSize;  /* Dimensione max della matrice (per resize) */
   	integer iMatSize;     /* ordine della matrice */
   	std::vector<integer> iRow;       /* array di interi con:
			       * tabella di SparseData/indici di riga
			       * di Y12LUSolver */
   	std::vector<integer> iCol;       /* array di interi con:
	                       * keys di SparseData/indici di colonna
			       * di Y12LUSolver */
   	std::vector<doublereal> dMat;    /* reali con la matrice */
   	std::vector<doublereal> dVec;    /* reali con residuo/soluzione */
   
	mutable SpMapMatrixHandler MH; /* sparse MatrixHandler */
   	VectorHandler* pVH;   /* puntatore a VectorHandler */
   	Y12LUSolver* pLU;     /* puntatore a Y12LUSolver */
   
   	bool bHasBeenReset;   /* flag di matrice resettata */

   	/* Prepara i vettori e la matrice per il solutore */
   	void PacVec(void); 
   
   	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
    	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
   	void Factor(void);

#ifdef DEBUG
   	/* Usata per il debug */
   	void IsValid(void) const;
#endif /* DEBUG */
 
public:
   	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
   	Y12SparseLUSolutionManager(integer iSize,
				   integer iWorkSpaceSize = 0,
				   const doublereal& dPivotFactor = 1.0);
   
   	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
   	~Y12SparseLUSolutionManager(void);

   	/* Inizializza il gestore delle matrici */
   	void MatrInit(const doublereal& dResetVal = 0.);
   
   	/* Risolve il sistema */
   	void Solve(void);
   
   	/* sposta il puntatore al vettore del residuo */
   	void ChangeResPoint(doublereal* pRes){
		pLU->pdRhs = pRes;
	};
   
   	/* sposta il puntatore al vettore della soluzione */
   	void ChangeSolPoint(doublereal* pSol){
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

/* Y12SparseLUSolutionManager - end */

#endif /* USE_Y12 */

#endif /* Y12WRAP_H */

