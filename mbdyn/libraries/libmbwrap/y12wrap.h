/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 * il Y12SparseSolutionManager rende disponibile lo spazio per la matrice, 
 * la soluzione ed il termine noto. Inoltre rende disponibili
 * metodi per l'handling della matrice e dei vettori 
 * (inserimento e lettura dei coefficienti). Infine consente la soluzione del 
 * problema lineare mediante fattorizzazione LU e successiva sostituzone.
 * 
 * Uso consigliato:
 * // creare un oggetto Y12SparseSolutionManager:
 *         Y12SparseSolutionManager SM(size, workspacesize, pivotfactor);
 * 
 * // ottenere l'handling alla matrice ed al termine noto
 *         MatrixHandler* pMH = SM.pMatHdl();
 *         VectorHandler* pRH = SM.pResHdl();
 * 
 * // Ogni volta che si desidera riassemblare la matrice:
 * // inizializzare il MatrixHandler
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
 * // Y12SparseSolutionManager::Solve()
 * // semplicemente esegue una nuova sostituzione all'indietro.
 *         SM.Solve();
 * 
 * // Se si desidera modificare la matrice per una nuova soluzione, occorre
 * // inizializzare di nuovo il MatrixHandler, con:
 *         SM.MatrReset();
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

#include <iostream>
#include <vector>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "submat.h"
#include "spmapmh.h"
#include "ls.h"

/* Y12Solver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

class Y12Solver : public LinearSolver  {
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
		ErrFactorization(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
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

	mutable integer iMaxSize;
	mutable integer iCurSize;

	mutable integer *piRow;
	mutable integer *piCol;
	mutable doublereal *pdMat;

	integer *pir, *pic;

	/*
	 * NOTE: Y12 alters the index arrays :(
	 */
	bool bDuplicateIndices;		/* true if need to duplicate indices */
	std::vector<integer> iRow;
	mutable std::vector<integer> iCol;

	mutable integer iN;		/* ordine della matrice */
	mutable integer iNonZeroes; 	/* coeff. non nulli */

	integer *piHA;    		/* vettore di lavoro */
	mutable integer iIFLAG[10];     /* vettore di lavoro */
	doublereal dAFLAG[8];		/* vettore di lavoro */
	doublereal* pdPIVOT;    	/* vettore di lavoro */

	void PutError(integer rc) const; /* scrive l'errore */

	/* Fattorizza la matrice */
	void Factor(void);

public:
	/* Costruttore: si limita ad allocare la memoria */
	Y12Solver(integer iMatOrd, integer iWorkSpaceSize,
			doublereal*  pdTmpRhs, 
			integer iPivotParam, bool bDupInd = false);
	/* Distruttore */
	~Y12Solver(void);

#ifdef DEBUG	
	void IsValid(void) const;
#endif /* DEBUG */

	/* Risolve */
	void Solve(void) const;

	/* Index Form */
	void MakeCompactForm(SparseMatrixHandler& mh,
			std::vector<doublereal>& Ax,
			std::vector<integer>& Ar,
			std::vector<integer>& Ac,
			std::vector<integer>& Ap) const;
};

/* Y12Solver - end */


/* Y12SparseSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class Y12SparseSolutionManager : public SolutionManager {
protected:
	integer iMatSize;		/* ordine della matrice */
	std::vector<integer> iRow;	/* array di interi con
					 * indici di riga di Y12Solver */
	std::vector<integer> iCol;	/* array di interi con
					 * indici di colonna di Y12Solver */
   	std::vector<integer> iColStart;	/* array di interi con
					 * indici di colonna CC */
	std::vector<doublereal> dMat;	/* reali con la matrice */
	std::vector<doublereal> dVec;	/* reali con residuo/soluzione */

	mutable SpMapMatrixHandler MH;	/* sparse MatrixHandler */
	mutable MyVectorHandler VH;	/* puntatore a VectorHandler */

	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
	void Factor(void);

#ifdef DEBUG
	/* Usata per il debug */
	void IsValid(void) const;
#endif /* DEBUG */

	virtual void MakeIndexForm(void);
public:
	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
	Y12SparseSolutionManager(integer iSize,
			integer iWorkSpaceSize = 0,
			const doublereal& dPivotFactor = 1.0,
			bool bDupInd = false);

	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
	~Y12SparseSolutionManager(void);

	/* Inizializza il gestore delle matrici */
	void MatrReset(void);

	/* Risolve il sistema */
	void Solve(void);

	/* Rende disponibile l'handler per la matrice */
	MatrixHandler* pMatHdl(void) const {
		return &MH;
	};

	/* Rende disponibile l'handler per il termine noto */
	VectorHandler* pResHdl(void) const {
		return &VH;
	};

	/* Rende disponibile l'handler per la soluzione (e' lo stesso 
	 * del termine noto, ma concettualmente sono separati) */
	VectorHandler* pSolHdl(void) const {
		return &VH;
	};
};

/* Y12SparseSolutionManager - end */

/* Y12SparseCCSolutionManager - begin */
template <class CC>
class Y12SparseCCSolutionManager: public Y12SparseSolutionManager {
protected:
	bool CCReady;
	CompactSparseMatrixHandler *Ac;

	virtual void MatrReset(void);
	virtual void MakeIndexForm(void);
	
public:
	Y12SparseCCSolutionManager(integer Dim, integer /* unused */ = 0, 
			doublereal dPivot = -1.);
	virtual ~Y12SparseCCSolutionManager(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
	
	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;
};

/* Y12SparseCCSolutionManager - end */

#endif /* USE_Y12 */

#endif /* Y12WRAP_H */

