/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2004
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
 *                          ParNaive C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/


#ifndef PARNAIVEWRAP_H
#define PARNAIVEWRAP_H

#include "ac/iostream"
#include "ac/pthread.h"
#include <vector>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "submat.h"
#include "naivemh.h"
#include "ls.h"

/* ParNaiveSolver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

struct ParNaiveSolverData;

class ParNaiveSolver : public LinearSolver  {
public:
private:
// 	mutable integer *Aip;
// 	mutable integer *App;
// 	mutable doublereal *Axp;

	integer iSize;
// 	mutable integer iN;     	/* ordine della matrice */
// 	mutable integer iNonZeroes;

	doublereal dMinPiv;
	mutable std::vector<integer> piv;
// 	doublereal dPivotFactor;

//	mutable bool bFirstSol;		/* true se prima backsubst */
//	mutable bool bRegenerateMatrix;	/* true se prima backsubst */

	NaiveMatrixHandler *const A;
	ParNaiveSolverData *sld;

	unsigned nThreads;

	enum Op {
		FACTOR,
		EXIT
	};

	struct thread_data_t {
		pthread_t		thread;
		ParNaiveSolver		*pSLUS;
		unsigned		threadNumber;
		sem_t			sem;
		void			*pdgstrf_threadarg;
	} *thread_data;

	Op		thread_operation;
	unsigned	thread_count;
	pthread_mutex_t	thread_mutex;
	pthread_cond_t	thread_cond;

	/* Thread process */
	static void *thread_op(void *arg);

	/* Fattorizza la matrice */
	void Factor(void);

	void EndOfOp(void);

public:
	/* Costruttore: si limita ad allocare la memoria */
	ParNaiveSolver(unsigned nt, const integer &size, 
		const doublereal& dMP,
		NaiveMatrixHandler *const a);

	/* Distruttore */
	~ParNaiveSolver(void);

#ifdef DEBUG	
	void IsValid(void) const;
#endif /* DEBUG */

	/* Risolve */
	void Solve(void) const;

// 	/* Index Form */
// 	void MakeCompactForm(SparseMatrixHandler& mh,
// 			std::vector<doublereal>& Ax,
// 			std::vector<integer>& Ar,
// 			std::vector<integer>& Ac,
// 			std::vector<integer>& Ap) const;
};

/* ParNaiveSolver - end */



/* ParNaiveSparseSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class ParNaiveSparseSolutionManager : public SolutionManager {
public: 
	class ErrGeneric {};

private:

protected:
	mutable NaiveMatrixHandler A;

	mutable MyVectorHandler VH;
	mutable MyVectorHandler XH;

	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
	void Factor(void);

#ifdef DEBUG
	/* Usata per il debug */
	void IsValid(void) const;
#endif /* DEBUG */

// 	virtual void MakeCompressedColumnForm(void);
public:
	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
	ParNaiveSparseSolutionManager(unsigned nt, 
			const integer Dim,
			const doublereal dMP = 1.e-9);

	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
	~ParNaiveSparseSolutionManager(void);

	/* Inizializza il gestore delle matrici */
	void MatrReset(void);

	/* Risolve il sistema */
	void Solve(void);

	/* Rende disponibile l'handler per la matrice */
	MatrixHandler* pMatHdl(void) const {
		return &A;
	};

	/* Rende disponibile l'handler per il termine noto */
	VectorHandler* pResHdl(void) const {
#ifdef DEBUG
		VH.IsValid();
#endif /* DEBUG */
		return &VH;
	};

	/* Rende disponibile l'handler per la soluzione (e' lo stesso 
	 * del termine noto, ma concettualmente sono separati) */
	VectorHandler* pSolHdl(void) const {
#ifdef DEBUG
		VH.IsValid();
#endif /* DEBUG */
		return &VH;
	};
};

/* ParNaiveSparseSolutionManager - end */



#endif /* PARNAIVEWRAP_H */

