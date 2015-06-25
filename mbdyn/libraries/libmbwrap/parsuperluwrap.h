/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2015
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
 *                          SuperLU C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/

/*
 * Wrapper for SuperLU sparse smp LU solution 
 * http://www.cs.berkeley.edu/~xiaoye/superlu_mt.html
 */

#ifndef PARSUPERLUWRAP_H
#define PARSUPERLUWRAP_H

#ifdef USE_SUPERLU_MT

#include <iostream>
#include "ac/pthread.h"
#include <vector>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "submat.h"
#include "spmapmh.h"
#include "ls.h"

/* ParSuperLUSolver - begin */

/*
 * Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate
 */

struct ParSuperLUSolverData;

class ParSuperLUSolver : public LinearSolver  {
public:
private:
	mutable integer *Aip;
	mutable integer *App;
	mutable doublereal *Axp;

	mutable integer iN;     	/* ordine della matrice */
	mutable integer iNonZeroes;

	doublereal dPivotFactor;

	mutable bool bFirstSol;		/* true se prima backsubst */
	mutable bool bRegenerateMatrix;	/* true se prima backsubst */

	ParSuperLUSolverData *sld;

	unsigned nThreads;

	enum Op {
		FACTOR,
		EXIT
	};

	struct thread_data_t {
		pthread_t		thread;
		ParSuperLUSolver	*pSLUS;
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
	ParSuperLUSolver(unsigned nt, integer iMatOrd, const doublereal &dPivot);

	/* Distruttore */
	~ParSuperLUSolver(void);

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

/* ParSuperLUSolver - end */



/* ParSuperLUSparseSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class ParSuperLUSparseSolutionManager : public SolutionManager {
protected:
	integer iMatSize;		/* ordine della matrice */
	std::vector<integer> Ai;	/* array di interi con
					 * indici di riga di ParSuperLUSolver */
   	std::vector<integer> Adummy;	/* dummy */
   	std::vector<integer> Ap;	/* array di interi con
						 * indici di colonna CC */
	std::vector<doublereal> Ax;	/* reali con la matrice */
	std::vector<doublereal> xb;	/* reali con residuo/soluzione */

	mutable SpMapMatrixHandler MH;	/* sparse MatrixHandler */
	mutable MyVectorHandler VH;	/* puntatore a VectorHandler */

	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
	void Factor(void);

#ifdef DEBUG
	/* Usata per il debug */
	void IsValid(void) const;
#endif /* DEBUG */

	virtual void MakeCompressedColumnForm(void);
public:
	/* Costruttore: usa e mette a disposizione matrici che gli sono date */
	ParSuperLUSparseSolutionManager(unsigned nt, integer iSize,
			const doublereal& dPivotFactor = 1.0);

	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
	~ParSuperLUSparseSolutionManager(void);

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

/* ParSuperLUSparseSolutionManager - end */

/* ParSuperLUSparseCCSolutionManager - begin */
template <class CC>
class ParSuperLUSparseCCSolutionManager: public ParSuperLUSparseSolutionManager {
protected:
	bool CCReady;
	CompactSparseMatrixHandler *Ac;

	virtual void MatrReset(void);
	virtual void MakeCompressedColumnForm(void);
	
public:
	ParSuperLUSparseCCSolutionManager(unsigned nt, integer Dim,
			const doublereal &dPivot = -1.);
	virtual ~ParSuperLUSparseCCSolutionManager(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
	
	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;
};

/* ParSuperLUSparseCCSolutionManager - end */

#endif /* USE_SUPERLU_MT */

#endif /* PARSUPERLUWRAP_H */

