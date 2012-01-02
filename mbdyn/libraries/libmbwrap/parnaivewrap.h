/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2012
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
 * Copyright (C) 2004-2012
 *
 * Marco Morandini
 *
 */

/*****************************************************************************
 *                                                                           *
 *                          ParNaive C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/


#ifndef PARNAIVEWRAP_H
#define PARNAIVEWRAP_H

#ifdef USE_NAIVE_MULTITHREAD

#include <atomic_ops.h>

#include <iostream>
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

	integer iSize;

	doublereal dMinPiv;
	mutable std::vector<integer> piv;
	integer **ppril;
	integer *pnril;
	mutable std::vector<doublereal> fwd;
	std::vector<integer>	todo; 
	mutable std::vector<AO_t> row_locks;
	mutable std::vector<AO_TS_t> col_locks;

	NaiveMatrixHandler *A;

	unsigned nThreads;

	enum Op {
		FACTOR,
		SOLVE,
		EXIT
	};

	struct thread_data_t {
		pthread_t		thread;
		ParNaiveSolver		*pSLUS;
		unsigned		threadNumber;
		sem_t			sem;
		int			retval;
	} *thread_data;

	mutable Op	thread_operation;
	mutable unsigned thread_count;
	mutable pthread_mutex_t	thread_mutex;
	mutable pthread_cond_t	thread_cond;

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

	void SetMat(NaiveMatrixHandler *const a);
};

/* ParNaiveSolver - end */



/* ParNaiveSparseSolutionManager - begin */

/*
 * Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione
 */

class ParNaiveSparseSolutionManager : public SolutionManager {
protected:
	mutable NaiveMatrixHandler *A;

	mutable MyVectorHandler VH;
	mutable MyVectorHandler XH;

	/* Fattorizza la matrice (non viene mai chiamato direttamente, 
	 * ma da Solve se la matrice ancora non e' stata fattorizzata) */
	void Factor(void);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const;
#endif /* DEBUG */

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
	MatrixHandler* pMatHdl(void) const;

	/* Rende disponibile l'handler per il termine noto */
	VectorHandler* pResHdl(void) const;

	/* Rende disponibile l'handler per la soluzione */
	VectorHandler* pSolHdl(void) const;
};

/* ParNaiveSparseSolutionManager - end */

/* ParNaiveSparsePermSolutionManager - begin */

class ParNaiveSparsePermSolutionManager: public ParNaiveSparseSolutionManager {
private:
	const doublereal dMinPiv;
	void ComputePermutation();
	void BackPerm();

protected:
	/* See NaiveSparsePermSolutionManager */
	enum {
		PERM_NO,
		PERM_INTERMEDIATE,
		PERM_READY
	} ePermState;
	
	mutable std::vector<integer> perm;
	mutable std::vector<integer> invperm;

	virtual void MatrReset(void);
	
public:
	ParNaiveSparsePermSolutionManager(unsigned nt,
		const integer Dim,
		const doublereal dMP = 1.e-9);
	virtual ~ParNaiveSparsePermSolutionManager(void);

	/* Risolve il sistema Backward Substitution; fattorizza se necessario */
	virtual void Solve(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
};

/* ParNaiveSparsePermSolutionManager - end */

#endif /* USE_NAIVE_MULTITHREAD */

#endif /* PARNAIVEWRAP_H */

