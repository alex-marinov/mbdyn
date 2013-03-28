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

/*
 * The Naive Solver is copyright (C) 2004 by
 * Paolo Mantegazza <mantegazza@aero.polimi.it>
 */

#ifndef NaiveSolutionManager_hh
#define NaiveSolutionManager_hh

#include <iostream>
#include <vector>

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "naivemh.h"

	
/* NaiveSolver - begin */

class NaiveSolver: public LinearSolver {
private:
	integer iSize;
	doublereal dMinPiv;
	mutable std::vector<integer> piv;
	NaiveMatrixHandler *A;

	void Factor(void) throw(LinearSolver::ErrFactor);

public:
	NaiveSolver(const integer &size, const doublereal &dMP,
			NaiveMatrixHandler *const a = 0);
	~NaiveSolver(void);

	void SetMat(NaiveMatrixHandler *const a);
	void Reset(void);
	void Solve(void) const;
};

/* NaiveSolver - end */

/* NaiveSparseSolutionManager - begin */

class NaiveSparseSolutionManager: public SolutionManager {
protected:
	mutable NaiveMatrixHandler *A;
	mutable MyVectorHandler VH;

	ScaleWhen ms;
	std::vector<doublereal> msr;
	std::vector<doublereal> msc;

	template <class MH>
	void ScaleMatrixAndRightHandSide(MH& mh);
	void ScaleSolution(void);

public:
	NaiveSparseSolutionManager(const integer Dim, const doublereal dMP = 1.e-9,
		SolutionManager::ScaleWhen ms = SolutionManager::NEVER);
	virtual ~NaiveSparseSolutionManager(void);
#ifdef DEBUG
	virtual void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

	/* Inizializzatore generico */
	virtual void MatrReset(void);
	
	/* Risolve il sistema Backward Substitution; fattorizza se necessario */
	virtual void Solve(void);

	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;

	/* Rende disponibile l'handler per il termine noto */
	virtual MyVectorHandler* pResHdl(void) const;

	/* Rende disponibile l'handler per la soluzione */
	virtual MyVectorHandler* pSolHdl(void) const;
};

/* NaiveSparseSolutionManager - end */


/* NaiveSparsePermSolutionManager - begin */

template<class T>
class NaiveSparsePermSolutionManager: public NaiveSparseSolutionManager {
private:
	const doublereal dMinPiv;
	mutable MyVectorHandler TmpH;

	void ComputePermutation(void);
	void BackPerm(void);

protected:
	enum {
		PERM_NO,
		PERM_INTERMEDIATE,
		PERM_READY
	} ePermState;
	
	mutable std::vector<integer> perm;
	mutable std::vector<integer> invperm;

	virtual void MatrReset(void);
	
public:
	NaiveSparsePermSolutionManager(const integer Dim, const doublereal dMP = 1.e-9,
		SolutionManager::ScaleWhen ms = SolutionManager::NEVER);
	virtual ~NaiveSparsePermSolutionManager(void);

	/* Risolve il sistema Backward Substitution; fattorizza se necessario */
	virtual void Solve(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
};

// class NaiveSparseCuthillMcKeePermSolutionManager: public NaiveSparseSolutionManager {
// private:
// 	const doublereal dMinPiv;
// 	mutable MyVectorHandler TmpH;
// 
// 	void ComputePermutation();
// 	void BackPerm();
// 
// protected:
// 	enum {
// 		PERM_NO,
// 		PERM_INTERMEDIATE,
// 		PERM_READY
// 	} ePermState;
// 	
// 	mutable std::vector<integer> perm;
// 	mutable std::vector<integer> invperm;
// 
// 	virtual void MatrReset(void);
// 	
// public:
// 	NaiveSparseCuthillMcKeePermSolutionManager(const integer Dim, const doublereal dMP = 1.e-9);
// 	virtual ~NaiveSparseCuthillMcKeePermSolutionManager(void);
// 
// 	/* Risolve il sistema Backward Substitution; fattorizza se necessario */
// 	virtual void Solve(void);
// 
// 	/* Inizializzatore "speciale" */
// 	virtual void MatrInitialize(void);
// };

/* NaiveSparsePermSolutionManager - end */


// supported permutations
class Colamd_ordering;
// #ifdef USE_BOOST
class rcmk_ordering;
class king_ordering;
class sloan_ordering;
class md_ordering;
// #endif USE_BOOST
// #ifdef USE_METIS
class metis_ordering;
// #endif USE_METIS
// #ifdef HAVE_UMFPACK
class amd_ordering;
// #endif


#endif // NaiveSolutionManager_hh

