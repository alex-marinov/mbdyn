/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2012
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* December 2001 
 * Modified to add a Sparse matrix in row form and to implement methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 2001-2012
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 * Umfpack is used by permission; please read its Copyright,
 * License and Availability note.
 */

#ifndef UmfpackSparseSolutionManager_hh
#define UmfpackSparseSolutionManager_hh

#ifdef USE_UMFPACK

#include <iostream>
#include <vector>

extern "C" {
#include <umfpack.h>
}

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"

	
/* UmfpackSolver - begin */

class UmfpackSolver: public LinearSolver {
public:
	enum Scale {
		SCALE_NONE = UMFPACK_SCALE_NONE,
		SCALE_MAX = UMFPACK_SCALE_MAX,
		SCALE_SUM = UMFPACK_SCALE_SUM,

		SCALE_UNDEF
	};

private:
	integer iSize;
	mutable doublereal *Axp;
	mutable integer *Aip;
	mutable integer *App;

	void *Symbolic;
	mutable doublereal Control[UMFPACK_CONTROL];
	mutable doublereal Info[UMFPACK_INFO];
	mutable void *Numeric;
	bool bHaveCond;

	bool bPrepareSymbolic(void);
	
	void Factor(void);
	void Solve(bool bTranspose) const;

public:
	UmfpackSolver(const integer &size,
		const doublereal &dPivot,
		const doublereal &dDropTolerance,
		const unsigned blockSize,
		Scale scale = SCALE_UNDEF);
	~UmfpackSolver(void);

	void Reset(void);
	void Solve(void) const;
	void SolveT(void) const;

	void MakeCompactForm(SparseMatrixHandler&,
			std::vector<doublereal>& Ax,
			std::vector<integer>& Ar,
			std::vector<integer>& Ac,
			std::vector<integer>& Ap) const;

	virtual bool bGetConditionNumber(doublereal& dCond);
};

/* UmfpackSolver - end */

/* UmfpackSparseSolutionManager - begin */

class UmfpackSparseSolutionManager: public SolutionManager {
protected:
	mutable SpMapMatrixHandler A;

	std::vector<doublereal> x;
	std::vector<doublereal> b;

	mutable MyVectorHandler xVH, bVH;

	std::vector<doublereal> Ax;
	std::vector<integer> Ai;
	std::vector<integer> Adummy;
	std::vector<integer> Ap;

	/* Passa in forma di Compressed Column (callback per solve,
	 * richiesto da SpMap e CC Matrix Handler) */
	virtual void MakeCompressedColumnForm(void);
	
	/* Backward Substitution */
	void BackSub(doublereal t_iniz = 0.);
   
public:
	UmfpackSparseSolutionManager(integer Dim,
		doublereal dPivot = -1.,
		doublereal dDropTolerance = 0.,
		const unsigned blockSize = 0,
		UmfpackSolver::Scale scale = UmfpackSolver::SCALE_UNDEF);
	virtual ~UmfpackSparseSolutionManager(void);
#ifdef DEBUG
	virtual void IsValid(void) const {
		NO_OP;
	};
#endif /* DEBUG */

	/* Inizializzatore generico */
	virtual void MatrReset(void);
	
	/* Risolve il sistema Backward Substitution; fattorizza se necessario */
	virtual void Solve(void);
	virtual void SolveT(void);

	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;

	/* Rende disponibile l'handler per il termine noto */
	virtual MyVectorHandler* pResHdl(void) const;

	/* Rende disponibile l'handler per la soluzione */
	virtual MyVectorHandler* pSolHdl(void) const;
};

/* UmfpackSparseSolutionManager - end */

/* UmfpackSparseCCSolutionManager - begin */

template <class CC>
class UmfpackSparseCCSolutionManager: public UmfpackSparseSolutionManager {
protected:
	bool CCReady;
	CompactSparseMatrixHandler *Ac;

	virtual void MatrReset(void);
	virtual void MakeCompressedColumnForm(void);
	
public:
	UmfpackSparseCCSolutionManager(integer Dim,
		doublereal dPivot = -1.,
		doublereal dDropTolerance = 0.,
		const unsigned& blockSize = 0,
		UmfpackSolver::Scale scale = UmfpackSolver::SCALE_UNDEF);
	virtual ~UmfpackSparseCCSolutionManager(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
	
	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;
};

/* UmfpackSparseCCSolutionManager - end */

#endif /* USE_UMFPACK */

#endif /* UmfpackSparseSolutionManager_hh */

