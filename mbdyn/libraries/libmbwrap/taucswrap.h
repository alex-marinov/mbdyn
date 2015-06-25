/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2015
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
 * Copyright (C) 2001-2015
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
 * Copyright (C) 1996-2015
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
 * Taucs is used by permission; please read its Copyright,
 * License and Availability note.
 */

#ifndef TaucsSparseSolutionManager_hh
#define TaucsSparseSolutionManager_hh

#ifdef USE_TAUCS

#include <iostream>
#include <vector>

/*
 * This is just to silence warnings. 
 * Do not use the taucs
 * complex solver!!!!"
*/
#include <complex>

extern "C" {
#include <taucs.h>
}

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"

	
/* TaucsSolver - begin */

class TaucsSolver: public LinearSolver {
private:
	integer iSize;
	mutable taucs_ccs_matrix A;
	mutable doublereal *Axp;
	mutable integer *Aip;
	mutable integer *App;
	mutable bool Symbolic;

	mutable void *Factorization;

	void Factor(void);

public:
	TaucsSolver(const integer &size);
	~TaucsSolver(void);

	void Reset(void);
	void Solve(void) const;

	void MakeCompactForm(SparseMatrixHandler&,
			std::vector<doublereal>& Ax,
			std::vector<integer>& Ar,
			std::vector<integer>& Ac,
			std::vector<integer>& Ap) const;
};

/* TaucsSolver - end */

/* TaucsSparseSolutionManager - begin */

class TaucsSparseSolutionManager: public SolutionManager {
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
	TaucsSparseSolutionManager(integer Dim);
	virtual ~TaucsSparseSolutionManager(void);
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

/* TaucsSparseSolutionManager - end */

/* TaucsSparseCCSolutionManager - begin */

template <class CC>
class TaucsSparseCCSolutionManager: public TaucsSparseSolutionManager {
protected:
	bool CCReady;
	CompactSparseMatrixHandler *Ac;

	virtual void MatrReset(void);
	virtual void MakeCompressedColumnForm(void);
	
public:
	TaucsSparseCCSolutionManager(integer Dim);
	virtual ~TaucsSparseCCSolutionManager(void);

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);
	
	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const;
};

/* TaucsSparseCCSolutionManager - end */

#endif /* USE_TAUCS */

#endif /* TaucsSparseSolutionManager_hh */

