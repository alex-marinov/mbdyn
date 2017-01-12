/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2017
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
 * Copyright (C) 2001-2017
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

/*
 * Taucs is used by permission; please read its Copyright,
 * License and Availability note.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_TAUCS
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "taucswrap.h"


/* TaucsSolver - begin */
	
TaucsSolver::TaucsSolver(const integer &size)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(false),
Factorization(0)
{
	NO_OP;
}

TaucsSolver::~TaucsSolver(void)
{
	if (Factorization) {
		/* de-allocate factorization */
		taucs_linsolve(NULL,&Factorization,0, NULL,NULL,NULL,NULL);
		Factorization = 0;
	}
}

void
TaucsSolver::Reset(void)
{
	if (Factorization) {
		/* de-allocate factorization */
		taucs_linsolve(NULL,&Factorization,0, NULL,NULL,NULL,NULL);
		Factorization = 0;
	}

	bHasBeenReset = true;
}

void
TaucsSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		const_cast<TaucsSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}
		
	int status;
// 	char* solve [] = {
// 		"taucs.factor=false",
// 		"taucs.factor.ordering=colamd",
// 		NULL};
	
	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
	status = taucs_linsolve(&A, &Factorization, 1, 
		LinearSolver::pdSol, 
		LinearSolver::pdRhs, 
		NULL, 
		NULL);
	if (status != TAUCS_SUCCESS) {
		silent_cerr("Taucs back-solve failed" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
TaucsSolver::Factor(void)
{
	int status;

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
// 	A.n = iSize;
// 	A.m = iSize;
// 	A.flags = TAUCS_DOUBLE;
// 	A.colptr = App;
// 	A.rowind = Aip;
// 	A.values.d = Axp;
	
// 	if (Factorization) {
// 		/* de-allocate factorization */
// 		taucs_linsolve(NULL,&Factorization,0, NULL,NULL,NULL,NULL);
// 	}
	//if (Symbolic == 0 && !bPrepareSymbolic()) {
	char * options_factor_prevordering[] = { 
		"taucs.factor.mf=true", 
		"taucs.factor.ordering=colamd",
		"taucs.factor.symbolic=false",
		0};
	char * options_factor_ordering[] = { 
		"taucs.factor.mf=true", 
		"taucs.factor.ordering=colamd",
		"taucs.factor.symbolic=true",
		0};
	char ** opts;
	if (Symbolic) {
		opts = options_factor_prevordering;
	} else {
		opts = options_factor_ordering;
	}
	/* factor */
	status = taucs_linsolve(&A, &Factorization, 0, NULL, NULL, opts, NULL);
	if (status != TAUCS_SUCCESS) {
		if (Symbolic == true) {
			/* try to re-do symbolic analisys */
			Symbolic = false;
			const_cast<TaucsSolver *>(this)->Factor();
		} else {
			/* bail out */
			silent_cerr("Taucs factorization failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	Symbolic = true;
}

void
TaucsSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Ai,
		std::vector<integer>& Ac,
		std::vector<integer>& Ap) const
{
	if (!bHasBeenReset) {
		return;
	}
	
	mh.MakeCompressedColumnForm(Ax, Ai, Ap, 0);

	Axp = &(Ax[0]);
	Aip = &(Ai[0]);
	App = &(Ap[0]);
	
	A.n = iSize;
	A.m = iSize;
	A.flags = TAUCS_DOUBLE;
	A.colptr = App;
	A.rowind = Aip;
	A.values.d = Axp;
}

/* TaucsSolver - end */

/* TaucsSparseSolutionManager - begin */

TaucsSparseSolutionManager::TaucsSparseSolutionManager(integer Dim)
: A(Dim),
x(Dim),
b(Dim),
xVH(Dim, &x[0]),
bVH(Dim, &b[0])
{
	SAFENEWWITHCONSTRUCTOR(pLS, TaucsSolver, TaucsSolver(Dim));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&x[0]);
	pLS->SetSolutionManager(this);
}

TaucsSparseSolutionManager::~TaucsSparseSolutionManager(void) 
{
	NO_OP;
}

void
TaucsSparseSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

void
TaucsSparseSolutionManager::MakeCompressedColumnForm(void)
{
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
void
TaucsSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();
	pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
TaucsSparseSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
TaucsSparseSolutionManager::pResHdl(void) const
{
	return &bVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
TaucsSparseSolutionManager::pSolHdl(void) const
{
	return &xVH;
}

/* TaucsSparseSolutionManager - end */

template <class CC>
TaucsSparseCCSolutionManager<CC>::TaucsSparseCCSolutionManager(integer Dim)
: TaucsSparseSolutionManager(Dim),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
TaucsSparseCCSolutionManager<CC>::~TaucsSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
TaucsSparseCCSolutionManager<CC>::MatrReset(void)
{
	/* FIXME */
	pLS->Reset();
}

template <class CC>
void
TaucsSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

		ASSERT(Ac == 0);

		SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
template <class CC>
void
TaucsSparseCCSolutionManager<CC>::MatrInitialize(void)
{
	SAFEDELETE(Ac);
	Ac = 0;

	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
TaucsSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class TaucsSparseCCSolutionManager<CColMatrixHandler<0> >;
template class TaucsSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* TaucsSparseCCSolutionManager - end */

#endif /* USE_TAUCS */

