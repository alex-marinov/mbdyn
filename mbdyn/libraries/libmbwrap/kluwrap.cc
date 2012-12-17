/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2006
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
 * Copyright (C) 2001-2005
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
 * Copyright (C) 1996-2005
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_KLU
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "kluwrap.h"


/* KLUSolver - begin */
	
KLUSolver::KLUSolver(const integer &size, const doublereal &dPivot)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(0),
Numeric(0)
{
	klu_defaults(&Control);

	if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
		/*
		 * 1.0: true partial pivoting
		 * 0.0: treated as 1.0
		 * 
		 * default: 0.1
		 */
		Control.tol = dPivot;
	}
	// Control.ordering = 1; /* colamd */
	Control.ordering = 0; /* amd */

}

KLUSolver::~KLUSolver(void)
{
	if (Symbolic) {
		klu_free_symbolic(&Symbolic, &Control);
	}
	ASSERT(Symbolic == 0);
	
	if (Numeric) {
		klu_free_numeric(&Numeric, &Control);
	}
	ASSERT(Numeric == 0);
}

void
KLUSolver::Reset(void)
{
	if (Numeric) {
		//klu_free_numeric(&Numeric, &Control);
		//ASSERT(Numeric == 0);
	}

	bHasBeenReset = true;
}

void
KLUSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		const_cast<KLUSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}
		

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	int ok = klu_solve(Symbolic, Numeric, iSize, 1, LinearSolver::pdRhs, 
			&Control);
	if (!ok || Control.status != KLU_OK) {
		silent_cerr("KLUWRAP_solve failed" << std::endl);
		
		/* de-allocate memory */
		if (Numeric) {
			klu_free_numeric(&Numeric, &Control);
		}
		ASSERT(Numeric == 0);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
}

void
KLUSolver::Factor(void)
{

	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	if (Symbolic == 0 && !bPrepareSymbolic()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (Numeric != 0) {
		int ok;
		ok = klu_refactor(App, Aip, Axp, Symbolic, 
			Numeric, &Control);
		if (!ok) {
			// FIXME: might be too late!
			klu_free_numeric(&Numeric, &Control);
			klu_free_symbolic(&Symbolic, &Control);
			if (!bPrepareSymbolic()) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}
	
	if (Numeric == 0) {
		Numeric = klu_factor(App, Aip, Axp, Symbolic, 
				&Control);
	}
			
	if (Control.status != KLU_OK) {
		silent_cerr("KLUWRAP_numeric failed" << std::endl);

		/* de-allocate memory */
		if (Symbolic) {
			klu_free_symbolic(&Symbolic, &Control);
		}
		if (Numeric) {
			klu_free_numeric(&Numeric, &Control);
		}
		ASSERT(Symbolic == 0);
		ASSERT(Numeric == 0);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
KLUSolver::MakeCompactForm(SparseMatrixHandler& mh,
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
}

bool 
KLUSolver::bPrepareSymbolic(void)
{
	Symbolic = klu_analyze(iSize, App, Aip, &Control);
	if (Control.status != KLU_OK) {
		silent_cerr("KLUWRAP_symbolic failed" << std::endl);

		/* de-allocate memory */
		if (Symbolic) {
			klu_free_symbolic(&Symbolic, &Control);
		}
		ASSERT(Symbolic == 0);

		return false;
	}

	return true;
}

bool KLUSolver::bGetConditionNumber(doublereal& dCond)
{
	if (Symbolic == 0 || Numeric == 0) { // FIXME: Is it really safe?
		return false;
	}

	const bool bOK = klu_rcond(Symbolic, Numeric, &Control);

	dCond = 1. / Control.rcond;

	return bOK;
}

/* KLUSolver - end */

/* KLUSparseSolutionManager - begin */

KLUSparseSolutionManager::KLUSparseSolutionManager(integer Dim,
		doublereal dPivot)
: A(Dim),
b(Dim),
bVH(Dim, &b[0])
{
	SAFENEWWITHCONSTRUCTOR(pLS, KLUSolver,
			KLUSolver(Dim, dPivot));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&b[0]);
	pLS->SetSolutionManager(this);
}

KLUSparseSolutionManager::~KLUSparseSolutionManager(void) 
{
	NO_OP;
}

void
KLUSparseSolutionManager::MatrReset(void)
{
	pLS->Reset();
}

void
KLUSparseSolutionManager::MakeCompressedColumnForm(void)
{
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
KLUSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();
	pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
KLUSparseSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
KLUSparseSolutionManager::pResHdl(void) const
{
	return &bVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
KLUSparseSolutionManager::pSolHdl(void) const
{
	return &bVH;
}

/* KLUSparseSolutionManager - end */

template <class CC>
KLUSparseCCSolutionManager<CC>::KLUSparseCCSolutionManager(integer Dim,
		doublereal dPivot)
: KLUSparseSolutionManager(Dim, dPivot),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
KLUSparseCCSolutionManager<CC>::~KLUSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
KLUSparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

template <class CC>
void
KLUSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

		if (Ac == 0) {
			SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
		}

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
template <class CC>
void
KLUSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
KLUSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class KLUSparseCCSolutionManager<CColMatrixHandler<0> >;
template class KLUSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* KLUSparseCCSolutionManager - end */

#endif /* USE_KLU */

