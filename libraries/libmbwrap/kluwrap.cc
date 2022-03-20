/* $Header$ */
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
#include <cassert>
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "kluwrap.h"
#include "dgeequ.h"
#include "cscmhtpl.h"

/* KLUSolver - begin */
	
KLUSolver::KLUSolver(const integer &size, const doublereal &dPivot, Scale scale)
: LinearSolver(0),
iSize(size),
Axp(0),
Aip(0),
App(0),
Symbolic(0),
Numeric(0),
iNumNonZeros(-1)
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

	switch (scale) {
	case SCALE_UNDEF:
		// use the default value provided by KLU
		break;
	default:
		Control.scale = scale;
	}
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
KLUSolver::ResetSymbolic(void) const
{
	if (Symbolic) {
		klu_free_symbolic(&Symbolic, &Control);
		ASSERT(Symbolic == 0);
	}

        ResetNumeric();
}

void
KLUSolver::ResetNumeric(void) const
{
	if (Numeric) {
		klu_free_numeric(&Numeric, &Control);
		ASSERT(Numeric == 0);
	}

	bHasBeenReset = true;
}

void
KLUSolver::Reset(void)
{
        ResetNumeric();
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

        if (iNumNonZeros != mh.Nz()) {
                ResetSymbolic();
        }
        
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

	const bool bOK = klu_condest(App, Axp, Symbolic, Numeric, &Control);

	dCond = Control.condest;

	return bOK;
}

/* KLUSolver - end */

/* KLUSparseSolutionManager - begin */
template <typename MatrixHandlerType>
KLUSparseSolutionManager<MatrixHandlerType>::KLUSparseSolutionManager(integer Dim,
								      doublereal dPivot,
								      const ScaleOpt& s)
     : A(Dim, Dim),
b(Dim),
bVH(Dim, &b[0]),
scale(s),
pMatScale(0)
{
	KLUSolver::Scale kscale = KLUSolver::SCALE_UNDEF;

	switch (scale.algorithm) {
	case SCALEA_UNDEF:
		scale.when = SCALEW_NEVER;
		break;

	case SCALEA_NONE:
		kscale = KLUSolver::SCALE_NONE;
		scale.when = SCALEW_NEVER;
		break;

	case SCALEA_ROW_MAX:
		kscale = KLUSolver::SCALE_MAX;
		scale.when = SCALEW_NEVER; // Do not scale twice! Use built in scaling from KLU
		break;

	case SCALEA_ROW_SUM:
		kscale = KLUSolver::SCALE_SUM;
		scale.when = SCALEW_NEVER; // Do not scale twice! Use built in scaling from KLU
		break;

	default:
		// Allocate MatrixScale<T> on demand
		kscale = KLUSolver::SCALE_NONE; // Do not scale twice!
	}

	SAFENEWWITHCONSTRUCTOR(pLS, KLUSolver,
			KLUSolver(Dim, dPivot, kscale));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&b[0]);
	pLS->SetSolutionManager(this);
}

template <typename MatrixHandlerType>
KLUSparseSolutionManager<MatrixHandlerType>::~KLUSparseSolutionManager(void)
{
	if (pMatScale) {
		SAFEDELETE(pMatScale);
	}
}

template <typename MatrixHandlerType>
void
KLUSparseSolutionManager<MatrixHandlerType>::MatrReset(void)
{
	pLS->Reset();
}

template <typename MatrixHandlerType>
void
KLUSparseSolutionManager<MatrixHandlerType>::MakeCompressedColumnForm(void)
{
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

        CSCMatrixHandlerTpl<doublereal, integer, 0> Acsc(&Ax.front(), &Ai.front(), &Ap.front(), A.iGetNumCols(), A.Nz());
        
        ScaleMatrixAndRightHandSide(Acsc);
}

template <typename MatrixHandlerType>
template <typename MH>
void KLUSparseSolutionManager<MatrixHandlerType>::ScaleMatrixAndRightHandSide(MH& mh)
{
	if (scale.when != SCALEW_NEVER) {
		MatrixScale<MH>& rMatScale = GetMatrixScale<MH>();

		if (pLS->bReset()) {
			if (!rMatScale.bGetInitialized()
				|| scale.when == SolutionManager::SCALEW_ALWAYS) {
				// (re)compute
				rMatScale.ComputeScaleFactors(mh);
			}
			// in any case scale matrix and right-hand-side
			rMatScale.ScaleMatrix(mh);

			if (silent_err) {
				rMatScale.Report(std::cerr);
			}
		}

		rMatScale.ScaleRightHandSide(bVH);
	}
}

template <typename MatrixHandlerType>
template <typename MH>
MatrixScale<MH>& KLUSparseSolutionManager<MatrixHandlerType>::GetMatrixScale()
{
	if (pMatScale == 0) {
		pMatScale = MatrixScale<MH>::Allocate(scale);
	}

	// Will throw std::bad_cast if the type does not match
	return dynamic_cast<MatrixScale<MH>&>(*pMatScale);
}

template <typename MatrixHandlerType>
void KLUSparseSolutionManager<MatrixHandlerType>::ScaleSolution(void)
{
	if (scale.when != SCALEW_NEVER) {
		ASSERT(pMatScale != 0);
		// scale solution
		pMatScale->ScaleSolution(bVH);
	}
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
template <typename MatrixHandlerType>
void
KLUSparseSolutionManager<MatrixHandlerType>::Solve(void)
{
	MakeCompressedColumnForm();

	pLS->Solve();

	ScaleSolution();
}

/* Rende disponibile l'handler per la matrice */
template <typename MatrixHandlerType>
MatrixHandler*
KLUSparseSolutionManager<MatrixHandlerType>::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
template <typename MatrixHandlerType>
MyVectorHandler*
KLUSparseSolutionManager<MatrixHandlerType>::pResHdl(void) const
{
	return &bVH;
}

/* Rende disponibile l'handler per la soluzione */
template <typename MatrixHandlerType>
MyVectorHandler*
KLUSparseSolutionManager<MatrixHandlerType>::pSolHdl(void) const
{
	return &bVH;
}

template class KLUSparseSolutionManager<SpMapMatrixHandler>;

#ifdef USE_SPARSE_AUTODIFF
template class KLUSparseSolutionManager<SpGradientSparseMatrixHandler>;
#endif

/* KLUSparseSolutionManager - end */

template <class CC>
KLUSparseCCSolutionManager<CC>::KLUSparseCCSolutionManager(integer Dim,
		doublereal dPivot,
		const ScaleOpt& scale)
     : KLUSparseSolutionManager<SpMapMatrixHandler>(Dim, dPivot, scale),
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

	ScaleMatrixAndRightHandSide(*Ac);

}

/* Inizializzatore "speciale" */
template <class CC>
void
KLUSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	if (Ac) {
		// If a DirCColMatrixHandler is in use and matrix scaling is enabled
		// an uncaught exception (MatrixHandler::ErrRebuildMatrix) will be thrown
		// if zero entries in the matrix become nonzero.
		// For that reason we have to reinitialize Ac!
		SAFEDELETE(Ac);
		Ac = 0;
	}

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

