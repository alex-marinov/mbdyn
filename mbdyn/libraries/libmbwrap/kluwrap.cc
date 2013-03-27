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
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "kluwrap.h"
#include "dgeequ.h"

/* KLUSolver - begin */
	
KLUSolver::KLUSolver(const integer &size, const doublereal &dPivot, Scale scale)
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

	const bool bOK = klu_condest(App, Axp, Symbolic, Numeric, &Control);

	dCond = Control.condest;

	return bOK;
}

/* KLUSolver - end */

/* KLUSparseSolutionManager - begin */

KLUSparseSolutionManager::KLUSparseSolutionManager(integer Dim,
		doublereal dPivot, KLUSolver::Scale scale, ScaleWhen ms)
: A(Dim),
b(Dim),
bVH(Dim, &b[0]),
ms(ms)
{
	SAFENEWWITHCONSTRUCTOR(pLS, KLUSolver,
			KLUSolver(Dim, dPivot, scale));

	(void)pLS->pdSetResVec(&b[0]);
	(void)pLS->pdSetSolVec(&b[0]);
	pLS->SetSolutionManager(this);

	if (scale != KLUSolver::SCALE_NONE && ms != NEVER) {
		ASSERT(0); // avoid double scaling
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
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
	ScaleMatrixAndRightHandSide<SpMapMatrixHandler>(A);

	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

template <class MH>
void KLUSparseSolutionManager::ScaleMatrixAndRightHandSide(MH& mh)
{
#if defined(DEBUG)
	static bool bPrintMat = true;
#endif

	if (ms != SolutionManager::NEVER) {
		if (pLS->bReset()) {
			if (msr.empty() || ms == SolutionManager::ALWAYS) {
#if defined(DEBUG)
				if (bPrintMat) {
					silent_cout("% matrix before scaling:\n");
					silent_cout("Apre=[");
					for (typename MH::const_iterator i = mh.begin(); i != mh.end(); ++i) {
						if (i->dCoef != 0.) {
							silent_cout(i->iRow + 1 << ", " << i->iCol + 1 << ", " << i->dCoef << ";\n");
						}
					}
					silent_cout("];\n");
				}
#endif
				// (re)compute
				doublereal rowcnd = -1., colcnd = -1., amax = -1.;
				dgeequ<MH>(mh, msr, msc, rowcnd, colcnd, amax);
#if defined(DEBUG)
				if (bPrintMat) {
					silent_cout("% scale factors:\n");
					silent_cout("msr=[");
					for (std::vector<doublereal>::const_iterator i = msr.begin(); i < msr.end(); ++i) {
						silent_cout(*i << ";\n");
					}

					silent_cout("];\n");
					silent_cout("msc=[");
					for (std::vector<doublereal>::const_iterator i = msc.begin(); i < msc.end(); ++i) {
						silent_cout(*i << ";\n");
					}
					silent_cout("];\n");
				}
#endif
			}
			// in any case scale matrix and right-hand-side
			dgeequ_scale<MH>(mh, msr, msc);

#if defined(DEBUG)
			if (bPrintMat) {
				silent_cout("% matrix after scaling:\n");
				silent_cout("Apost=[");
				for (typename MH::const_iterator i = mh.begin(); i != mh.end(); ++i) {
					if (i->dCoef != 0.) {
						silent_cout(i->iRow + 1 << ", " << i->iCol + 1 << ", " << i->dCoef << ";\n");
					}
				}
				silent_cout("];\n");
			}

			bPrintMat = false;
#endif
		}
		dgeequ_scale(bVH, &msr[0]);
	}
}

void KLUSparseSolutionManager::ScaleSolution(void)
{
	if (ms != SolutionManager::NEVER) {
		// scale solution
		dgeequ_scale(bVH, &msc[0]);
	}
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
void
KLUSparseSolutionManager::Solve(void)
{
	MakeCompressedColumnForm();

	pLS->Solve();

	ScaleSolution();
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
		doublereal dPivot, KLUSolver::Scale scale, ScaleWhen ms)
: KLUSparseSolutionManager(Dim, dPivot, scale, ms),
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

	ScaleMatrixAndRightHandSide<CC>(*dynamic_cast<CC *>(Ac));

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

