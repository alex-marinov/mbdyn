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
 * Umfpack is used by permission; please read its Copyright,
 * License and Availability note.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_UMFPACK
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "umfpackwrap.h"
#include "dgeequ.h"
#include <cstring>

/* in some cases, int and int32_t differ */

#ifdef USE_UMFPACK_LONG

#define UMFPACKWRAP_defaults 		umfpack_dl_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_dl_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_dl_free_numeric

#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_dl_symbolic(size, size, app, aip, axp, sym, ctrl, info)

#define UMFPACKWRAP_report_info 	umfpack_dl_report_info
#define UMFPACKWRAP_report_status 	umfpack_dl_report_status
#define UMFPACKWRAP_numeric 		umfpack_dl_numeric
#define UMFPACKWRAP_solve 		umfpack_dl_solve
#define UMFPACKWRAP_report_numeric      umfpack_dl_report_numeric
#define UMFPACKWRAP_report_symbolic     umfpack_dl_report_symbolic
#else // ! USE_UMFPACK_LONG

#define UMFPACKWRAP_defaults 		umfpack_di_defaults
#define UMFPACKWRAP_free_symbolic 	umfpack_di_free_symbolic
#define UMFPACKWRAP_free_numeric 	umfpack_di_free_numeric

#define UMFPACKWRAP_symbolic(size, app, aip, axp, sym, ctrl, info) \
	umfpack_di_symbolic(size, size, app, aip, axp, sym, ctrl, info)

#define UMFPACKWRAP_report_info 	umfpack_di_report_info
#define UMFPACKWRAP_report_status 	umfpack_di_report_status
#define UMFPACKWRAP_numeric 		umfpack_di_numeric
#define UMFPACKWRAP_solve 		umfpack_di_solve
#define UMFPACKWRAP_report_numeric      umfpack_di_report_numeric
#define UMFPACKWRAP_report_symbolic     umfpack_di_report_symbolic
#endif // ! USE_UMFPACK_LONG

/* required factorization type (A * x = b) */
#define SYS_VALUE 			UMFPACK_A
#define SYS_VALUET 			UMFPACK_Aat

#ifdef MBDYN_ENABLE_PROFILE
#include <chrono>
#endif
#include "cscmhtpl.h"

/* UmfpackSolver - begin */
	
UmfpackSolver::UmfpackSolver(const integer &size,
			     const doublereal &dPivot,
			     const doublereal &dDropTolerance,
			     const unsigned blockSize,
			     Scale scale,
			     integer iMaxIter,
			     integer iVerbose)
: LinearSolver(nullptr),
iSize(size),
Axp(nullptr),
Aip(nullptr),
App(nullptr),
iNumNonZeros(-1),
Symbolic(nullptr),
Numeric(nullptr),
bHaveCond(false)
{
	// silence static analyzers
	memset(&Info[0], 0, sizeof(Info));
	UMFPACKWRAP_defaults(Control);

	Control[UMFPACK_PRL] = iVerbose;
        
	if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
		/*
		 * 1.0: true partial pivoting
		 * 0.0: treated as 1.0
		 * 
		 * default: 0.1
		 */
		Control[UMFPACK_PIVOT_TOLERANCE] = dPivot;
	}

	if (dDropTolerance != 0.) {
#ifdef UMFPACK_DROPTOL
		ASSERT(dDropTolerance > 0.);
		Control[UMFPACK_DROPTOL] = dDropTolerance;
#endif
	}

	if (blockSize > 0) {
		Control[UMFPACK_BLOCK_SIZE] = blockSize;
	}

	if (scale != SCALE_UNDEF) {
		Control[UMFPACK_SCALE] = scale;
	}

	if (iMaxIter >= 0) {
		Control[UMFPACK_IRSTEP] = iMaxIter;
	}
}

UmfpackSolver::~UmfpackSolver(void)
{
	UMFPACKWRAP_free_symbolic(&Symbolic);
	ASSERT(Symbolic == nullptr);

	UMFPACKWRAP_free_numeric(&Numeric);
	ASSERT(Numeric == nullptr);
}

void UmfpackSolver::ResetNumeric(void) const
{
	bHaveCond = false;

	if (Numeric) {
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == nullptr);
	}

	bHasBeenReset = true;
}

void UmfpackSolver::ResetSymbolic(void) const
{
        if (Symbolic) {
                UMFPACKWRAP_free_symbolic(&Symbolic);
                ASSERT(Symbolic == nullptr);
        }

        ResetNumeric();
}

void
UmfpackSolver::Reset(void)
{
        ResetNumeric();
}

void
UmfpackSolver::Solve(void) const
{
	Solve(false);
}

void
UmfpackSolver::SolveT(void) const
{
	Solve(true);
}

void
UmfpackSolver::Solve(bool bTranspose) const
{
	if (bHasBeenReset) {
      		const_cast<UmfpackSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}
		
	int status;

#ifdef MBDYN_ENABLE_PROFILE
	using namespace std::chrono;
        auto start = high_resolution_clock::now();
#endif
	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
	status = UMFPACKWRAP_solve((bTranspose ? SYS_VALUET : SYS_VALUE),
			App, Aip, Axp,
			LinearSolver::pdSol, LinearSolver::pdRhs, 
			Numeric, Control, Info);

	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info);
		UMFPACKWRAP_report_status(Control, status) ;
		silent_cerr("UMFPACKWRAP_solve failed" << std::endl);
		
		/* de-allocate memory */
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == nullptr);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (Control[UMFPACK_IRSTEP] > 0 && Info[UMFPACK_IR_TAKEN] >= Control[UMFPACK_IRSTEP]) {
		silent_cerr("Umfpack warning: UMFPACK_IR_TAKEN = " << Info[UMFPACK_IR_TAKEN]
		             << " >= UMFPACK_IRSTEP = " << Control[UMFPACK_IRSTEP] <<  std::endl
		             << "Umfpack warning: The flag \"max iterations\" of the linear solver is too small or the condition number of the Jacobian matrix is too high" << std::endl);
	}

	UMFPACKWRAP_report_info(Control, Info);

#ifdef MBDYN_ENABLE_PROFILE	
        auto dt = high_resolution_clock::now() - start;
        silent_cout("UMFPACK: solve takes " << dt.count() << "ns\n");
#endif
}

void
UmfpackSolver::Factor(void)
{
	int status;
#ifdef MBDYN_ENABLE_PROFILE
	using namespace std::chrono;
        auto start = high_resolution_clock::now();
#endif
	/*
	 * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
	 */
		
	if (Symbolic == nullptr && !bPrepareSymbolic()) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

        UMFPACKWRAP_report_symbolic (Symbolic, Control) ;
	UMFPACKWRAP_report_info(Control, Info);

	status = UMFPACKWRAP_numeric(App, Aip, Axp, Symbolic, 
			&Numeric, Control, Info);
	if (status == UMFPACK_ERROR_different_pattern) {
		UMFPACKWRAP_free_symbolic(&Symbolic);
		if (!bPrepareSymbolic()) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		status = UMFPACKWRAP_numeric(App, Aip, Axp, Symbolic, 
				&Numeric, Control, Info);
	}

        if (status == UMFPACK_WARNING_singular_matrix) {
            throw ErrNoPivot(-1, MBDYN_EXCEPT_ARGS);
        }
        
	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info);
		UMFPACKWRAP_report_status(Control, status);
		silent_cerr("UMFPACKWRAP_numeric failed" << std::endl);

		/* de-allocate memory */
		UMFPACKWRAP_free_symbolic(&Symbolic);
		UMFPACKWRAP_free_numeric(&Numeric);
		ASSERT(Numeric == nullptr);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	} else {
		bHaveCond = true;
	}
		
        UMFPACKWRAP_report_numeric (Numeric, Control);
	UMFPACKWRAP_report_info(Control, Info);

#ifdef MBDYN_ENABLE_PROFILE
        auto dt = high_resolution_clock::now() - start;
        silent_cout("UMFPACK: factorization takes " << dt.count() << "ns\n");
#endif
}

void
UmfpackSolver::MakeCompactForm(SparseMatrixHandler& mh,
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

        if (mh.Nz() != iNumNonZeros) {
                ResetSymbolic();
        }

        iNumNonZeros = mh.Nz();
}

bool 
UmfpackSolver::bPrepareSymbolic(void)
{
	int status;

	status = UMFPACKWRAP_symbolic(iSize, App, Aip, Axp,
			&Symbolic, Control, Info);
	if (status != UMFPACK_OK) {
		UMFPACKWRAP_report_info(Control, Info) ;
		UMFPACKWRAP_report_status(Control, status);
		silent_cerr("UMFPACKWRAP_symbolic failed" << std::endl);

		/* de-allocate memory */
		UMFPACKWRAP_free_symbolic(&Symbolic);
		ASSERT(Symbolic == nullptr);

		return false;
	}

	return true;
}

bool UmfpackSolver::bGetConditionNumber(doublereal& dCond)
{
	if (!bHaveCond) {
		return false;
	}

	dCond = 1. / Info[UMFPACK_RCOND];

	return true;
}

/* UmfpackSolver - end */

/* UmfpackSparseSolutionManager - begin */
template <typename SparseMatrixHandlerType>
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::UmfpackSparseSolutionManager(integer Dim,
										    doublereal dPivot,
										    doublereal dDropTolerance,
										    const unsigned blockSize,
										    const ScaleOpt& s,
										    integer iMaxIter,
										    integer iVerbose)
: A(Dim, Dim),
x(Dim),
b(Dim),
xVH(Dim, &x[0]),
bVH(Dim, &b[0]),
scale(s),
pMatScale(nullptr)
{
	UmfpackSolver::Scale uscale = UmfpackSolver::SCALE_UNDEF;

	switch (scale.algorithm) {
	case SCALEA_UNDEF:
		scale.when = SCALEW_NEVER; // Do not scale twice
		break;

	case SCALEA_NONE:
		uscale = UmfpackSolver::SCALE_NONE;
		scale.when = SCALEW_NEVER;
		break;

	case SCALEA_ROW_MAX:
		uscale = UmfpackSolver::SCALE_MAX;
		scale.when = SCALEW_NEVER; // Do not scale twice! Use built in scaling from Umfpack
		break;

	case SCALEA_ROW_SUM:
		uscale = UmfpackSolver::SCALE_SUM;
		scale.when = SCALEW_NEVER; // Do not scale twice! Use built in scaling from Umfpack
		break;

	default:
		// Allocate MatrixScale<T> on demand
		uscale = UmfpackSolver::SCALE_NONE; // Do not scale twice!
	}

	SAFENEWWITHCONSTRUCTOR(pLS, UmfpackSolver,
			       UmfpackSolver(Dim, dPivot, dDropTolerance, blockSize, uscale, iMaxIter, iVerbose));

	pLS->pdSetResVec(&b[0]);
	pLS->pdSetSolVec(&x[0]);
	pLS->SetSolutionManager(this);
}

template <typename SparseMatrixHandlerType>
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::~UmfpackSparseSolutionManager(void) 
{
	if (pMatScale) {
		SAFEDELETE(pMatScale);
		pMatScale = nullptr;
	}
}

template <typename SparseMatrixHandlerType>
void UmfpackSparseSolutionManager<SparseMatrixHandlerType>::MatrReset(void)
{
	pLS->Reset();
}

template <typename SparseMatrixHandlerType>
void UmfpackSparseSolutionManager<SparseMatrixHandlerType>::MatrInitialize(void)
{
        pGetSolver()->ResetSymbolic();
}

template <typename SparseMatrixHandlerType>
void
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::MakeCompressedColumnForm(void)
{

#ifdef MBDYN_ENABLE_PROFILE
	using namespace std::chrono;
	auto start = high_resolution_clock::now();
#endif
	pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
	
#ifdef MBDYN_ENABLE_PROFILE
	auto dt = high_resolution_clock::now() - start;
	silent_cout("UMFPACK: making compact form takes " << dt.count() << "ns\n");
	start = high_resolution_clock::now();
#endif
        CSCMatrixHandlerTpl<doublereal, integer, 0> Acsc(&Ax.front(), &Ai.front(), &Ap.front(), A.iGetNumCols(), A.Nz());
	ScaleMatrixAndRightHandSide(Acsc);
        
#ifdef MBDYN_ENABLE_PROFILE
        dt = high_resolution_clock::now() - start;
        silent_cout("UMFPACK: scaling A takes " << dt.count() << "ns\n");
#endif
}

template <typename SparseMatrixHandlerType>
template <typename MH>
void UmfpackSparseSolutionManager<SparseMatrixHandlerType>::ScaleMatrixAndRightHandSide(MH& mh)
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

template <typename SparseMatrixHandlerType>
template <typename MH>
MatrixScale<MH>& UmfpackSparseSolutionManager<SparseMatrixHandlerType>::GetMatrixScale()
{
	if (pMatScale == nullptr) {
		pMatScale = MatrixScale<MH>::Allocate(scale);
	}

	// Will throw std::bad_cast if the type does not match
	return dynamic_cast<MatrixScale<MH>&>(*pMatScale);
}

template <typename SparseMatrixHandlerType>
void UmfpackSparseSolutionManager<SparseMatrixHandlerType>::ScaleSolution(void)
{
	if (scale.when != SCALEW_NEVER) {
		ASSERT(pMatScale != nullptr);
		// scale solution
		pMatScale->ScaleSolution(xVH);
	}
}

/* Risolve il sistema Fattorizzazione + Backward Substitution */
template <typename SparseMatrixHandlerType>
void
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::Solve(void)
{
	MakeCompressedColumnForm();

	pLS->Solve();

	ScaleSolution();
}

/* Risolve il sistema (trasposto) Fattorizzazione + Backward Substitution */
template <typename SparseMatrixHandlerType>
void
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::SolveT(void)
{
	MakeCompressedColumnForm();

	pLS->SolveT();

	ScaleSolution();
}

/* Rende disponibile l'handler per la matrice */
template <typename SparseMatrixHandlerType>
MatrixHandler*
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
template <typename SparseMatrixHandlerType>
MyVectorHandler*
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::pResHdl(void) const
{
	return &bVH;
}

/* Rende disponibile l'handler per la soluzione */
template <typename SparseMatrixHandlerType>
MyVectorHandler*
UmfpackSparseSolutionManager<SparseMatrixHandlerType>::pSolHdl(void) const
{
	return &xVH;
}

template class UmfpackSparseSolutionManager<SpMapMatrixHandler>;

#ifdef USE_SPARSE_AUTODIFF
template class UmfpackSparseSolutionManager<SpGradientSparseMatrixHandler>;
#endif

/* UmfpackSparseSolutionManager - end */

template <class CC>
UmfpackSparseCCSolutionManager<CC>::UmfpackSparseCCSolutionManager(integer Dim,
								   doublereal dPivot,
								   doublereal dDropTolerance,
								   const unsigned& blockSize,
								   const ScaleOpt& scale,
								   integer iMaxIter,
								   integer iVerbose)
: UmfpackSparseSolutionManager(Dim, dPivot, dDropTolerance, blockSize, scale, iMaxIter, iVerbose),
CCReady(false),
Ac(nullptr)
{
	NO_OP;
}

template <class CC>
UmfpackSparseCCSolutionManager<CC>::~UmfpackSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

		if (Ac == nullptr) {
			SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
		}

		CCReady = true;
	}

	ScaleMatrixAndRightHandSide(*Ac);
}

/* Inizializzatore "speciale" */
template <class CC>
void
UmfpackSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	if (Ac) {
		// If a DirCColMatrixHandler is in use and matrix scaling is enabled
		// an uncaught exception (MatrixHandler::ErrRebuildMatrix) will be thrown
		// if zero entries in the matrix become nonzero.
		// For that reason we have to reinitialize Ac!
		SAFEDELETE(Ac);
		Ac = nullptr;
	}

        UmfpackSparseSolutionManager<SpMapMatrixHandler>::MatrInitialize();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
UmfpackSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &A;
	}

	ASSERT(Ac != nullptr);
	return Ac;
}

template class UmfpackSparseCCSolutionManager<CColMatrixHandler<0> >;
template class UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* UmfpackSparseCCSolutionManager - end */

#endif /* USE_UMFPACK */

