/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "myassert.h"
#include "ac/sys_sysinfo.h"

#include "parser.h"

#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "harwrap.h"
#include "mschwrap.h"
#include "y12wrap.h"
#include "umfpackwrap.h"
#include "kluwrap.h"
#include "parsuperluwrap.h"
#include "superluwrap.h"
#include "lapackwrap.h"
#include "taucswrap.h"
#include "naivewrap.h"
#include "parnaivewrap.h"


#include "linsol.h"

/* solver data */
const LinSol::solver_t solver[] = {
	{ "Empty", NULL,
		LinSol::EMPTY_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE,
       		-1., -1. },
	{ "Harwell", NULL,
		LinSol::HARWELL_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE,
		-1., -1. },
	{ "Lapack", NULL,
		LinSol::LAPACK_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE,
		-1., -1. },
	{ "Meschach", NULL,
		LinSol::MESCHACH_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE,
		-1., -1. },
	{ "Naive", NULL,
		LinSol::NAIVE_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE |
			LinSol::SOLVER_FLAGS_ALLOWS_COLAMD |
			LinSol::SOLVER_FLAGS_ALLOWS_MMDATA |
			LinSol::SOLVER_FLAGS_ALLOWS_MDAPLUSAT |
			LinSol::SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE |
			LinSol::SOLVER_FLAGS_ALLOWS_KING |
			LinSol::SOLVER_FLAGS_ALLOWS_SLOAN |
			LinSol::SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION |
			LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
#ifdef USE_NAIVE_MULTITHREAD
			LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT |
#endif /* USE_NAIVE_MULTITHREAD */
			0,
		LinSol::SOLVER_FLAGS_NONE,
		1.e-8, -1. },
	{ "SuperLU", NULL, 
		LinSol::SUPERLU_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_COLAMD |
			LinSol::SOLVER_FLAGS_ALLOWS_MMDATA |
			LinSol::SOLVER_FLAGS_ALLOWS_MAP |
			LinSol::SOLVER_FLAGS_ALLOWS_CC |
			LinSol::SOLVER_FLAGS_ALLOWS_DIR |
#ifdef USE_SUPERLU_MT
			LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT |
#endif /* USE_SUPERLU_MT */
			0,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_COLAMD,
		1., -1. },
	{ "Taucs", NULL, 
		LinSol::TAUCS_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP,
       		-1., -1. },
	{ "Umfpack", "umfpack3", 
		LinSol::UMFPACK_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		.1,
#ifdef UMFPACK_DROPTOL
		0.
#else // ! UMFPACK_DROPTOL
		-1.
#endif // ! UMFPACK_DROPTOL
		},
	{ "KLU", NULL, 
		LinSol::KLU_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		.1 },
	{ "Y12", NULL,
		LinSol::Y12_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		-1., -1. },
	{ NULL, NULL, 
		LinSol::EMPTY_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE,
		-1., -1. }
};

/*
 * Default solver
 */
LinSol::SolverType LinSol::defaultSolver = 
#if defined(USE_UMFPACK)
	LinSol::UMFPACK_SOLVER
#elif /* !USE_UMFPACK */ defined(USE_KLU)
	LinSol::KLU_SOLVER
#else /* !USE_KLU */
	// Naive always present
	LinSol::NAIVE_SOLVER
#endif
	;

LinSol::LinSol(void)
: currSolver(LinSol::defaultSolver),
solverFlags(0),
nThreads(1),
iWorkSpaceSize(0),
blockSize(0),
dPivotFactor(-1.),
dDropTolerance(0.),
ms(SolutionManager::NEVER),
scale(SCALE_UNDEF)
{
	NO_OP;
}

LinSol::~LinSol(void)
{
	NO_OP;
}

LinSol::SolverType
LinSol::GetSolver(void) const
{
	return currSolver;
}

bool
LinSol::SetSolver(LinSol::SolverType t, unsigned f)
{
	/* check flags */
	if ((::solver[t].s_flags & f) != f) {
		return false;
	}

	solverFlags = f;
	
	switch (t) {
	case LinSol::UMFPACK_SOLVER:
#ifdef USE_UMFPACK
		currSolver = t;
		return true;
#endif /* USE_UMFPACK */

	case LinSol::KLU_SOLVER:
#ifdef USE_KLU
		currSolver = t;
		return true;
#endif /* USE_UMFPACK */

	case LinSol::SUPERLU_SOLVER:
#ifdef USE_SUPERLU
		currSolver = t;
		return true;
#endif /* USE_SUPERLU */

	case LinSol::LAPACK_SOLVER:
#ifdef USE_LAPACK
		currSolver = t;
		return true;
#endif /* USE_LAPACK */

	case LinSol::TAUCS_SOLVER:
#ifdef USE_TAUCS
		currSolver = t;
		return true;
#endif /* USE_TAUCS */

	case LinSol::Y12_SOLVER:
#ifdef USE_Y12
		currSolver = t;
		return true;
#endif /* USE_Y12 */

	case LinSol::HARWELL_SOLVER:
#ifdef USE_HARWELL
		currSolver = t;
		return true;
#endif /* USE_HARWELL */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		currSolver = t;
		return true;
#endif /* USE_MESCHACH */

		/* else */
		silent_cerr(::solver[t].s_name << " unavailable" << std::endl);
		return false;

	case LinSol::NAIVE_SOLVER:
		currSolver = t;
		return true;

	case LinSol::EMPTY_SOLVER:
		currSolver = t;
		return true;

	default:
		return false;
	}
}

unsigned
LinSol::GetSolverFlags(void) const
{
	return solverFlags;
}

unsigned
LinSol::GetSolverFlags(SolverType t) const
{
	return ::solver[t].s_flags;
}

const char *const
LinSol::GetSolverName(void) const
{
	return ::solver[currSolver].s_name;
}

const char *const
LinSol::GetSolverName(SolverType t) const
{
	return ::solver[t].s_name;
}

bool
LinSol::SetSolverFlags(unsigned f)
{
	if ((::solver[currSolver].s_flags & f) == f) {
		solverFlags = f;
		return true;
	}

	return false;
}

bool
LinSol::AddSolverFlags(unsigned f)
{
	if ((::solver[currSolver].s_flags & f) == f) {
		solverFlags |= f;
		return true;
	}

	return false;
}

bool
LinSol::MaskSolverFlags(unsigned f)
{
	if ((::solver[currSolver].s_flags & f) == f) {
		solverFlags &= ~f;
		return true;
	}

	return false;
}

bool
LinSol::SetNumThreads(unsigned nt)
{
	if (GetSolverFlags(currSolver) & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) {
		if (nt == 0) {
			solverFlags &= ~LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT;

		} else {
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT;
		}

		nThreads = nt;

		return true;
	}

	return false;
}

unsigned
LinSol::GetNumThreads(void) const
{
	return nThreads;
}

integer
LinSol::iGetWorkSpaceSize(void) const
{
	return iWorkSpaceSize;
}

const doublereal&
LinSol::dGetPivotFactor(void) const
{
	return dPivotFactor;
}

const doublereal&
LinSol::dGetDropTolerance(void) const
{
	return dDropTolerance;
}

bool
LinSol::SetWorkSpaceSize(integer i)
{
	switch (currSolver) {
	case LinSol::Y12_SOLVER:
		iWorkSpaceSize = i;
		break;

	default:
		return false;
	}

	return true;
}

bool
LinSol::SetPivotFactor(const doublereal& d)
{
	if (::solver[currSolver].s_pivot_factor == -1.) {
		return false;
	}

	dPivotFactor = d;

	return true;
}

bool
LinSol::SetDropTolerance(const doublereal& d)
{
	if (::solver[currSolver].s_drop_tolerance == -1.) {
		return false;
	}

	dDropTolerance = d;

	return true;
}

unsigned
LinSol::GetBlockSize(void) const
{
	return blockSize;
}

bool
LinSol::SetBlockSize(unsigned bs)
{
	switch (currSolver) {
	case LinSol::UMFPACK_SOLVER:
		blockSize = bs;
		break;

	default:
		return false;
	}

	return true;
}

SolutionManager::ScaleWhen
LinSol::GetScaleWhen(void) const
{
	return ms;
}

bool
LinSol::SetScaleWhen(SolutionManager::ScaleWhen ms)
{
	switch (currSolver) {
	case LinSol::NAIVE_SOLVER:
		this->ms = ms;
		break;

	default:
		return false;
	}

	return true;
}

LinSol::Scale
LinSol::GetScale(void) const
{
	return scale;
}

bool
LinSol::SetScale(Scale s)
{
	switch (currSolver) {
	case LinSol::UMFPACK_SOLVER:
		this->scale = s;
		break;

	default:
		return false;
	}

	return true;
}

SolutionManager *const
LinSol::GetSolutionManager(integer iNLD, integer iLWS) const
{
	SolutionManager *pCurrSM = NULL;
	const unsigned type = (solverFlags & LinSol::SOLVER_FLAGS_TYPE_MASK);
	const unsigned perm = (solverFlags & LinSol::SOLVER_FLAGS_PERM_MASK);
	const bool mt = (solverFlags & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT);

	/* silence warning */
	if (mt) {
		NO_OP;
	}

	ASSERT((::solver[currSolver].s_flags & solverFlags) == solverFlags);

	if (iLWS == 0) {
		iLWS = iWorkSpaceSize;
	}

   	switch (currSolver) {
     	case LinSol::Y12_SOLVER: 
#ifdef USE_Y12
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, iLWS, dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, iLWS, dPivotFactor));
			break;
		}

		default:
      			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				Y12SparseSolutionManager,
				Y12SparseSolutionManager(iNLD, iLWS,
					dPivotFactor));
			break;
		}
      		break;
#else /* !USE_Y12 */
      		silent_cerr("Configure with --with-y12 "
			"to enable Y12 solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_Y12 */

     	case LinSol::SUPERLU_SOLVER: 
#ifdef USE_SUPERLU
#ifdef USE_SUPERLU_MT
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef ParSuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(nThreads, iNLD, dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef ParSuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(nThreads, iNLD, dPivotFactor));
			break;
		}

		default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				ParSuperLUSparseSolutionManager,
				ParSuperLUSparseSolutionManager(nThreads, iNLD,
					dPivotFactor));
			break;
		}
		break;
		
		if (nThreads == 1) {
			silent_cerr("warning, using multithread SuperLU with only one thread; "
				<< std::endl);
		}
#else /* !USE_SUPERLU_MT */
		if (nThreads > 1 && mt) {
			silent_cerr("multithread SuperLU solver support not compiled; "
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor));
			break;
		}

		default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				SuperLUSparseSolutionManager,
				SuperLUSparseSolutionManager(iNLD,
					dPivotFactor));
			break;
		}
#endif /* !USE_SUPERLU_MT */
		break;
#else /* !USE_SUPERLU */
      		silent_cerr("Configure with --with-superlu "
			"to enable superlu solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_SUPERLU */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			MeschachSparseSolutionManager,
			MeschachSparseSolutionManager(iNLD, iLWS,
				dPivotFactor));
		break;
#else /* !USE_MESCHACH */
		silent_cerr("Configure with --with-meschach "
			"to enable Meschach solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_MESCHACH */

	case LinSol::LAPACK_SOLVER:
#ifdef USE_LAPACK
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			LapackSolutionManager,
			LapackSolutionManager(iNLD, dPivotFactor));
		break;
#else /* !USE_LAPACK */
		silent_cerr("Configure with --with-lapack "
			"to enable Lapack dense solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_LAPACK */

     	case LinSol::TAUCS_SOLVER: 
#ifdef USE_TAUCS
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef TaucsSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM, CCSM(iNLD));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef TaucsSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM, CCSM(iNLD));
			break;
		}

		default:
      			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				TaucsSparseSolutionManager,
				TaucsSparseSolutionManager(iNLD));
			break;
		}
      		break;
#else /* !USE_TAUCS */
		silent_cerr("Configure with --with-taucs "
			"to enable Taucs sparse solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_TAUCS */

 	case LinSol::HARWELL_SOLVER:
#ifdef USE_HARWELL
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			HarwellSparseSolutionManager,
			HarwellSparseSolutionManager(iNLD, iLWS,
				dPivotFactor));
      		break;
#else /* !USE_HARWELL */
      		silent_cerr("Configure with --with-harwell "
			"to enable Harwell solver" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_HARWELL */

	case LinSol::UMFPACK_SOLVER: 
#ifdef USE_UMFPACK
		{
		UmfpackSolver::Scale s = UmfpackSolver::SCALE_UNDEF;
		switch (scale) {
		case SCALE_NONE:
			s = UmfpackSolver::SCALE_NONE;
			break;

		case SCALE_MAX:
			s = UmfpackSolver::SCALE_MAX;
			break;

		case SCALE_SUM:
			s = UmfpackSolver::SCALE_SUM;
			break;

		default:
			// don't do anything
			break;
		}

		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor, dDropTolerance, blockSize, s));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor, dDropTolerance, blockSize, s));
			break;
		}

		default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				UmfpackSparseSolutionManager,
				UmfpackSparseSolutionManager(iNLD,
					dPivotFactor, dDropTolerance, blockSize, s));
			break;
		}
		} break;
#else /* !USE_UMFPACK */
      		silent_cerr("Configure with --with-umfpack "
			"to enable Umfpack solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_UMFPACK */

	case LinSol::KLU_SOLVER:
#ifdef USE_KLU
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef KLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef KLUSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor));
			break;
		}

		default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				KLUSparseSolutionManager,
				KLUSparseSolutionManager(iNLD, dPivotFactor));
			break;
		}
      		break;
#else /* !USE_KLU */
      		silent_cerr("Configure with --with-klu "
			"to enable KLU solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_KLU */

	case LinSol::NAIVE_SOLVER:
		if (perm == LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<Colamd_ordering>,
					NaiveSparsePermSolutionManager<Colamd_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#ifdef USE_BOOST
#ifdef HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP
		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<rcmk_ordering>,
					NaiveSparsePermSolutionManager<rcmk_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
				silent_cerr("multithread naive solver with"
					"reverse Cuthill-McKee permutation not"
					"available yet. Patches welcome"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#endif /* HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP */

#if 0 /* ?!? */
 		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_MMDATA) {
 			if (nThreads == 1) {
 				SAFENEWWITHCONSTRUCTOR(pCurrSM,
 					NaiveSparsePermSolutionManager<amd_ordering>,
 					NaiveSparsePermSolutionManager<amd_ordering>(iNLD, dPivotFactor, ms));
 			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
 				silent_cerr("multithread naive solver with"
 					"approximate minimum degree permutation not"
 					"available yet. Patches welcome"
 					<< std::endl);
 				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else
 				silent_cerr("multithread naive solver support not compiled; "
 					"you can configure --enable-multithread-naive "
 					"on a linux ix86 to get it"
 					<< std::endl);
 				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* USE_NAIVE_MULTITHREAD */
 			}
#endif

#ifdef HAVE_BOOST_GRAPH_KING_ORDERING_HPP
		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_KING) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<king_ordering>,
					NaiveSparsePermSolutionManager<king_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
				silent_cerr("multithread naive solver with"
					"king permutation not"
					"available yet. Patches welcome"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#endif /* HAVE_BOOST_GRAPH_KING_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP
		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_SLOAN) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<sloan_ordering>,
					NaiveSparsePermSolutionManager<sloan_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
				silent_cerr("multithread naive solver with"
					"sloan permutation not"
					"available yet. Patches welcome"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#endif /* HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP
		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_MDAPLUSAT) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<md_ordering>,
					NaiveSparsePermSolutionManager<md_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
				silent_cerr("multithread naive solver with"
					"md permutation not"
					"available yet. Patches welcome"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#endif /* HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP */
#endif /* USE_BOOST */

		} else if (perm == LinSol::SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION) {
#ifdef USE_METIS
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<metis_ordering>,
					NaiveSparsePermSolutionManager<metis_ordering>(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
#if 0
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nThreads, iNLD, dPivotFactor));
#endif
				silent_cerr("multithread naive solver with"
					"nested dissection permutation not"
					"available yet. Patches welcome"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* ! USE_NAIVE_MULTITHREAD */
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_NAIVE_MULTITHREAD */
			}
#else /* ! USE_METIS */
			silent_cerr("you should not get here("<< __FILE__ << ":" <<
				__LINE__ << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_METIS */

		} else {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparseSolutionManager,
					NaiveSparseSolutionManager(iNLD, dPivotFactor, ms));
			} else {
#ifdef USE_NAIVE_MULTITHREAD
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					ParNaiveSparseSolutionManager,
					ParNaiveSparseSolutionManager(nThreads, iNLD, dPivotFactor));
#else
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* USE_NAIVE_MULTITHREAD */
			}
		}
		break;

	case LinSol::EMPTY_SOLVER:
		break;
		
   	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	}

	return pCurrSM;
}

