/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "myassert.h"
#include "ac/sys_sysinfo.h"

#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "harwrap.h"
#include "mschwrap.h"
#include "y12wrap.h"
#include "umfpackwrap.h"
#include "superluwrap.h"

#include "mbpar.h"

#include "linearsolver.h"

/* private data */
static struct solver_t {
	const char *const	s_name;
	const char *const	s_alias;
	enum LinSol::SolverType	s_type;
	unsigned		s_flags;
	unsigned		s_default_flags;
} solver[] = {
	{ "Harwell", NULL,
		LinSol::HARWELL_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE },
	{ "Meschach", NULL,
		LinSol::MESCHACH_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE },
	{ "Y12", NULL,
		LinSol::Y12_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP },
	{ "Umfpack", "umfpack3", 
		LinSol::UMFPACK_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP },
	{ "SuperLU", NULL, 
		LinSol::SUPERLU_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT },
	{ "Empty", NULL,
		LinSol::EMPTY_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE },
	{ NULL, NULL, 
		LinSol::EMPTY_SOLVER,
		LinSol::SOLVER_FLAGS_NONE,
		LinSol::SOLVER_FLAGS_NONE },
};

/*
 * Default solver
 */
LinSol::SolverType LinSol::defaultSolver = 
#if defined(USE_UMFPACK)
	LinSol::UMFPACK_SOLVER
#elif /* !USE_UMFPACK */ defined(USE_Y12)
	LinSol::Y12_SOLVER
#elif /* !USE_Y12 */ defined(USE_SUPERLU)
	LinSol::SUPERLU_SOLVER
#elif /* !USE_SUPERLU */ defined(USE_HARWELL)
	LinSol::HARWELL_SOLVER
#elif /* !USE_HARWELL */ defined(USE_MESCHACH)
	LinSol::MESCHACH_SOLVER
#else /* !USE_MESCHACH */
	LinSol::EMPTY_SOLVER
/* FIXME: remove this error if no solver becomes acceptable :) */
#error "need a solver!"
#endif /* !USE_MESCHACH */
	;

LinSol::LinSol(void)
: CurrSolver(LinSol::defaultSolver),
nThreads(1),
iWorkSpaceSize(0),
dPivotFactor(1.)
{
	NO_OP;
}

LinSol::~LinSol(void)
{
	NO_OP;
}

void
LinSol::Read(HighParser &HP, bool bAllowEmpty)
{
   	/* parole chiave */
   	const char* sKeyWords[] = { 
		::solver[LinSol::HARWELL_SOLVER].s_name,
		::solver[LinSol::MESCHACH_SOLVER].s_name,
		::solver[LinSol::Y12_SOLVER].s_name,
		::solver[LinSol::UMFPACK_SOLVER].s_name,
		::solver[LinSol::UMFPACK_SOLVER].s_alias,
		::solver[LinSol::SUPERLU_SOLVER].s_name,
		::solver[LinSol::EMPTY_SOLVER].s_name,
		NULL
	};

	enum KeyWords {
		HARWELL,
		MESCHACH,
		Y12,
		UMFPACK,
		UMFPACK3,
		SUPERLU,
		EMPTY,

		LASTKEYWORD
	};

   	/* tabella delle parole chiave */
   	KeyTable K(HP, sKeyWords);
   
	switch(KeyWords(HP.GetWord())) {
	case MESCHACH:
#ifdef USE_MESCHACH
		CurrSolver = LinSol::MESCHACH_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using meschach sparse LU solver"
				<< std::endl);
#endif /* USE_MESCHACH */
		break;

	case Y12:
#ifdef USE_Y12
		/*
		 * FIXME: use CC as default???
		 */
		CurrSolver = LinSol::Y12_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using y12 sparse LU solver" << std::endl);
#endif /* USE_Y12 */
		break;

	case SUPERLU:
#ifdef USE_SUPERLU
		/*
		 * FIXME: use CC as default???
		 */
		CurrSolver = LinSol::SUPERLU_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using SuperLU sparse LU solver" << std::endl);
#endif /* USE_SUPERLU */
		break;

	case UMFPACK3:
		pedantic_cerr("\"umfpack3\" is deprecated; "
				"use \"umfpack\" instead" << std::endl);
	case UMFPACK:
#ifdef USE_UMFPACK
		/*
		 * FIXME: use CC as default???
		 */
		CurrSolver = LinSol::UMFPACK_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using umfpack sparse LU solver" << std::endl);
#endif /* USE_UMFPACK */
		break;

	case HARWELL:
#ifdef USE_HARWELL
		CurrSolver = LinSol::HARWELL_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using harwell sparse LU solver" << std::endl);
#endif /* USE_HARWELL */
		break;

	case EMPTY:
		if (!bAllowEmpty) {
			std::cerr << "empty solver not allowed at line " << HP.GetLineData() << std::endl;
			THROW(ErrGeneric());
		}

		CurrSolver = LinSol::EMPTY_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"No LU solver" << std::endl);
		break;

	default:
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Unknown solver; switching to default" << std::endl);
		break;
	}

	solverFlags |= ::solver[CurrSolver].s_default_flags;

	/* map? */
	if (HP.IsKeyWord("map")) {
		if (::solver[CurrSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MAP) {
			solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MAP;
			pedantic_cout("using map matrix handling for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("map is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
		}

	/* CC? */
	} else if (HP.IsKeyWord("column" "compressed") || HP.IsKeyWord("cc")) {
		if (::solver[CurrSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_CC) {
			solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_CC;
			pedantic_cout("using column compressed matrix handling for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("column compressed is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
		}

	/* direct? */
	} else if (HP.IsKeyWord("direct" "access") || HP.IsKeyWord("dir")) {
		if (::solver[CurrSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_DIR) {
			solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_DIR;
			pedantic_cout("using direct access matrix handling for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("direct is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
		}
	}

	/* mutithread? */
	if (HP.IsKeyWord("multi" "thread") || HP.IsKeyWord("mt")) {
		nThreads = HP.GetInt();

		if (::solver[CurrSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT) {
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MT;
			if (nThreads < 1) {
				silent_cerr("illegal thread number, using 1" << std::endl);
				nThreads = 1;
			}

		} else {
			pedantic_cerr("multithread is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
			nThreads = 1;
		}
	} else {
		if (::solver[CurrSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT) {
			int n = get_nprocs();

			if (n > 1) {
				silent_cout("no multithread requested "
						"with a potential of " << n
						<< " CPUs" << std::endl);
				nThreads = n;

			} else {
				nThreads = 1;
			}
		}
	}

	if (HP.IsKeyWord("workspace" "size")) {
		iWorkSpaceSize = HP.GetInt();
		if (iWorkSpaceSize < 0) {
			iWorkSpaceSize = 0;
		}

		switch (CurrSolver) {
		case LinSol::EMPTY_SOLVER:
		case LinSol::MESCHACH_SOLVER:
		case LinSol::UMFPACK_SOLVER:
			pedantic_cerr("workspace size is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
			break;

		default:
			break;
		}
	}

	if (HP.IsKeyWord("pivot" "factor")) {
		dPivotFactor = HP.GetReal();
		if (dPivotFactor <= 0. || dPivotFactor > 1.) {
			dPivotFactor = 1.;
		}

		switch (CurrSolver) {
		case LinSol::EMPTY_SOLVER:
			pedantic_cerr("pivot factor is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
			break;

		default:
			break;
		}
	}

	if (HP.IsKeyWord("block" "size")) {
		integer i = HP.GetInt();
		if (i < 1) {
			silent_cerr("illegal negative block size; "
					"using default" << std::endl);
			blockSize = 0;
		} else {
			blockSize = (unsigned)i;
		}

		switch (CurrSolver) {
		case LinSol::UMFPACK_SOLVER:
			break;

		default:
			pedantic_cerr("workspace size is meaningless for "
					<< ::solver[CurrSolver].s_name
					<< " solver" << std::endl);
			break;
		}
	}
}

LinSol::SolverType
LinSol::GetSolver(void) const
{
	return CurrSolver;
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
		CurrSolver = t;
		return true;
#endif /* USE_UMFPACK */

	case LinSol::SUPERLU_SOLVER:
#ifdef USE_SUPERLU
		CurrSolver = t;
		return true;
#endif /* USE_SUPERLU */

	case LinSol::Y12_SOLVER:
#ifdef USE_Y12
		CurrSolver = t;
		return true;
#endif /* USE_Y12 */

	case LinSol::HARWELL_SOLVER:
#ifdef USE_HARWELL
		CurrSolver = t;
		return true;
#endif /* USE_HARWELL */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		CurrSolver = t;
		return true;
#endif /* USE_MESCHACH */

		/* else */
		silent_cerr(::solver[t].s_name << " unavailable" << std::endl);
		return false;

	case LinSol::EMPTY_SOLVER:
		CurrSolver = t;
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
	return ::solver[CurrSolver].s_name;
}

const char *const
LinSol::GetSolverName(SolverType t) const
{
	return ::solver[t].s_name;
}

bool
LinSol::SetSolverFlags(unsigned f)
{
	if ((::solver[CurrSolver].s_flags & f) == f) {
		solverFlags = f;
		return true;
	}

	return false;
}

bool
LinSol::AddSolverFlags(unsigned f)
{
	if ((::solver[CurrSolver].s_flags & f) == f) {
		solverFlags |= f;
		return true;
	}

	return false;
}

bool
LinSol::MaskSolverFlags(unsigned f)
{
	if ((::solver[CurrSolver].s_flags & f) == f) {
		solverFlags &= ~f;
		return true;
	}

	return false;
}

bool
LinSol::SetNumThreads(unsigned nt)
{
	if (GetSolverFlags(CurrSolver) & LinSol::SOLVER_FLAGS_ALLOWS_MT) {
		if (nt == 0) {
			solverFlags &= ~LinSol::SOLVER_FLAGS_ALLOWS_MT;

		} else {
			solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MT;
		}

		nThreads = nt;

		return true;
	}

	return false;
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

SolutionManager *const
LinSol::GetSolutionManager(integer iNLD, integer iLWS) const
{
	SolutionManager *pCurrSM = NULL;
	unsigned type = (solverFlags & LinSol::SOLVER_FLAGS_TYPE_MASK);
	bool mt = (solverFlags & LinSol::SOLVER_FLAGS_ALLOWS_MT);

	ASSERT((::solver[CurrSolver].s_flags & solverFlags) == solverFlags);

	if (iLWS == 0) {
		iLWS = iWorkSpaceSize;
	}

   	switch (CurrSolver) {
     	case LinSol::Y12_SOLVER: 
#ifdef USE_Y12
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, iLWS,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, iLWS,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}

		default:
      			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				Y12SparseSolutionManager,
				Y12SparseSolutionManager(iNLD, iLWS,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

     	case LinSol::SUPERLU_SOLVER: 
#ifdef USE_SUPERLU
		if (!mt) {
			silent_cerr("warning: SuperLU supperted only in multithread form" << std::endl);
		}

		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(nThreads, iNLD,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(nThreads, iNLD,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}

		default:
      			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				SuperLUSparseSolutionManager,
				SuperLUSparseSolutionManager(nThreads, iNLD,
					dPivotFactor == -1. ? 1. : dPivotFactor));
			break;
		}
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			MeschachSparseSolutionManager,
			MeschachSparseSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
		break;
#else /* !USE_MESCHACH */
		std::cerr << "Configure with --with-meschach "
			"to enable Meschach solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_MESCHACH */

 	case LinSol::HARWELL_SOLVER:
#ifdef USE_HARWELL
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			HarwellSparseSolutionManager,
			HarwellSparseSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_HARWELL */
      		std::cerr << "Configure with --with-harwell "
			"to enable Harwell solver" << std::endl;
		THROW(ErrGeneric());
#endif /* !USE_HARWELL */

	case LinSol::UMFPACK_SOLVER:
#ifdef USE_UMFPACK
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor, blockSize));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					CCSM(iNLD, dPivotFactor, blockSize));
			break;
		}

		default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
				UmfpackSparseSolutionManager,
				UmfpackSparseSolutionManager(iNLD,
					dPivotFactor, blockSize));
			break;
		}
      		break;
#else /* !USE_UMFPACK */
      		std::cerr << "Configure with --with-umfpack "
			"to enable Umfpack solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_UMFPACK */

	case LinSol::EMPTY_SOLVER:
		break;
		
   	default:
		ASSERT(0);
		THROW(ErrGeneric());

	}

	return pCurrSM;
}

