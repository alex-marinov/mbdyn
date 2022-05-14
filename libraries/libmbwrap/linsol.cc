/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "myassert.h"
#include "ac/sys_sysinfo.h"

#include "parser.h"

#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "harwrap.h"
#include "y12wrap.h"
#include "umfpackwrap.h"
#include "kluwrap.h"
#include "parsuperluwrap.h"
#include "superluwrap.h"
#include "lapackwrap.h"
#include "taucswrap.h"
#include "naivewrap.h"
#include "parnaivewrap.h"
#include "pardisowrap.h"
#include "pastixwrap.h"
#include "qrwrap.h"
#include "strumpackwrap.h"
#include "wsmpwrap.h"
#ifdef USE_TRILINOS
#include "aztecoowrap.h"
#endif
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
		1.e-5, -1. },
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
	        LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS|LinSol::SOLVER_FLAGS_ALLOWS_GRAD,
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
	        LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_GRAD|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		.1 },
	{ "Y12", NULL,
		LinSol::Y12_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_CC|LinSol::SOLVER_FLAGS_ALLOWS_DIR|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP|LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
		-1., -1. },
        { "Pardiso", NULL,
		LinSol::PARDISO_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	        LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP,
		1e-13, -1. },
        { "Pardiso_64", NULL,
		LinSol::PARDISO_64_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	        LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP,
		1e-13, -1. },
        { "Pastix", NULL,
		LinSol::PASTIX_SOLVER,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	        LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT |
	        LinSol::SOLVER_FLAGS_ALLOWS_SCOTCH |
	        LinSol::SOLVER_FLAGS_ALLOWS_METIS |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_SVD |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_PQRCP | 
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRCP |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_TQRCP |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRRT,
		LinSol::SOLVER_FLAGS_ALLOWS_MAP,
		-1., -1. },
        { "QR", NULL,
                LinSol::QR_SOLVER,
                LinSol::SOLVER_FLAGS_NONE,
                LinSol::SOLVER_FLAGS_NONE,
                -1., -1. },
        { "SPQR", NULL,
                LinSol::SPQR_SOLVER,
                LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	        LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
                LinSol::SOLVER_FLAGS_ALLOWS_AMD |
                LinSol::SOLVER_FLAGS_ALLOWS_METIS |
                LinSol::SOLVER_FLAGS_ALLOWS_GIVEN,
                LinSol::SOLVER_FLAGS_ALLOWS_MAP,
                -1., -1. },
        { "STRUMPACK", NULL,
                LinSol::STRUMPACK_SOLVER,
                LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	        LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS |
	        LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT |
	        LinSol::SOLVER_FLAGS_ALLOWS_METIS |
	        LinSol::SOLVER_FLAGS_ALLOWS_SCOTCH |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_HSS |
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_BLR | 
	        LinSol::SOLVER_FLAGS_ALLOWS_COMPRESSION_HODLR,
                LinSol::SOLVER_FLAGS_ALLOWS_MAP,
                -1., -1. },	
	{ "Watson", NULL,
	  LinSol::WATSON_SOLVER,
	  LinSol::SOLVER_FLAGS_ALLOWS_MAP |
	  LinSol::SOLVER_FLAGS_ALLOWS_CC |
	  LinSol::SOLVER_FLAGS_ALLOWS_DIR |
	  LinSol::SOLVER_FLAGS_ALLOWS_GRAD |
#ifdef WSMP_SOLVER_ENABLE_MULTITHREAD_FACTOR // Disable by default because valgrind reports a race condition in WSMP Version 20.05.20
	  LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT | 
#endif
	  LinSol::SOLVER_FLAGS_ALLOWS_MT_ASS,
	  LinSol::SOLVER_FLAGS_ALLOWS_MAP,
	  -1., -1. },
        {"AztecOO", NULL,
         LinSol::AZTECOO_SOLVER,
         LinSol::SOLVER_FLAGS_PRECOND_MASK,
         LinSol::SOLVER_FLAGS_ALLOWS_PRECOND_UMFPACK,
         -1, -1},
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
dLowRankCompressTol(0.01),
dLowRankCompressMinRatio(1.),
iMaxIter(0), // Restore the original behavior by default
dTolRes(1e-10),
iVerbose(0)
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

        case LinSol::PARDISO_SOLVER:
        case LinSol::PARDISO_64_SOLVER:
#ifdef USE_PARDISO
             currSolver = t;
             return true;
#endif
        case LinSol::PASTIX_SOLVER:
#ifdef USE_PASTIX
                currSolver = t;
                return true;
#endif
        case LinSol::QR_SOLVER:
                currSolver = t;
                return true;
#ifdef USE_SUITESPARSE_QR
        case LinSol::SPQR_SOLVER:
                currSolver = t;
                return true;
#endif
#ifdef USE_STRUMPACK
	case LinSol::STRUMPACK_SOLVER:
	     currSolver = t;
	     return true;
#endif
#ifdef USE_WSMP
	case LinSol::WATSON_SOLVER:
	     currSolver = t;
	     return true;
#endif
#ifdef USE_TRILINOS
        case LinSol::AZTECOO_SOLVER:
             currSolver = t;
             return true;
#endif
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
LinSol::AddSolverFlags(unsigned mask, unsigned flag)
{
        ASSERT((mask & flag) == flag);
	
        if ((::solver[currSolver].s_flags & flag) == flag) {
	        solverFlags &= ~mask;
		solverFlags |= flag;
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

bool
LinSol::SetLowRankCompressTol(const doublereal& d)
{
     if (0 == (::solver[currSolver].s_flags & LinSol::SOLVER_FLAGS_COMPRESSION_MASK)) {
	  return false;
     }

     dLowRankCompressTol = d;

     return true;
}

bool LinSol::SetLowRankCompressMinRatio(const doublereal& d)
{
     if (0 == (::solver[currSolver].s_flags & LinSol::SOLVER_FLAGS_COMPRESSION_MASK)) {
	  return false;
     }

     dLowRankCompressMinRatio = d;

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

bool
LinSol::SetScale(const SolutionManager::ScaleOpt& s)
{
	switch (currSolver) {
	case LinSol::NAIVE_SOLVER:
	case LinSol::UMFPACK_SOLVER:
	case LinSol::KLU_SOLVER:
        case LinSol::PASTIX_SOLVER:
	case LinSol::WATSON_SOLVER:
		scale = s;
		break;

	default:
		return false;
	}

	return true;
}

integer
LinSol::GetMaxIterations(void) const
{
	return iMaxIter;
}

bool
LinSol::SetMaxIterations(integer iMaxIterations)
{
	switch (currSolver) {
	case LinSol::UMFPACK_SOLVER:
        case LinSol::PARDISO_SOLVER:
        case LinSol::PARDISO_64_SOLVER:
        case LinSol::PASTIX_SOLVER:
	case LinSol::STRUMPACK_SOLVER:
	case LinSol::WATSON_SOLVER:
        case LinSol::AZTECOO_SOLVER:
		iMaxIter = iMaxIterations;
		break;

	default:
		return false;
	}
        
	return true;
}

bool LinSol::SetTolerance(doublereal dToleranceRes)
{
        switch (currSolver) {
        case LinSol::AZTECOO_SOLVER:
                dTolRes = dToleranceRes;
                break;
        default:
                return false;
        }

        return true;               
}

bool LinSol::SetVerbose(integer iVerb)
{
     switch (currSolver) {
     case LinSol::UMFPACK_SOLVER:
     case LinSol::PARDISO_SOLVER:
     case LinSol::PARDISO_64_SOLVER:
     case LinSol::PASTIX_SOLVER:
     case LinSol::STRUMPACK_SOLVER:
     case LinSol::AZTECOO_SOLVER:
	  iVerbose = iVerb;
	  break;
	  
     default:
	  return false;
     }
     
     return true;
}

SolutionManager *const
LinSol::GetSolutionManager(integer iNLD,
#ifdef USE_MPI
                           MPI::Intracomm& oComm,
#endif
                           integer iLWS) const
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
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					       CCSM(iNLD, dPivotFactor, dDropTolerance, blockSize, scale, iMaxIter, iVerbose));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					       CCSM(iNLD, dPivotFactor, dDropTolerance, blockSize, scale, iMaxIter, iVerbose));
			break;
		}
		case LinSol::SOLVER_FLAGS_ALLOWS_GRAD: {
		        typedef UmfpackSparseSolutionManager<SpGradientSparseMatrixHandler> UMFSM;
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
					       UMFSM,
					       UMFSM(iNLD, dPivotFactor, dDropTolerance, blockSize, scale, iMaxIter, iVerbose));
			break;		     
		}
		default:
		        typedef UmfpackSparseSolutionManager<SpMapMatrixHandler> UMFSM;
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
					       UMFSM,
					       UMFSM(iNLD, dPivotFactor, dDropTolerance, blockSize, scale, iMaxIter, iVerbose));
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
		{
		switch (type) {
		case LinSol::SOLVER_FLAGS_ALLOWS_DIR: {
			typedef KLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					       CCSM(iNLD, dPivotFactor, scale));
			break;
		}

		case LinSol::SOLVER_FLAGS_ALLOWS_CC: {
			typedef KLUSparseCCSolutionManager<CColMatrixHandler<0> > CCSM;
	      		SAFENEWWITHCONSTRUCTOR(pCurrSM, CCSM,
					       CCSM(iNLD, dPivotFactor, scale));
			break;
		}
		     
		case LinSol::SOLVER_FLAGS_ALLOWS_GRAD: {
		     typedef KLUSparseSolutionManager<SpGradientSparseMatrixHandler> KLUSM;
		     SAFENEWWITHCONSTRUCTOR(pCurrSM,
					    KLUSM,
					    KLUSM(iNLD, dPivotFactor, scale));
		     break;
		}
                     
		default:
		     typedef KLUSparseSolutionManager<SpMapMatrixHandler> KLUSM;
			  SAFENEWWITHCONSTRUCTOR(pCurrSM,
						 KLUSM,
						 KLUSM(iNLD, dPivotFactor, scale));
			break;
		}

		} break;
#else /* !USE_KLU */
      		silent_cerr("Configure with --with-klu "
			"to enable KLU solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_KLU */
        case LinSol::PARDISO_SOLVER:
        case LinSol::PARDISO_64_SOLVER:
#ifdef USE_PARDISO
                {
                        typedef PardisoSolutionManager<SpGradientSparseMatrixHandler, MKL_INT> PARDISO_SM_GRAD;
                        typedef PardisoSolutionManager<SpGradientSparseMatrixHandler, long long> PARDISO_64_SM_GRAD;
                        typedef PardisoSolutionManager<SpMapMatrixHandler, MKL_INT> PARDISO_SM_MAP;
                        typedef PardisoSolutionManager<SpMapMatrixHandler, long long> PARDISO_64_SM_MAP;
                        
                        switch (type) {
                        case LinSol::SOLVER_FLAGS_ALLOWS_GRAD: {
                                if (currSolver == LinSol::PARDISO_64_SOLVER) {
                                        SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                                               PARDISO_64_SM_GRAD,
                                                               PARDISO_64_SM_GRAD(iNLD, dPivotFactor, nThreads, iMaxIter, iVerbose));
                                } else {
                                        SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                                               PARDISO_SM_GRAD,
                                                               PARDISO_SM_GRAD(iNLD, dPivotFactor, nThreads, iMaxIter, iVerbose));
                                }
                        } break;
                        default:
                                if (currSolver == LinSol::PARDISO_64_SOLVER) {
                                        SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                                               PARDISO_64_SM_MAP,
                                                               PARDISO_64_SM_MAP(iNLD, dPivotFactor, nThreads, iMaxIter, iVerbose));
                                } else {
                                        SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                                               PARDISO_SM_MAP,
                                                               PARDISO_SM_MAP(iNLD, dPivotFactor, nThreads, iMaxIter, iVerbose));
                                }
                        }
		} break;
#else /* !USE_PARDISO */
      		silent_cerr("Configure with --with-pardiso "
                            "to enable Pardiso solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_PARDISO */
	case LinSol::PASTIX_SOLVER: 
#ifdef USE_PASTIX
		{
                    switch (type) {
		    case LinSol::SOLVER_FLAGS_ALLOWS_GRAD: {
                        SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                               PastixSolutionManager<SpGradientSparseMatrixHandler>,
                                               PastixSolutionManager<SpGradientSparseMatrixHandler>(iNLD, nThreads, iMaxIter, scale, solverFlags, dLowRankCompressTol, dLowRankCompressMinRatio, iVerbose));
                    } break;
                    default:
			SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                               PastixSolutionManager<SpMapMatrixHandler>,
                                               PastixSolutionManager<SpMapMatrixHandler>(iNLD, nThreads, iMaxIter, scale, solverFlags, dLowRankCompressTol, dLowRankCompressMinRatio, iVerbose));
                    }
		} break;
#else /* !USE_PASTIX */
      		silent_cerr("Configure with --with-pastix "
			"to enable Pastix solver" << std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* !USE_PASTIX */
        case LinSol::QR_SOLVER:
                SAFENEWWITHCONSTRUCTOR(pCurrSM,
                                       QrDenseSolutionManager,
                                       QrDenseSolutionManager(iNLD));
                break;
#ifdef USE_SUITESPARSE_QR
        case LinSol::SPQR_SOLVER:
	     switch (type) {
	     case LinSol::SOLVER_FLAGS_ALLOWS_GRAD:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 QrSparseSolutionManager<SpGradientSparseMatrixHandler>,
					 QrSparseSolutionManager<SpGradientSparseMatrixHandler>(iNLD, solverFlags));
		  break;
	     default:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 QrSparseSolutionManager<SpMapMatrixHandler>,
					 QrSparseSolutionManager<SpMapMatrixHandler>(iNLD, solverFlags));
		  break;
	     }
	     break;
#endif
#ifdef USE_STRUMPACK
	case LinSol::STRUMPACK_SOLVER:
	     switch (type) {
	     case LinSol::SOLVER_FLAGS_ALLOWS_GRAD:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 StrumpackSolutionManager<SpGradientSparseMatrixHandler>,
					 StrumpackSolutionManager<SpGradientSparseMatrixHandler>(iNLD, nThreads, iMaxIter, solverFlags, iVerbose));
		  break;
	     default:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 StrumpackSolutionManager<SpMapMatrixHandler>,
					 StrumpackSolutionManager<SpMapMatrixHandler>(iNLD, nThreads, iMaxIter, solverFlags, iVerbose));
		  break;
	     }
	     break;
#endif
#ifdef USE_WSMP
	case LinSol::WATSON_SOLVER:
	     switch (type) {
	     case LinSol::SOLVER_FLAGS_ALLOWS_CC:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 WsmpSparseCCSolutionManager<CColMatrixHandler<0> >,
					 WsmpSparseCCSolutionManager<CColMatrixHandler<0> >(iNLD, dPivotFactor, blockSize, nThreads, scale, iMaxIter));
		  break;
	     case LinSol::SOLVER_FLAGS_ALLOWS_DIR:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 WsmpSparseCCSolutionManager<DirCColMatrixHandler<0> >,
					 WsmpSparseCCSolutionManager<DirCColMatrixHandler<0> >(iNLD, dPivotFactor, blockSize, nThreads, scale, iMaxIter));
		  break;
	     case LinSol::SOLVER_FLAGS_ALLOWS_GRAD:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 WsmpSparseSolutionManager<SpGradientSparseMatrixHandler>,
					 WsmpSparseSolutionManager<SpGradientSparseMatrixHandler>(iNLD, dPivotFactor, blockSize, nThreads, scale, iMaxIter));
		  break;
	     default:
		  SAFENEWWITHCONSTRUCTOR(pCurrSM,
					 WsmpSparseSolutionManager<SpMapMatrixHandler>,
					 WsmpSparseSolutionManager<SpMapMatrixHandler>(iNLD, dPivotFactor, blockSize, nThreads, scale, iMaxIter));
	     }
	     break;
#endif
#ifdef USE_TRILINOS
        case LinSol::AZTECOO_SOLVER:
             pCurrSM = pAllocateAztecOOSolutionManager(
#ifdef USE_MPI
                  oComm,
#endif
                  iNLD,
                  iMaxIter,
                  dTolRes,
                  iVerbose,
                  solverFlags);
             break;
#endif
	case LinSol::NAIVE_SOLVER:
		if (perm == LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) {
			if (nThreads == 1) {
				SAFENEWWITHCONSTRUCTOR(pCurrSM,
					NaiveSparsePermSolutionManager<Colamd_ordering>,
					NaiveSparsePermSolutionManager<Colamd_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparsePermSolutionManager<rcmk_ordering>(iNLD, dPivotFactor, scale));
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
 					NaiveSparsePermSolutionManager<amd_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparsePermSolutionManager<king_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparsePermSolutionManager<sloan_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparsePermSolutionManager<md_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparsePermSolutionManager<metis_ordering>(iNLD, dPivotFactor, scale));
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
					NaiveSparseSolutionManager(iNLD, dPivotFactor, scale));
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

