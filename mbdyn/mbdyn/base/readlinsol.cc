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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "ac/sys_sysinfo.h"
#include "dataman.h"
#include "readlinsol.h"

void
ReadLinSol(LinSol& cs, HighParser &HP, bool bAllowEmpty)
{
   	/* parole chiave */
   	const char* sKeyWords[] = { 
		::solver[LinSol::EMPTY_SOLVER].s_name,
		::solver[LinSol::HARWELL_SOLVER].s_name,
		::solver[LinSol::LAPACK_SOLVER].s_name,
		::solver[LinSol::MESCHACH_SOLVER].s_name,
		::solver[LinSol::NAIVE_SOLVER].s_name,
		::solver[LinSol::SUPERLU_SOLVER].s_name,
		::solver[LinSol::TAUCS_SOLVER].s_name,
		::solver[LinSol::UMFPACK_SOLVER].s_name,
		::solver[LinSol::UMFPACK_SOLVER].s_alias,
		::solver[LinSol::Y12_SOLVER].s_name,
		NULL
	};

	enum KeyWords {
		EMPTY,
		HARWELL,
		LAPACK,
		MESCHACH,
		NAIVE,
		SUPERLU,
		TAUCS,
		UMFPACK,
		UMFPACK3,
		Y12,

		LASTKEYWORD
	};

   	/* tabella delle parole chiave */
   	KeyTable K(HP, sKeyWords);

	bool bGotIt = false;	
	switch (KeyWords(HP.GetWord())) {
	case EMPTY:
		if (!bAllowEmpty) {
			silent_cerr("empty solver not allowed at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		cs.currSolver = LinSol::EMPTY_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"No LU solver" << std::endl);
		bGotIt = true;
		break;

	case HARWELL:
		cs.currSolver = LinSol::HARWELL_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using harwell sparse LU solver" << std::endl);
#ifdef USE_HARWELL
		bGotIt = true;
#endif /* USE_HARWELL */
		break;

	case LAPACK:
		cs.currSolver = LinSol::LAPACK_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using Lapack dense LU solver" << std::endl);
#ifdef USE_LAPACK
		bGotIt = true;
#endif /* USE_LAPACK */
		break;

	case MESCHACH:
		cs.currSolver = LinSol::MESCHACH_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using meschach sparse LU solver"
				<< std::endl);
#ifdef USE_MESCHACH
		bGotIt = true;
#endif /* USE_MESCHACH */
		break;

	case NAIVE:
		cs.currSolver = LinSol::NAIVE_SOLVER;
		bGotIt = true;
		break;

	case SUPERLU:
		/*
		 * FIXME: use CC as default???
		 */
		cs.currSolver = LinSol::SUPERLU_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using SuperLU sparse LU solver" << std::endl);
#ifdef USE_SUPERLU
		bGotIt = true;
#endif /* USE_SUPERLU */
		break;

	case TAUCS:
		/*
		 * FIXME: use CC as default???
		 */
		cs.currSolver = LinSol::TAUCS_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using Taucs sparse solver" << std::endl);
#ifdef USE_TAUCS
		bGotIt = true;
#endif /* USE_TAUCS */
		break;

	case UMFPACK3:
		pedantic_cerr("\"umfpack3\" is deprecated; "
				"use \"umfpack\" instead" << std::endl);
	case UMFPACK:
		/*
		 * FIXME: use CC as default???
		 */
		cs.currSolver = LinSol::UMFPACK_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using umfpack sparse LU solver" << std::endl);
#ifdef USE_UMFPACK
		bGotIt = true;
#endif /* USE_UMFPACK */
		break;

	case Y12:
		/*
		 * FIXME: use CC as default???
		 */
		cs.currSolver = LinSol::Y12_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using y12 sparse LU solver" << std::endl);
#ifdef USE_Y12
		bGotIt = true;
#endif /* USE_Y12 */
		break;

	default:
		silent_cerr("unknown solver" << std::endl);
		throw ErrGeneric();
	}

	if (!bGotIt) {
		silent_cerr(::solver[cs.currSolver].s_name << " solver "
			"not available at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	cs.solverFlags = ::solver[cs.currSolver].s_default_flags;

	/* map? */
	if (HP.IsKeyWord("map")) {
		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MAP) {
			cs.solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			cs.solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MAP;
			pedantic_cout("using map matrix handling for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("map is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
		}

	/* CC? */
	} else if (HP.IsKeyWord("column" "compressed") || HP.IsKeyWord("cc")) {
		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_CC) {
			cs.solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			cs.solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_CC;
			pedantic_cout("using column compressed matrix handling for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("column compressed is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
		}

	/* direct? */
	} else if (HP.IsKeyWord("direct" "access") || HP.IsKeyWord("dir")) {
		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_DIR) {
			cs.solverFlags &= ~LinSol::SOLVER_FLAGS_TYPE_MASK;
			cs.solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_DIR;
			pedantic_cout("using direct access matrix handling for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("direct is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
		}
	}

	/* colamd? */
	if (HP.IsKeyWord("colamd")) {
		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) {
			cs.solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_COLAMD;
			pedantic_cout("using colamd preordering for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("colamd preordering is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
		}
	}

	/* multithread? */
	if (HP.IsKeyWord("multi" "thread") || HP.IsKeyWord("mt")) {
		cs.nThreads = HP.GetInt();

		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) {
			cs.solverFlags |= LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT;
			if (cs.nThreads < 1) {
				silent_cerr("illegal thread number, using 1" << std::endl);
				cs.nThreads = 1;
			}

		} else {
			pedantic_cerr("multithread is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
			cs.nThreads = 1;
		}
	} else {
		if (::solver[cs.currSolver].s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) {
			int n = get_nprocs();

			if (n > 1) {
				silent_cout("no multithread requested "
						"with a potential of " << n
						<< " CPUs" << std::endl);
				cs.nThreads = n;

			} else {
				cs.nThreads = 1;
			}
		}
	}

	if (HP.IsKeyWord("workspace" "size")) {
		cs.iWorkSpaceSize = HP.GetInt();
		if (cs.iWorkSpaceSize < 0) {
			cs.iWorkSpaceSize = 0;
		}

		switch (cs.currSolver) {
		case LinSol::Y12_SOLVER:
			break;

		default:
			pedantic_cerr("workspace size is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
			break;
		}
	}

	if (HP.IsKeyWord("pivot" "factor")) {
		cs.dPivotFactor = HP.GetReal();

		if (::solver[cs.currSolver].s_pivot_factor == -1.) {
			pedantic_cerr("pivot factor is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
			cs.dPivotFactor = -1.;

		} else if (cs.dPivotFactor <= 0. || cs.dPivotFactor > 1.) {
			cs.dPivotFactor = ::solver[cs.currSolver].s_pivot_factor;
		}

	} else {
		if (::solver[cs.currSolver].s_pivot_factor != -1.) {
			cs.dPivotFactor = ::solver[cs.currSolver].s_pivot_factor;
		}
	}

	if (HP.IsKeyWord("block" "size")) {
		integer i = HP.GetInt();
		if (i < 1) {
			silent_cerr("illegal negative block size; "
					"using default" << std::endl);
			cs.blockSize = 0;

		} else {
			cs.blockSize = (unsigned)i;
		}

		switch (cs.currSolver) {
		case LinSol::UMFPACK_SOLVER:
			break;

		default:
			pedantic_cerr("block size is meaningless for "
					<< ::solver[cs.currSolver].s_name
					<< " solver" << std::endl);
			break;
		}
	}
}

