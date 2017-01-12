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
		::solver[LinSol::KLU_SOLVER].s_name,
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
		KLU,
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
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		cs.SetSolver(LinSol::EMPTY_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"No LU solver" << std::endl);
		bGotIt = true;
		break;

	case HARWELL:
		cs.SetSolver(LinSol::HARWELL_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using harwell sparse LU solver" << std::endl);
#ifdef USE_HARWELL
		bGotIt = true;
#endif /* USE_HARWELL */
		break;

	case LAPACK:
		cs.SetSolver(LinSol::LAPACK_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using Lapack dense LU solver" << std::endl);
#ifdef USE_LAPACK
		bGotIt = true;
#endif /* USE_LAPACK */
		break;

	case MESCHACH:
		cs.SetSolver(LinSol::MESCHACH_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using meschach sparse LU solver"
				<< std::endl);
#ifdef USE_MESCHACH
		bGotIt = true;
#endif /* USE_MESCHACH */
		break;

	case NAIVE:
		cs.SetSolver(LinSol::NAIVE_SOLVER);
		bGotIt = true;
		break;

	case SUPERLU:
		/*
		 * FIXME: use CC as default???
		 */
		cs.SetSolver(LinSol::SUPERLU_SOLVER);
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
		cs.SetSolver(LinSol::TAUCS_SOLVER);
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
		cs.SetSolver(LinSol::UMFPACK_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using umfpack sparse LU solver" << std::endl);
#ifdef USE_UMFPACK
		bGotIt = true;
#endif /* USE_UMFPACK */
		break;

	case KLU:
		/*
		 * FIXME: use CC as default???
		 */
		cs.SetSolver(LinSol::KLU_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using KLU sparse LU solver" << std::endl);
#ifdef USE_KLU
		bGotIt = true;
#endif /* USE_KLU */
		break;

	case Y12:
		/*
		 * FIXME: use CC as default???
		 */
		cs.SetSolver(LinSol::Y12_SOLVER);
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using y12 sparse LU solver" << std::endl);
#ifdef USE_Y12
		bGotIt = true;
#endif /* USE_Y12 */
		break;

	default:
		silent_cerr("unknown solver" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const LinSol::solver_t	currSolver = ::solver[cs.GetSolver()];

	if (!bGotIt) {
		silent_cerr(currSolver.s_name << " solver "
			"not available; requested at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	cs.SetSolverFlags(currSolver.s_default_flags);

	/* map? */
	if (HP.IsKeyWord("map")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MAP) {
			cs.MaskSolverFlags(LinSol::SOLVER_FLAGS_TYPE_MASK);
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_MAP);
			pedantic_cout("using map matrix handling for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		} else {
			pedantic_cerr("map is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}

	/* CC? */
	} else if (HP.IsKeyWord("column" "compressed") || HP.IsKeyWord("cc")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_CC) {
			cs.MaskSolverFlags(LinSol::SOLVER_FLAGS_TYPE_MASK);
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_CC);
			pedantic_cout("using column compressed matrix handling for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("column compressed is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}

	/* direct? */
	} else if (HP.IsKeyWord("direct" "access") || HP.IsKeyWord("dir")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_DIR) {
			cs.MaskSolverFlags(LinSol::SOLVER_FLAGS_TYPE_MASK);
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_DIR);
			pedantic_cout("using direct access matrix handling for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		} else {
			pedantic_cerr("direct is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	}

	/* colamd? */
	if (HP.IsKeyWord("colamd")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_COLAMD);
			pedantic_cout("using colamd symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		} else {
			pedantic_cerr("colamd preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* amd ata? */
	} else if (HP.IsKeyWord("mmdata")) {
		silent_cerr("approximate minimum degree solver support is still TODO"
			"task: detect (or import) the MD library;" 
			"uncomment the relevant bits in naivewrap;"
			"remove this check (readlinsol.cc)."
			"Patches welcome"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MMDATA) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_MMDATA);
			pedantic_cout("using mmd symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		} else {
			pedantic_cerr("mmdata preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* minimum degree ?*/
	} else if (HP.IsKeyWord("minimum" "degree")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MDAPLUSAT) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_MDAPLUSAT);
			pedantic_cout("using minimum degree symmetric preordering of A+A^T for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		} else {
			pedantic_cerr("md preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* Reverse Kuthill McKee? */
	} else if (HP.IsKeyWord("rcmk")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE);
			pedantic_cout("using rcmk symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("rcmk preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* king ?*/
	} else if (HP.IsKeyWord("king")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_KING) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_KING);
			pedantic_cout("using king symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("king preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* sloan ? */
	} else if (HP.IsKeyWord("sloan")) {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_KING) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_KING);
			pedantic_cout("using sloan symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("sloan preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
	/* nested dissection ? */
	} else if (HP.IsKeyWord("nested" "dissection")) {
#ifdef USE_METIS
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION);
			pedantic_cout("using nested dissection symmetric preordering for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			pedantic_cerr("nested dissection preordering is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}
#else //!USE_METIS
		silent_cerr("nested dissection permutation not built in;"
			"please configure --with-metis to get it"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif //USE_METIS
	}

	/* multithread? */
	if (HP.IsKeyWord("multi" "thread") || HP.IsKeyWord("mt")) {
		int nThreads = HP.GetInt();

		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) {
			cs.AddSolverFlags(LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT);
			if (nThreads < 1) {
				silent_cerr("illegal thread number, using 1" << std::endl);
				nThreads = 1;
			}
			cs.SetNumThreads(nThreads);

		} else if (nThreads != 1) {
			pedantic_cerr("multithread is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
		}


	} else {
		if (currSolver.s_flags & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) {
			int nThreads = get_nprocs();

			if (nThreads > 1) {
				silent_cout("no multithread requested "
						"with a potential "
						"of " << nThreads << " CPUs"
						<< std::endl);
			}

			cs.SetNumThreads(nThreads);
		}
	}

	if (HP.IsKeyWord("workspace" "size")) {
		integer iWorkSpaceSize = HP.GetInt();
		if (iWorkSpaceSize < 0) {
			iWorkSpaceSize = 0;
		}

		switch (cs.GetSolver()) {
		case LinSol::Y12_SOLVER:
			cs.SetWorkSpaceSize(iWorkSpaceSize);
			break;

		default:
			pedantic_cerr("workspace size is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
			break;
		}
	}

	if (HP.IsKeyWord("pivot" "factor")) {
		doublereal dPivotFactor = HP.GetReal();

		if (currSolver.s_pivot_factor == -1.) {
			pedantic_cerr("pivot factor is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			if (dPivotFactor <= 0. || dPivotFactor > 1.) {
				silent_cerr("pivot factor " << dPivotFactor
						<< " is out of bounds; "
						"using default "
						"(" << currSolver.s_pivot_factor << ")"
						<< std::endl);
				dPivotFactor = currSolver.s_pivot_factor;
			}
			cs.SetPivotFactor(dPivotFactor);
		}

	} else {
		if (currSolver.s_pivot_factor != -1.) {
			cs.SetPivotFactor(currSolver.s_pivot_factor);
		}
	}

	if (HP.IsKeyWord("drop" "tolerance")) {
		doublereal dDropTolerance = HP.GetReal();

		if (currSolver.s_drop_tolerance == -1.) {
			pedantic_cerr("\"drop tolerance\" is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);

		} else {
			if (dDropTolerance < 0.) {
				silent_cerr("drop tolerance " << dDropTolerance
						<< " is out of bounds; "
						"using default "
						"(" << currSolver.s_drop_tolerance << ")"
						<< std::endl);
				dDropTolerance = currSolver.s_drop_tolerance;
			}
			cs.SetDropTolerance(dDropTolerance);
		}

	} else {
		if (currSolver.s_drop_tolerance != -1.) {
			cs.SetDropTolerance(currSolver.s_drop_tolerance);
		}
	}

	if (HP.IsKeyWord("block" "size")) {
		integer blockSize = HP.GetInt();
		if (blockSize < 1) {
			silent_cerr("illegal negative block size; "
					"using default" << std::endl);
			blockSize = 0;
		}

		switch (cs.GetSolver()) {
		case LinSol::UMFPACK_SOLVER:
			cs.SetBlockSize(blockSize);
			break;

		default:
			pedantic_cerr("block size is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
			break;
		}
	}

	if (HP.IsKeyWord("scale")) {
		switch (cs.GetSolver()) {
		case LinSol::NAIVE_SOLVER:
		case LinSol::KLU_SOLVER:
		case LinSol::UMFPACK_SOLVER: {
			SolutionManager::ScaleOpt scale;

			if (HP.IsKeyWord("no")) {
				scale.algorithm = SolutionManager::SCALEA_NONE;

			} else if (HP.IsKeyWord("max") || HP.IsKeyWord("row" "max")) {
				scale.algorithm = SolutionManager::SCALEA_ROW_MAX;

			} else if (HP.IsKeyWord("sum") || HP.IsKeyWord("row" "sum")) {
				scale.algorithm = SolutionManager::SCALEA_ROW_SUM;
			} else if (HP.IsKeyWord("column" "max")) {
				scale.algorithm = SolutionManager::SCALEA_COL_MAX;
			} else if (HP.IsKeyWord("column" "sum")) {
				scale.algorithm = SolutionManager::SCALEA_COL_SUM;
			} else if (HP.IsKeyWord("lapack")) {
				scale.algorithm = SolutionManager::SCALEA_LAPACK;
			} else if (HP.IsKeyWord("iterative")) {
				scale.algorithm = SolutionManager::SCALEA_ITERATIVE;
			} else if (HP.IsKeyWord("row" "max" "column" "max") || HP.IsKeyWord("sinkhorn" "knopp")) {
				scale.algorithm = SolutionManager::SCALEA_ROW_MAX_COL_MAX;
			}

			if (HP.IsKeyWord("scale" "tolerance")) {
				scale.dTol = HP.GetReal();
			}

			if (HP.IsKeyWord("scale" "iterations")) {
				scale.iMaxIter = HP.GetInt();
			}

			if (HP.IsKeyWord("once")) {
				scale.when = SolutionManager::SCALEW_ONCE;

			} else if (HP.IsKeyWord("always")) {
				scale.when = SolutionManager::SCALEW_ALWAYS;

			} else if (HP.IsKeyWord("never")) {
				scale.when = SolutionManager::SCALEW_NEVER;
			}

			switch (scale.when) {
			case SolutionManager::SCALEW_ONCE:
			case SolutionManager::SCALEW_ALWAYS:
				switch (scale.algorithm) {
				case SolutionManager::SCALEA_UNDEF:
					scale.algorithm = SolutionManager::SCALEA_LAPACK; // Restore the original behavior for Naive
					break;

				default:
					// Use the value provided by the input file
					;
				}
				break;

			case SolutionManager::SCALEW_NEVER:
				scale.algorithm = SolutionManager::SCALEA_NONE;
				break;
			}

			if (HP.IsKeyWord("verbose") && HP.GetYesNoOrBool()) {
				scale.uFlags |= SolutionManager::SCALEF_VERBOSE;
			}

			if (HP.IsKeyWord("warnings") && HP.GetYesNoOrBool()) {
				scale.uFlags |= SolutionManager::SCALEF_WARN;
			}

			unsigned uCondFlag = 0;

			if (HP.IsKeyWord("print" "condition" "number")) {
				if (HP.IsKeyWord("norm")) {
					if (HP.IsKeyWord("inf")) {
						uCondFlag = SolutionManager::SCALEF_COND_NUM_INF;
					} else {
						const doublereal dNorm = HP.GetReal();

						if (dNorm != 1.) {
							silent_cerr("Only one norm or infinity norm are supported for condition numbers at line " << HP.GetLineData() << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						uCondFlag = SolutionManager::SCALEF_COND_NUM_1;
					}
				} else {
					uCondFlag = SolutionManager::SCALEF_COND_NUM_1;
				}
			}

			if (uCondFlag != 0 && HP.GetYesNoOrBool()) {
				scale.uFlags |= uCondFlag;
			}

			if (!cs.SetScale(scale)) {
				silent_cerr("Warning: Scale options are not available for "
						    << cs.GetSolverName()
						    << " at line "
						    << HP.GetLineData() << std::endl);
			}

			} break;

		default:
			pedantic_cerr("scale is meaningless for "
					<< currSolver.s_name
					<< " solver" << std::endl);
			break;
		}
	}

	if (HP.IsKeyWord("max" "iterations")) {
		if (!cs.SetMaxIterations(HP.GetInt())) {
			silent_cerr("Warning: iterative refinement is not supported by " << cs.GetSolverName() << " at line " << HP.GetLineData() << std::endl);
		}
	}

	switch (cs.GetSolver()) {
	case LinSol::NAIVE_SOLVER:
		if (!(cs.GetSolverFlags() & LinSol::SOLVER_FLAGS_ALLOWS_COLAMD)) {
			silent_cout("warning: \"naive\" solver should be used with \"colamd\"" << std::endl);
		}
		break;

	// add more warnings...

	default:
		break;
	}
}

std::ostream & RestartLinSol(std::ostream& out, const LinSol& cs)
{
	out << cs.GetSolverName();
	
	const LinSol::solver_t	currSolver = ::solver[cs.GetSolver()];
	
	if (cs.GetSolverFlags() != currSolver.s_default_flags) {
		unsigned f = cs.GetSolverFlags();
		if((f & LinSol::SOLVER_FLAGS_ALLOWS_MAP) == 
			LinSol::SOLVER_FLAGS_ALLOWS_MAP) {
			/*Map*/
			out << ", map ";
		}
		if((f & LinSol::SOLVER_FLAGS_ALLOWS_CC) == 
			LinSol::SOLVER_FLAGS_ALLOWS_CC) {
			/*column compressed*/
			out << ", cc ";
		}
		if((f & LinSol::SOLVER_FLAGS_ALLOWS_DIR) == 
			LinSol::SOLVER_FLAGS_ALLOWS_DIR) {
			/*direct access*/
			out << ", dir ";
		}
		if((f & LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) == 
			LinSol::SOLVER_FLAGS_ALLOWS_COLAMD) {
			/*colamd*/
			out << ", colamd ";
		}
		unsigned nt = cs.GetNumThreads();
		if(((f & LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) == 
			 LinSol::SOLVER_FLAGS_ALLOWS_MT_FCT) && 
			(nt != 1)) {
			/*multi thread*/
			out << ", mt , " << nt;
		}
	}
	integer ws = cs.iGetWorkSpaceSize();
	if(ws > 0) {
		/*workspace size*/
		out << ", workspace size, " << ws;
	}
	doublereal pf = cs.dGetPivotFactor();
	if(pf != -1.) {
		/*pivot factor*/
		out << ", pivot factor, " << pf;
	}
	unsigned bs = cs.GetBlockSize();
	if(bs != 0) {
		/*block size*/
		out << ", block size, " << bs ;
	}
	out << ";" << std::endl;
	return out;
}
