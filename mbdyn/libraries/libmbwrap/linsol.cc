/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#include <myassert.h>

#include <harwrap.h>
#include <mschwrap.h>
#include <y12wrap.h>
#include <umfpackwrap.h>

#include <mbpar.h>

#include <integr.h>

/*
 * Default solver
 */
LinSol::SolverType LinSol::defaultSolver =
#if defined(USE_UMFPACK)
	LinSol::UMFPACK_SOLVER
#elif /* !USE_UMFPACK */ defined(USE_Y12)
	LinSol::Y12_SOLVER
#elif /* !USE_Y12 */ defined(USE_HARWELL)
	LinSol::HARWELL_SOLVER
#elif /* !USE_HARWELL */ defined(USE_MESCHACH)
	LinSol::MESCHACH_SOLVER
#else /* !USE_MESCHACH */
	LinSol::EMPTY_SOLVER
/* FIXME: remove this error if no solver becomes acceptaable :) */
#error "need a solver!"
#endif /* !USE_MESCHACH */
	;

const char *psSolverNames[] = {
	"Harwell",
	"Meschach",
	"Y12",
	"Umfpack",
	"Empty",
	NULL
};

LinSol::LinSol(void)
: CurrSolver(LinSol::defaultSolver),
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
LinSol::Read(MBDynParser &HP, bool bAllowEmpty)
{
   	/* parole chiave */
   	const char* sKeyWords[] = { 
		"harwell",
		"meschach",
		"y12",
		"umfpack",
		"umfpack3",
		"empty",
		NULL
	};

	enum KeyWords {
		HARWELL,
		MESCHACH,
		Y12,
		UMFPACK,
		UMFPACK3,
		EMPTY,

		LASTKEYWORD
	};

   	/* tabella delle parole chiave */
   	KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   	/* cambia la tabella del parser */
   	HP.PutKeyTable(K);

	switch(KeyWords(HP.GetWord())) {
	case MESCHACH:
#ifdef USE_MESCHACH
		CurrSolver = LinSol::MESCHACH_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using meschach sparse LU solver"
				<< std::endl);
		break;
#endif /* USE_MESCHACH */

	case Y12:
#ifdef USE_Y12
		CurrSolver = LinSol::Y12_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using y12 sparse LU solver" << std::endl);
		break;
#endif /* USE_Y12 */

	case UMFPACK3:
		pedantic_cerr("\"umfpack3\" is deprecated; "
				"use \"umfpack\" instead" << std::endl);
	case UMFPACK:
#ifdef USE_UMFPACK
		CurrSolver = LinSol::UMFPACK_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using umfpack sparse LU solver" << std::endl);
		break;
#endif /* USE_UMFPACK */

	case HARWELL:
#ifdef USE_HARWELL
		CurrSolver = LinSol::HARWELL_SOLVER;
		DEBUGLCOUT(MYDEBUG_INPUT,
				"Using harwell sparse LU solver" << std::endl);
		break;
#endif /* USE_HARWELL */

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
					<< psSolverNames[CurrSolver] << " solver" 
					<< std::endl);
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
					<< psSolverNames[CurrSolver] << " solver" 
					<< std::endl);
			break;

		default:
			break;
		}
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "Workspace size: " << iWorkSpaceSize
			<< ", pivor factor: " << dPivotFactor << std::endl);
}

LinSol::SolverType
LinSol::GetSolver(void) const
{
	return CurrSolver;
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

	if (iLWS == 0) {
		iLWS = iWorkSpaceSize;
	}

   	switch (CurrSolver) {
     	case LinSol::Y12_SOLVER: 
#ifdef USE_Y12
      		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			Y12SparseLUSolutionManager,
			Y12SparseLUSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_Y12 */
      		std::cerr << "Configure with --with-y12 "
			"to enable Y12 solver" << std::endl;
      		THROW(ErrGeneric());
#endif /* !USE_Y12 */

	case LinSol::MESCHACH_SOLVER:
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			MeschachSparseLUSolutionManager,
			MeschachSparseLUSolutionManager(iNLD, iLWS,
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
			HarwellSparseLUSolutionManager,
			HarwellSparseLUSolutionManager(iNLD, iLWS,
				dPivotFactor == -1. ? 1. : dPivotFactor));
      		break;
#else /* !USE_HARWELL */
      		std::cerr << "Configure with --with-harwell "
			"to enable Harwell solver" << std::endl;
		THROW(ErrGeneric());
#endif /* !USE_HARWELL */

	case LinSol::UMFPACK_SOLVER:
#ifdef USE_UMFPACK
		SAFENEWWITHCONSTRUCTOR(pCurrSM,
			UmfpackSparseLUSolutionManager,
			UmfpackSparseLUSolutionManager(iNLD, 
				0, dPivotFactor));
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

