/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/* Continua il DataManager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include "ac/float.h"

#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */

#include "dataman.h"
#include "dataman_.h"
#include "modules.h"
#include "readlinsol.h"

#include "drive.h"
#include "drive_.h"
#include "filedrv.h"
#include "presnode.h"
#include "j2p.h"
#include "sah.h"

class NotAllowed {};

/* Legge i dati di controllo */

void
DataManager::ReadControl(MBDynParser& HP,
	const char* sOutputFileName,
	const char* sInputFileName)
{
	DEBUGCOUTFNAME("DataManager::ReadControl");

	/* attach self to parser ... */
	HP.SetDataManager(this);

	/* parole chiave del blocco di controllo */
	const char* sKeyWords[] = {
		"end",
		"control" "data",

		psReadControlNodes[Node::STRUCTURAL],
		psReadControlNodes[Node::ELECTRIC],
		psReadControlNodes[Node::ABSTRACT],
		psReadControlNodes[Node::PARAMETER],
		psReadControlNodes[Node::HYDRAULIC],

		psReadControlElems[Elem::AUTOMATICSTRUCTURAL],
		psReadControlElems[Elem::GRAVITY],
		psReadControlElems[Elem::BODY],
		psReadControlElems[Elem::JOINT],
		psReadControlElems[Elem::BEAM],
		psReadControlElems[Elem::PLATE],
		psReadControlElems[Elem::AIRPROPERTIES],
		psReadControlElems[Elem::ROTOR],
		psReadControlElems[Elem::AEROMODAL],
		psReadControlElems[Elem::AERODYNAMIC],
		psReadControlElems[Elem::FORCE],
		psReadControlElems[Elem::GENEL],
		psReadControlElems[Elem::ELECTRICBULK],
		psReadControlElems[Elem::ELECTRIC],
		psReadControlElems[Elem::HYDRAULIC],
		psReadControlElems[Elem::BULK],
		psReadControlElems[Elem::LOADABLE],
		psReadControlElems[Elem::EXTERNAL],
		psReadControlElems[Elem::SOCKETSTREAM_OUTPUT],
		"RTAI" "output",

		psReadControlDrivers[Drive::FILEDRIVE],

		"loadable" "path",

		"skip" "initial" "joint" "assembly",
		"use",
		"in" "assembly",
		"initial" "stiffness",
		"stiffness",
		"omega" "rotates",
		"initial" "tolerance",
		"tolerance",
		"max" "initial" "iterations",
		"max" "iterations",
		"epsilon",

		"solver",		/* deprecated */
		"linear" "solver",

		"print",

		"title",
		"make" "restart" "file",
		"output" "file" "name",
		"output" "precision",
		"output" "frequency", /* deprecated */
		"output" "meter",
		"output" "results",
		"default" "output",
		"all",
		"none",
		"reference" "frames",

		"default" "scale",

		"read" "solution" "array",

		"select" "timeout",
		"model",

		0
	};


	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		END = 0,
		CONTROLDATA,
		STRUCTURALNODES,
		ELECTRICNODES,
		ABSTRACTNODES,
		PARAMETERNODES,
		HYDRAULICNODES,

		AUTOMATICSTRUCTURAL,
		GRAVITY,
		RIGIDBODIES,
		JOINTS,
		BEAMS,
		PLATES,
		AIRPROPERTIES,
		ROTORS,
		AEROMODALS,
		AERODYNAMICELEMENTS,
		FORCES,
		GENELS,
		ELECTRICBULKELEMENTS,
		ELECTRICELEMENTS,
		HYDRAULICELEMENTS,
		BULKELEMENTS,
		LOADABLEELEMENTS,
		EXTERNALELEMENTS,
		SOCKETSTREAMOUTPUTELEMENTS,
		RTAIOUTPUTELEMENTS,

		FILEDRIVERS,

		LOADABLEPATH,

		SKIPINITIALJOINTASSEMBLY,
		USE,
		INASSEMBLY,
		INITIALSTIFFNESS,
		STIFFNESS,
		OMEGAROTATES,
		INITIALTOLERANCE,
		TOLERANCE,
		MAXINITIALITERATIONS,
		MAXITERATIONS,
		EPSILON,

		SOLVER,	/* deprecated */
		LINEARSOLVER,

		PRINT,

		TITLE,
		MAKERESTARTFILE,
		OUTPUTFILENAME,
		OUTPUTPRECISION,
		OUTPUTFREQUENCY,
		OUTPUTMETER,

		OUTPUTRESULTS,
		DEFAULTOUTPUT,
		ALL,
		NONE,
		REFERENCEFRAMES,

		DEFAULTSCALE,

		READSOLUTIONARRAY,

		SELECTTIMEOUT,
		MODEL,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
 	KeyTable K(HP, sKeyWords);

	KeyWords CurrDesc;
	while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
		switch (CurrDesc) {

		/******** Nodes *********/

		/* Numero di nodi strutturali attesi */
		case STRUCTURALNODES: {
#if defined(USE_STRUCT_NODES)
			int iDmy = HP.GetInt();
			NodeData[Node::STRUCTURAL].iNum = iDmy;
			DofData[DofOwner::STRUCTURALNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Structural nodes: " << iDmy << std::endl);
#else /* USE_STRUCT_NODES */
			silent_cerr("you're not allowed to use structural nodes"
				<< std::endl);
			throw ErrGeneric();
#endif /* USE_STRUCT_NODES */
		} break;

		/* Numero di nodi elettrici attesi */
		case ELECTRICNODES: {
#if defined(USE_ELECTRIC_NODES)
			int iDmy = HP.GetInt();
			NodeData[Node::ELECTRIC].iNum = iDmy;
			DofData[DofOwner::ELECTRICNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Electric nodes: " << iDmy << std::endl);
#else /* USE_ELECTRIC_NODES */
			silent_cerr("you're not allowed to use electric nodes" << std::endl);
			throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
		} break;

		/* Numero di nodi astratti attesi */
		case ABSTRACTNODES: {
#if defined(USE_ELECTRIC_NODES)
			int iDmy = HP.GetInt();
			NodeData[Node::ABSTRACT].iNum = iDmy;
			DofData[DofOwner::ABSTRACTNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Abstract nodes: " << iDmy << std::endl);
#else /* USE_ELECTRIC_NODES */
			silent_cerr("you're not allowed to use abstract nodes" << std::endl);
			throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
		} break;

		/* Numero di nodi astratti attesi */
		case PARAMETERNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::PARAMETER].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Parameter nodes: " << iDmy << std::endl);
 		} break;

		/* Numero di nodi idraulici attesi */
		case HYDRAULICNODES: {
#if defined(USE_HYDRAULIC_NODES)
			int iDmy = HP.GetInt();
			NodeData[Node::HYDRAULIC].iNum = iDmy;
			DofData[DofOwner::HYDRAULICNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic nodes: " << iDmy << std::endl);
#else /* defined(USE_HYDRAULIC_NODES) */
			silent_cerr("you're not allowed to use hydraulic nodes" << std::endl);
			throw ErrGeneric();
#endif /* defined(USE_HYDRAULIC_NODES) */
		} break;

		/******** Elements *********/

		/* Numero di corpi rigidi attesi */
#if defined(USE_STRUCT_NODES)
		case AUTOMATICSTRUCTURAL: {
			int iDmy = HP.GetInt();
#ifdef DEBUG
#if 0
	  ElemData[Elem::AUTOMATICSTRUCTURAL].iExpectedNum = iDmy;
#endif /* 0 */
			DEBUGLCOUT(MYDEBUG_INPUT, "Automatic structural elements expected: "
				<< iDmy << std::endl);
#else
			iDmy = 0;
#endif
		} break;

		/* Accelerazione di gravita' */
		case GRAVITY: {
			if (ElemData[Elem::GRAVITY].iExpectedNum > 0) {
				silent_cerr("warning: gravity acceleration already defined;" << std::endl
					<< "only one definition will be considered" << std::endl);
			}
			ElemData[Elem::GRAVITY].iExpectedNum = 1;
			DEBUGLCOUT(MYDEBUG_INPUT, "Gravity acceleration expected in elements data" << std::endl);
		} break;

		/* Numero di corpi rigidi attesi */
		case RIGIDBODIES: {
			int iDmy = HP.GetInt();
			ElemData[Elem::BODY].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Rigid bodies: " << iDmy << std::endl);
		} break;

		/* Numero di vincoli attesi */
		case JOINTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::JOINT].iExpectedNum = iDmy;
			DofData[DofOwner::JOINT].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Joints: " << iDmy << std::endl);
			if (iDmy > 0 ) {
				bInitialJointAssemblyToBeDone = true;
			}
		} break;

		/* Numero di travi attese */
		case BEAMS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::BEAM].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Beams: " << iDmy << std::endl);
			if (iDmy > 0 ) {
				bInitialJointAssemblyToBeDone = true;
			}
		} break;

#if defined(USE_AERODYNAMIC_ELEMS)
		/* Elementi aerodinamici: proprieta' dell'aria */
		case AIRPROPERTIES: {
			if (ElemData[Elem::AIRPROPERTIES].iExpectedNum > 0) {
				silent_cerr("warning: air properties already defined;" << std::endl
					<< "only one definition will be considered" << std::endl);
			}
			ElemData[Elem::AIRPROPERTIES].iExpectedNum = 1;
			DEBUGLCOUT(MYDEBUG_INPUT, "Air properties expected in elements data" << std::endl);
		} break;

		/* Elementi aerodinamici: rotori */
		case ROTORS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::ROTOR].iExpectedNum = iDmy;
			DofData[DofOwner::ROTOR].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Rotors: " << iDmy << std::endl);
		} break;

		case AEROMODALS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::AEROMODAL].iExpectedNum = iDmy;
			DofData[DofOwner::AEROMODAL].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Aeromodals: " << iDmy << std::endl);
		} break;

		/* Elementi aerodinamici: vari elementi aerodinamici senza dof */
		case AERODYNAMICELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::AERODYNAMIC].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements: " << iDmy << std::endl);
		} break;
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

		/* Numero di forze e coppie attese */
		case FORCES: {
			int iDmy = HP.GetInt();
			ElemData[Elem::FORCE].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Forces: " << iDmy << std::endl);
		} break;

#if defined(USE_ELECTRIC_NODES)
		/* Numero di vincoli attesi */
		case GENELS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::GENEL].iExpectedNum = iDmy;
			DofData[DofOwner::GENEL].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Genels: " << iDmy << std::endl);
		} break;

		/* Numero di elementi elettrici attesi */
		case ELECTRICELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::ELECTRIC].iExpectedNum = iDmy;
			DofData[DofOwner::ELECTRIC].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Electric elements: " << iDmy << std::endl);
		} break;
#endif /* USE_ELECTRIC_NODES */

		/* Numero di elementi idraulici attesi */
		case HYDRAULICELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::HYDRAULIC].iExpectedNum = iDmy;
			DofData[DofOwner::HYDRAULIC].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic elements: " << iDmy << std::endl);
		} break;

		/* Numero di elementi elettrici attesi */
		case BULKELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::BULK].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Bulk elements: " << iDmy << std::endl);
		} break;

		/* Numero di elementi caricabili attesi */
		case LOADABLEELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::LOADABLE].iExpectedNum = iDmy;
			DofData[DofOwner::LOADABLE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Loadable elements: " << iDmy << std::endl);
		} break;

#if defined(HAVE_RUNTIME_LOADING)
		case LOADABLEPATH: {
#if defined(HAVE_LTDL_H)
			bool add(false);

			if (!loadableElemInitialized) {
				module_initialize();
				loadableElemInitialized = true;
			}

			if (HP.IsKeyWord("set")) {
				add = false;
			} else if (HP.IsKeyWord("add")) {
				add = true;
			}

			const char *s = HP.GetFileName();
			if (s == NULL) {
				silent_cerr("missing path "
					"in \"loadable path\" statement "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			if (add) {
				if (lt_dladdsearchdir(s) != 0) {
					silent_cerr("unable to add path "
						"\"" << s << "\" "
						"in \"loadable path\" "
						"statement at line "
						<< HP.GetLineData()
						<< std::endl);
					throw ErrGeneric();
				}
			} else {
				if (lt_dlsetsearchpath(s) != 0) {
					silent_cerr("unable to set path "
						"\"" << s << "\" "
						"in \"loadable path\" "
						"statement at line "
						<< HP.GetLineData()
						<< std::endl);
					throw ErrGeneric();
				}
			}

#else /* HAVE_LTDL_H */
			silent_cerr("loadable path allowed "
				"only in presence of libltdl (ignored)"
				<< std::endl);
				(void)HP.GetStringWithDelims();
#endif /* HAVE_LTDL_H */
		} break;
#endif /* defined(HAVE_RUNTIME_LOADING) */

		case EXTERNALELEMENTS: {
#ifdef USE_EXTERNAL
			int iDmy = HP.GetInt();
			ElemData[Elem::EXTERNAL].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "External elements: " << iDmy
				<< std::endl);
#else /* USE_EXTERNAL */
			silent_cerr("cannot use external elements "
				"when not compiled with -DUSE_EXTERNAL"
				<< std::endl);
#endif /* USE_EXTERNAL */
		} break;

		case SOCKETSTREAMOUTPUTELEMENTS:
		case RTAIOUTPUTELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::SOCKETSTREAM_OUTPUT].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "RTAI output elements: " << iDmy
				<< std::endl);
#ifndef USE_RTAI
			silent_cerr("cannot use RTAI output elements "
				"when not configured --with-rtai"
				<< std::endl);
#endif /* ! USE_RTAI */
		} break;

		/* Numero di drivers attesi */
		case FILEDRIVERS: {
			int iDmy = HP.GetInt();
			DriveData[Drive::FILEDRIVE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "File drivers: " << iDmy << std::endl);
		} break;

		/********* Miscellaneous *********/

#if defined(USE_STRUCT_NODES)
		/* Spegne il flag di assemblaggio iniziale;
		 * di default viene eseguito solo se sono definiti vincoli */
		case SKIPINITIALJOINTASSEMBLY: {
			bSkipInitialJointAssembly = true;
			DEBUGLCOUT(MYDEBUG_INPUT, "Skipping initial joint assembly" << std::endl);
		} break;

		/* Uso di diversi tipi di elementi nell'assemblaggio iniziale */
		case USE:
			while (true) {
				switch (KeyWords(HP.GetWord())) {
				/* Esce dal ciclo */
				case INASSEMBLY:
					goto EndOfUse;

				case RIGIDBODIES:
					ElemData[Elem::BODY].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Rigid bodies will be used "
						"in initial joint assembly"
						<< std::endl);
					break;

				case GRAVITY:
					ElemData[Elem::GRAVITY].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Gravity will be used "
						"in initial joint assembly"
						<< std::endl);
					break;

				case FORCES:
					ElemData[Elem::FORCE].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Forces will be used "
						"in initial joint assembly"
						<< std::endl);
					break;

				/* Lo lascio per backwards compatibility */
				case BEAMS:
#if 0
					ElemData[Elem::BEAM].ToBeUsedInAssembly(true);
#endif /* 0 */
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Beams are used by default "
						"in initial joint assembly"
						<< std::endl);
					break;

#if defined(USE_AERODYNAMIC_ELEMS)
				case AERODYNAMICELEMENTS:
					ElemData[Elem::AERODYNAMIC].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Aerodynamic Elements will be used "
						"in initial joint assembly"
						<< std::endl);

					if (!ElemData[Elem::AIRPROPERTIES].bToBeUsedInAssembly()) {
						ElemData[Elem::AIRPROPERTIES].ToBeUsedInAssembly(true);
					}

					break;
#endif /* USE_AERODYNAMIC_ELEMS */

				case LOADABLEELEMENTS:
					ElemData[Elem::LOADABLE].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Loadable Elements will be used "
						"in initial joint assembly"
						<< std::endl);
					break;

				/* Elemento non autorizzato */
				default:
					silent_cerr("Element type at line "
						<< HP.GetLineData()
						<< " is not allowed; aborting ..."
						<< std::endl);

					throw DataManager::ErrElemNotAllowedInAssembly();

				/* Errore */
				case UNKNOWN:
					silent_cerr("Unknown element type "
						"at line " << HP.GetLineData() << "; "
						"aborting ..." << std::endl);

					throw DataManager::ErrUnknownElem();
				}
			}

EndOfUse:
			break;

		/* Rigidezza delle molle fittizie usate
		 * nell'assemblaggio iniziale;
		 * se viene fornito un solo valore, viene usato sia per posizione
		 * che per velocita';
		 * se ne vengono forniti due, il primo vale per la posizione ed
		 * il secondo per la velocita'
		 */
		case INITIALSTIFFNESS:
			pedantic_cout("\"initial stiffness\" deprecated at line "
				<< HP.GetLineData() << "; use \"stiffness\""
			<< std::endl);

		case STIFFNESS:
			dInitialPositionStiffness =
				HP.GetReal(dDefaultInitialStiffness);

			if (HP.IsArg()) {
				dInitialVelocityStiffness =
					HP.GetReal(dDefaultInitialStiffness);
			} else {
				dInitialVelocityStiffness =
					dInitialPositionStiffness;
			}

			DEBUGLCOUT(MYDEBUG_INPUT, "Initial position stiffness: "
				<< dInitialPositionStiffness << std::endl
				<< "Initial velocity stiffness: "
				<< dInitialVelocityStiffness << std::endl);

			break;

		/* Omega solidale con il nodo o con il rif. globale */
		case OMEGAROTATES:
			if (HP.IsKeyWord("yes")) {
				bOmegaRotates = true;
			} else if (HP.IsKeyWord("no")) {
				bOmegaRotates = false;
			} else {
				silent_cerr("Invalid option at line "
					<< HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric();
			}

			break;

		/* Tolleranza nell'assemblaggio iniziale; viene calcolata come:
		 * sqrt(sum(res^2)/(1.+sum(sol^2)) */
		case INITIALTOLERANCE:
			pedantic_cout("\"initial tolerance\" deprecated at line "
				<< HP.GetLineData() << "; use \"tolerance\""
				<< std::endl);

		case TOLERANCE:
			dInitialAssemblyTol =
				HP.GetReal(dDefaultInitialAssemblyTol);
			DEBUGLCOUT(MYDEBUG_INPUT, "Initial assembly tolerance: "
				<< dInitialAssemblyTol << std::endl);
			break;

		/* Numero massimo di iterazioni nell'assemblaggio iniziale;
		 * di default ne e' consentita solo una, indice di condizioni
		 * iniziali corrette */
		case MAXINITIALITERATIONS:
			pedantic_cout("\"max initial iterations\" deprecated at line "
				<< HP.GetLineData() << "; use \"max iterations\""
				<< std::endl);
		case MAXITERATIONS:
			iMaxInitialIterations =
				HP.GetInt(iDefaultMaxInitialIterations);
			DEBUGLCOUT(MYDEBUG_INPUT, "Max initial iterations: "
				<< iMaxInitialIterations << std::endl);
			break;
#endif /* USE_STRUCT_NODES */

		case EPSILON:
			dEpsilon = HP.GetReal();
			if (dEpsilon <= 0.) {
				silent_cerr("illegal \"epsilon\"=" << dEpsilon
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}
			break;

		case PRINT:
			while (HP.IsArg()) {
				if (HP.IsKeyWord("dof" "stats")) {
					uPrintFlags |= PRINT_DOFSTATS;

				} else if (HP.IsKeyWord("dof" "description")) {
					uPrintFlags |= (PRINT_DOFSTATS | PRINT_DOFDESCRIPTION);

				} else if (HP.IsKeyWord("equation" "description")) {
					uPrintFlags |= (PRINT_DOFSTATS | PRINT_EQDESCRIPTION);

				} else if (HP.IsKeyWord("all")) {
					uPrintFlags = ~PRINT_NONE;

				} else if (HP.IsKeyWord("none")) {
					uPrintFlags = PRINT_NONE;

				} else {
					silent_cerr("unknown print flag at line "
						<< HP.GetLineData() << std::endl);
					throw ErrGeneric();
				}
			}
			break;

		case SOLVER:
			silent_cerr("\"solver\" keyword at line "
				<< HP.GetLineData() << " is deprecated; "
				"use \"linear solver\" instead" << std::endl);
		case LINEARSOLVER:
			ReadLinSol(CurrSolver, HP);
			break;

		/* Titolo */
		case TITLE: {
			ASSERT(sSimulationTitle == NULL);
			if (sSimulationTitle != NULL) {
				SAFEDELETEARR(sSimulationTitle);
			}
			const char* sTmp(HP.GetStringWithDelims());
			SAFESTRDUP(sSimulationTitle, sTmp);
			DEBUGLCOUT(MYDEBUG_INPUT, "Simulation title: "
				"\"" << sSimulationTitle << '"' << std::endl);
		} break;

		/* Crea il file di restart */
		case MAKERESTARTFILE:
			DEBUGLCOUT(MYDEBUG_INPUT, "Restart file will be generated " << std::endl);
			if (HP.IsArg()) {
				if (HP.IsKeyWord("iterations")) {
					RestartEvery = ITERATIONS;
					iRestartIterations = HP.GetInt();
					DEBUGLCOUT(MYDEBUG_INPUT,
						"every " << iRestartIterations
						<< " iterations" << std::endl);
				} else if (HP.IsKeyWord("time")) {
					RestartEvery = TIME;
					dRestartTime = HP.GetReal();
					DEBUGLCOUT(MYDEBUG_INPUT,
						"every " << pdRestartTime[0]
						<< " time units" << std::endl);
				} else if (HP.IsKeyWord("times")) {
					RestartEvery = TIMES;
					iNumRestartTimes = HP.GetInt();
					if (iNumRestartTimes < 1) {
						silent_cerr("illegal number of restart times "
							<< iNumRestartTimes << std::endl);
						throw ErrGeneric();
					}
					SAFENEWARR(pdRestartTimes, doublereal, iNumRestartTimes);
					for (integer i = 0; i < iNumRestartTimes; i++) {
						pdRestartTimes[i] = HP.GetReal();
						DEBUGLCOUT(MYDEBUG_INPUT,
							"    at time "
							<< pdRestartTimes[0]
							<< std::endl);
					}
				} else {
					silent_cerr("Error: unrecognized restart option at line "
						<< HP.GetLineData() << std::endl);

					throw DataManager::ErrGeneric();
				}
			} else {
				RestartEvery = ATEND;
			}

			if (HP.IsKeyWord("with" "solution" "array")) {
				saveXSol = true;
			}
			break;

		case OUTPUTFILENAME:
			silent_cerr("\"output file name\" no longer supported at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();

		case OUTPUTPRECISION: {
			int iPrec = HP.GetInt();
			OutHdl.SetPrecision(iPrec);
		} break;

		case OUTPUTFREQUENCY: {
			integer iFreq = HP.GetInt();
			if (iFreq < 1) {
				silent_cerr("Illegal output frequency " << iFreq
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}
			SAFENEWWITHCONSTRUCTOR(pOutputMeter,
				MeterDriveCaller,
				MeterDriveCaller(&DrvHdl, 
					std::numeric_limits<double>::min(),
					std::numeric_limits<double>::max(), 
					iFreq)
				);
		} break;

		case OUTPUTMETER:
			if (pOutputMeter != 0) {
				silent_cerr("Output meter already defined" << std::endl);
				throw ErrGeneric();
			}
			pOutputMeter = ReadDriveData(this, HP, false);
			break;

		case OUTPUTRESULTS:
			while (HP.IsArg()) {
				/* require support for ADAMS/View .res output */
				if (HP.IsKeyWord("adams")) {
#if defined USE_ADAMS
					ResMode |= RES_ADAMS;

					if (HP.IsArg()) {
						if (HP.IsKeyWord("model" "name")) {
							if (sAdamsModelName != 0) {
								pedantic_cerr("line " << HP.GetLineData()
									<< ": ADAMS output model name "
									"already defined; replacing..."
									<< std::endl);
								SAFEDELETEARR(sAdamsModelName);
								sAdamsModelName = 0;
							}

							const char *tmp = HP.GetStringWithDelims();
							SAFESTRDUP(sAdamsModelName, tmp);
						}

						/* default; conservative: output is very verbose */
						if ( HP.IsKeyWord("velocity")) {
							if (HP.IsKeyWord("yes")) {
								bAdamsVelocity = true;

							} else if (HP.IsKeyWord("no")) {
								bAdamsVelocity = false;

							} else {
								silent_cerr("unknown value "
									"for \"velocity\" flag at line "
									<< HP.GetLineData() << std::endl);
								throw ErrGeneric();
							}
						}

						if (HP.IsKeyWord("acceleration")) {
							if (HP.IsKeyWord("yes")) {
								bAdamsAcceleration = true;

							} else if (HP.IsKeyWord("no")) {
								bAdamsAcceleration = false;

							} else {
								silent_cerr("unknown value "
									"for \"acceleration\" flag at line "
									<< HP.GetLineData() << std::endl);
								throw ErrGeneric();
							}
						}
					}

					if (sAdamsModelName == 0) {
						SAFESTRDUP(sAdamsModelName, "mbdyn");
					}
#else /* !USE_ADAMS */
					silent_cerr("Please rebuild with ADAMS output enabled"
						<< std::endl);
					throw ErrGeneric();
#endif /* USE_ADAMS */
				/* require support for MotionView output */
				} else if (HP.IsKeyWord("motion" "view")) {
#if defined USE_MOTIONVIEW
					ResMode |= RES_MOTIONVIEW;

					/*
					 * add output info
					 */
#else /* !USE_MOTIONVIEW */
					silent_cerr("Please rebuild with MotionView output enabled"
						<< std::endl);
					throw ErrGeneric();
#endif /* USE_MOTIONVIEW */

				} else if (HP.IsKeyWord("netcdf")) {
#ifdef USE_NETCDF
					ResMode |= RES_NETCDF;
					if (HP.IsKeyWord("sync")) {
						bNetCDFsync = true;
					}
#else /* ! USE_NETCDF */
					silent_cerr("Please rebuild with NetCDF output enabled"
						<< std::endl);
					throw ErrGeneric();
#endif /* ! USE_NETCDF */

				} else {
					silent_cerr("unknown \"output results\" "
						"mode at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric();
				}
			}
			break;

		case DEFAULTOUTPUT:
			while (HP.IsArg()) {
				KeyWords CurrDefOut(KeyWords(HP.GetWord()));
				switch (CurrDefOut) {
				case ALL:
					for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
						ElemData[iCnt].DefaultOut(true);
					}
					for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
						NodeData[iCnt].DefaultOut(true);
					}
					break;

				case NONE:
					for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
						ElemData[iCnt].DefaultOut(false);
					}
					for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
						NodeData[iCnt].DefaultOut(false);
					}
					break;

#if defined(USE_STRUCT_NODES)
				case REFERENCEFRAMES:
					bOutputFrames = true;
					break;
#endif /* USE_STRUCT_NODES */

#if defined(USE_STRUCT_NODES)
				case STRUCTURALNODES:
					NodeData[Node::STRUCTURAL].DefaultOut(true);
					break;
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
				case ELECTRICNODES:
					NodeData[Node::ELECTRIC].DefaultOut(true);
					break;

				case ABSTRACTNODES:
					NodeData[Node::ABSTRACT].DefaultOut(true);
					break;
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
				case HYDRAULICNODES:
					NodeData[Node::HYDRAULIC].DefaultOut(true);
					break;
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_STRUCT_NODES)
				case GRAVITY:
					ElemData[Elem::GRAVITY].DefaultOut(true);
					break;

				case RIGIDBODIES:
					ElemData[Elem::BODY].DefaultOut(true);
					break;

				case JOINTS:
					ElemData[Elem::JOINT].DefaultOut(true);
					break;

				case BEAMS:
					ElemData[Elem::BEAM].DefaultOut(true);
					break;

				case PLATES:
					ElemData[Elem::PLATE].DefaultOut(true);
					break;

#if defined(USE_AERODYNAMIC_ELEMS)
				case AIRPROPERTIES:
					ElemData[Elem::AIRPROPERTIES].DefaultOut(true);
					break;

				case ROTORS:
					ElemData[Elem::ROTOR].DefaultOut(true);
					break;

				case AEROMODALS:
					ElemData[Elem::AEROMODAL].DefaultOut(true);
					break;

				case AERODYNAMICELEMENTS:
					ElemData[Elem::AERODYNAMIC].DefaultOut(true);
					break;
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

				case FORCES:
					ElemData[Elem::FORCE].DefaultOut(true);
					break;

#if defined(USE_ELECTRIC_NODES)
				case GENELS:
					ElemData[Elem::GENEL].DefaultOut(true);
					break;

				case ELECTRICBULKELEMENTS:
					ElemData[Elem::ELECTRICBULK].DefaultOut(true);
					break;

				case ELECTRICELEMENTS:
					ElemData[Elem::ELECTRIC].DefaultOut(true);
					break;
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC)
				case HYDRAULICELEMENTS:
					ElemData[Elem::HYDRAULIC].DefaultOut(true);
					break;
#endif /* USE_HYDRAULIC */
				case LOADABLEELEMENTS:
					ElemData[Elem::LOADABLE].DefaultOut(true);
					break;
#ifdef USE_EXTERNAL
				case EXTERNALELEMENTS:
					ElemData[Elem::EXTERNAL].DefaultOut(true);
					break;
#endif /* USE_EXTERNAL */

				case UNKNOWN:
					silent_cerr("warning: unknown output case at line "
						<< HP.GetLineData() << std::endl);
					ASSERT(0);
					break;

				default:
					silent_cerr("case " << sKeyWords[CurrDesc] << " at line "
						<< HP.GetLineData() << " is not allowed" << std::endl);
					ASSERT(0);
					break;
				}
			}
			break;

		case DEFAULTSCALE:
			while (HP.IsArg()) {
				KeyWords CurrDefOut(KeyWords(HP.GetWord()));
				doublereal dScale = HP.GetReal(1.);

				switch (CurrDefOut) {
				case ALL:
					for (int iCnt = 0; iCnt < DofOwner::LASTDOFTYPE; iCnt++) {
						DofData[iCnt].dDefScale= dScale;
					}
					break;

#if defined(USE_STRUCT_NODES)
				case STRUCTURALNODES:
					DofData[DofOwner::STRUCTURALNODE].dDefScale = dScale;
					break;
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
				case ELECTRICNODES:
					DofData[DofOwner::ELECTRICNODE].dDefScale = dScale;
					break;

				case ABSTRACTNODES:
					DofData[DofOwner::ABSTRACTNODE].dDefScale = dScale;
					break;
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
				case HYDRAULICNODES:
					DofData[DofOwner::HYDRAULICNODE].dDefScale = dScale;
					break;
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_STRUCT_NODES)
				case JOINTS:
					DofData[DofOwner::JOINT].dDefScale = dScale;
					break;

#if defined(USE_AERODYNAMIC_ELEMS)
				case ROTORS:
					DofData[DofOwner::ROTOR].dDefScale = dScale;
					break;

				case AEROMODALS:
					DofData[DofOwner::AEROMODAL].dDefScale = dScale;
					break;
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

#if defined(USE_ELECTRIC_NODES)
				case GENELS:
					DofData[DofOwner::GENEL].dDefScale = dScale;
					break;

				case ELECTRICBULKELEMENTS:
					DofData[DofOwner::ELECTRICBULK].dDefScale = dScale;
					break;

				case ELECTRICELEMENTS:
					DofData[DofOwner::ELECTRIC].dDefScale = dScale;
					break;
#endif /* USE_ELECTRIC_NODES */

#if defined(USE_HYDRAULIC_NODES)
				case HYDRAULICELEMENTS:
					DofData[DofOwner::HYDRAULIC].dDefScale = dScale;
					break;
#endif /* USE_HYDRAULIC */

				case LOADABLEELEMENTS:
					DofData[DofOwner::LOADABLE].dDefScale = dScale;
					break;

				case UNKNOWN:
					silent_cerr("warning: unknown output case at line "
						<< HP.GetLineData() << std::endl);
					ASSERT(0);
					break;

				default:
					silent_cerr("case " << sKeyWords[CurrDesc] << " at line "
						<< HP.GetLineData() << " is not allowed" << std::endl);
					ASSERT(0);
					break;
				}
			}
			break;

			/* add more entries ... */
		case READSOLUTIONARRAY:{
			int len = strlen(sInputFileName) + sizeof(".X");
			SAFENEWARR(solArrFileName, char, len);
			snprintf(solArrFileName, len, "%s.X", sInputFileName);
		} break;

		case SELECTTIMEOUT: {
			int timeout = HP.GetInt();
			if (timeout <= 0) {
				silent_cerr("illegal select timeout " << timeout
					<< " at line " << HP.GetLineData()
					<< std::endl);
			} else {
				SocketUsersTimeout = 60*timeout;
			}
		} break;

		case MODEL:
			if (HP.IsKeyWord("static")) {
				bStaticModel = true;
			}
			break;

		case UNKNOWN:
			/*
			 * If description is not in key table the parser
			 * returns UNKNONW, so "default" can be used to
			 * intercept control cases that are not allowed.
			 */
			DEBUGCERR("");
			silent_cerr("unknown description at line "
				<< HP.GetLineData() << std::endl);
			ASSERT(0);
			break;

		default:
			silent_cerr("case " << sKeyWords[CurrDesc] << " at line "
				<< HP.GetLineData() << " is not allowed" << std::endl);
			throw DataManager::ErrGeneric();
		}
	}

	if (KeyWords(HP.GetWord()) != CONTROLDATA) {
		DEBUGCERR("");
		silent_cerr("\"end: control data;\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	/* inizializza l'output handler */
	OutHdl.Init(sOutputFileName, 0);

	/* FIXME: from now on, NetCDF is enabled */
	// OutHdl.ClearText();
	if (bOutput(RES_NETCDF)) {
		OutHdl.SetNetCDF(OutputHandler::NETCDF);
		OutHdl.SetNetCDF(OutputHandler::STRNODES);
	}

	integer iOutputFrequency = 0;
	if (pOutputMeter == 0) {
		SAFENEW(pOutputMeter, OneDriveCaller);
		iOutputFrequency = 1;

	} else {
		MeterDriveCaller *pMDC = dynamic_cast<MeterDriveCaller *>(pOutputMeter);
		if (pMDC != 0) {
			iOutputFrequency = pMDC->iGetSteps();
		}
	}

	if (iOutputFrequency == 0) {
		OutHdl.Log() << "output frequency: custom" << std::endl;

	} else {
		OutHdl.Log() << "output frequency: " << iOutputFrequency << std::endl;
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "End of control data" << std::endl);
} /* End of DataManager::ReadControl() */


/* Legge se un item deve scrivere sull'output handler */

flag
DataManager::fReadOutput(MBDynParser& HP, enum Elem::Type t)
{
	flag fDef = fGetDefaultOutputFlag(t);
	if (!HP.IsKeyWord("output")) {
		return fDef;
	}

	if (HP.IsKeyWord("no")) {
		return flag(0);
	}

	if (HP.IsKeyWord("yes")) {
		return flag(1);
	}

	if (HP.IsKeyWord("default")) {
		return fDef;
	}

	silent_cerr("Unknown output flag for element \""
		<< psElemNames[t] << "\" at line " << HP.GetLineData()
		<< std::endl);
	throw DataManager::ErrGeneric();
} /* End of DataManager::fReadOutput */

flag
DataManager::fReadOutput(MBDynParser& HP, enum Node::Type t)
{
	flag fDef = fGetDefaultOutputFlag(t);
	if (!HP.IsKeyWord("output")) {
		return fDef;
	}

	if (HP.IsKeyWord("no")) {
		return flag(0);
	}

	if (HP.IsKeyWord("yes")) {
		return flag(1);
	}

	if (HP.IsKeyWord("default")) {
		return fDef;
	}

	silent_cerr("Unknown output flag for node \""
		<< psNodeNames[t] << "\" at line " << HP.GetLineData()
		<< std::endl);
	throw DataManager::ErrGeneric();
} /* End of DataManager::fReadOutput */

doublereal
DataManager::dReadScale(MBDynParser& HP, enum DofOwner::Type t)
{
	doublereal d = dGetDefaultScale(t);

	if (!HP.IsKeyWord("scale")) {
		return d;
	}

	if (!HP.IsKeyWord("default")) {
		d = HP.GetReal(d);
	}

	return d;
}


/* legge i nodi e li costruisce */

int
DataManager::ReadScalarAlgebraicNode(MBDynParser& HP,
	unsigned int uLabel, Node::Type type,
	doublereal& dX)
{
	/* verifica di esistenza del nodo */
	if (pFindNode(type, uLabel) != NULL) {
		silent_cerr("line " << HP.GetLineData()
      			<< ": " << psNodeNames[type] << "(" << uLabel
      			<< ") already defined" << std::endl);
		throw DataManager::ErrGeneric();
	}

	if (HP.IsArg()) {
		/* eat keyword "value" */
		if (!HP.IsKeyWord("value")) {
			pedantic_cerr(psNodeNames[type] << "(" << uLabel
     				<< "): initial value specified without "
     				"\"value\" keyword (deprecated)" << std::endl);
		}
		dX = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Initial value x = " << dX
			<< " is supplied" << std::endl);
		return 1;
	}

	return 0;
}

int
DataManager::ReadScalarDifferentialNode(MBDynParser& HP,
	unsigned int uLabel, Node::Type type,
	doublereal& dX, doublereal& dXP)
{
	if (ReadScalarAlgebraicNode(HP, uLabel, type, dX) == 1) {
		if (HP.IsKeyWord("derivative")) {
			dX = HP.GetReal();
			DEBUGLCOUT(MYDEBUG_INPUT,
				"Initial derivative value xp = "
				<< dXP << " is supplied" << std::endl);
			return 1;
		}
	}

	return 0;
}

void
DataManager::ReadNodes(MBDynParser& HP)
{
	DEBUGCOUTFNAME("DataManager::ReadNodes");

	/* parole chiave del blocco di controllo */
	const char* sKeyWords[] = {
		"end",
		"nodes",

		psReadNodesNodes[Node::STRUCTURAL],
		psReadNodesNodes[Node::ELECTRIC],
		psReadNodesNodes[Node::ABSTRACT],
		psReadNodesNodes[Node::PARAMETER],
		psReadNodesNodes[Node::HYDRAULIC],

		"output",
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		END = 0,
		NODES,

		STRUCTURAL,
		ELECTRIC,
		ABSTRACT,
		PARAMETER,
		HYDRAULIC,

		OUTPUT,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* struttura di servizio che conta i nodi tipo per tipo */
	int iNumTypes[Node::LASTNODETYPE];
	for (int i = 0; i < Node::LASTNODETYPE; i++) {
		iNumTypes[i] = NodeData[i].iNum;
	}

	int iMissingNodes = iTotNodes;
	DEBUGLCOUT(MYDEBUG_INPUT, "Expected nodes: " << iMissingNodes << std::endl);

	KeyWords CurrDesc;
	while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
		if (CurrDesc == OUTPUT) {
			DEBUGLCOUT(MYDEBUG_INPUT, "nodes to be output: ");

			Node::Type Typ;
			flag fOutput = 1;
			switch (KeyWords(HP.GetWord())) {

			case STRUCTURAL:
#if defined(USE_STRUCT_NODES)
				DEBUGLCOUT(MYDEBUG_INPUT, "structural" << std::endl);
				Typ = Node::STRUCTURAL;
				if (HP.IsKeyWord("accelerations")) {
					fOutput |= 2;
				}
#else /* USE_STRUCT_NODES */
				silent_cerr("you're not allowed to use structural nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_STRUCT_NODES */
				break;

			case ELECTRIC:
#if defined(USE_ELECTRIC_NODES)
				DEBUGLCOUT(MYDEBUG_INPUT, "electric" << std::endl);
				Typ = Node::ELECTRIC;
#else /* USE_ELECTRIC_NODES */
				silent_cerr("you're not allowed to use electric nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
				break;

			case ABSTRACT:
#if defined(USE_ELECTRIC_NODES)
				DEBUGLCOUT(MYDEBUG_INPUT, "abstract" << std::endl);
				Typ = Node::ABSTRACT;
#else /* USE_ELECTRIC_NODES */
				silent_cerr("you're not allowed to use abstract nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
				break;

			case PARAMETER:
				DEBUGLCOUT(MYDEBUG_INPUT, "parameter" << std::endl);
				Typ = Node::PARAMETER;
				break;

			case HYDRAULIC:
#if defined (USE_HYDRAULIC_NODES)
				DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic" << std::endl);
				Typ = Node::HYDRAULIC;
#else /* defined (USE_HYDRAULIC_NODES) */
				silent_cerr("you're not allowed to use hydraulic nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* defined (USE_HYDRAULIC_NODES) */
				break;

			default:
				silent_cerr("Error: unknown node type, cannot modify output"
					<< std::endl);

				throw DataManager::ErrUnknownNode();
			}

			while (HP.IsArg()) {
				if (HP.IsKeyWord("range")) {
					unsigned int uL = (unsigned int)HP.GetInt();
					unsigned int uEndL = (unsigned int)HP.GetInt();

					if (uEndL < uL) {
						silent_cerr("End label " << uEndL
							<< " must be larger "
							"than or equal to start label " << uL
							<< std::endl);
						throw ErrGeneric();
					}

					for (; uL <= uEndL; uL++) {
						Node* pN = pFindNode(Typ, uL);
						if (pN != NULL) {
							DEBUGLCOUT(MYDEBUG_INPUT,
								"node " << uL << std::endl);
							pN->SetOutputFlag(fOutput);
						}
					}

				} else {
					unsigned int uL = (unsigned int)HP.GetInt();

					Node* pN = pFindNode(Typ, uL);
					if (pN == 0) {
						silent_cerr(psNodeNames[Typ]
							<< "(" << uL << ") "
							"is not defined; "
							"output cannot be modified"
							<< std::endl);

					} else {
						DEBUGLCOUT(MYDEBUG_INPUT, "node " << uL << std::endl);
						pN->SetOutputFlag(fOutput);
					}
				}
			}
		} else {
			/* puntatore al puntatore al nodo */
			Node** ppN = NULL;

			/* legge la label */
			unsigned int uLabel = 0;
			if (CurrDesc >= 0) {
				uLabel = unsigned(HP.GetInt());
			}

			/* Nome del nodo */
			const char *sName = NULL;
			if (HP.IsKeyWord("name")) {
				const char *sTmp = HP.GetStringWithDelims();

				SAFESTRDUP(sName, sTmp);
			}

			/* in base al tipo, avviene l'allocazione */
			switch (CurrDesc) {
			/* Struttura del blocco di lettura dei nodi:
			 *
			 * - test sul numero di nodi letti rispetto a quelli dichiarati
			 *   nel blocco di controllo:
			 *       if(iNumTypes[Node::??]-- <= 0)
			 *         < gestione dell'errore >
			 *
			 * - verifica di esistenza del nodo
			 *       if(pFindNode(Node::??, uLabel) != NULL)
			 *         < gestione dell'errore >
			 *
			 * - lettura dati specifici
			 *
			 * - allocazione e costruzione del nodo:
			 *   - assegnazione del puntatore:
			 *       ppN = NodeData[Node::??].ppFirstNode+
			 *         iNumTypes[Node::??];
			 *
			 *   - allocazione e costruzione:
			 *       SAFENEW((??Node*)*ppN, ??Node(uLabel));
			 *
			 * - correzione del DofOwner relativo al nodo:
			 *       (DofData[DofOwner::??].pFirstDofOwner+
			 *         iNumTypes[Node::??])->iNumDofs =
			 *         DofData[DofOwner::??].iSize;
			 *
			 * - scrittura dei dati specifici dell'oggetto creato.
			 *   In alternativa i dati possono essere passato tutti
			 *   attraverso il costruttore.
			 */


			/* Nodi strutturali */
			case STRUCTURAL: {
#if defined(USE_STRUCT_NODES)
				silent_cout("Reading "
					<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ")"
					<< std::endl);

				/* verifica che non siano gia' stati letti tutti
				 * quelli previsti */
				if (iNumTypes[Node::STRUCTURAL]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ") "
						"exceeds structural nodes number"
						<< std::endl);

					throw DataManager::ErrGeneric();
				}

				/* verifica di esistenza del nodo */
				if (pFindStructNode(uLabel) != NULL) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ") "
						"already defined" << std::endl);

					throw DataManager::ErrGeneric();
				}


				/* lettura dei dati specifici */

				/* allocazione e creazione */
				int i = NodeData[Node::STRUCTURAL].iNum
					- iNumTypes[Node::STRUCTURAL] - 1;
				ppN = NodeData[Node::STRUCTURAL].ppFirstNode + i;
				DofOwner* pDO = DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner + i;

				*ppN = ReadStructNode(this, HP, pDO, uLabel);
#else /* USE_STRUCT_NODES */
				silent_cerr("you're not allowed to use structural nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_STRUCT_NODES */
			} break;

			/* nodi elettrici */
			case ELECTRIC: {
#if defined(USE_ELECTRIC_NODES)
				silent_cout("Reading "
					<< psNodeNames[Node::ELECTRIC] << "(" << uLabel << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::ELECTRIC]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::ELECTRIC] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric();
				}

				/* Initial values */
				doublereal dx(0.);
				doublereal dxp(0.);

				ReadScalarDifferentialNode(HP, uLabel, Node::ELECTRIC, dx, dxp);
				doublereal dScale = dReadScale(HP, DofOwner::ELECTRICNODE);
				flag fOut = fReadOutput(HP, Node::ELECTRIC);

				/* allocazione e creazione */
				int i = NodeData[Node::ELECTRIC].iNum
					- iNumTypes[Node::ELECTRIC] - 1;
				ppN = NodeData[Node::ELECTRIC].ppFirstNode + i;
				DofOwner* pDO = DofData[DofOwner::ELECTRICNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				SAFENEWWITHCONSTRUCTOR(*ppN,
					ElectricNode,
					ElectricNode(uLabel, pDO, dx, dxp, fOut));

#else /* USE_ELECTRIC_NODES */
				silent_cerr("you're not allowed to use electric nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
			} break;

			/* nodi astratti */
			case ABSTRACT: {
#if defined(USE_ELECTRIC_NODES)
				silent_cout("Reading "
					<< psNodeNames[Node::ABSTRACT] << "(" << uLabel << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::ABSTRACT]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::ABSTRACT] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric();
				}

				/* lettura dei dati specifici */
				doublereal dx(0.);
				doublereal dxp(0.);
				ReadScalarDifferentialNode(HP, uLabel, Node::ABSTRACT, dx, dxp);
				doublereal dScale = dReadScale(HP, DofOwner::ABSTRACTNODE);
				flag fOut = fReadOutput(HP, Node::ABSTRACT);

				/* allocazione e creazione */
				int i = NodeData[Node::ABSTRACT].iNum
					- iNumTypes[Node::ABSTRACT] - 1;
				ppN = NodeData[Node::ABSTRACT].ppFirstNode + i;
				DofOwner* pDO = DofData[DofOwner::ABSTRACTNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				SAFENEWWITHCONSTRUCTOR(*ppN,
					AbstractNode,
					AbstractNode(uLabel, pDO, dx, dxp, fOut));

#else /* USE_ELECTRIC_NODES */
				silent_cerr("you're not allowed to use abstract nodes"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_ELECTRIC_NODES */
			} break;

			/* parametri */
			case PARAMETER:
				silent_cout("Reading "
					<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ")"
					<< std::endl);

				/* verifica che non siano gia' stati letti tutti
				 * quelli previsti */
				if (iNumTypes[Node::PARAMETER]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"exceeds parameters number" << std::endl);

					throw DataManager::ErrGeneric();
				}

				/* verifica di esistenza del nodo */
				if (pFindNode(Node::PARAMETER, uLabel) != NULL) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"already defined" << std::endl);

					throw DataManager::ErrGeneric();
				}

				/* bound a elemento */
				if (HP.IsKeyWord("element")) {
					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"will be linked to an element" << std::endl);
					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					int i = NodeData[Node::PARAMETER].iNum
						- iNumTypes[Node::PARAMETER] - 1;
					ppN = NodeData[Node::PARAMETER].ppFirstNode + i;

					SAFENEWWITHCONSTRUCTOR(*ppN,
						Elem2Param,
						Elem2Param(uLabel, &DummyDofOwner, fOut));

				} else if (HP.IsKeyWord("sample" "and" "hold") ||
					HP.IsKeyWord("sample'n'hold"))
				{

					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"is a sample-and-hold" << std::endl);

					ScalarDof SD(ReadScalarDof(this, HP, 0));

					DriveCaller *pDC = NULL;
					SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller,
						TimeDriveCaller(&DrvHdl));

					doublereal dSP = HP.GetReal();
					if (dSP <= 0.) {
						silent_cerr("illegal sample period "
							"for SampleAndHold(" << uLabel << ") "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw ErrGeneric();
					}

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					int i = NodeData[Node::PARAMETER].iNum
						- iNumTypes[Node::PARAMETER] - 1;
					ppN = NodeData[Node::PARAMETER].ppFirstNode + i;

					SAFENEWWITHCONSTRUCTOR(*ppN,
						SampleAndHold,
						SampleAndHold(uLabel, &DummyDofOwner,
							SD.pNode, pDC, dSP, fOut));

				/* strain gage */
				} else if (HP.IsKeyWord("strain" "gage")) {
					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"is a strain gage" << std::endl);

					doublereal dY = HP.GetReal();
					doublereal dZ = HP.GetReal();

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					int i = NodeData[Node::PARAMETER].iNum
						- iNumTypes[Node::PARAMETER] - 1;
					ppN = NodeData[Node::PARAMETER].ppFirstNode + i;

					SAFENEWWITHCONSTRUCTOR(*ppN,
						StrainGageParam,
						StrainGageParam(uLabel, &DummyDofOwner,
							dY, dZ, fOut));

				/* parametro generico */
				} else {

					/* lettura dei dati specifici */
					doublereal dX(0.);
					if (HP.IsArg()) {
						dX = HP.GetReal();
						DEBUGLCOUT(MYDEBUG_INPUT,
							"Initial value x = " << dX
							<< " is supplied" << std::endl);
					}

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					int i = NodeData[Node::PARAMETER].iNum
						- iNumTypes[Node::PARAMETER] - 1;
					ppN = NodeData[Node::PARAMETER].ppFirstNode + i;

					SAFENEWWITHCONSTRUCTOR(*ppN,
						ParameterNode,
						ParameterNode(uLabel, &DummyDofOwner,
							dX, fOut));
				}

				break;

#if defined(USE_HYDRAULIC_NODES)
			/* nodi idraulici */
			case HYDRAULIC: {
				silent_cout("Reading "
					<< psNodeNames[Node::HYDRAULIC] << "(" << uLabel << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::HYDRAULIC]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::HYDRAULIC] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric();
				}

				/* lettura dei dati specifici */
				doublereal dx(0.);
				ReadScalarAlgebraicNode(HP, uLabel, Node::ABSTRACT, dx);
				doublereal dScale = dReadScale(HP, DofOwner::HYDRAULICNODE);
				flag fOut = fReadOutput(HP, Node::HYDRAULIC);

				/* allocazione e creazione */
				int i = NodeData[Node::HYDRAULIC].iNum
					- iNumTypes[Node::HYDRAULIC] - 1;
				ppN = NodeData[Node::HYDRAULIC].ppFirstNode + i;
				DofOwner* pDO = DofData[DofOwner::HYDRAULICNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				SAFENEWWITHCONSTRUCTOR(*ppN,
					PressureNode,
					PressureNode(uLabel, pDO, dx, fOut));

			} break;
#endif /* USE_HYDRAULIC_NODES */


			/* aggiungere eventuali nuovi tipi di nodo */


			/* in caso di tipo errato */
			case UNKNOWN:
				DEBUGCERR("");
				silent_cerr("unknown node type at line "
					<< HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric();

			default:
				DEBUGCERR("");
				silent_cerr("node type " << sKeyWords[CurrDesc]
					<< " at line " << HP.GetLineData()
					<< " is not allowed" << std::endl);

				throw DataManager::ErrGeneric();
			}

			/* verifica dell'allocazione - comune a tutti i casi */
			if (*ppN == NULL) {
				DEBUGCERR("");
				silent_cerr("error in allocation "
					"of " << psNodeNames[CurrDesc] << "(" << uLabel << ")"
					<< std::endl);

				throw ErrMemory();
			}

			if (sName != NULL) {
				(*ppN)->PutName(sName);
				SAFEDELETEARR(sName);
			}

			/* Decrementa i nodi attesi */
			iMissingNodes--;
		}
	}

	if (KeyWords(HP.GetWord()) != NODES) {
		DEBUGCERR("");
		silent_cerr("\"end: nodes;\" expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric();
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric();
	}

	if (iMissingNodes > 0) {
		DEBUGCERR("");
		silent_cerr("warning: " << iMissingNodes
			<< " nodes are missing;" << std::endl);
		for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
			if (iNumTypes[iCnt] > 0) {
				silent_cerr("  " << iNumTypes[iCnt]
					<< ' ' << psNodeNames[iCnt] << std::endl);
			}
		}

		throw DataManager::ErrMissingNodes();
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "End of nodes data" << std::endl);
} /* End of DataManager::ReadNodes() */


void
DataManager::ReadDrivers(MBDynParser& HP)
{
	DEBUGCOUTFNAME("DataManager::ReadDrivers");

	/* parole chiave del blocco di controllo */
	const char* sKeyWords[] = {
		"end",
		"drivers",
		"file",
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		END = 0,
		DRIVERS,
		FILEDRIVE,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* strutture di conteggio dei drivers letti */
	int iNumTypes[Drive::LASTDRIVETYPE];
	for (int i = 0; i < Drive::LASTDRIVETYPE; i++) {
		iNumTypes[i] = DriveData[i].iNum;
	}

	int iMissingDrivers = iTotDrive;
	DEBUGLCOUT(MYDEBUG_INPUT, "Expected drivers: " << iMissingDrivers
		<< std::endl);

	KeyWords CurrDesc;
	while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
		/* puntatore al puntatore al driver */
		Drive** ppD = NULL;

		/* legge la label */
		unsigned int uLabel = 0;
		if (CurrDesc >= 0) {
			uLabel = unsigned(HP.GetInt());
		}

		/* in base al tipo, avviene l'allocazione */
		switch (CurrDesc) {
		/* drivers */
		case FILEDRIVE: {
			silent_cout("Reading "
				<< psDriveNames[Drive::FILEDRIVE] << "(" << uLabel << ")"
				<< std::endl);

			if (iNumTypes[Drive::FILEDRIVE]-- <= 0) {
				DEBUGCERR("");
				silent_cerr("line " << HP.GetLineData() << ": "
					<< psDriveNames[Drive::FILEDRIVE] << "(" << uLabel << ") "
					"exceeds file drivers number" << std::endl);
				throw DataManager::ErrGeneric();
			}

			/* allocazione e creazione */
			int i = DriveData[Drive::FILEDRIVE].iNum
				- iNumTypes[Drive::FILEDRIVE] - 1;
			ppD = DriveData[Drive::FILEDRIVE].ppFirstDrive + i;

			*ppD = ReadFileDriver(this, HP, uLabel);
		} break;

		/* aggiungere qui i nuovi tipi */

		/* in caso di tipo sconosciuto */
		default:
			DEBUGCERR("");
			silent_cerr("unknown drive type at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();

		}

		/* verifica dell'allocazione */
		if (*ppD == NULL) {
			DEBUGCERR("");
			silent_cerr("error in allocation "
				"of " << psDriveNames[Drive::FILEDRIVE] << "(" << uLabel << ")"
				<< std::endl);
			throw ErrMemory();
		}

		/* decrementa il totale degli elementi mancanti */
		iMissingDrivers--;
	}

	if (KeyWords(HP.GetWord()) != DRIVERS) {
		DEBUGCERR("");
		silent_cerr("\"end: drivers;\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	if (iMissingDrivers > 0) {
		silent_cerr("warning, " << iMissingDrivers
			<< " drivers are missing" << std::endl);
		throw DataManager::ErrGeneric();
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "End of drivers data" << std::endl);
} /* End of ReadDrivers */


/* Legge un legame costitutivo monodimensionale */
ConstitutiveLaw1D*
DataManager::ReadConstLaw1D(MBDynParser& HP, ConstLawType::Type& T)
{
	return ReadCL1D(this, HP, T);
}


/* Legge un legame costitutivo tridimensionale */
ConstitutiveLaw3D*
DataManager::ReadConstLaw3D(MBDynParser& HP, ConstLawType::Type& T)
{
	return ReadCL3D(this, HP, T);
}


/* Legge un legame costitutivo esadimensionale */
ConstitutiveLaw6D*
DataManager::ReadConstLaw6D(MBDynParser& HP, ConstLawType::Type& T)
{
	return ReadCL6D(this, HP, T);
}

/* DataManager - end */

Node*
DataManager::ReadNode(MBDynParser& HP, Node::Type type)
{
	unsigned int uNode = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << std::endl);

	/* verifica di esistenza del nodo */
	Node* pNode;
	if ((pNode = pFindNode(type, uNode)) == NULL) {
		silent_cerr(": " << psNodeNames[type] << " node " << uNode
			<< " not defined at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	return pNode;
}

Elem*
DataManager::ReadElem(MBDynParser& HP, Elem::Type type)
{
	unsigned int uElem = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Element " << uElem << std::endl);

	/* verifica di esistenza dell'elemento */
	Elem* pElem;

	if ((pElem = (Elem*)pFindElem(type, uElem)) == NULL) {
		silent_cerr(": " << psElemNames[type] << uElem
			<< " not defined at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	return pElem;
}


static int
GetDofOrder(MBDynParser& HP, Node* pNode, int iIndex)
{
	/*
	 * Ordine del grado di liberta' da considerare
	 * (stato, oppure derivata se esiste)
	 */
	int iOrder = 0;

	if (HP.IsKeyWord("differential")) {
		iOrder = 1;

	} else if (HP.IsKeyWord("order")) {
		iOrder = HP.GetInt();
		if (iOrder < 0 || iOrder > 2) {
			silent_cerr(psNodeNames[pNode->GetNodeType()]
				<< "(" << pNode->GetLabel() << "): "
				"illegal order " << iOrder
				<< " at line " << HP.GetLineData()
				<< std::endl);
	  		throw DataManager::ErrGeneric();

		} else if (iOrder == 2) {
			DynamicStructNode *pStrNode = dynamic_cast<DynamicStructNode *>(pNode);
			if (pStrNode == 0) {
				silent_cerr(psNodeNames[pNode->GetNodeType()]
					<< "(" << pNode->GetLabel() << "): "
					"order " << iOrder << " not allowed "
					"at line " << HP.GetLineData()
					<< std::endl);
	  			throw DataManager::ErrGeneric();
			}
			pStrNode->ComputeAccelerations(true);
		}
	}

	if (iOrder > 0) {
		if (pNode->GetDofType(iIndex-1) != DofOrder::DIFFERENTIAL) {
			silent_cerr(psNodeNames[pNode->GetNodeType()]
				<< "(" << pNode->GetLabel()
				<< "): invalid order for index " << iIndex
				<< " variable at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric();
		}
		return iOrder;
	}

	if (HP.IsKeyWord("algebraic")) {
		return 0;
	} /* else */

	silent_cerr(psNodeNames[pNode->GetNodeType()]
		<< "(" << pNode->GetLabel()
		<< "): unknown or illegal order for index " << iIndex << std::endl
      		<< "(hint: you may need to specify "
		"\"differential\" or \"algebraic\" or \"order <n>\" when referencing " << std::endl
		<< "a generic degree of freedom at line "
		<< HP.GetLineData() << ")" << std::endl);

	throw DataManager::ErrGeneric();
}

ScalarDof
ReadScalarDof(const DataManager* pDM, MBDynParser& HP, flag fOrder)
{
	/* tabella delle parole chiave */
	KeyTable KDof(HP, psReadNodesNodes);

	/* Label del nodo */
	unsigned int uNode = (unsigned int)HP.GetInt();
	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << std::endl);

	/* Tipo del nodo */
	Node::Type Type = Node::Type(HP.GetWord());
	if (Type == Node::UNKNOWN) {
		silent_cerr("line " << HP.GetLineData() << ": "
			"unknown node type" << std::endl);
		throw ErrGeneric();
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "Node type: " << psNodeNames[Type] << std::endl);

	/* verifica di esistenza del nodo */
	Node* pNode = pDM->pFindNode(Type, uNode);
	if (pNode == 0) {
		silent_cerr(psNodeNames[Type] << "(" << uNode << ") not defined"
			" at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
	}

	/* si procura il numero di dof del nodo */
	unsigned int iMaxIndex = pNode->iGetNumDof();
	DEBUGLCOUT(MYDEBUG_INPUT, "max index: " << iMaxIndex << std::endl);

	/* se il nodo ha piu' di un dof, chiede quale dof si desidera */
	unsigned int iIndex = 1;
	if (iMaxIndex > 1) {
		iIndex = HP.GetInt();
		if (iIndex > iMaxIndex) {
			silent_cerr("Illegal index " << iIndex << ", "
				<< psNodeNames[Type] << "(" << uNode << ") "
				"has only " << iMaxIndex << " dofs "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
		DEBUGLCOUT(MYDEBUG_INPUT, "index: " << iIndex << std::endl);
	}

	/* se e' richiesto l'order del dof e se il dof e' differenziale ... */
	int iOrder = 0;
	if (fOrder && pNode->iGetNumDof() > 0
		&& pNode->GetDofType(iIndex-1) == DofOrder::DIFFERENTIAL)
	{
		iOrder = GetDofOrder(HP, pNode, iIndex);
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "order: " << iOrder << std::endl);

	/* se il nodo non e' scalare, alloca un Node2Scalar che lo wrappa */
	if (iMaxIndex > 1) {
		NodeDof nd(pNode->GetLabel(), iIndex-1, pNode);

		pNode = NULL;
		/* Chi dealloca questa memoria? ci vorrebbe l'handle */
		SAFENEWWITHCONSTRUCTOR(pNode, Node2Scalar, Node2Scalar(nd));
		pedantic_cerr(psNodeNames[Type] << "(" << uNode << "): "
			"warning, possibly allocating a NodeDof "
			"that nobody will delete "
			"at line " << HP.GetLineData() << std::endl);
	}

	return ScalarDof(dynamic_cast<ScalarNode *>(pNode), iOrder, iIndex);
}

/* Legge una shape1D;
 * NOTA: il proprietario del puntatore alla Shape la deve distruggere */

#if (defined(USE_STRUCT_NODES) && defined(USE_AERODYNAMIC_ELEMS))
Shape*
ReadShape(MBDynParser& HP)
{
	DEBUGCOUTFNAME("ReadShape");

	const char* sKeyWords[] = {
		"const",
		"linear",
		"piecewise" "linear",
		"parabolic",
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		SHAPECONST = 0,
		LINEAR,
		PIECEWISELINEAR,
		PARABOLIC,
		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di drive */
	KeyWords CurrKeyWord;
	if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
		CurrKeyWord = SHAPECONST;
	}

#ifdef DEBUG
	if (CurrKeyWord >= 0) {
		std::cout << "shape type: " << sKeyWords[CurrKeyWord] << std::endl;
	}
#endif /* DEBUG */

	Shape* pS = NULL;

	switch (CurrKeyWord) {
	/* forma costante */
	case SHAPECONST: {
		/* lettura dei dati specifici */
		doublereal dConst = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Const value: " << dConst << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			ConstShape1D,
			ConstShape1D(dConst));
	} break;

	/* forma lineare */
	case LINEAR: {
		/* lettura dei dati specifici */
		doublereal da0;
		doublereal da1;
		if (HP.IsKeyWord("coefficients")) {
			da0 = HP.GetReal();
			da1 = HP.GetReal();
		} else {
			doublereal dm = HP.GetReal();
			doublereal dp = HP.GetReal();
			da0 = (dp + dm)/2.;
			da1 = (dp - dm)/2;
		}

		DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: "
			<< da0 << ", " << da1 << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			LinearShape1D,
			LinearShape1D(da0, da1));
	} break;

	/* forma lineare a tratti (costante al di fuori del dominio definito) */
	case PIECEWISELINEAR: {
		int np = HP.GetInt();
		if (np <= 0) {
			silent_cerr("Illegal number of points " << np
       				<< " for piecewise linear shape at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		doublereal *px = NULL;
		doublereal *pv = NULL;

		SAFENEWARR(px, doublereal, np);
		SAFENEWARR(pv, doublereal, np);

		px[0] = HP.GetReal();
		if (px[0] < -1. || px[0] > 1.) {
			silent_cerr("Illegal value " << px[0]
				<< " for first point abscissa (must be -1. < x < 1.) "
				"in piecewise linear shape at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
		pv[0] = HP.GetReal();

		for (int i = 1; i < np; i++) {
			px[i] = HP.GetReal();
			if (px[i] <= px[i-1] || px[i] > 1.) {
				silent_cerr("Illegal value " << px[i]
					<< " for point " << i + 1 << " abscissa "
					"(must be " << px[i - 1] << " < x < 1.) "
					"in piecewise linear shape at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}
			pv[i] = HP.GetReal();
		}

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			PiecewiseLinearShape1D,
			PiecewiseLinearShape1D(np, px, pv));
	} break;

	/* forma lineare */
	case PARABOLIC: {
		/* lettura dei dati specifici */
		doublereal dm = HP.GetReal();
		doublereal da0 = HP.GetReal();
		doublereal dp = HP.GetReal();
		doublereal da1 = (dp - dm)/2.;
		doublereal da2 = (dp + dm)/2. - da0;
		DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: "
			<< da0 << ", " << da1 << ", " << da2 << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			ParabolicShape1D,
			ParabolicShape1D(da0, da1, da2));
	} break;

	/* Non c'e' default in quanto all'inizio il default e' stato messo
	 * pari a SHAPECONST */
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw ErrGeneric();
	}

	ASSERT(pS != NULL);
	return pS;
} /* ReadShape */

#endif /* STRUCT && AERODYNAMIC */

