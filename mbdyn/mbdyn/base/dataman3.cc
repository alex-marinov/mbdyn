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

/* Continua il DataManager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <unistd.h>
#include <cfloat>

#if defined(USE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* USE_RUNTIME_LOADING && HAVE_LTDL_H */

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

#include "thermalnode.h"

#include "aeroelem.h"
#include "beam.h"

#include "rbk_impl.h"

class NotAllowed {};

/* Legge i dati di controllo */

void
DataManager::ReadControl(MBDynParser& HP,
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
		psReadControlNodes[Node::THERMAL],
		psReadControlNodes[Node::ABSTRACT],
		psReadControlNodes[Node::PARAMETER],
		psReadControlNodes[Node::HYDRAULIC],

		psReadControlElems[Elem::AUTOMATICSTRUCTURAL],
		psReadControlElems[Elem::GRAVITY],
		psReadControlElems[Elem::BODY],
		psReadControlElems[Elem::JOINT],
		psReadControlElems[Elem::JOINT_REGULARIZATION],
		psReadControlElems[Elem::BEAM],
		psReadControlElems[Elem::PLATE],
		psReadControlElems[Elem::AIRPROPERTIES],
		psReadControlElems[Elem::INDUCEDVELOCITY],
		psReadControlElems[Elem::AEROMODAL],
		psReadControlElems[Elem::AERODYNAMIC],
		psReadControlElems[Elem::FORCE],
		psReadControlElems[Elem::INERTIA],
		psReadControlElems[Elem::GENEL],
		psReadControlElems[Elem::ELECTRICBULK],
		psReadControlElems[Elem::ELECTRIC],
		psReadControlElems[Elem::THERMAL],
		psReadControlElems[Elem::HYDRAULIC],
		psReadControlElems[Elem::BULK],
		psReadControlElems[Elem::LOADABLE],
		psReadControlElems[Elem::EXTERNAL],
		psReadControlElems[Elem::SOCKETSTREAM_OUTPUT],
			"RTAI" "output",	// deprecated
			"rotors",		// deprecated

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
			"accelerations",

		"default" "orientation",
		"default" "beam" "output",
		"default" "aerodynamic" "output",
		"default" "scale",

		"finite" "difference" "jacobian" "meter",

		"read" "solution" "array",

		"select" "timeout",
		"model",

		"rigid" "body" "kinematics",

		0
	};


	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		END = 0,
		CONTROLDATA,
		STRUCTURALNODES,
		ELECTRICNODES,
		THERMALNODES,
		ABSTRACTNODES,
		PARAMETERNODES,
		HYDRAULICNODES,

		AUTOMATICSTRUCTURAL,
		GRAVITY,
		RIGIDBODIES,
		JOINTS,
			JOINT_REGULARIZATIONS,
		BEAMS,
		PLATES,
		AIRPROPERTIES,
		INDUCEDVELOCITYELEMENTS,
		AEROMODALS,
		AERODYNAMICELEMENTS,
		FORCES,
		INERTIA,
		GENELS,
		ELECTRICBULKELEMENTS,
		ELECTRICELEMENTS,
		THERMALELEMENTS,
		HYDRAULICELEMENTS,
		BULKELEMENTS,
		LOADABLEELEMENTS,
		EXTERNALELEMENTS,
		SOCKETSTREAMOUTPUTELEMENTS,
			RTAIOUTPUTELEMENTS,	// deprecated
			ROTORS,			// deprecated

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
			ACCELERATIONS,

		DEFAULTORIENTATION,
		DEFAULTBEAMOUTPUT,
		DEFAULTAERODYNAMICOUTPUT,
		DEFAULTSCALE,

		FDJAC_METER,

		READSOLUTIONARRAY,

		SELECTTIMEOUT,
		MODEL,
		RIGIDBODYKINEMATICS,

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
			int iDmy = HP.GetInt();
			NodeData[Node::STRUCTURAL].iExpectedNum = iDmy;
			DofData[DofOwner::STRUCTURALNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Structural nodes: " << iDmy << std::endl);
		} break;

		/* Numero di nodi elettrici attesi */
		case ELECTRICNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::ELECTRIC].iExpectedNum = iDmy;
			DofData[DofOwner::ELECTRICNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Electric nodes: " << iDmy << std::endl);
		} break;

		/* Numero di nodi termici attesi */
		case THERMALNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::THERMAL].iExpectedNum = iDmy;
			DofData[DofOwner::THERMALNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Thermal nodes: " << iDmy << std::endl);
		} break;

		/* Numero di nodi astratti attesi */
		case ABSTRACTNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::ABSTRACT].iExpectedNum = iDmy;
			DofData[DofOwner::ABSTRACTNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Abstract nodes: " << iDmy << std::endl);
		} break;

		/* Numero di nodi astratti attesi */
		case PARAMETERNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::PARAMETER].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Parameter nodes: " << iDmy << std::endl);
 		} break;

		/* Numero di nodi idraulici attesi */
		case HYDRAULICNODES: {
			int iDmy = HP.GetInt();
			NodeData[Node::HYDRAULIC].iExpectedNum = iDmy;
			DofData[DofOwner::HYDRAULICNODE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Hydraulic nodes: " << iDmy << std::endl);
		} break;

		/******** Elements *********/

		/* Numero di corpi rigidi attesi */
		case AUTOMATICSTRUCTURAL: {
#ifdef DEBUG
			int iDmy =
#endif // DEBUG
			HP.GetInt();
#ifdef DEBUG
#if 0
	  ElemData[Elem::AUTOMATICSTRUCTURAL].iExpectedNum = iDmy;
#endif // 0
			DEBUGLCOUT(MYDEBUG_INPUT, "Automatic structural elements expected: "
				<< iDmy << std::endl);
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

		/* Numero di inerzie attese */
		case INERTIA: {
			int iDmy = HP.GetInt();
			ElemData[Elem::INERTIA].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Inertia: " << iDmy << std::endl);
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

		/* Numero di regolarizzazion vincoli attesi */
		case JOINT_REGULARIZATIONS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::JOINT_REGULARIZATION].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Joint regularizations: " << iDmy << std::endl);
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

		/* Numero di piastre attese */
		case PLATES: {
			int iDmy = HP.GetInt();
			ElemData[Elem::PLATE].iExpectedNum = iDmy;
			DofData[DofOwner::PLATE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Plates: " << iDmy << std::endl);
			if (iDmy > 0 ) {
				bInitialJointAssemblyToBeDone = true;
			}
		} break;

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
		case ROTORS:
			silent_cerr("deprecated \"rotors\", use \"induced velocity elements\" instead at line " << HP.GetLineData() << std::endl);
			// fallthru
		case INDUCEDVELOCITYELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::INDUCEDVELOCITY].iExpectedNum = iDmy;
			DofData[DofOwner::INDUCEDVELOCITY].iNum = iDmy;
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
			DofData[DofOwner::AERODYNAMIC].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Aerodynamic Elements: " << iDmy << std::endl);
		} break;

		/* Numero di forze e coppie attese */
		case FORCES: {
			int iDmy = HP.GetInt();
			ElemData[Elem::FORCE].iExpectedNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Forces: " << iDmy << std::endl);
		} break;

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

		/* Numero di elementi termici attesi */
		case THERMALELEMENTS: {
			int iDmy = HP.GetInt();
			ElemData[Elem::THERMAL].iExpectedNum = iDmy;
			DofData[DofOwner::THERMAL].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "Thermal elements: " << iDmy << std::endl);
		} break;

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

		case LOADABLEPATH: {
#if defined(USE_RUNTIME_LOADING)
			bool add(false);

			if (!moduleInitialized) {
				module_initialize();
				moduleInitialized = true;
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
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (add) {
				if (lt_dladdsearchdir(s) != 0) {
					silent_cerr("unable to add path "
						"\"" << s << "\" "
						"in \"loadable path\" "
						"statement at line "
						<< HP.GetLineData()
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			} else {
				if (lt_dlsetsearchpath(s) != 0) {
					silent_cerr("unable to set path "
						"\"" << s << "\" "
						"in \"loadable path\" "
						"statement at line "
						<< HP.GetLineData()
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

#else // ! USE_RUNTIME_LOADING
			silent_cerr("loadable path allowed "
				"only in presence of libltdl (ignored)"
				<< std::endl);
				(void)HP.GetStringWithDelims();
#endif // ! USE_RUNTIME_LOADING
		} break;

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
			if (CurrDesc == RTAIOUTPUTELEMENTS) {
				silent_cerr("cannot use "
					"RTAI output elements "
					"when not configured "
					"--with-rtai; "
					"using \"socket stream output\" "
					"instead"
					<< std::endl);
			}
#endif /* ! USE_RTAI */
		} break;

		/* Numero di drivers attesi */
		case FILEDRIVERS: {
			int iDmy = HP.GetInt();
			DriveData[Drive::FILEDRIVE].iNum = iDmy;
			DEBUGLCOUT(MYDEBUG_INPUT, "File drivers: " << iDmy << std::endl);
		} break;

		/********* Miscellaneous *********/

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

				case PLATES:
					ElemData[Elem::PLATE].ToBeUsedInAssembly(true);
					DEBUGLCOUT(MYDEBUG_INPUT,
						"Plates will be used "
						"in initial joint assembly"
						<< std::endl);
					break;

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
						<< " is not allowed; aborting..."
						<< std::endl);

					throw DataManager::ErrElemNotAllowedInAssembly(MBDYN_EXCEPT_ARGS);

				/* Errore */
				case UNKNOWN:
					silent_cerr("Unknown element type "
						"at line " << HP.GetLineData() << "; "
						"aborting..." << std::endl);

					throw DataManager::ErrUnknownElem(MBDYN_EXCEPT_ARGS);
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
			if (!HP.GetYesNo(bOmegaRotates)) {
				silent_cerr("Invalid option at line "
					<< HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

		case EPSILON:
			dEpsilon = HP.GetReal();
			if (dEpsilon <= 0.) {
				silent_cerr("illegal \"epsilon\"=" << dEpsilon
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case PRINT:
			while (HP.IsArg()) {
				if (HP.IsKeyWord("dof" "stats")) {
					uPrintFlags |= PRINT_DOF_STATS;

				} else if (HP.IsKeyWord("dof" "description")) {
					uPrintFlags |= (PRINT_DOF_STATS | PRINT_DOF_DESCRIPTION);

				} else if (HP.IsKeyWord("equation" "description")) {
					uPrintFlags |= (PRINT_DOF_STATS | PRINT_EQ_DESCRIPTION);

				} else if (HP.IsKeyWord("description")) {
					uPrintFlags |= (PRINT_DOF_STATS | PRINT_DESCRIPTION);

				} else if (HP.IsKeyWord("element" "connection")) {
					uPrintFlags |= PRINT_EL_CONNECTION;

				} else if (HP.IsKeyWord("node" "connection")) {
					uPrintFlags |= PRINT_NODE_CONNECTION;

				} else if (HP.IsKeyWord("connection")) {
					uPrintFlags |= PRINT_CONNECTION;

				} else if (HP.IsKeyWord("all")) {
					uPrintFlags = ~PRINT_TO_FILE;

				} else if (HP.IsKeyWord("none")) {
					uPrintFlags = PRINT_NONE;

				} else if (HP.IsKeyWord("to" "file")) {
					uPrintFlags |= PRINT_TO_FILE;

				} else {
					silent_cerr("unknown print flag at line "
						<< HP.GetLineData() << std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
						"every " << dRestartTime
						<< " time units" << std::endl);
				} else if (HP.IsKeyWord("times")) {
					RestartEvery = TIMES;
					iNumRestartTimes = HP.GetInt();
					if (iNumRestartTimes < 1) {
						silent_cerr("illegal number of restart times "
							<< iNumRestartTimes << std::endl);
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

		case OUTPUTPRECISION: {
			int iPrec = HP.GetInt();
			OutHdl.SetPrecision(iPrec);
		} break;

		case OUTPUTFREQUENCY: {
			if (pOutputMeter != 0) {
				silent_cerr("Output meter/frequency already defined" << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			integer iFreq = HP.GetInt();
			if (iFreq < 1) {
				silent_cerr("Illegal output frequency " << iFreq
					<< " at line " << HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			SAFENEWWITHCONSTRUCTOR(pOutputMeter,
				MeterDriveCaller,
				MeterDriveCaller(&DrvHdl, 
					-std::numeric_limits<double>::max(),
					std::numeric_limits<double>::max(), 
					iFreq)
				);
		} break;

		case OUTPUTMETER:
			if (pOutputMeter != 0) {
				silent_cerr("Output meter/frequency already defined" << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pOutputMeter = HP.GetDriveCaller(false);
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
						if (HP.IsKeyWord("velocity")) {
							if (!HP.GetYesNo(bAdamsVelocity)) {
								silent_cerr("unknown value "
									"for \"velocity\" flag at line "
									<< HP.GetLineData() << std::endl);
								throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
							}
						}

						if (HP.IsKeyWord("acceleration")) {
							if (!HP.GetYesNo(bAdamsAcceleration)) {
								silent_cerr("unknown value "
									"for \"acceleration\" flag at line "
									<< HP.GetLineData() << std::endl);
								throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
							}
						}
					}

					if (sAdamsModelName == 0) {
						SAFESTRDUP(sAdamsModelName, "mbdyn");
					}
#else /* !USE_ADAMS */
					silent_cerr("Please rebuild with ADAMS output enabled"
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* USE_MOTIONVIEW */

				} else if (HP.IsKeyWord("netcdf")) {
					ResMode |= RES_NETCDF;
					if (HP.IsKeyWord("sync")) {
#ifdef USE_NETCDF
						bNetCDFsync = true;
#endif // USE_NETCDF
					}
					if (HP.IsKeyWord("no" "text")) {
#ifdef USE_NETCDF
						bNetCDFnoText = true;
#endif // USE_NETCDF
					}
#ifndef USE_NETCDF
					silent_cerr("\"netcdf\" ignored; please rebuild with NetCDF output enabled"
						<< std::endl);
#endif /* ! USE_NETCDF */

				} else {
					silent_cerr("unknown \"output results\" "
						"mode at line " << HP.GetLineData()
						<< std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

				case REFERENCEFRAMES:
					bOutputFrames = true;
					break;

				case ACCELERATIONS:
					bOutputAccels = true;
					break;

				case STRUCTURALNODES:
					NodeData[Node::STRUCTURAL].DefaultOut(true);
					break;

				case ELECTRICNODES:
					NodeData[Node::ELECTRIC].DefaultOut(true);
					break;

				case THERMALNODES:
					NodeData[Node::THERMAL].DefaultOut(true);
					break;

				case ABSTRACTNODES:
					NodeData[Node::ABSTRACT].DefaultOut(true);
					break;

				case HYDRAULICNODES:
					NodeData[Node::HYDRAULIC].DefaultOut(true);
					break;

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

				case AIRPROPERTIES:
					ElemData[Elem::AIRPROPERTIES].DefaultOut(true);
					break;

				case ROTORS:
				case INDUCEDVELOCITYELEMENTS:
					ElemData[Elem::INDUCEDVELOCITY].DefaultOut(true);
					break;

				case AEROMODALS:
					ElemData[Elem::AEROMODAL].DefaultOut(true);
					break;

				case AERODYNAMICELEMENTS:
					ElemData[Elem::AERODYNAMIC].DefaultOut(true);
					break;

				case FORCES:
					ElemData[Elem::FORCE].DefaultOut(true);
					break;

				case GENELS:
					ElemData[Elem::GENEL].DefaultOut(true);
					break;

				case ELECTRICBULKELEMENTS:
					ElemData[Elem::ELECTRICBULK].DefaultOut(true);
					break;

				case ELECTRICELEMENTS:
					ElemData[Elem::ELECTRIC].DefaultOut(true);
					break;

				case THERMALELEMENTS:
					ElemData[Elem::THERMAL].DefaultOut(true);
					break;

				case HYDRAULICELEMENTS:
					ElemData[Elem::HYDRAULIC].DefaultOut(true);
					break;
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

		case DEFAULTORIENTATION:
			od = ReadOrientationDescription(HP);
			break;

		case DEFAULTBEAMOUTPUT:
			ReadBeamCustomOutput(this, HP, unsigned(-1), Beam::VISCOELASTIC,
				ElemData[Elem::BEAM].uOutputFlags, ElemData[Elem::BEAM].od);
			break;

		case DEFAULTAERODYNAMICOUTPUT:
			ReadAerodynamicCustomOutput(this, HP, unsigned(-1),
				ElemData[Elem::AERODYNAMIC].uOutputFlags, ElemData[Elem::AERODYNAMIC].od);
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

				case STRUCTURALNODES:
					DofData[DofOwner::STRUCTURALNODE].dDefScale = dScale;
					break;

				case ELECTRICNODES:
					DofData[DofOwner::ELECTRICNODE].dDefScale = dScale;
					break;

				case THERMALNODES:
					DofData[DofOwner::THERMALNODE].dDefScale = dScale;
					break;

				case ABSTRACTNODES:
					DofData[DofOwner::ABSTRACTNODE].dDefScale = dScale;
					break;

				case HYDRAULICNODES:
					DofData[DofOwner::HYDRAULICNODE].dDefScale = dScale;
					break;

				case JOINTS:
					DofData[DofOwner::JOINT].dDefScale = dScale;
					break;

				case ROTORS:
				case INDUCEDVELOCITYELEMENTS:
					DofData[DofOwner::INDUCEDVELOCITY].dDefScale = dScale;
					break;

				case AEROMODALS:
					DofData[DofOwner::AEROMODAL].dDefScale = dScale;
					break;

				case GENELS:
					DofData[DofOwner::GENEL].dDefScale = dScale;
					break;

				case ELECTRICBULKELEMENTS:
					DofData[DofOwner::ELECTRICBULK].dDefScale = dScale;
					break;

				case ELECTRICELEMENTS:
					DofData[DofOwner::ELECTRIC].dDefScale = dScale;
					break;

				case THERMALELEMENTS:
					DofData[DofOwner::THERMAL].dDefScale = dScale;
					break;

				case HYDRAULICELEMENTS:
					DofData[DofOwner::HYDRAULIC].dDefScale = dScale;
					break;

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

		case FDJAC_METER: {
#ifdef MBDYN_FDJAC
			if (pFDJacMeter != 0) {
				silent_cerr("\"finite difference jacobian meter\" already defined" << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
#endif // MBDYN_FDJAC
			DriveCaller *pTmp = HP.GetDriveCaller(false);
#ifdef MBDYN_FDJAC
			pFDJacMeter = pTmp;
#else // !MBDYN_FDJAC
			silent_cerr("warning, \"finite difference jacobian meter\" not supported (ignored)" << std::endl);
			SAFEDELETE(pTmp);
#endif // !MBDYN_FDJAC
		} break;

		case READSOLUTIONARRAY:{
			int len = strlen(sInputFileName) + sizeof(".X");
			SAFENEWARR(solArrFileName, char, len);
			snprintf(solArrFileName, len, "%s.X", sInputFileName);
		} break;

		case SELECTTIMEOUT:
#ifdef USE_SOCKET
			if (HP.IsKeyWord("forever")) {
				SocketUsersTimeout = 0;
			} else {
				int timeout = HP.GetInt();
				if (timeout <= 0) {
					silent_cerr("warning: illegal select timeout " << timeout
						<< " (ignored) at line " << HP.GetLineData()
						<< std::endl);
				} else {
					SocketUsersTimeout = 60*timeout;
				}
			}
#else // ! USE_SOCKET
			silent_cerr("warning: \"select timeout\" not allowed (ignored) "
				"at line " << HP.GetLineData()
				<< " because the current architecture "
				"apparently does not support sockets"
				<< std::endl);
#endif // ! USE_SOCKET
			break;

		case MODEL:
			if (HP.IsKeyWord("static")) {
				bStaticModel = true;
			}
			break;

		case RIGIDBODYKINEMATICS: {
			if (HP.IsKeyWord("const")) {
				Vec3 X(Zero3);
				Mat3x3 R(Eye3);
				Vec3 V(Zero3);
				Vec3 W(Zero3);
				Vec3 XPP(Zero3);
				Vec3 WP(Zero3);

				bool bGot(false);

				if (HP.IsKeyWord("position")) {
					X = HP.GetPosAbs(AbsRefFrame);
					if (!X.IsNull()) {
						bGot = true;
					}
				}

				if (HP.IsKeyWord("orientation")) {
					R = HP.GetRotAbs(AbsRefFrame);
					if (!R.IsExactlySame(Eye3)) {
						bGot = true;
					}
				}

				if (HP.IsKeyWord("velocity")) {
					V = HP.GetVecAbs(AbsRefFrame);
					if (!V.IsNull()) {
						bGot = true;
					}
				}

				if (HP.IsKeyWord("angular" "velocity")) {
					W = HP.GetVecAbs(AbsRefFrame);
					if (!W.IsNull()) {
						bGot = true;
					}
				}

				if (HP.IsKeyWord("acceleration")) {
					XPP = HP.GetVecAbs(AbsRefFrame);
					if (!XPP.IsNull()) {
						bGot = true;
					}
				}

				if (HP.IsKeyWord("angular" "acceleration")) {
					WP = HP.GetVecAbs(AbsRefFrame);
					if (!WP.IsNull()) {
						bGot = true;
					}
				}

				if (!bGot) {
					silent_cerr("null rigid body kinematics "
						"at line " << HP.GetLineData()
						<< std::endl);
					break;
				}

				SAFENEWWITHCONSTRUCTOR(pRBK,
					ConstRigidBodyKinematics,
					ConstRigidBodyKinematics(X, R, V, W, XPP, WP));

			} else if (HP.IsKeyWord("drive")) {
				TplDriveCaller<Vec3> *pXDrv(0);
				TplDriveCaller<Vec3> *pThetaDrv(0);
				TplDriveCaller<Vec3> *pVDrv(0);
				TplDriveCaller<Vec3> *pWDrv(0);
				TplDriveCaller<Vec3> *pXPPDrv(0);
				TplDriveCaller<Vec3> *pWPDrv(0);

				bool bGot(false);

				if (HP.IsKeyWord("position")) {
					pXDrv = ReadDCVecRel(this, HP, AbsRefFrame);
					bGot = true;
				}

				if (HP.IsKeyWord("orientation")) {
					pThetaDrv = ReadDCVecRel(this, HP, AbsRefFrame);
					bGot = true;
				}

				if (HP.IsKeyWord("velocity")) {
					pVDrv = ReadDCVecRel(this, HP, AbsRefFrame);
					bGot = true;
				}

				if (HP.IsKeyWord("angular" "velocity")) {
					pWDrv = ReadDCVecRel(this, HP, AbsRefFrame);
					bGot = true;
				}

				if (HP.IsKeyWord("acceleration")) {
					pXPPDrv = ReadDC3D(this, HP);
					bGot = true;
				}

				if (HP.IsKeyWord("angular" "acceleration")) {
					pWPDrv = ReadDC3D(this, HP);
					bGot = true;
				}

				if (!bGot) {
					silent_cerr("null rigid body kinematics "
						"at line " << HP.GetLineData()
						<< std::endl);
					break;
				}

				SAFENEWWITHCONSTRUCTOR(pRBK,
					DriveRigidBodyKinematics,
					DriveRigidBodyKinematics(pXDrv,
						pThetaDrv, pVDrv, pWDrv,
						pXPPDrv, pWPDrv));

			} else {
				silent_cerr("unknown rigid body kinematics "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (HP.IsArg()) {
				silent_cerr("Semicolon expected at line "
					<< HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		} break;

		case UNKNOWN:
			/*
			 * If description is not in key table the parser
			 * returns UNKNONW, so "default" can be used to
			 * intercept control cases that are not allowed.
			 */
			DEBUGCERR("");
			silent_cerr("unknown description at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			ASSERT(0);
			break;

		default:
			silent_cerr("case " << sKeyWords[CurrDesc] << " at line "
				<< HP.GetLineData() << " is not allowed" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (KeyWords(HP.GetWord()) != CONTROLDATA) {
		DEBUGCERR("");
		silent_cerr("\"end: control data;\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef USE_NETCDF
	/* FIXME: from now on, NetCDF is enabled */
	if (bNetCDFnoText) {
		// enables or disables text output
		OutHdl.ClearText();
	}
#endif // USE_NETCDF

	if (bOutput(RES_NETCDF)) {
		OutHdl.SetNetCDF(OutputHandler::NETCDF);
		OutHdl.SetNetCDF(OutputHandler::STRNODES);
		OutHdl.SetNetCDF(OutputHandler::INERTIA);
		OutHdl.SetNetCDF(OutputHandler::JOINTS);
		OutHdl.SetNetCDF(OutputHandler::BEAMS);
		OutHdl.SetNetCDF(OutputHandler::AERODYNAMIC);
		OutHdl.SetNetCDF(OutputHandler::LOADABLE);
		OutHdl.SetNetCDF(OutputHandler::FORCES);
		// OutHdl.SetNetCDF(OutputHandler::PLATES);
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

doublereal
DataManager::dReadScale(MBDynParser& HP, enum DofOwner::Type t) const
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
	doublereal& dX) const
{
	/* verifica di esistenza del nodo */
	if (pFindNode(type, uLabel) != NULL) {
		silent_cerr("line " << HP.GetLineData()
      			<< ": " << psNodeNames[type] << "(" << uLabel
      			<< ") already defined" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	doublereal& dX, doublereal& dXP) const
{
	if (ReadScalarAlgebraicNode(HP, uLabel, type, dX) == 1) {
		if (HP.IsKeyWord("derivative")) {
			dXP = HP.GetReal();
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
		psReadNodesNodes[Node::THERMAL],
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
		THERMAL,
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
		iNumTypes[i] = NodeData[i].iExpectedNum;
	}

	int iMissingNodes = iTotNodes;
	DEBUGLCOUT(MYDEBUG_INPUT, "Expected nodes: " << iMissingNodes << std::endl);

	NodeVecType::iterator ni = Nodes.begin();

	KeyWords CurrDesc;
	while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
		if (CurrDesc == OUTPUT) {
			DEBUGLCOUT(MYDEBUG_INPUT, "nodes to be output: ");

			Node::Type Typ;
			flag fOutput = 1;
			switch (KeyWords(HP.GetWord())) {

			case STRUCTURAL:
				DEBUGLCOUT(MYDEBUG_INPUT, "structural" << std::endl);
				Typ = Node::STRUCTURAL;
				if (HP.IsKeyWord("accelerations")) {
					fOutput |= 2;
				}
				break;

			case ELECTRIC:
				DEBUGLCOUT(MYDEBUG_INPUT, "electric" << std::endl);
				Typ = Node::ELECTRIC;
				break;

			case THERMAL:
				DEBUGLCOUT(MYDEBUG_INPUT, "thermal" << std::endl);
				Typ = Node::THERMAL;
				break;

			case ABSTRACT:
				DEBUGLCOUT(MYDEBUG_INPUT, "abstract" << std::endl);
				Typ = Node::ABSTRACT;
				break;

			case PARAMETER:
				DEBUGLCOUT(MYDEBUG_INPUT, "parameter" << std::endl);
				Typ = Node::PARAMETER;
				break;

			case HYDRAULIC:
				DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic" << std::endl);
				Typ = Node::HYDRAULIC;
				break;

			default:
				silent_cerr("Error: unknown node type, cannot modify output"
					<< std::endl);

				throw DataManager::ErrUnknownNode(MBDYN_EXCEPT_ARGS);
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
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
			std::string sName;
			if (HP.IsKeyWord("name")) {
				const char *sTmp = HP.GetStringWithDelims();
				sName = sTmp;
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
				silent_cout("Reading "
					<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				/* verifica che non siano gia' stati letti tutti
				 * quelli previsti */
				if (iNumTypes[Node::STRUCTURAL]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ") "
						"exceeds structural nodes number"
						<< std::endl);

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* verifica di esistenza del nodo */
				if (pFindNode(Node::STRUCTURAL, uLabel) != NULL) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::STRUCTURAL] << "(" << uLabel << ") "
						"already defined" << std::endl);

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}


				/* lettura dei dati specifici */

				/* allocazione e creazione */
				int i = NodeData[Node::STRUCTURAL].iExpectedNum
					- iNumTypes[Node::STRUCTURAL] - 1;
				DofOwner* pDO = DofData[DofOwner::STRUCTURALNODE].pFirstDofOwner + i;

				Node *pN = ReadStructNode(this, HP, pDO, uLabel);
				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::STRUCTURAL], uLabel, pN);
				}
			} break;

			/* nodi elettrici */
			case ELECTRIC: {
				silent_cout("Reading "
					<< psNodeNames[Node::ELECTRIC] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::ELECTRIC]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::ELECTRIC] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* Initial values */
				doublereal dx(0.);
				doublereal dxp(0.);

				ReadScalarDifferentialNode(HP, uLabel, Node::ELECTRIC, dx, dxp);
				doublereal dScale = dReadScale(HP, DofOwner::ELECTRICNODE);
				flag fOut = fReadOutput(HP, Node::ELECTRIC);

				/* allocazione e creazione */
				int i = NodeData[Node::ELECTRIC].iExpectedNum
					- iNumTypes[Node::ELECTRIC] - 1;
				DofOwner* pDO = DofData[DofOwner::ELECTRICNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				Node *pN = 0;
				SAFENEWWITHCONSTRUCTOR(pN,
					ElectricNode,
					ElectricNode(uLabel, pDO, dx, dxp, fOut));
				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::ELECTRIC], uLabel, pN);
				}

			} break;

			case THERMAL: {
				silent_cout("Reading "
					<< psNodeNames[Node::THERMAL] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::THERMAL]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::THERMAL] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* Initial values */
				doublereal dx(0.);
				doublereal dxp(0.);

				ReadScalarDifferentialNode(HP, uLabel, Node::THERMAL, dx, dxp);
				doublereal dScale = dReadScale(HP, DofOwner::THERMALNODE);
				flag fOut = fReadOutput(HP, Node::THERMAL);

				/* allocazione e creazione */
				int i = NodeData[Node::THERMAL].iExpectedNum
					- iNumTypes[Node::THERMAL] - 1;
				DofOwner* pDO = DofData[DofOwner::THERMALNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				Node *pN = 0;
				SAFENEWWITHCONSTRUCTOR(pN,
					ThermalNode,
					ThermalNode(uLabel, pDO, dx, dxp, fOut));
				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::THERMAL], uLabel, pN);
				}

			} break;

			/* nodi astratti */
			case ABSTRACT: {
				silent_cout("Reading "
					<< psNodeNames[Node::ABSTRACT] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				bool bAlgebraic(false);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::ABSTRACT]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::ABSTRACT] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (HP.IsKeyWord("algebraic")) {
					bAlgebraic = true;

				} else if (!HP.IsKeyWord("differential")) {
					pedantic_cout("unspecified "
						<< psNodeNames[Node::ABSTRACT] << "(" << uLabel << ") "
						"at line " << HP.GetLineData() << "; "
						"assuming \"differential\"" << std::endl);
				}

				/* lettura dei dati specifici */
				doublereal dx(0.);
				doublereal dxp(0.);

				if (bAlgebraic) {
					ReadScalarAlgebraicNode(HP, uLabel, Node::ABSTRACT, dx);

				} else {
					ReadScalarDifferentialNode(HP, uLabel, Node::ABSTRACT, dx, dxp);
				}
				doublereal dScale = dReadScale(HP, DofOwner::ABSTRACTNODE);
				flag fOut = fReadOutput(HP, Node::ABSTRACT);

				/* allocazione e creazione */
				int i = NodeData[Node::ABSTRACT].iExpectedNum
					- iNumTypes[Node::ABSTRACT] - 1;
				DofOwner* pDO = DofData[DofOwner::ABSTRACTNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				Node *pN = 0;
				if (bAlgebraic) {
					SAFENEWWITHCONSTRUCTOR(pN,
						ScalarAlgebraicNode,
						ScalarAlgebraicNode(uLabel, pDO, dx, fOut));

				} else {
					SAFENEWWITHCONSTRUCTOR(pN,
						ScalarDifferentialNode,
						ScalarDifferentialNode(uLabel, pDO, dx, dxp, fOut));
				}

				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::ABSTRACT], uLabel, pN);
				}
			} break;

			/* parametri */
			case PARAMETER: {
				silent_cout("Reading "
					<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				/* verifica che non siano gia' stati letti tutti
				 * quelli previsti */
				if (iNumTypes[Node::PARAMETER]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"exceeds parameters number" << std::endl);

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* verifica di esistenza del nodo */
				if (pFindNode(Node::PARAMETER, uLabel) != NULL) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"already defined" << std::endl);

					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				Node *pN = 0;

				/* bound a elemento */
				if (HP.IsKeyWord("element")) {
					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"will be linked to an element" << std::endl);
					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					SAFENEWWITHCONSTRUCTOR(pN,
						Elem2Param,
						Elem2Param(uLabel, &DummyDofOwner, fOut));

				} else if (HP.IsKeyWord("sample" "and" "hold") ||
					HP.IsKeyWord("sample'n'hold"))
				{

					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"is a sample-and-hold" << std::endl);

					ScalarDof SD(ReadScalarDof(this, HP, false, false));

					DriveCaller *pDC = NULL;
					SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller,
						TimeDriveCaller(&DrvHdl));

					doublereal dSP = HP.GetReal();
					if (dSP <= 0.) {
						silent_cerr("illegal sample period "
							"for SampleAndHold(" << uLabel << ") "
							"at line " << HP.GetLineData()
							<< std::endl);
						throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					SAFENEWWITHCONSTRUCTOR(pN,
						SampleAndHold,
						SampleAndHold(uLabel, &DummyDofOwner,
							SD.pNode, pDC, dSP, fOut));

				/* strain gage */
				} else if (HP.IsKeyWord("strain" "gage")	/* deprecated */
					|| HP.IsKeyWord("beam" "strain" "gage"))
				{
					DEBUGLCOUT(MYDEBUG_INPUT,
						psNodeNames[Node::PARAMETER] << "(" << uLabel << ") "
						"is a strain gage" << std::endl);

					doublereal dY = HP.GetReal();
					doublereal dZ = HP.GetReal();

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					SAFENEWWITHCONSTRUCTOR(pN,
						StrainGageParam,
						StrainGageParam(uLabel, &DummyDofOwner,
							dY, dZ, fOut));

				/* parametro generico */
				} else {

					/* lettura dei dati specifici */
					doublereal dX(0.);
					if (HP.IsArg()) {
						/* eat keyword "value" */
						if (!HP.IsKeyWord("value")) {
							pedantic_cerr("ParameterNode(" << uLabel << "): "
								"initial value specified without "
    				 				"\"value\" keyword (deprecated)" << std::endl);
						}
						dX = HP.GetReal();
						DEBUGLCOUT(MYDEBUG_INPUT,
							"Initial value x = " << dX
							<< " is supplied" << std::endl);
					}

					flag fOut = fReadOutput(HP, Node::PARAMETER);

					/* allocazione e creazione */
					SAFENEWWITHCONSTRUCTOR(pN,
						ParameterNode,
						ParameterNode(uLabel, &DummyDofOwner,
							dX, fOut));
				}

				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::PARAMETER], uLabel, pN);
				}
				} break;

			/* nodi idraulici */
			case HYDRAULIC: {
				silent_cout("Reading "
					<< psNodeNames[Node::HYDRAULIC] << "(" << uLabel << ( sName.empty() ? "" : ( std::string(", \"") + sName + "\"" ) ) << ")"
					<< std::endl);

				/*
				 * verifica che non siano gia' stati letti tutti
				 * quelli previsti
				 */
				if (iNumTypes[Node::HYDRAULIC]-- <= 0) {
					silent_cerr("line " << HP.GetLineData() << ": "
						<< psNodeNames[Node::HYDRAULIC] << "(" << uLabel << ") "
						"exceeds declared number" << std::endl);
					throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* lettura dei dati specifici */
				doublereal dx(0.);
				ReadScalarAlgebraicNode(HP, uLabel, Node::ABSTRACT, dx);
				doublereal dScale = dReadScale(HP, DofOwner::HYDRAULICNODE);
				flag fOut = fReadOutput(HP, Node::HYDRAULIC);

				/* allocazione e creazione */
				int i = NodeData[Node::HYDRAULIC].iExpectedNum
					- iNumTypes[Node::HYDRAULIC] - 1;
				DofOwner* pDO = DofData[DofOwner::HYDRAULICNODE].pFirstDofOwner + i;
				pDO->SetScale(dScale);

				Node *pN = 0;
				SAFENEWWITHCONSTRUCTOR(pN,
					PressureNode,
					PressureNode(uLabel, pDO, dx, fOut));

				if (pN != 0) {
					ppN = InsertNode(NodeData[Node::HYDRAULIC], uLabel, pN);
				}
			} break;


			/* aggiungere eventuali nuovi tipi di nodo */


			/* in caso di tipo errato */
			case UNKNOWN:
				DEBUGCERR("");
				silent_cerr("unknown node type at line "
					<< HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

			default:
				DEBUGCERR("");
				silent_cerr("node type " << sKeyWords[CurrDesc]
					<< " at line " << HP.GetLineData()
					<< " is not allowed" << std::endl);

				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			/* verifica dell'allocazione - comune a tutti i casi */
			if (*ppN == NULL) {
				DEBUGCERR("");
				silent_cerr("error in allocation "
					"of " << psNodeNames[CurrDesc] << "(" << uLabel << ")"
					<< std::endl);

				throw ErrMemory(MBDYN_EXCEPT_ARGS);
			}

			if (!sName.empty()) {
				(*ppN)->PutName(sName);
			}

			*ni = *ppN;
			ni++;

			/* Decrementa i nodi attesi */
			iMissingNodes--;
		}
	}

	if (KeyWords(HP.GetWord()) != NODES) {
		DEBUGCERR("");
		silent_cerr("\"end: nodes;\" expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

		throw DataManager::ErrMissingNodes(MBDYN_EXCEPT_ARGS);
	}

	/* count & initialize node array */
	unsigned iNumNodes = 0;
	for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		iNumNodes += NodeData[iCnt].NodeContainer.size();
	}

	ASSERT(ni == Nodes.end());
	ASSERT(iNumNodes == Nodes.size());

#if 0
	Nodes.resize(iNumNodes);
	for (int iCnt = 0, iNode = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		for (NodeContainerType::const_iterator p = NodeData[iCnt].NodeContainer.begin();
			p != NodeData[iCnt].NodeContainer.end();
			++p, ++iNode)
		{
			Nodes[iNode] = p->second;
		}
	}
#endif

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
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

		}

		/* verifica dell'allocazione */
		if (*ppD == NULL) {
			DEBUGCERR("");
			silent_cerr("error in allocation "
				"of " << psDriveNames[Drive::FILEDRIVE] << "(" << uLabel << ")"
				<< std::endl);
			throw ErrMemory(MBDYN_EXCEPT_ARGS);
		}

		/* decrementa il totale degli elementi mancanti */
		iMissingDrivers--;
	}

	if (KeyWords(HP.GetWord()) != DRIVERS) {
		DEBUGCERR("");
		silent_cerr("\"end: drivers;\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iMissingDrivers > 0) {
		silent_cerr("warning, " << iMissingDrivers
			<< " drivers are missing" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "End of drivers data" << std::endl);
} /* End of ReadDrivers */

/* DataManager - end */

Node*
DataManager::ReadNode(MBDynParser& HP, Node::Type type) const
{
	integer iNode = HP.GetInt();
	if (iNode < 0) {
		silent_cerr("DataManager::ReadNode: invalid node label " << iNode
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	unsigned int uNode = (unsigned int)iNode;

	DEBUGLCOUT(MYDEBUG_INPUT, "DataManager::ReadNode: " << psNodeNames[type] << "(" << uNode << ")" << std::endl);

	/* verifica di esistenza del nodo */
	Node* pNode = pFindNode(type, uNode);
	if (pNode == 0) {
		silent_cerr("DataManager::ReadNode: " << psNodeNames[type] << "(" << uNode << ")"
			" not defined at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pNode;
}

Elem*
DataManager::ReadElem(MBDynParser& HP, Elem::Type type) const
{
	integer iElem = HP.GetInt();
	if (iElem < 0) {
		silent_cerr("DataManager::ReadElem: invalid node label " << iElem
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	unsigned int uElem = (unsigned int)iElem;

	DEBUGLCOUT(MYDEBUG_INPUT, "DataManager::ReadElem: " << psNodeNames[type] << "(" << uElem << ")" << std::endl);

	/* verifica di esistenza dell'elemento */
	Elem* pElem = dynamic_cast<Elem *>(pFindElem(type, uElem));
	if (pElem == 0) {
		silent_cerr("DataManager::ReadElem: " << psElemNames[type] << "(" << uElem << ")"
			" not defined at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	  		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (iOrder == 2) {
			DynamicStructNode *pStrNode = dynamic_cast<DynamicStructNode *>(pNode);
			if (pStrNode == 0) {
				silent_cerr(psNodeNames[pNode->GetNodeType()]
					<< "(" << pNode->GetLabel() << "): "
					"order " << iOrder << " not allowed "
					"at line " << HP.GetLineData()
					<< std::endl);
	  			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}

ScalarDof
ReadScalarDof(const DataManager* pDM, MBDynParser& HP, bool bDof, bool bOrder)
{
	/* tabella delle parole chiave */
	KeyTable KDof(HP, psReadNodesNodes);

	/* Label del nodo */
	int iNode = HP.GetInt();
	if (iNode < 0) {
		silent_cerr("ReadScalarDof: invalid node label " << iNode
			<< " at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	unsigned uNode = unsigned(iNode);
	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << std::endl);

	/* Tipo del nodo */
	Node::Type Type = Node::Type(HP.GetWord());
	if (Type == Node::UNKNOWN) {
		silent_cerr("ReadScalarDof: unknown node type "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "Node type: " << psNodeNames[Type] << std::endl);

	/* verifica di esistenza del nodo */
	Node* pNode = pDM->pFindNode(Type, uNode);
	if (pNode == 0) {
		silent_cerr(psNodeNames[Type] << "(" << uNode << ") not defined"
			" at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* si procura il numero di dof del nodo */
	unsigned int iMaxIndex = pNode->iGetNumDof();
	DEBUGLCOUT(MYDEBUG_INPUT, "max index: " << iMaxIndex << std::endl);

	if (bDof && iMaxIndex == 0) {
		silent_cerr(psNodeNames[Type] << "(" << uNode << ") must have dofs "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* se il nodo ha piu' di un dof, chiede quale dof si desidera */
	unsigned int iIndex = 1;
	if (iMaxIndex > 1) {
		iIndex = HP.GetInt();
		if (iIndex > iMaxIndex) {
			silent_cerr("Illegal index " << iIndex << ", "
				<< psNodeNames[Type] << "(" << uNode << ") "
				"has only " << iMaxIndex << " dofs "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		DEBUGLCOUT(MYDEBUG_INPUT, "index: " << iIndex << std::endl);
	}

	/* se e' richiesto l'order del dof e se il dof e' differenziale ... */
	int iOrder = 0;
	if (bOrder && pNode->iGetNumDof() > 0
		&& pNode->GetDofType(iIndex-1) == DofOrder::DIFFERENTIAL)
	{
		iOrder = GetDofOrder(HP, pNode, iIndex);
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "order: " << iOrder << std::endl);

	/* se il nodo non e' scalare, alloca un Node2Scalar che lo wrappa */
	if (iMaxIndex > 1) {
		NodeDof nd(iIndex - 1, pNode);

		pNode = 0;
		/* Chi dealloca questa memoria? ci vorrebbe l'handle */
		SAFENEWWITHCONSTRUCTOR(pNode, Node2Scalar, Node2Scalar(nd));
	}

	return ScalarDof(dynamic_cast<ScalarNode *>(pNode), iOrder, iIndex);
}

OrientationDescription
ReadOrientationDescription(MBDynParser& HP)
{
	OrientationDescription dod = UNKNOWN_ORIENTATION_DESCRIPTION;

	if (HP.IsKeyWord("euler" "123")) {
		dod = EULER_123;

	} else if (HP.IsKeyWord("euler" "313")) {
		dod = EULER_313;

	} else if (HP.IsKeyWord("euler" "321")) {
		dod = EULER_321;

	} else if (HP.IsKeyWord("orientation" "vector")) {
		dod = ORIENTATION_VECTOR;

	} else if (HP.IsKeyWord("orientation" "matrix")) {
		dod = ORIENTATION_MATRIX;

	} else {
		silent_cerr("Unknown orientation description "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return dod;
}

OrientationDescription
ReadOptionalOrientationDescription(DataManager *pDM, MBDynParser& HP)
{
	OrientationDescription dod = UNKNOWN_ORIENTATION_DESCRIPTION;

	if (HP.IsKeyWord("orientation" "description")) {
		dod = ReadOrientationDescription(HP);

	} else if (pDM != 0) {
		/* get a sane default */
		dod = pDM->GetOrientationDescription();

	} else {
		dod = EULER_123;
		silent_cerr("Warning, data manager not defined yet, "
			"using default orientation (\"euler 123\") "
			"at line " << HP.GetLineData() << std::endl);
	}

	return dod;
}


