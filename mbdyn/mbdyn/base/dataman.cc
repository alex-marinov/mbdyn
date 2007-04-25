/* $Header$ */
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

/* datamanager */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "friction.h"

extern "C" {
#include <strings.h>
#include <time.h>
}

#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */

#include "solver.h"
#include "invsolver.h"
#include "constltp.h"
#include "dataman_.h"
/* add-ons for math parser */
#include "dofpgin.h"
#include "privpgin.h"
#include "dummypgin.h"
#ifdef USE_TCL
#include "tclpgin.h"
#endif /* USE_TCL */
#include "modelns.h"

/* temporary? */
#include "beam.h"

/* To allow direct loading of modules */
#include <modules.h>

/* To handle  of Elem2Param */
#include "j2p.h"
/* DataManager - begin */

/* linka i singoli DriveCaller al DriveHandler posseduto dal DataManager */
extern void SetDrvHdl(DriveHandler*);

const bool bDefaultInitialJointAssemblyToBeMade(false);
const bool bDefaultSkipInitialJointAssembly(false);
const doublereal dDefaultInitialStiffness = 1.;
const bool bDefaultOmegaRotates(false);
const doublereal dDefaultInitialAssemblyTol = 1.e-6;
const integer iDefaultMaxInitialIterations = 1;

/*
 * costruttore: inizializza l'oggetto, legge i dati e crea le strutture di
 * gestione di Dof, nodi, elementi e drivers.
 */

DataManager::DataManager(MBDynParser& HP,
	unsigned OF,
	Solver* pS,
	doublereal dInitialTime,
	const char* sOutputFileName,
	const char* sInputFileName,
	bool bAbortAfterInput)
:
SolverDiagnostics(OF),
#ifdef USE_MULTITHREAD
nThreads(0),
#endif /* USE_MULTITHREAD */
MathPar(HP.GetMathParser()),
pSolver(pS),
DrvHdl(HP.GetMathParser()),
OutHdl(sOutputFileName),        /* ..takes output filename as argument */
pTime(0),
pXCurr(0), pXPrimeCurr(0), 
/* Inverse Dynamics: */
pXPrimePrimeCurr(0),
#if defined(USE_STRUCT_NODES)
bInitialJointAssemblyToBeDone(bDefaultInitialJointAssemblyToBeMade),
bSkipInitialJointAssembly(bDefaultSkipInitialJointAssembly),
bOutputFrames(false),
bOutputAccels(false),
dInitialPositionStiffness(dDefaultInitialStiffness),
dInitialVelocityStiffness(dDefaultInitialStiffness),
bOmegaRotates(bDefaultOmegaRotates),
dInitialAssemblyTol(dDefaultInitialAssemblyTol),
iMaxInitialIterations(iDefaultMaxInitialIterations),
dEpsilon(1.),
CurrSolver(pS->GetLinearSolver()),
bStaticModel(false),
/* auto-detect if running inverse dynamics */
bInverseDynamics(dynamic_cast<InverseSolver *>(pS) != 0),
#endif /* USE_STRUCT_NODES */
#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
loadableElemInitialized(false),
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */
uPrintFlags(PRINT_NONE),		/* Morandini, 2003-11-17 */
sSimulationTitle(0),
RestartEvery(NEVER),
iRestartIterations(0),
dRestartTime(0.),
pdRestartTimes(0),
iNumRestartTimes(0),
iCurrRestartTime(0),
iCurrRestartIter(0),
dLastRestartTime(dInitialTime),
saveXSol(false),
solArrFileName(0),
pOutputMeter(0),
iOutputCount(0),
ResMode(RES_TEXT),
#ifdef USE_NETCDF
/* NetCDF stuff */
bNetCDFsync(false),
Var_Step(0),
Var_Time(0),
Var_TimeStep(0),
#endif /* USE_NETCDF */
od(EULER_123),
#if defined(USE_ADAMS) || defined(USE_MOTIONVIEW)
iOutputBlock(1),
#endif /* defined(USE_ADAMS) || defined(USE_MOTIONVIEW) */
#ifdef USE_ADAMS
sAdamsModelName(0),
bAdamsVelocity(false),
bAdamsAcceleration(false),
iAdamsOutputNodes(0),
iAdamsOutputParts(0),
adamsNoab(0),
#endif /* USE_ADAMS */

/* ElemManager */
ElemIter(),
ppDrive(0),
iTotDrive(0),
iMaxWorkNumRows(0),
iMaxWorkNumCols(0),
iWorkIntSize(0),
iWorkDoubleSize(0),
pWorkMatA(0),
pWorkMatB(0),
pWorkMat(0),
pWorkVec(0),

/* NodeManager */
NodeIter(),
iTotNodes(0),
ppNodes(0),

/* DofManager */
iTotDofOwners(0),
pDofOwners(0),
iTotDofs(0),
pDofs(0),
DofIter(),
SocketUsersTimeout(0)
{
	DEBUGCOUTFNAME("DataManager::DataManager");

	/* pseudocostruttori */
	ElemManager();
	NodeManager();
	DofManager();

	/* registra il plugin per i dofs */
	HP.GetMathParser().RegisterPlugIn("dof", dof_plugin, this);

	/* registra il plugin per i dati privati dei nodi */
	HP.GetMathParser().RegisterPlugIn("node", node_priv_plugin, this);

	/* registra il plugin per i dati privati degli elementi */
	HP.GetMathParser().RegisterPlugIn("element", elem_priv_plugin, this);

#ifdef USE_TCL
	/* registra il plugin per il tcl */
	HP.GetMathParser().RegisterPlugIn("tcl", tcl_plugin, 0);
#else /* !USE_TCL */
	HP.GetMathParser().RegisterPlugIn("tcl", dummy_plugin,
		(void *)"configure with --with-tcl to use tcl plugin");
#endif /* USE_TCL */

	/* registra il namespace del modello */
	HP.GetMathParser().RegisterNameSpace(new ModelNameSpace(this));

	/* Setta il tempo al valore iniziale */
	pTime = (Var *)MathPar.GetSymbolTable().Get("Time");
	ASSERT(pTime != 0);
	if (pTime == 0) {
		DEBUGCERR("");
		silent_cerr("error in setting Time symbol" << std::endl);

		throw DataManager::ErrGeneric();
	}
	SetTime(dInitialTime, false);

	DEBUGLCOUT(MYDEBUG_INIT, "Global symbol table:"
		<< MathPar.GetSymbolTable() << std::endl);

	/*
	 * Possiede MathParser, con relativa SymbolTable.
	 * Crea ExternKeyTable, MBDynParser,
	 * e legge i dati esterni. Quindi, quando trova
	 * i dati di controllo, chiama la relativa
	 * funzione di lettura (distinta per comodita')
	 */

	/* initialize pre-loaded loadable elements... */
#ifdef STATIC_MODULES
	SetLoadableElemModule(module_wheel2_lc.name, &module_wheel2_lc);
#endif /* STATIC_MODULES */

	/* parole chiave */
	const char* sKeyWords[] = {
		"begin",
		"controldata",
		"scalarfunction",
		"nodes",
		"elements",
		"drivers",
		"output",
		0
	};

	/* enum delle parole chiave */
	enum KeyWords {
		BEGIN = 0,
		CONTROLDATA,
		SCALARFUNCTION,
		NODES,
		ELEMENTS,
		DRIVERS,
		OUTPUT,
		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	KeyWords CurrDesc = KeyWords(HP.GetDescription());
	/* legge i dati di controllo */
	if (CurrDesc != BEGIN) {
		DEBUGCERR("");
		silent_cerr("<begin> expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric();
	}

	if (KeyWords(HP.GetWord()) != CONTROLDATA) {
		DEBUGCERR("");
		silent_cerr("<begin: control data;> expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric();
	}

	ReadControl(HP, sInputFileName);
	try {
		CurrDesc = KeyWords(HP.GetDescription());
	} catch (EndOfFile) {
		NO_OP;
	}

	/* fine lettura dati di controllo */



	/*
	 * a questo punto ElemManager contiene i numeri totali di elementi;
	 * NodeManager contiene i numeri totali di nodi;
	 * DofManager contiene i numeri totali di potenziali possessori di gradi
	 * di liberta'; quindi si puo' cominciare ad allocare matrici
	 */



	/* Costruzione struttura DofData e creazione array DofOwner */
	DofDataInit();

	/* Costruzione struttura NodeData e creazione array Node* */
	NodeDataInit();

	/*
	 * legge i nodi, crea gli item della struttura ppNodes
	 * e contemporaneamente aggiorna i dof
	 */
	if (iTotNodes > 0) {
		if (CurrDesc != BEGIN) {
			DEBUGCERR("");
			silent_cerr("<begin> expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		if (KeyWords(HP.GetWord()) != NODES) {
			DEBUGCERR("");
			silent_cerr("<begin: nodes;> expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		ReadNodes(HP);
		try {
			CurrDesc = KeyWords(HP.GetDescription());
		} catch (EndOfFile) {}
	} else {
		DEBUGCERR("");
		silent_cerr("warning, no nodes are defined" << std::endl);
	}
	/* fine lettura nodi */

	/*
	 * Costruzione struttura ElemData, DriveData
	 * e creazione array Elem* e Drive*
	 */
	ElemDataInit();

	/* legge i drivers, crea la struttura ppDrive */
	if (iTotDrive > 0) {
		if (CurrDesc != BEGIN) {
			DEBUGCERR("");
			silent_cerr("\"begin\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		if (KeyWords(HP.GetWord()) != DRIVERS) {
			DEBUGCERR("");
			silent_cerr("\"begin: drivers;\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		ReadDrivers(HP);
		try {
			CurrDesc = KeyWords(HP.GetDescription());
		} catch (EndOfFile) {}
	} else {
		DEBUGCERR("warning, no drivers are defined" << std::endl);
	}

	/* fine lettura drivers */


	/*
	 * legge gli elementi, crea la struttura ppElems
	 * e contemporaneamente aggiorna i dof
	 */
	if (!Elems.empty()) {
		if (CurrDesc != BEGIN) {
			DEBUGCERR("");
			silent_cerr("\"begin\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		if (KeyWords(HP.GetWord()) != ELEMENTS) {
			DEBUGCERR("");
			silent_cerr("\"begin: elements;\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric();
		}

		ReadElems(HP);
#if 0
		/* FIXME: we don't check extra statements after "end: elements;"
		 * because there might be more... */
		try {
			CurrDesc = KeyWords(HP.GetDescription());
		} catch (EndOfFile) {}
#endif
	} else {
		DEBUGCERR("");
		silent_cerr("warning, no elements are defined" << std::endl);
	}

#ifdef USE_STRUCT_NODES
	if (bOutputFrames) {
		OutHdl.Open(OutputHandler::REFERENCEFRAMES);
		HP.OutputFrames(OutHdl.ReferenceFrames());
	}

#ifdef USE_AERODYNAMIC_ELEMS
	if (!ElemData[Elem::AIRPROPERTIES].ElemMap.empty()) {
		OutHdl.Open(OutputHandler::AIRPROPS);
	}
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

	/* fine lettura elementi */

	for (int i = 0; i < Node::LASTNODETYPE; i++) {
		if (NodeData[i].iNum > 0
			&& NodeData[i].OutFile != OutputHandler::UNKNOWN)
		{
			OutHdl.Open(NodeData[i].OutFile);
		}
	}

	for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
		if (!ElemData[i].ElemMap.empty()
			&& ElemData[i].OutFile != OutputHandler::UNKNOWN)
		{
			OutHdl.Open(ElemData[i].OutFile);
		}
	}

	/* Inizializza il drive handler */
	DrvHdl.iRandInit(0);
	DrvHdl.iMeterInit(0);

	/* Verifica dei dati di controllo */
#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
		/* mostra in modo succinto il numero di DofOwner per tipo */
		for (int i = 0; i < DofOwner::LASTDOFTYPE; i++) {
			std::cout << "DofType " << i << " (" << psDofOwnerNames[i]
				<< "), n. of owners: " << DofData[i].iNum << std::endl;
		}

		/* mostra in modo succinto il numero di nodi per tipo */
		for (int i = 0; i < Node::LASTNODETYPE; i++) {
			std::cout << "NodeType " << i << " (" << psNodeNames[i]
				<< "), n. of nodes: " << NodeData[i].iNum << std::endl;
		}

		/* mostra in modo succinto il numero di elementi per tipo */
		for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
			std::cout << "Element Type " << i << " (" << psElemNames[i]
				<< "), n. of elems: " << ElemData[i].ElemMap.size() << std::endl;
		}

		/* mostra in modo succinto il numero di drivers per tipo */
		for (int i = 0; i < Drive::LASTDRIVETYPE; i++) {
			std::cout << "DriveType " << i << " (" << psDriveNames[i]
				<< "), n. of drivers: " << DriveData[i].iNum << std::endl;
		}
	}
#endif /* DEBUG */

#ifdef USE_ADAMS
	/* Se richiesto, inizializza il file di output AdamsRes */
	if (bAdamsOutput()) {
		AdamsResOutputInit();
		AdamsResOutput(iOutputBlock, "INPUT", "MBDyn");
	}
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW
	if (bMotionViewOutput()) {
		MotionViewResOutputInit(sOutputFileName);
		MotionViewResOutput(1, "INPUT", "MBDyn");
	}
#endif /* USE_MOTIONVIEW */

	if (bAbortAfterInput) {
		silent_cout("Only input is required" << std::endl);
		return;
	}

	/* waits for all pending sockets to connect */
	WaitSocketUsers();

	/* Qui intercetto la struttura dei Dof prima che venga costruita e modifico
	 * in modo che sia adatta all'assemblaggio dei vincoli; quindi la resetto
	 * e la lascio generare correttamente.
	 * L'assemblaggio iniziale viene gestito dal DataManager perche'
	 * concettualmente la consistenza dei dati deve essere assicurata da lui,
	 * che conosce gli aspetti fisici del problema
	 */

#if defined(USE_STRUCT_NODES)
	if (bInitialJointAssemblyToBeDone) {
		if (!bSkipInitialJointAssembly && !bInverseDynamics) {
			InitialJointAssembly();
		} else {
			silent_cout("Skipping initial joints assembly" << std::endl);
		}
	} else {
		silent_cout("No initial assembly is required since no joints are defined"
			<< std::endl);
	}
#endif /* USE_STRUCT_NODES */

	/* Costruzione dei dati dei Dof definitivi da usare nella simulazione */

	bool bIsSquare;

	/* Aggiornamento DofOwner degli elementi (e dei nodi) */
	if(bInverseDynamics)	{
		bIsSquare = InverseDofOwnerSet();
	} else	{
		DofOwnerSet();
	}

	/* Verifica dei dati di controllo */
#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
		/* mostra in modo succinto i DofOwners */
		int k = 0;
		for (int i = 0; i < DofOwner::LASTDOFTYPE; i++) {
			std::cout << "DofType " << i << ':' << std::endl;
			for (int j = 0; j < DofData[i].iNum; j++) {
				std::cout << "DofOwner " << j << ", n. of dofs: "
					<< pDofOwners[k++].iNumDofs << std::endl;
			}
		}
	}
#endif /* DEBUG */

	/* a questo punto i dof sono completamente scritti in quanto a numero.
	 * Rimane da sistemare il vettore con gli indici "interni" ed il tipo */

	/* Costruzione array DofOwner e Dof */
	if(bInverseDynamics)	{
		InverseDofInit(bIsSquare);
	} else	{
		DofInit();
	}

	/* Aggiornamento Dof di proprieta' degli elementi e dei nodi */
//	if(bInverseDynamics)	{
//		InverseDofOwnerInit();
//	} else	{
		DofOwnerInit();
//	}

	/* Creazione strutture di lavoro per assemblaggio e residuo */
	ElemAssInit();

#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
		for (int iCnt = 0; iCnt < iTotDofs; iCnt++) {
			std::cout << "Dof " << std::setw(4) << iCnt+1 << ": order "
				<< pDofs[iCnt].Order << std::endl;
		}
	}
#endif /* DEBUG */

	/* Se richiesto, esegue l'output delle condizioni iniziali */
#if defined(USE_ADAMS) || defined(USE_MOTIONVIEW)
	iOutputBlock++;
#endif /* defined(USE_ADAMS) || defined(USE_MOTIONVIEW) */

#ifdef USE_ADAMS
	if (bAdamsOutput()) {
		AdamsResOutput(iOutputBlock, "INITIAL CONDITIONS", "MBDyn");
	}
#endif /* USE_ADAMS */

	/* Se richiesto, esegue l'output delle condizioni iniziali */
#ifdef USE_MOTIONVIEW
	if (bMotionViewOutput()) {
		MotionViewResOutput(iOutputBlock, "INITIAL CONDITIONS", "MBDyn");
	}
#endif /* USE_ADAMS */

} /* End of DataManager::DataManager() */


/* Distruttore: se richiesto, crea il file di restart; quindi
 * chiama i distruttori degli oggetti propri e libera la memoria di lavoro */

DataManager::~DataManager(void)
{
	/* Se e' richiesto il file di restart, il distruttore del DataManager
	 * crea il file e forza gli oggetti a scrivere il loro contributo nel modo
	 * opportuno
	 */
	if (RestartEvery == ATEND) {
		MakeRestart();
	}

#ifdef USE_ADAMS
	if (bAdamsOutput()) {
		AdamsResOutputFini();
	}
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW
	if (bMotionViewOutput()) {
		MotionViewResOutputFini();
	}
#endif /* USE_MOTIONVIEW */

	if (sSimulationTitle != 0) {
		SAFEDELETEARR(sSimulationTitle);
	}

	if (pOutputMeter) {
		SAFEDELETE(pOutputMeter);
	}

	ElemManagerDestructor();
	NodeManagerDestructor();
	DofManagerDestructor();

#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
	if (loadableElemInitialized) {
		if (lt_dlexit()) {
			std::cerr << "lt_dlexit failed" << std::endl;
			throw ErrGeneric();
		}
	}
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */
} /* End of DataManager::DataManager() */


void
DataManager::OutputOpen(const OutputHandler::OutFiles o)
{
	if (!OutHdl.IsOpen(o)) {
		OutHdl.Open(o);
	}
}


bool
DataManager::bOutput(ResType t) const
{
	return (ResMode & t) ? true : false;
}

void DataManager::MakeRestart(void)
{
	silent_cout("Making restart file ..." << std::endl);
	OutHdl.RestartOpen(saveXSol);
	/* Inizializzazione del file di restart */
	time_t tCurrTime(time(0));
	OutHdl.Restart() << "# Restart file prepared by Mbdyn, "
		<< ctime(&tCurrTime) << std::endl << std::endl;
	/* Dati iniziali */
	OutHdl.Restart() << "begin: data;" << std::endl
		<< "# uncomment this line to use the default integrator" << std::endl
		<< "  integrator: multistep;" << std::endl
		<< "end: data;" << std::endl << std::endl
		<< "# the following block contains data for the multistep integrator"
		<< std::endl;
	pSolver->Restart(OutHdl.Restart(), RestartEvery);

	/* Dati di controllo */
	OutHdl.Restart() << "begin: control data;" << std::endl;

	/* Nodi */
	for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		if (NodeData[iCnt].iNum > 0) {
			OutHdl.Restart() << "  " << psReadControlNodes[iCnt] << ": "
				<< NodeData[iCnt].iNum << ';' << std::endl;
		}
	}

	/* Drivers */
	for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE; iCnt++) {
		if (DriveData[iCnt].iNum > 0) {
			OutHdl.Restart() << "  "
				<< psReadControlDrivers[iCnt] << ": "
				<< DriveData[iCnt].iNum << ';' << std::endl;
		}
	}

	/* Elementi */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		if (!ElemData[iCnt].ElemMap.empty()) {
			if (ElemData[iCnt].bIsUnique()) {
				OutHdl.Restart() << "  " << psReadControlElems[iCnt]
					<< ';' << std::endl;
			} else {
				OutHdl.Restart() << "  " << psReadControlElems[iCnt] << ": "
					<< ElemData[iCnt].ElemMap.size() << ';' << std::endl;
			}
		}
	}

	if (sSimulationTitle != 0) {
		OutHdl.Restart() << "  title: \""
			<< sSimulationTitle << "\";" << std::endl;
	}

#if defined(USE_STRUCT_NODES)
	OutHdl.Restart() << std::endl
		<< "# comment this line if the model is to be modified!" << std::endl
		<< "  skip initial joint assembly;" << std::endl
		<< "# uncomment the following lines to improve the satisfaction of constraints"
		<< std::endl
		<< "  # initial stiffness: " << dInitialPositionStiffness << ", "
		<< dInitialVelocityStiffness << ';' << std::endl
		<< "  # initial tolerance: " << dInitialAssemblyTol << ';' << std::endl
		<< "  # max initial iterations: " << iMaxInitialIterations
		<< ';' << std::endl;
#endif // USE_STRUCT_NODES
	OutHdl.Restart() << "# uncomment this line if restart file is to be made again"
		<< std::endl
		<< "  # make restart file;" << std::endl
		<< "# remember: it will replace the present file if you don't change its name"
		<< std::endl
		<< "  default output: none";
	for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {
		if (NodeData[iCnt].bDefaultOut()) {
			OutHdl.Restart() << ", " << psReadControlNodes[iCnt];
		}
	}
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		if (ElemData[iCnt].bDefaultOut()) {
			OutHdl.Restart() << ", " << psReadControlElems[iCnt];
		}
	}
	OutHdl.Restart() << "; " << std::endl;
	if (saveXSol) {
		OutHdl.Restart() << "  read solution array;" << std::endl;
	}
	OutHdl.Restart() << "end: control data;" << std::endl << std::endl;

	/* Dati dei nodi */
	OutHdl.Restart() << "begin: nodes;" << std::endl;
	for (Node** ppTmpNode = ppNodes;
		ppTmpNode < ppNodes+iTotNodes;
		ppTmpNode++)
	{
		(*ppTmpNode)->Restart(OutHdl.Restart());
	}
	OutHdl.Restart() << "end: nodes;" << std::endl << std::endl;

	/* Dati dei driver */
	if (iTotDrive > 0) {
		OutHdl.Restart() << "begin: drivers;" << std::endl;
		for (Drive** ppTmpDrv = ppDrive;
			ppTmpDrv < ppDrive+iTotDrive;
			ppTmpDrv++)
		{
			(*ppTmpDrv)->Restart(OutHdl.Restart());
#if 0
			OutHdl.Restart()
				<< "  # file driver " << (*ppTmpDrv)->GetLabel()
				<< " is required" << std::endl;
#endif
		}
		OutHdl.Restart() << "end: drivers;" << std::endl << std::endl;
	}

	/* Dati degli elementi */
	OutHdl.Restart() << "begin: elements;" << std::endl;
	for (ElemVecType::const_iterator pTmpEl = Elems.begin();
		pTmpEl != Elems.end();
		pTmpEl++)
	{
		(*pTmpEl)->Restart(OutHdl.Restart());
	}

	Node** pLastNP = NodeData[Node::PARAMETER].ppFirstNode
		+ NodeData[Node::PARAMETER].iNum;
	for (Node** pTmpNP = NodeData[Node::PARAMETER].ppFirstNode; pTmpNP < pLastNP; pTmpNP++) {
		((Elem2Param* )*pTmpNP) ->RestartBind(OutHdl.Restart());
	}

	OutHdl.Restart() << "end: elements;" << std::endl;

	if (saveXSol) {
		OutHdl.RestartXSol().write((char*)(pXCurr->pdGetVec()),
			(pXCurr->iGetSize())*sizeof(double));
		OutHdl.RestartXSol().write((char*)(pXPrimeCurr->pdGetVec()),
			(pXPrimeCurr->iGetSize())*sizeof(double));
		OutHdl.Close(OutputHandler::RESTARTXSOL);
	}

	OutHdl.Close(OutputHandler::RESTART);
}

NamedValue *
DataManager::InsertSym(const char* const s, const Real& v, int redefine)
{
	return MathPar.InsertSym(s, v, redefine);
}

NamedValue *
DataManager::InsertSym(const char* const s, const Int& v, int redefine)
{
	return MathPar.InsertSym(s, v, redefine);
}

/* default orientation description */
void
DataManager::SetOrientationDescription(OrientationDescription od)
{
	this->od = od;
}

OrientationDescription
DataManager::GetOrientationDescription(void) const
{
	return this->od;
}

bool
DataManager::bOutputAccelerations(void) const
{
	return bOutputAccels;
}

