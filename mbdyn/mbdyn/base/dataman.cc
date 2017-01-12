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

/* datamanager */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

extern "C" {
#include <strings.h>
#include <time.h>
}

#include "dataman.h"
#include "friction.h"

#if defined(USE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif // USE_RUNTIME_LOADING && HAVE_LTDL_H

#include "solver.h"
#include "invsolver.h"
#include "constltp.h"
#include "dataman_.h"
/* add-ons for math parser */
#include "dofpgin.h"
#include "privpgin.h"
#include "dummypgin.h"
#include "modelns.h"

/* To allow direct loading of modules */
#include "modules.h"

/* To handle Elem2Param */
#include "j2p.h"

// To handle user-defined elements
#include "userelem.h"

// To handle gusts
#include "gust.h"

// To handle cleanups
#include "cleanup.h"

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
MBPar(HP),
MathPar(HP.GetMathParser()),
pSolver(pS),
DrvHdl(HP.GetMathParser()),
OutHdl(sOutputFileName, 0),
pXCurr(0), pXPrimeCurr(0), 
/* Inverse Dynamics: */
pXPrimePrimeCurr(0),
pLambdaCurr(0),
bInitialJointAssemblyToBeDone(bDefaultInitialJointAssemblyToBeMade),
bSkipInitialJointAssembly(bDefaultSkipInitialJointAssembly),
bOutputFrames(false),
bOutputAccels(false),
bOutputDriveCaller(false),
dInitialPositionStiffness(dDefaultInitialStiffness),
dInitialVelocityStiffness(dDefaultInitialStiffness),
bOmegaRotates(bDefaultOmegaRotates),
dInitialAssemblyTol(dDefaultInitialAssemblyTol),
iMaxInitialIterations(iDefaultMaxInitialIterations),
dEpsilon(1.),
CurrSolver(pS->GetLinearSolver()),
pRBK(0),
bStaticModel(false),
/* auto-detect if running inverse dynamics */
bInverseDynamics(dynamic_cast<InverseSolver *>(pS) != 0),
iIDNodeTotNumDofs(0),
iIDJointTotNumDofs(0),
#if defined(USE_RUNTIME_LOADING)
moduleInitialized(false),
#endif // USE_RUNTIME_LOADING
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
#ifdef MBDYN_FDJAC
pFDJacMeter(0),
#endif // MBDYN_FDJAC
ResMode(RES_TEXT),
#ifdef USE_NETCDF
// NetCDF stuff
bNetCDFsync(false),
bNetCDFnoText(false),
Var_Step(0),
Var_Time(0),
Var_TimeStep(0),
Var_Eig_lStep(0),
Var_Eig_dTime(0),
Var_Eig_dCoef(0),
Var_Eig_dAplus(0),
Var_Eig_dAminus(0),
Var_Eig_dAlpha(0),
Var_Eig_Idx(0),
Var_Eig_dVR(0),
Var_Eig_dVL(0),
#endif // USE_NETCDF
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

#ifdef USE_SOCKET
SocketUsersTimeout(0),
#endif // USE_SOCKET

/* ElemManager */
ElemIter(),
ppDrive(0),
iTotDrive(0),
iMaxWorkNumRowsRes(1), // Allocating a work space size >= 1 will be safe in any case
iMaxWorkNumRowsJac(1),
iMaxWorkNumColsJac(1),
iMaxWorkNumItemsJac(1),
pWorkMatA(0),
pWorkMatB(0),
pWorkMat(0),
pWorkVec(0),

/* NodeManager */
iTotNodes(0),

/* DofManager */
iTotDofOwners(0),
DofOwners(),
iTotDofs(0),
Dofs()
{
	DEBUGCOUTFNAME("DataManager::DataManager");

	mbdyn_cleanup_register(datamanager_cleanup, &ppCleanupData);
	*ppCleanupData = (void *)this;

	OutHdl.SetExceptions(std::ios::badbit); // terminate if disk is full

	/* pseudocostruttori */
	ElemManager();
	NodeManager();
	DofManager();

	InitDriveData();
	InitUDE();
	InitGustData();

	/* registra il plugin per i dofs */
	HP.GetMathParser().RegisterPlugIn("dof", dof_plugin, this);

	/* registra il plugin per i dati privati dei nodi */
	HP.GetMathParser().RegisterPlugIn("node", node_priv_plugin, this);

	/* registra il plugin per i dati privati degli elementi */
	HP.GetMathParser().RegisterPlugIn("element", elem_priv_plugin, this);

	/* registra il namespace del modello */
	
	HP.GetMathParser().RegisterNameSpace(new ModelNameSpace(this));

	/* Setta il tempo al valore iniziale */
	SetTime(dInitialTime, 0., 0, false);

	DEBUGLCOUT(MYDEBUG_INIT, "Global symbol table:"
		<< MathPar.GetSymbolTable() << std::endl);

	/*
	 * Possiede MathParser, con relativa SymbolTable.
	 * Crea ExternKeyTable, MBDynParser,
	 * e legge i dati esterni. Quindi, quando trova
	 * i dati di controllo, chiama la relativa
	 * funzione di lettura (distinta per comodita')
	 */

	/* parole chiave */
	const char* sKeyWords[] = {
		"begin",
		"control" "data",
		"scalar" "function",
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

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (KeyWords(HP.GetWord()) != CONTROLDATA) {
		DEBUGCERR("");
		silent_cerr("<begin: control data;> expected at line "
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (KeyWords(HP.GetWord()) != NODES) {
			DEBUGCERR("");
			silent_cerr("<begin: nodes;> expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
	bool bGotDrivers = false;
	if (iTotDrive > 0) {
		if (CurrDesc != BEGIN) {
			DEBUGCERR("");
			silent_cerr("\"begin\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (KeyWords(HP.GetWord()) != DRIVERS) {
			DEBUGCERR("");
			silent_cerr("\"begin: drivers;\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		bGotDrivers = true;

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

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		switch (KeyWords(HP.GetWord())) {
		case ELEMENTS:
			break;

		case DRIVERS:
			if (!bGotDrivers) {
				silent_cerr("got unexpected \"begin: drivers;\" "
					"at line " << HP.GetLineData() << std::endl
					<< "(hint: define \"file drivers\" in \"control data\" block)"
					<< std::endl);
			}
			// fallthru

		default:
			DEBUGCERR("");
			silent_cerr("\"begin: elements;\" expected at line "
				<< HP.GetLineData() << std::endl);

			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	if (bOutputFrames) {
		OutHdl.Open(OutputHandler::REFERENCEFRAMES);
		HP.OutputFrames(OutHdl.ReferenceFrames());
	}

	/* fine lettura elementi */

	// if output is defined and at least one node wants to output,
	// open file
	for (int i = 0; i < Node::LASTNODETYPE; i++) {
		if (NodeData[i].OutFile != OutputHandler::UNKNOWN)
		{
			for (NodeContainerType::const_iterator n = NodeData[i].NodeContainer.begin();
				n != NodeData[i].NodeContainer.end(); ++n)
			{
				if (n->second->bToBeOutput()) {
					OutHdl.Open(NodeData[i].OutFile);
					break;
				}
			}
		}
	}

	// if output is defined and at least one element wants to output,
	// open file
	for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
		if (ElemData[i].OutFile != OutputHandler::UNKNOWN)
		{
			for (ElemContainerType::const_iterator e = ElemData[i].ElemContainer.begin();
				e != ElemData[i].ElemContainer.end(); ++e)
			{
				if (e->second->bToBeOutput()) {
					OutHdl.Open(ElemData[i].OutFile);
					break;
				}
			}
		}
	}

	// open drive output & trave files only if needed
	const MBDynParser::DCType& DC = MBPar.GetDriveCallerContainer();
	for (MBDynParser::DCType::const_iterator i = DC.begin(); i != DC.end(); ++i) {
		if (i->second->fToBeTraced()) {
			OutHdl.Open(OutputHandler::TRACES);
			break;
		}
	}

	for (MBDynParser::DCType::const_iterator i = DC.begin(); i != DC.end(); ++i) {
		if (i->second->bToBeOutput()) {
			OutHdl.Open(OutputHandler::DRIVECALLERS);
			break;
		}
	}

	/* Verifica dei dati di controllo */
#ifdef DEBUG
	if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
		/* mostra in modo succinto il numero di DofOwner per tipo */
		for (int i = 0; i < DofOwner::LASTDOFTYPE; i++) {
			std::cout << "DofType " << i << " "
				"(" << psDofOwnerNames[i] << "), "
				"n. of owners: " << DofData[i].iNum
				<< std::endl;
		}

		/* mostra in modo succinto il numero di nodi per tipo */
		for (int i = 0; i < Node::LASTNODETYPE; i++) {
			std::cout << "NodeType " << i << " "
				"(" << psNodeNames[i] << "), "
				"n. of nodes: " << NodeData[i].NodeContainer.size()
				<< std::endl;
		}

		/* mostra in modo succinto il numero di elementi per tipo */
		for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
			std::cout << "Element Type " << i << " "
				"(" << psElemNames[i] << "), "
				"n. of elems: " << ElemData[i].ElemContainer.size()
				<< std::endl;
		}

		/* mostra in modo succinto il numero di drivers per tipo */
		for (int i = 0; i < Drive::LASTDRIVETYPE; i++) {
			std::cout << "DriveType " << i << " "
				"(" << psDriveNames[i] << "), "
				"n. of drivers: " << DriveData[i].iNum
				<< std::endl;
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

#ifdef USE_SOCKET
	/* waits for all pending sockets to connect */
	WaitSocketUsers();
#endif // USE_SOCKET

	/* Qui intercetto la struttura dei Dof prima che venga costruita e modifico
	 * in modo che sia adatta all'assemblaggio dei vincoli; quindi la resetto
	 * e la lascio generare correttamente.
	 * L'assemblaggio iniziale viene gestito dal DataManager perche'
	 * concettualmente la consistenza dei dati deve essere assicurata da lui,
	 * che conosce gli aspetti fisici del problema
	 */

	if (!bInitialJointAssemblyToBeDone) {
		for (int i = 0; i < Elem::LASTELEMTYPE; ++i) {
			if (ElemData[i].bToBeUsedInAssembly() && !ElemData[i].ElemContainer.empty()) {
				bInitialJointAssemblyToBeDone = true;
				break;
			}
		}
	}

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

	/* Costruzione dei dati dei Dof definitivi da usare nella simulazione */

	/* Aggiornamento DofOwner degli elementi (e dei nodi) */
	if (bInverseDynamics) {
		IDDofOwnerSet();
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
					<< DofOwners[k++].iNumDofs << std::endl;
			}
		}
	}
#endif /* DEBUG */

	/* a questo punto i dof sono completamente scritti in quanto a numero.
	 * Rimane da sistemare il vettore con gli indici "interni" ed il tipo */

	/* Costruzione array DofOwner e Dof */
	if(bInverseDynamics)	{
		IDDofInit();
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
				<< Dofs[iCnt].Order << std::endl;
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
	*ppCleanupData = 0;

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
		sSimulationTitle = 0;
	}

	if (pOutputMeter) {
		SAFEDELETE(pOutputMeter);
		pOutputMeter = 0;
	}

#ifdef MBDYN_FDJAC
	if (pFDJacMeter) {
		SAFEDELETE(pFDJacMeter);
		pFDJacMeter = 0;
	}
#endif // MBDYN_FDJAC

	if (pRBK) {
		SAFEDELETE(pRBK);
		pRBK = 0;
	}

	ElemManagerDestructor();
	NodeManagerDestructor();
	DofManagerDestructor();

	DestroyDriveData();
	DestroyUDE();
	DestroyGustData();

#if defined(USE_RUNTIME_LOADING)
	if (moduleInitialized) {
		module_finalize();
	}
#endif // USE_RUNTIME_LOADING
} /* End of DataManager::DataManager() */

int
DataManager::Cleanup(void)
{
#ifdef USE_SOCKET
	DeleteSocketUsers();
#endif // USE_SOCKET
	return 0;
}

void
DataManager::OutputOpen(const OutputHandler::OutFiles o)
{
	if (!OutHdl.IsOpen(o)) {
		OutHdl.Open(o);
	}
}

void
DataManager::OutputEigOpen(const unsigned uCurrEigSol)
{
	if(!OutHdl.IsOpen(OutputHandler::EIGENANALYSIS)) {
		OutHdl.Open(OutputHandler::EIGENANALYSIS, uCurrEigSol);
		ASSERT(OutHdl.IsOpen(OutputHandler::EIGENANALYSIS));
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
		if (!NodeData[iCnt].NodeContainer.empty()) {
			OutHdl.Restart() << "  " << psReadControlNodes[iCnt] << ": "
				<< NodeData[iCnt].NodeContainer.size() << ';' << std::endl;
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
		if (!ElemData[iCnt].ElemContainer.empty()) {
			if (ElemData[iCnt].bIsUnique()) {
				OutHdl.Restart() << "  " << psReadControlElems[iCnt]
					<< ';' << std::endl;
			} else {
				OutHdl.Restart() << "  " << psReadControlElems[iCnt] << ": "
					<< ElemData[iCnt].ElemContainer.size() << ';' << std::endl;
			}
		}
	}

	if (sSimulationTitle != 0) {
		OutHdl.Restart() << "  title: \""
			<< sSimulationTitle << "\";" << std::endl;
	}

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
	for (NodeVecType::const_iterator n = Nodes.begin(); n != Nodes.end(); ++n) {
		(*n)->Restart(OutHdl.Restart());
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
	for (ElemVecType::const_iterator e = Elems.begin();
		e != Elems.end(); ++e)
	{
		(*e)->Restart(OutHdl.Restart());
	}

	for (NodeContainerType::const_iterator n = NodeData[Node::PARAMETER].NodeContainer.begin();
		n != NodeData[Node::PARAMETER].NodeContainer.end(); ++n)
	{
		dynamic_cast<const Elem2Param *>(n->second)->RestartBind(OutHdl.Restart());
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

/* default output */
void
DataManager::SetOutput(Elem::Type t, unsigned flags, OrientationDescription od)
{
	ElemData[t].uOutputFlags = flags;
	ElemData[t].od = od;
}

void
DataManager::GetOutput(Elem::Type t, unsigned& flags, OrientationDescription& od) const
{
	flags = ElemData[t].uOutputFlags;
	od = ElemData[t].od;
}

bool
DataManager::bOutputAccelerations(void) const
{
	return bOutputAccels;
}

bool
DataManager::bOutputDriveCallers(void) const
{
	return bOutputDriveCaller;
}

DataManager::NodeContainerType::const_iterator
DataManager::begin(Node::Type t) const
{
	return NodeData[t].NodeContainer.begin();
}

DataManager::NodeContainerType::const_iterator
DataManager::end(Node::Type t) const
{
	return NodeData[t].NodeContainer.end();
}

DataManager::ElemContainerType::const_iterator
DataManager::begin(Elem::Type t) const
{
	return ElemData[t].ElemContainer.begin();
}

DataManager::ElemContainerType::const_iterator
DataManager::end(Elem::Type t) const
{
	return ElemData[t].ElemContainer.end();
}

extern "C" int
datamanager_cleanup(void *data)
{
	int rc = 0;

	if (data != 0) {
		DataManager *pDM = (DataManager *)data;

		silent_cerr("DataManager cleanup" << std::endl);

		rc = pDM->Cleanup();
	}

	return rc;
}

bool
DataManager::PushCurrData(const std::string& name, const TypedValue& value)
{
	ModelNameSpace *pMNS = dynamic_cast<ModelNameSpace *>(MathPar.GetNameSpace("model"));
	if (pMNS == 0) {
		return false;
	}

	return pMNS->PushCurrData(name, value);
}

bool
DataManager::PopCurrData(const std::string& name)
{
	ModelNameSpace *pMNS = dynamic_cast<ModelNameSpace *>(MathPar.GetNameSpace("model"));
	if (pMNS == 0) {
		return false;
	}

	return pMNS->PopCurrData(name);
}

