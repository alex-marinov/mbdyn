/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* gestore dei dati */

#ifndef DATAMAN_H
#define DATAMAN_H

#include <iostream>
#include <list>
#include <map>
#include <string>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "demangle.h"

#include "mbpar.h"
#include "constltp.h"
#include "shape.h"

/* da elman.h */
#include "solman.h"
#include "submat.h"
#include "veciter.h"

#include "elem.h"      /* Classe di base di tutti gli elementi */
#include "driven.h"
#include "output.h"

#include "drive.h"     /* Drive vari */
#include "tpldrive.h"  /* Drive vari */

/* da nodeman.h */
#include "node.h"
#include "strnode.h"
#include "elecnode.h"

#include "solverdiagnostics.h"
#include "nonlin.h"
#include "linsol.h"
#include "converged.h"
#include "invdyn.h"

#ifdef USE_SOCKET
#include "usesock.h"
#endif // USE_SOCKET

struct LoadableCalls;
class Solver;

#include "datamanforward.h"

/* DataManager - begin */

class DataManager : public SolutionDataManager, public SolverDiagnostics,
	public DataManagerErrors {

protected:
	void **ppCleanupData;

#ifdef USE_MULTITHREAD
	/* from input file, or auto-detected */
	unsigned int nThreads;
#endif /* USE_MULTITHREAD */

	/* Handler vari */
	MBDynParser& MBPar;	/* Received from Solver */
	MathParser& MathPar;	/* Received from Solver */
	Solver* pSolver;

	/* loadable elements */
	std::map<std::string, const LoadableCalls *> MapOfLoadableElemHandlers;

	DriveHandler DrvHdl;
	mutable OutputHandler OutHdl;

	/* Puntatori ai vettori soluzione durante il passo */
	mutable VectorHandler* pXCurr;
	mutable VectorHandler* pXPrimeCurr;
	
	/* used by inverse dynamics: */
	mutable VectorHandler* pXPrimePrimeCurr;
	mutable VectorHandler* pLambdaCurr;

	/* Parametri usati durante l'assemblaggio iniziale */
	bool bInitialJointAssemblyToBeDone;
	bool bSkipInitialJointAssembly;
	bool bOutputFrames;
	bool bOutputAccels;
	bool bOutputDriveCaller;
	doublereal dInitialPositionStiffness;
	doublereal dInitialVelocityStiffness;
	bool bOmegaRotates;
	doublereal dInitialAssemblyTol;
	integer iMaxInitialIterations;
	doublereal dEpsilon;
	LinSol CurrSolver;

	RigidBodyKinematics *pRBK;
	bool bStaticModel;

	// inverse dynamics
	bool bInverseDynamics;
	int iIDNodeTotNumDofs;
	int iIDJointTotNumDofs;

#ifdef USE_RUNTIME_LOADING
	bool moduleInitialized;
#endif // USE_RUNTIME_LOADING

	enum PrintFlags {
		PRINT_NONE		= 0x00U,

		PRINT_DOF_STATS		= 0x01U,

		PRINT_DOF_DESCRIPTION	= 0x02U,
		PRINT_EQ_DESCRIPTION	= 0x04U,
		PRINT_DESCRIPTION	= (PRINT_DOF_DESCRIPTION|PRINT_EQ_DESCRIPTION),

		PRINT_NODE_CONNECTION	= 0x10U,
		PRINT_EL_CONNECTION	= 0x20U,
		PRINT_CONNECTION	= (PRINT_NODE_CONNECTION|PRINT_EL_CONNECTION),

		PRINT_TO_FILE 		= 0x1000U
	};
	unsigned uPrintFlags;
	/* Parametri vari */
	char* sSimulationTitle;

public:
	enum eRestart { NEVER, ATEND, ITERATIONS, TIME, TIMES };
protected:
	/* soft-restart stuff */
	eRestart RestartEvery;
	integer iRestartIterations;
	doublereal dRestartTime;

	doublereal *pdRestartTimes;
	integer iNumRestartTimes;
	mutable integer iCurrRestartTime;

	mutable integer iCurrRestartIter;
	mutable doublereal dLastRestartTime;

	bool saveXSol;
	char * solArrFileName;

	/* raw output stuff */
	DriveCaller *pOutputMeter;
	mutable integer iOutputCount;

#ifdef MBDYN_FDJAC
protected:
	DriveCaller *pFDJacMeter;

public:
	bool bFDJac(void) const;
#endif // MBDYN_FDJAC

	/* specialized output stuff */
public:
	enum ResType {
		RES_NONE	= 0x00,
		RES_TEXT	= 0x01,
		RES_NETCDF	= 0x02,
		RES_ADAMS	= 0x04,
		RES_MOTIONVIEW	= 0x08
	};

	bool bOutput(ResType t) const;

protected:
	int ResMode;

#ifdef USE_NETCDF
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* NetCDF stuff */
	bool bNetCDFsync;
	bool bNetCDFnoText;
	NcVar *Var_Step;
	NcVar *Var_Time;
	NcVar *Var_TimeStep;
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#endif /* USE_NETCDF */
	OrientationDescription od;

#if defined(USE_ADAMS) || defined(USE_MOTIONVIEW)
	mutable integer iOutputBlock;
#endif /* defined(USE_ADAMS) || defined(USE_MOTIONVIEW) */

#ifdef USE_ADAMS
	char *sAdamsModelName;
	bool bAdamsVelocity, bAdamsAcceleration;
	unsigned int iAdamsOutputNodes;
	unsigned int iAdamsOutputParts;
	std::streambuf::pos_type adamsNoab;
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW
	/* nothing yet */
#endif /* USE_MOTIONVIEW */

protected:
	/* chiamate dal costruttore per leggere i relativi articoli */
	void ReadControl(MBDynParser& HP, const char* sInputFileName);
	void ReadNodes(MBDynParser& HP);
	void ReadDrivers(MBDynParser& HP);
	void ReadElems(MBDynParser& HP);

	template <class T> T* Cast(Elem *pEl, bool bActive = false);

public:
	template <class T>
	flag fReadOutput(MBDynParser& HP, const T& t) const;
	doublereal dReadScale(MBDynParser& HP, enum DofOwner::Type t) const;

	bool bOutputAccelerations(void) const;
	bool bOutputDriveCallers(void) const;
	const doublereal& dGetInitialPositionStiffness(void) const;
	const doublereal& dGetInitialVelocityStiffness(void) const;
	bool bDoesOmegaRotate(void) const;

	void IncElemCount(Elem::Type type);

	/* additional CPU time, if any */
	virtual clock_t GetCPUTime(void) const {
		return 0;
	};

protected:
	/* chiamate a funzioni di inizializzazione */
	void InitialJointAssembly(void);

	void DofOwnerSet(void);
	void DofOwnerInit(void);
	
	/* inverse dynamics: */
	void IDDofOwnerSet(void);

public:
	/* costruttore - legge i dati e costruisce le relative strutture */
	DataManager(MBDynParser& HP,
		unsigned OF,
		Solver* pS,
		doublereal dInitialTime,
		const char* sOutputFileName,
		const char* sInputFileName,
		bool bAbortAfterInput);

	/* distruttore */
	virtual ~DataManager(void);

	int Cleanup(void);

	/* helpers */
	int ReadScalarAlgebraicNode(MBDynParser& HP, unsigned int uLabel,
			Node::Type type, doublereal& dX) const;
	int ReadScalarDifferentialNode(MBDynParser& HP, unsigned int uLabel,
			Node::Type type, doublereal& dX, doublereal& dXP) const;
	Node* ReadNode(MBDynParser& HP, Node::Type type) const;
	Elem* ReadElem(MBDynParser& HP, Elem::Type type) const;

	template <class Tbase, Node::Type type>
	Tbase *ReadNode(MBDynParser& HP) const;
	template <class Tder, class Tbase, Node::Type type>
	Tder *ReadNode(MBDynParser& HP) const;
	template <class Tbase, Elem::Type type>
	Tbase *ReadElem(MBDynParser& HP) const;
	template <class Tder, class Tbase, Elem::Type type>
	Tder *ReadElem(MBDynParser& HP) const;

	/* Funzioni usate dal metodo di integrazione */

	/* Setta il valore della variabile tempo nella tabella dei simboli
	 * del DataManager e nel DriveHandler */
	void SetTime(const doublereal& dTime, const doublereal& dTimeStep = -1.,
		const integer& iStep = -1, bool bServePending = true);
	doublereal dGetTime(void) const;

	NamedValue *InsertSym(const char* const s, const Real& v,
			int redefine = 0);
	NamedValue *InsertSym(const char* const s, const Int& v,
			int redefine = 0);

	/* Collega il DataManager ed il DriveHandler ai vettori soluzione */
	void LinkToSolution(VectorHandler& XCurr,
		VectorHandler& XPrimeCurr);

	/* inverse dynamics: */
	void LinkToSolution(VectorHandler& XCurr,
		VectorHandler& XPrimeCurr,
		VectorHandler& XPrimePrimeCurr,
		VectorHandler& LambdaCurr);
	
	/* Restituisce l'ostream al file di output,
	 * usato dai vari metodi per scrivere il log del calcolo */
	std::ostream& GetOutFile(void) const { return OutHdl.Output(); };
	std::ostream& GetLogFile(void) const { return OutHdl.Log(); };

	/* required for binary NetCDF output access */
	const OutputHandler* pGetOutHdl(void) const { return &OutHdl; };

	/* default orientation description */
	void SetOrientationDescription(OrientationDescription);
	OrientationDescription GetOrientationDescription(void) const;

	/* default beam output */
	void SetOutput(Elem::Type t, unsigned, OrientationDescription);
	void GetOutput(Elem::Type t, unsigned&, OrientationDescription&) const;

	/* Restituisce il DriveHandler */
	const DriveHandler* pGetDrvHdl(void) const { return &DrvHdl; };
	MathParser& GetMathParser(void) const { return MathPar; };
	MBDynParser& GetMBDynParser(void) const { return MBPar; };
	const Solver *GetSolver(void) const { return pSolver; };

	bool PushCurrData(const std::string& name, const TypedValue& value);
	bool PopCurrData(const std::string& name);

	/* Assembla lo jacobiano */
	virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef);

	/* Assembla le matrici per gli autovalori */
	virtual void AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl);

	/* Assembla il residuo */
	virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef)
		throw(ChangedEquationStructure);

	// inverse dynamics
	/* Constraints residual, switch iOrder*/
	virtual void AssConstrRes(VectorHandler& ResHdl,
		InverseDynamics::Order iOrder) 
		throw(ChangedEquationStructure) ;

	/* Elem residual, equilibrium with no constraints */
	virtual void AssRes(VectorHandler& ResHdl)
		throw(ChangedEquationStructure);
	
	/* Constraint Jacobian matrix*/	
	virtual void AssConstrJac(MatrixHandler& JacHdl);
	// end of inverse dynamics

protected:
	/* specialized functions, called by above general helpers */
	virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef,
			VecIter<Elem *> &Iter,
			VariableSubMatrixHandler& WorkMat);
	virtual void AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl,
			VecIter<Elem *> &Iter,
			VariableSubMatrixHandler& WorkMatA,
			VariableSubMatrixHandler& WorkMatB);
	virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef,
			VecIter<Elem *> &Iter,
			SubVectorHandler& WorkVec)
		throw(ChangedEquationStructure);

	// inverse dynamics
	void AssConstrJac(MatrixHandler& JacHdl,
		VecIter<Elem *> &Iter,
		VariableSubMatrixHandler& WorkMat);

	void AssConstrRes(VectorHandler& ResHdl,
		VecIter<Elem *> &Iter,
		SubVectorHandler& WorkVec,
		InverseDynamics::Order iOrder)
	throw(ChangedEquationStructure);
	
	void AssRes(VectorHandler& ResHdl,
		VecIter<Elem *> &Iter,
		SubVectorHandler& WorkVec)
	throw(ChangedEquationStructure);
	// end of inverse dynamics

protected:
	typedef std::vector<Converged::State> Converged_t;
	mutable Converged_t m_IsConverged;

public:
	// returns an idx to a newly created slot for convergence
	unsigned ConvergedRegister(void);
	// set the value of a slot; elements that register
	// using ConvergedRegister() should set to false
	// at first iteration, and set to true when convergence
	// is allowed
	void ConvergedSet(unsigned idx, Converged::State s);
	// returns true only if all slots are true
	bool IsConverged(void) const;
	bool EndOfSimulation(void) const;

public:
	virtual void OutputPrepare(void);

	/* stampa i risultati */
	virtual bool
	Output(long lStep, const doublereal& dTime,
		const doublereal& dTimeStep, bool force = false) const;
	virtual void
	Output(const VectorHandler& X, const VectorHandler& XP) const;

	void OutputOpen(const OutputHandler::OutFiles out);

#ifdef USE_ADAMS
	/* MSC's ADAMS/View .res output */
	bool bAdamsOutput(void) const;
	void AdamsResOutputInit(void);
	void AdamsResOutput(integer iBlock, const char *type,
			const char *id) const;
	void AdamsResOutputFini(void) const;
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW
	/* Altair's Motion View output */
	bool bMotionViewOutput(void) const;
	void MotionViewResOutputInit(const char *sOutputFileName);
	void MotionViewResOutput(integer iBlock, const char *type,
			const char *id) const;
	void MotionViewResOutputFini(void) const;
#endif /* USE_MOTIONVIEW */

	/* Prepara la soluzione con i valori iniziali */
	void SetValue(VectorHandler& X, VectorHandler& XP);

	/* Funzioni di aggiornamento dati durante la simulazione */
	virtual void MakeRestart(void);
	virtual void DerivativesUpdate(void) const;
	virtual void BeforePredict(VectorHandler& X, VectorHandler& XP,
			VectorHandler& XPrev, VectorHandler& XPPrev) const;
	virtual void AfterPredict(void) const;
	virtual void Update(void) const;
	virtual void AfterConvergence(void) const;
	
	/* Inverse Dynamics: */
	virtual void Update(InverseDynamics::Order iOrder) const;
	virtual void IDAfterConvergence(void) const;
	virtual void IDSetTest(NonlinearSolverTestRange *pResTest, NonlinearSolverTestRange *pSolTest);

	void bSetStaticModel(bool b) {
		bStaticModel = b;
	};
	bool bIsStaticModel(void) const {
		return bStaticModel;
	};
	const RigidBodyKinematics *pGetRBK(void) const {
		return pRBK;
	};

	/* Inverse Dynamics: */
	void bSetInverseDynamics(bool b) {
		bInverseDynamics = b;
	};
	bool bIsInverseDynamics(void) const {
		return bInverseDynamics;
	};

	/* socket select stuff */
#ifdef USE_SOCKET
protected:
	std::map<int, UseSocket *> SocketUsers;
	time_t SocketUsersTimeout;

	void WaitSocketUsers(void);
	void DeleteSocketUsers(void);

public:
	void RegisterSocketUser(UseSocket *pUS);
#endif // USE_SOCKET

	/* da ElemManager */
	friend class InitialAssemblyIterator;

	enum ModuleInsertMode {
		MIM_FAIL,
		MIM_REPLACE,
		MIM_IGNORE
	};

	/* loadable elements */
	const LoadableCalls *GetLoadableElemModule(std::string) const;
	void SetLoadableElemModule(std::string, const LoadableCalls *,
			ModuleInsertMode = MIM_FAIL);

public:
	/* FIXME: will be eliminated */
	enum DerivationTable {
		ELEM			= 0x0U,  // pleonastico
		DOFOWNER		= 0x1U,
		GRAVITYOWNER		= 0x2U,
		AIRPROPOWNER		= 0x4U,
		INITIALASSEMBLY		= 0x8U
	};
	/* end of FIXME: will be eliminated */

	enum DataFlags {
		NONE			= 0x00U,
		ISUNIQUE		= 0x01U,
		TOBEUSEDINASSEMBLY	= 0x02U,
		GENERATESINERTIAFORCES	= 0x04U,
		USESAIRPROPERTIES	= 0x08U,
		DEFAULTOUT		= 0x10U
	};

	/* element read functional object prototype */
	struct ElemRead {
		virtual ~ElemRead( void ) { NO_OP; };
		virtual Elem *
		Read(const DataManager *pDM, MBDynParser& HP,
			unsigned int uLabel, int CurrType) const = 0;
	};

	typedef std::map<std::string, DataManager::ElemRead *, ltstrcase> ElemReadType;
	typedef std::pair<unsigned, Elem*> KeyElemPair;
	typedef std::list<KeyElemPair> ElemContainerType;
	typedef std::map<unsigned, ElemContainerType::iterator> ElemMapToListType;

protected:

	/* struttura dei dati fondamentali degli elementi */
	struct ElemDataStructure {
		unsigned int iExpectedNum;	// numero di elementi del tipo
		const char *Desc;
		const char *ShortDesc;

		DofOwner::Type DofOwnerType;	// Tipo di DofOwner
		unsigned int iDerivation;	// Tabella delle derivazioni

		OutputHandler::OutFiles OutFile;	// Tipo di file in output

		unsigned uFlags;		// flags

		unsigned uOutputFlags;
		OrientationDescription od;


		/* helpers */
		void IsUnique(bool b) { if (b) { uFlags |= ISUNIQUE; } else { uFlags &= ~ISUNIQUE; } };
		void ToBeUsedInAssembly(bool b) { if (b) { uFlags |= TOBEUSEDINASSEMBLY; } else { uFlags &= ~TOBEUSEDINASSEMBLY; } };
		void GeneratesInertiaForces(bool b) { if (b) { uFlags |= GENERATESINERTIAFORCES; } else { uFlags &= ~GENERATESINERTIAFORCES; } };
		void UsesAirProperties(bool b) { if (b) { uFlags |= USESAIRPROPERTIES; } else { uFlags &= ~USESAIRPROPERTIES; } };
		void DefaultOut(bool b) { if (b) { uFlags |= DEFAULTOUT; } else { uFlags &= ~DEFAULTOUT; } };

		bool bIsUnique(void) const { return (uFlags & ISUNIQUE) == ISUNIQUE; };
		bool bToBeUsedInAssembly(void) const { return (uFlags & TOBEUSEDINASSEMBLY) == TOBEUSEDINASSEMBLY; };
		bool bGeneratesInertiaForces(void) const { return (uFlags & GENERATESINERTIAFORCES) == GENERATESINERTIAFORCES; };
		bool bUsesAirProperties(void) const { return (uFlags & USESAIRPROPERTIES) == USESAIRPROPERTIES; };
		bool bDefaultOut(void) const { return (uFlags & DEFAULTOUT) == DEFAULTOUT; };

		/* element read map */
		ElemReadType ElemRead;
		ElemContainerType ElemContainer;
		ElemMapToListType ElemMapToList;
	} ElemData[Elem::LASTELEMTYPE];

	Elem ** InsertElem(ElemDataStructure& eldata, unsigned int uLabel, Elem * pE) {
		eldata.ElemContainer.push_back(ElemContainerType::value_type(uLabel, pE));
		eldata.ElemMapToList[uLabel] = --eldata.ElemContainer.end();
		return &eldata.ElemContainer.back().second;
	};

#if 0
	/* element type map; will replace ElemData */
	typedef std::map<std::string, ElemDataStructure *, ltstrcase> ElemDataMapType;
	ElemDataMapType ElemDataMap;
#endif

	/* array of elements */
	typedef std::vector<Elem *> ElemVecType;
	ElemVecType Elems;

	/* NOTE: will be removed? */
	mutable VecIter<Elem *> ElemIter;
	/* end of NOTE: will be removed? */

	/* struttura dei drivers */
	struct {
		Drive** ppFirstDrive;
		unsigned int iNum;
	} DriveData[Drive::LASTDRIVETYPE];

	Drive** ppDrive;         /* puntatore ai drivers */
	unsigned int iTotDrive;  /* numero totale dei drivers */

	/* dati di lavoro */
	integer iMaxWorkNumRowsRes;	// the maximum number of rows of the residual subvector
	integer iMaxWorkNumRowsJac; // the maximum number of rows in a full Jacobian submatrix
	integer iMaxWorkNumColsJac;	// the maximum number of columns in a full Jacobian submatrix
	integer iMaxWorkNumItemsJac; // the maximum number of items in a sparse Jacobian submatrix

	VariableSubMatrixHandler *pWorkMatA;  /* SubMatrix di lavoro */
	VariableSubMatrixHandler *pWorkMatB;
	VariableSubMatrixHandler *pWorkMat;
	MySubVectorHandler *pWorkVec;

	/* ricerca elementi */
	Elem* pFindElem(Elem::Type Typ, unsigned int uElem,
			unsigned int iDeriv) const;
	Elem* pChooseElem(Elem* p, unsigned int iDeriv) const;

	Elem** ppFindElem(Elem::Type Typ, unsigned int uElem) const;

	flag fGetDefaultOutputFlag(const Elem::Type& t) const;

	Elem** ReadOneElem(MBDynParser& HP,
			unsigned int uLabel,
			const std::string& sName,
			int CurrType);

public:
	/* ricerca drives */
	Drive* pFindDrive(Drive::Type Typ, unsigned int uL) const;

	/* ricerca elementi*/
	Elem* pFindElem(Elem::Type Typ, unsigned int uElem = unsigned(-1)) const;

	template <class Tbase, Elem::Type type>
	Tbase *pFindElem(unsigned int uElem = unsigned(-1)) const;
	template <class Tder, class Tbase, Elem::Type type>
	Tder *pFindElem(unsigned int uElem = unsigned(-1)) const;

	const DataManager::ElemDataStructure& GetElemDataStructure(Elem::Type Typ) const { return ElemData[Typ]; };

	/* pseudocostruttore */
	void ElemManager(void);
	void ElemManagerDestructor(void);

	/* Funzioni di inizializzazione */

	/* Inizializzatore */
	void ElemDataInit(void);

	/* Preassemblaggio */
	void ElemAssInit(void);

	/* Funzioni di routine */

	/* Scrive i risultati */
	void ElemOutputPrepare(OutputHandler& OH);
	void ElemOutput(OutputHandler& OH) const;
	void ElemOutput(OutputHandler& OH,
			const VectorHandler& X, const VectorHandler& XP) const;
	void DriveOutput(OutputHandler& OH) const;
	void DriveTrace(OutputHandler& OH) const;
	DataManager::ElemContainerType::const_iterator begin(Elem::Type t) const;
	DataManager::ElemContainerType::const_iterator end(Elem::Type t) const;

	/* da NodeManager */
public:
	/* element read functional object prototype */
	struct NodeRead {
		virtual ~NodeRead(void) { NO_OP; };
		virtual Elem *
		Read(const DataManager *pDM, MBDynParser& HP,
			unsigned int uLabel, int CurrType) const = 0;
	};

	typedef std::map<std::string, DataManager::NodeRead *, ltstrcase> NodeReadType;
	typedef std::pair<unsigned, Node*> KeyNodePair;
	typedef std::list<KeyNodePair> NodeContainerType;
	typedef std::map<unsigned, NodeContainerType::iterator> NodeMapToListType;

protected:

	/* struttura dei dati dei nodi. Per ogni tipo:
	 * puntatore al puntatore al primo dato, numero degli item per tipo */
	struct NodeDataStructure {
		unsigned int iExpectedNum;	// numero di nodi del tipo
		unsigned uFlags;		// flags
		const char *Desc;
		const char *ShortDesc;

		/* helpers */
		void DefaultOut(bool b) { if (b) { uFlags |= DEFAULTOUT; } else { uFlags &= ~DEFAULTOUT; } };

		bool bDefaultOut(void) const { return (uFlags & DEFAULTOUT) == DEFAULTOUT; };

		OutputHandler::OutFiles OutFile; /* Tipo di file in output */

		/* element read map */
		NodeReadType NodeRead;
		NodeContainerType NodeContainer;
		NodeMapToListType NodeMapToList;
	} NodeData[Node::LASTNODETYPE];
	
	Node ** InsertNode(NodeDataStructure& nodedata, unsigned int uLabel, Node * pN) {
		nodedata.NodeContainer.push_back(NodeContainerType::value_type(uLabel, pN));
		nodedata.NodeMapToList[uLabel] = --nodedata.NodeContainer.end();
		return &nodedata.NodeContainer.back().second;
	};

	/* array of nodes */
	typedef std::vector<Node *> NodeVecType;
	NodeVecType Nodes;

	/* dati dei nodi: numero totale e puntatore all'array dei dati
	 * (ogni nodo ha il suo formato caratteristico, comunque derivato
	 * dalla classe Node) */
	unsigned int iTotNodes;

public:
	Node** ppFindNode(Node::Type Typ, unsigned int uNode) const;
	/* ricerca di nodi */
	Node* pFindNode(Node::Type Typ, unsigned int uNode) const;

	template <class Tbase, Node::Type type>
	Tbase *pFindNode(unsigned int uNode) const;
	template <class Tder, class Tbase, Node::Type type>
	Tder *pFindNode(unsigned int uNode) const;

protected:
	flag fGetDefaultOutputFlag(const Node::Type& t) const;

public:
	/* Pseudocostruttore */
	void NodeManager(void);
	void NodeManagerDestructor(void);

	/* inizializza le matrici ed alloca memoria */
	void NodeDataInit(void);

	DataManager::NodeContainerType::const_iterator begin(Node::Type t) const;
	DataManager::NodeContainerType::const_iterator end(Node::Type t) const;

	/* scrive i dati dei nodi */
	void NodeOutputPrepare(OutputHandler& OH);
	void NodeOutput(OutputHandler& OH) const;
	void NodeOutput(OutputHandler& OH,
			const VectorHandler& X, const VectorHandler& XP) const;

	/* da DofManager */
protected:

	/* struttura dei dati generali dei dof: numero totale per tipo,
	 * dimensione caratteristica (se esiste), puntatore al primo del tipo */
	struct {
		DofOwner* pFirstDofOwner;     /* punt. al 1o di ogni tipo */
		integer iNum;                 /* n. DofOwners per ogni tipo */
		integer iSize;                /* n. Dof (se fisso, es. nodi) */
		doublereal dDefScale;
	} DofData[DofOwner::LASTDOFTYPE];

	/* struttura dei dati dei dof di ogni ente possessore:
	 * totale dei possessori; per ognuno: indice del primo dof,
	 * numero di dof posseduti */
	integer iTotDofOwners;           /* numero totale di DofOwners */
	DofOwner* pDofOwners;            /* puntatore all'array dei DofOwner */

public:
	/* struttura dei dati di ogni singolo dof: totale dei dof; per ognuno:
	 * indice e tipo (algebrico o differenziale) */
	typedef std::vector<Dof> DofVecType;
	typedef DofVecType::const_iterator DofIterator_const;
	typedef DofVecType::iterator DofIterator;

	/* Restituisce il puntatore alla struttura dei dof */
	const DofVecType& GetDofs(void) const { return Dofs; };

	/* Restituisce il numero di dof per la costruzione delle matrici ecc. */
	integer iGetNumDofs(void) const { return iTotDofs; };

protected:
	DofVecType Dofs;
	integer iTotDofs;                /* numero totale di Dof */

	DofOwner DummyDofOwner; /* Per quelli che non hanno dof */

	doublereal dGetDefaultScale(DofOwner::Type t) const;

public:
	/* pseudocostruttore */
	void DofManager(void);
	void DofManagerDestructor(void);

	/* funzioni di inizializzazione */
	void DofDataInit(void);
	void DofInit(void);

	/* Inverse Dynamics: */
	void IDDofInit(void);

	int iIDGetNodeTotNumDofs(void) const;
	int iIDGetJointTotNumDofs(void) const;
	int iIDGetTotNumDofs(void) const;

	void SetScale(VectorHandler& XScale) const;

#if 0
	/* DataOut: entita' che richiedono solo l'output */
protected:
	struct {
		DataOut *pFirstDataOut;
		integer iNum;
	} OutData[OutData::LASTDATAOUTTYPE];
	integer iTotDataOut;
	DataOut **pDataOut;

public:
	void OutManager(void);
	void OutManagerDestructor(void);
#endif /* 0 */

public:
	const VectorHandler* GetpXCurr(void) const {
		return pXCurr;
	};

	const VectorHandler* GetpXPCurr(void) const {
		return pXPrimeCurr;
	}

public:
	virtual void PrintResidual(const VectorHandler& Res, integer iIterCnt) const;
	virtual void PrintSolution(const VectorHandler& Sol, integer iIterCnt) const;

	virtual const std::string& GetDofDescription(int i) const;
	virtual const std::string& GetEqDescription(int i) const;
	virtual DofOrder::Order GetDofType(int i) const;
	virtual DofOrder::Order GetEqType(int i) const;
};

// if bActive is true, the cast only succeeds when driven element is active
// otherwise it always succeeds
template <class T>
T*
DataManager::Cast(Elem *pEl, bool bActive)
{
	ASSERT(pEl != NULL);

	T *pT = dynamic_cast<T *>(pEl);

	if (pT == 0) {
		DrivenElem *pDE = dynamic_cast<DrivenElem *>(pEl);
		if (pDE == 0) {
			silent_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") as \"" << mbdyn_demangle<T>() << "\" (not driven)" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (bActive && !pDE->bIsActive()) {
			pedantic_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") as \"" << mbdyn_demangle<T>() << "\""
				" (driven but currently inactive)" << std::endl);
			return 0;
		}

		pT = dynamic_cast<T *>(pDE->pGetElem());
		if (pT == 0) {
			pedantic_cerr("unable to cast "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << ") as \"" << mbdyn_demangle<T>() << "\""
				" (driven but cast failed)" << std::endl);
		}
	}

	return pT;
}

template <class Tbase, Node::Type type>
Tbase *
DataManager::pFindNode(unsigned int uNode) const
{
	/* verifica di esistenza del nodo */
	Node* pNode = pFindNode(type, uNode);
	if (pNode == 0) {
		silent_cerr("DataManager::pFindNode: " << psNodeNames[type] << "(" << uNode << ") not found" << std::endl);
		return 0;
	}

	Tbase *pNodeBase = dynamic_cast<Tbase *>(pNode);
	if (pNodeBase == 0) {
		silent_cerr("DataManager::pFindNode: unable to cast " << psNodeNames[type] << "(" << pNode->GetLabel() << ") "
			"to \"" << mbdyn_demangle<Tbase>() << "\"" << std::endl);
		return 0;
	}

	return pNodeBase;
}

template <class Tder, class Tbase, Node::Type type>
Tder *
DataManager::pFindNode(unsigned int uNode) const
{
	Tbase *pNodeBase = pFindNode<Tbase, type>(uNode);
	if (pNodeBase == 0) {
		return 0;
	}

	Tder *pNodeDer = dynamic_cast<Tder *>(pNodeBase);
	if (pNodeDer == 0) {
		silent_cerr("DataManager::pFindNode: unable to cast " << psNodeNames[type] << "(" << pNodeBase->GetLabel() << ") "
			"from \"" << mbdyn_demangle<Tbase>() << "\" "
			"to \"" << mbdyn_demangle<Tder>() << "\"" << std::endl);
		return 0;
	}

	return pNodeDer;
}

template <class Tbase, Elem::Type type>
Tbase *
DataManager::pFindElem(unsigned int uElem) const
{
	/* verifica di esistenza dell'elemento */
	Elem* pElem = pFindElem(type, uElem);
	if (pElem == 0) {
		silent_cerr("DataManager::pFindElem: " << psElemNames[type] << "(" << uElem << ") not found" << std::endl);
		return 0;
	}

	Tbase *pElemBase = dynamic_cast<Tbase *>(pElem);
	if (pElemBase == 0) {
		silent_cerr("DataManager::pFindElem: unable to cast " << psElemNames[type] << "(" << pElem->GetLabel() << ") "
			"to \"" << mbdyn_demangle<Tbase>() << "\"" << std::endl);
		return 0;
	}

	return pElemBase;
}

template <class Tder, class Tbase, Elem::Type type>
Tder *
DataManager::pFindElem(unsigned int uElem) const
{
	Tbase *pElemBase = pFindElem<Tbase, type>(uElem);
	if (pElemBase == 0) {
		return 0;
	}

	Tder *pElemDer = dynamic_cast<Tder *>(pElemBase);
	if (pElemDer == 0) {
		silent_cerr("DataManager::pFindElem: unable to cast " << psElemNames[type] << "(" << pElemBase->GetLabel() << ") "
			"from \"" << mbdyn_demangle<Tbase>() << "\" "
			"to \"" << mbdyn_demangle<Tder>() << "\"" << std::endl);
		return 0;
	}

	return pElemDer;
}

template <class Tbase, Node::Type type>
Tbase *
DataManager::ReadNode(MBDynParser& HP) const
{
	Node *pNode = ReadNode(HP, type);
	ASSERT(pNode != 0);

	Tbase *pNodeBase = dynamic_cast<Tbase *>(pNode);
	if (pNodeBase == 0) {
		silent_cerr("DataManager::ReadNode: unable to cast " << psNodeNames[type] << "(" << pNode->GetLabel() << ") "
			"to \"" << mbdyn_demangle<Tbase>() << "\" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pNodeBase;
}

template <class Tder, class Tbase, Node::Type type>
Tder *
DataManager::ReadNode(MBDynParser& HP) const
{
	Tbase *pNodeBase = ReadNode<Tbase, type>(HP);
	ASSERT(pNodeBase != 0);

	Tder *pNodeDer = dynamic_cast<Tder *>(pNodeBase);
	if (pNodeDer == 0) {
		silent_cerr("DataManager::ReadNode: unable to cast " << psNodeNames[type] << "(" << pNodeBase->GetLabel() << ") "
			"from \"" << mbdyn_demangle<Tbase>() << "\" "
			"to \"" << mbdyn_demangle<Tder>() << "\" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pNodeDer;
}

template <class Tbase, Elem::Type type>
Tbase *
DataManager::ReadElem(MBDynParser& HP) const
{
	Elem *pElem = ReadElem(HP, type);
	ASSERT(pElem != 0);

	Tbase *pElemBase = dynamic_cast<Tbase *>(pElem);
	if (pElemBase == 0) {
		silent_cerr("DataManager::ReadElem: unable to cast " << psElemNames[type] << "(" << pElem->GetLabel() << ") "
			"to \"" << mbdyn_demangle<Tbase>() << "\" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pElemBase;
}

template <class Tder, class Tbase, Elem::Type type>
Tder *
DataManager::ReadElem(MBDynParser& HP) const
{
	Tbase *pElemBase = ReadElem<Tbase, type>(HP);
	ASSERT(pElemBase != 0);

	Tder *pElemDer = dynamic_cast<Tder *>(pElemBase);
	if (pElemBase == 0) {
		silent_cerr("DataManager::ReadElem: unable to cast " << psElemNames[type] << "(" << pElemBase->GetLabel() << ") "
			"from \"" << mbdyn_demangle<Tbase>() << "\" "
			"to \"" << mbdyn_demangle<Tder>() << "\" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pElemDer;
}

template <class T>
flag
DataManager::fReadOutput(MBDynParser& HP, const T& t) const
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

	return HP.GetBool(fDef);
}

/* DataManager - end */


/* Usato per iterare sugli elementi che partecipano
 * all'assemblaggio iniziale dei vincoli */

/* InitialAssemblyIterator - begin */

class InitialAssemblyIterator {
private:
	const DataManager::ElemDataStructure (*pElemData)[Elem::LASTELEMTYPE];
	mutable Elem::Type m_FirstType;
	mutable DataManager::ElemContainerType::const_iterator m_CurrElem;
	mutable Elem::Type m_CurrType;

public:
	InitialAssemblyIterator(const DataManager::ElemDataStructure
			(*pED)[Elem::LASTELEMTYPE]);
	InitialAssemblyElem* GetFirst(void) const;
	InitialAssemblyElem* GetNext(void) const;
};

/* InitialAssemblyIterator - end */

extern "C" int
datamanager_cleanup(void *);

extern ScalarDof
ReadScalarDof(const DataManager* pDM, MBDynParser& HP, bool bDof, bool bOrder);

extern OrientationDescription
ReadOrientationDescription(MBDynParser& HP);

extern OrientationDescription
ReadOptionalOrientationDescription(DataManager *pDM, MBDynParser& HP);

#endif /* DATAMAN_H */

