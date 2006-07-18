/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#include "ac/iostream"
#include <list>
#include <map>
#include <string>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "mbpar.h"
#include "constltp.h"
#include "shape.h"

/* da elman.h */
#include "solman.h"
#include "submat.h"
#include "veciter.h"

#include "elem.h"      /* Classe di base di tutti gli elementi */
#include "output.h"

#include "drive.h"     /* Drive vari */

/* da nodeman.h */
#include "node.h"
#include "strnode.h"
#include "elecnode.h"

#include "solverdiagnostics.h"
#include "linsol.h"

struct LoadableCalls;


#include "usesock.h"

class Solver;

/* DataManager - begin */

class DataManager : public SolutionDataManager, public SolverDiagnostics {
public:
	class ErrGeneric {};
	class ErrAssemblyDiverged {};
	class ErrAssemblyMaxIters {};
	class ErrElemNotAllowedInAssembly {};
	class ErrUnknownElem {};
	class ErrUnknownFunction {};
	class ErrUnknownNode {};
	class ErrMissingNodes {};
	class ErrNeedDataManager {};

private:
#ifdef USE_MULTITHREAD
	/* from input file, or auto-detected */
	unsigned int nThreads;
#endif /* USE_MULTITHREAD */
	
	/* Handler vari */
	MathParser& MathPar;      /* Received from MultiStepIntegrator */
	Solver* pSolver;

	/* loadable elements */
	std::map<std::string, const LoadableCalls *> MapOfLoadableElemHandlers;

protected:
	DriveHandler DrvHdl;
	mutable OutputHandler OutHdl;

	/* Puntatore alla variabile Time nella symbol table di MathPar */
	Var* pTime;

	/* Puntatori ai vettori soluzione durante il passo */
	const VectorHandler* pXCurr;
	const VectorHandler* pXPrimeCurr;

private:
	/* Parametri usati durante l'assemblaggio iniziale */
#if defined(USE_STRUCT_NODES)
	flag fInitialJointAssemblyToBeDone;
	flag fSkipInitialJointAssembly;
	flag fOutputFrames;
	doublereal dInitialPositionStiffness;
	doublereal dInitialVelocityStiffness;
	flag fOmegaRotates;
	doublereal dInitialAssemblyTol;
	integer iMaxInitialIterations;
	doublereal dEpsilon;
	LinSol CurrSolver;

	bool bStaticModel;
#endif /* USE_STRUCT_NODES */

#if defined(HAVE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
	bool loadableElemInitialized;
#endif /* HAVE_RUNTIME_LOADING && HAVE_LTDL_H */

	enum PrintFlags {
		PRINT_NONE		= 0x00U,
		PRINT_DOFSTATS		= 0x01U,
		PRINT_DOFDESCRIPTION	= 0x02U,
		PRINT_EQDESCRIPTION	= 0x04U
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

	/* specialized output stuff */
public:
	enum ResType {
		RES_NONE	= 0x00,
		RES_NATIVE	= 0x01,
		RES_ADAMS	= 0x02,
		RES_MOTIONVIEW	= 0x04
	};

	bool bOutput(ResType t) const;

protected:
	int ResMode;

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

private:
	/* chiamate dal costruttore per leggere i relativi articoli */
	void ReadControl(MBDynParser& HP, const char* sOutputFileName, 
				const char* sInputFileName);
	void ReadNodes(MBDynParser& HP);
	void ReadDrivers(MBDynParser& HP);
	void ReadElems(MBDynParser& HP);

public:
	flag fReadOutput(MBDynParser& HP, enum Elem::Type t);
	flag fReadOutput(MBDynParser& HP, enum Node::Type t);

	doublereal dReadScale(MBDynParser& HP, enum DofOwner::Type t);

	const doublereal& dGetInitialPositionStiffness(void) const;
	const doublereal& dGetInitialVelocityStiffness(void) const;
	flag fDoesOmegaRotate(void) const;

	void IncElemCount(Elem::Type type);

	/* legge i legami costitutivi */
	ConstitutiveLaw1D* ReadConstLaw1D(MBDynParser& HP,
			ConstLawType::Type& T);
	ConstitutiveLaw3D* ReadConstLaw3D(MBDynParser& HP,
			ConstLawType::Type& T);
	ConstitutiveLaw6D* ReadConstLaw6D(MBDynParser& HP,
			ConstLawType::Type& T);

	enum ConstLawDim {
		SCALAR = 0,
		DIM3,
		DIM6,
		LASTDIM
	};

	/* additional CPU time, if any */
	virtual clock_t GetCPUTime(void) const {
		return 0;
	};

private:
	/* chiamate a funzioni di inizializzazione */
#if defined(USE_STRUCT_NODES)
	void InitialJointAssembly(void);
#endif /* USE_STRUCT_NODES */

	void DofOwnerSet(void);
	void DofOwnerInit(void);
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

	/* helpers */
	int ReadScalarAlgebraicNode(MBDynParser& HP, unsigned int uLabel,
			Node::Type type, doublereal& dX);
	int ReadScalarDifferentialNode(MBDynParser& HP, unsigned int uLabel,
			Node::Type type, doublereal& dX, doublereal& dXP);
	Node* ReadNode(MBDynParser& HP, Node::Type type);
	Elem* ReadElem(MBDynParser& HP, Elem::Type type);

	/* Funzioni usate dal metodo di integrazione */

	/* Setta il valore della variabile tempo nella tabella dei simboli
	 * del DataManager e nel DriveHandler */
	void SetTime(doublereal dTime, bool bServePending = true);
	doublereal dGetTime() const;

	NamedValue *InsertSym(const char* const s, const Real& v,
			int redefine = 0);
	NamedValue *InsertSym(const char* const s, const Int& v,
			int redefine = 0);

	/* Collega il DataManager ed il DriveHandler ai vettori soluzione */
	void LinkToSolution(const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/* Restituisce il numero di dof per la costruzione delle matrici ecc. */
	integer iGetNumDofs(void) const { return iTotDofs; };

	/* Restituisce il puntatore alla struttura dei dof */
	VecIter<Dof>& GetDofIterator(void) /* const */ { return DofIter; };

	/* Restituisce l'ostream al file di output,
	 * usato dai vari metodi per scrivere il log del calcolo */
	std::ostream& GetOutFile(void) { return OutHdl.Output(); };
	std::ostream& GetLogFile(void) { return OutHdl.Log(); };

	/* Restituisce il DriveHandler */
	const DriveHandler* pGetDrvHdl(void) const { return &DrvHdl; };
	MathParser& GetMathParser(void) const { return MathPar; };

	/* Assembla lo jacobiano */
	virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef);

	/* Assembla le matrici per gli autovalori */
	virtual void AssMats(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl);

	/* Assembla il residuo */
	virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef)
		throw(ChangedEquationStructure);

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
			SubVectorHandler& WorkVec);

public:
	/* stampa i risultati */
	virtual void Output(bool force = false) const;
	virtual void Output(const VectorHandler& X,
			const VectorHandler& XP) const;
#if 0
	virtual void Output_pch(std::ostream& pch) const;
	virtual void Output_f06(std::ostream& f06,
			const VectorHandler& X) const;
	virtual void Output_f06(std::ostream& f06, const VectorHandler& Xr,
			const VectorHandler& Xi) const;
#endif
#if 0
	virtual void Output_OpenDX(std::ostream& dx, const VectorHandler& Xr,
			const VectorHandler& Xi) const;
#endif

	/* Aggiungere qui le funzioni che aprono i singoli stream */
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

	bool bIsStaticModel(void) const {
		return bStaticModel;
	};

	/* da ElemManager */
	friend class InitialAssemblyIterator;

	enum ModuleInsertMode {
		FAIL,
		REPLACE,
		IGNORE
	};

	/* loadable elements */
	const LoadableCalls *GetLoadableElemModule(std::string) const;
	void SetLoadableElemModule(std::string, const LoadableCalls *,
			ModuleInsertMode = FAIL);

public:
	/* FIXME: will be eliminated */
	enum DerivationTable {
		ELEM			= 0,  // pleonastico
		DOFOWNER		= 1,
		GRAVITYOWNER		= 2,
		AIRPROPOWNER		= 4,
		INITIALASSEMBLY		= 8
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

protected:
	/* struttura dei dati fondamentali degli elementi */
	struct ElemDataStructure {
		Elem** ppFirstElem;		// punt. al punt. al primo el. del tipo
		unsigned int iNum;		// numero di elementi del tipo

		DofOwner::Type DofOwnerType;	// Tipo di DofOwner
		unsigned int iDerivation;	// Tabella delle derivazioni

		OutputHandler::OutFiles OutFile;	// Tipo di file in output

		unsigned uFlags;		// flags

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
		typedef std::map<std::string, DataManager::ElemRead *> ElemReadType;
		ElemReadType ElemRead;

	} ElemData[Elem::LASTELEMTYPE];

#if 0
	/* element type map; will replace ElemData */
	typedef std::map<std::string, ElemDataStructure *> ElemDataMapType;
	ElemDataStrMapType ElemDataMap;
#endif

	mutable VecIter<Elem *> ElemIter;

	Elem** ppElems;         /* puntatore all'array di puntatori agli el. */
	unsigned int iTotElem;  /* numero totale di el. definiti */

	/* struttura dei drivers */
	struct {
		Drive** ppFirstDrive;
		unsigned int iNum;
	} DriveData[Drive::LASTDRIVETYPE];

	Drive** ppDrive;         /* puntatore ai drivers */
	unsigned int iTotDrive;  /* numero totale dei drivers */

	/* dati di lavoro */
	integer iMaxWorkNumRows; /* dimensioni max della matrice di lavoro */
	integer iMaxWorkNumCols; /*    ''         ''   */
	integer iWorkIntSize;    /* dimensioni degli spazi di lavoro */
	integer iWorkDoubleSize;

	VariableSubMatrixHandler *pWorkMatA;  /* SubMatrix di lavoro */
	VariableSubMatrixHandler *pWorkMatB;
	VariableSubMatrixHandler *pWorkMat;
	MySubVectorHandler *pWorkVec;

	/* ricerca elementi*/
	void* pFindElem(Elem::Type Typ, unsigned int uL,
			unsigned int iDeriv) const;
	void* pChooseElem(Elem* p, unsigned int iDeriv) const;

	Elem** ppFindElem(Elem::Type Typ, unsigned int uL) const;

	flag fGetDefaultOutputFlag(const Elem::Type& t) const;
	Elem** ReadOneElem(MBDynParser& HP,
			unsigned int uLabel,
			int CurrType);

public:
	/* ricerca drives */
	void* pFindDrive(Drive::Type Typ, unsigned int uL) const;

	/* ricerca elementi*/
	void* pFindElem(Elem::Type Typ, unsigned int uL) const;

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
	void ElemOutput(OutputHandler& OH) const;
	void ElemOutput(OutputHandler& OH,
			const VectorHandler& X, const VectorHandler& XP) const;
#if 0
	void ElemOutput_pch(std::ostream& pch) const;
	void ElemOutput_f06(std::ostream& f06, const VectorHandler& X) const;
	void ElemOutput_f06(std::ostream& f06, const VectorHandler& Xr,
			const VectorHandler& Xi) const;
#endif

	/* da NodeManager */
protected:

	/* struttura dei dati dei nodi. Per ogni tipo:
	 * puntatore al puntatore al primo dato, numero degli item per tipo */
	struct {
		Node** ppFirstNode;
		unsigned int iNum;
		unsigned uFlags;		// flags

		/* helpers */
		void DefaultOut(bool b) { if (b) { uFlags |= DEFAULTOUT; } else { uFlags &= ~DEFAULTOUT; } };

		bool bDefaultOut(void) const { return (uFlags & DEFAULTOUT) == DEFAULTOUT; };

		OutputHandler::OutFiles OutFile; /* Tipo di file in output */
	} NodeData[Node::LASTNODETYPE];

	VecIter<Node*> NodeIter;

	/* dati dei nodi: numero totale e puntatore all'array dei dati
	 * (ogni nodo ha il suo formato caratteristico, comunque derivato
	 * dalla classe Node) */
	unsigned int iTotNodes;
	Node** ppNodes;

public:
	Node** ppFindNode(Node::Type Typ, unsigned int uL) const;
	/* ricerca di nodi */
	Node* pFindNode(Node::Type Typ, unsigned int uL) const;
#if defined(USE_STRUCT_NODES)
	StructNode* pFindStructNode(unsigned int uL) const;
#endif // USE_STRUCT_NODES

#if defined(USE_ELECTRIC_NODES)
	ElectricNode* pFindElectricNode(unsigned int uL) const;
#endif // USE_ELECTRIC_NODES

protected:
	flag fGetDefaultOutputFlag(const Node::Type& t) const;

public:
	/* Pseudocostruttore */
	void NodeManager(void);
	void NodeManagerDestructor(void);

	/* inizializza le matrici ed alloca memoria */
	void NodeDataInit(void);

	/* scrive i dati dei nodi */
	void NodeOutput(OutputHandler& OH) const;
	void NodeOutput(OutputHandler& OH,
			const VectorHandler& X, const VectorHandler& XP) const;
#if 0
	void NodeOutput_pch(std::ostream& pch) const;
	void NodeOutput_f06(std::ostream& f06, const VectorHandler& X) const;
	void NodeOutput_f06(std::ostream& f06, const VectorHandler& Xr,
			const VectorHandler& Xi) const;
#endif

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

	/* struttura dei dati di ogni singolo dof: totale dei dof; per ognuno:
	 * indice e tipo (algebrico o differenziale) */
	integer iTotDofs;                /* numero totale di Dof */
	Dof* pDofs;                      /* puntatore all'array dei Dof */
	VecIter<Dof> DofIter;            /* Iteratore dei Dof */

	DofOwner DummyDofOwner; /* Per quelli che non hanno dof */

	doublereal dGetDefaultScale(DofOwner::Type t) const;
public:
	/* pseudocostruttore */
	void DofManager(void);
	void DofManagerDestructor(void);

	/* funzioni di inizializzazione */
	void DofDataInit(void);
	void DofInit(void);

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

	/* socket select stuff */
private:
	std::map<int, UseSocket *> SocketUsers;
	time_t SocketUsersTimeout;

	void WaitSocketUsers(void);

public:
	void RegisterSocketUser(UseSocket *pUS);
};

/* DataManager - end */


/* Usato per iterare sugli elementi che partecipano
 * all'assemblaggio iniziale dei vincoli */

/* InitialAssemblyIterator - begin */

class InitialAssemblyIterator {
private:
	const DataManager::ElemDataStructure (*pElemData)[Elem::LASTELEMTYPE];
	const Elem::Type FirstType;
	const Elem** ppFirst;
	mutable Elem::Type CurrType;
	mutable Elem** ppCurr;

public:
	InitialAssemblyIterator(const DataManager::ElemDataStructure
			(*pED)[Elem::LASTELEMTYPE]);
	InitialAssemblyElem* GetFirst(void) const;
	InitialAssemblyElem* GetNext(void) const;
};

/* InitialAssemblyIterator - end */


extern ScalarDof
ReadScalarDof(const DataManager* pDM, MBDynParser& HP, flag fOrder);

#if (defined(USE_STRUCT_NODES) && defined(USE_AERODYNAMIC_ELEMS))
extern Shape* ReadShape(MBDynParser& HP);
#endif /* STRUCT && AERODYNAMIC */

/* used by maps to compare strings case-insensitive */
struct ltstrcase {
	/* case-insensitive string comparison */
	bool operator()(const std::string& s1, const std::string& s2) const {
		return strcasecmp(s1.c_str(), s2.c_str()) < 0;
	};
};

#endif /* DATAMAN_H */

