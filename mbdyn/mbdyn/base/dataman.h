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

/* gestore dei dati */

#ifndef DATAMAN_H
#define DATAMAN_H

#include <ac/iostream>

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

#include <mbpar.h>
#include <constltp.h>
#include <shape.h>

/* da elman.h */
#include <solman.h>
#include <submat.h>
#include <veciter.h>

#include <elem.h>      /* Classe di base di tutti gli elementi */
#include <output.h>

#include <drive.h>     /* Drive vari */

/* da nodeman.h */
#include <node.h>
#include <strnode.h>
#include <elecnode.h>

/* DataManager - begin */

const int iGlobalSymbolTableInitialSize = 21;

class DataManager : public SolutionDataManager {
   friend DriveCaller* ReadDriveData(const DataManager*, MBDynParser&, const DriveHandler*);   
   friend ScalarDof ReadScalarDof(const DataManager*, MBDynParser&, flag);
   
 public: 
   class ErrGeneric {};
   class ErrAssemblyDiverged {};
   class ErrAssemblyMaxIters {};
   class ErrElemNotAllowedInAssembly {};
   class ErrUnknownElem {};
   class ErrUnknownNode {};
   class ErrMissingNodes {};
   
 private:
   /* Handler vari */
   MathParser& MathPar;      /* Received from MultiStepIntegrator */
   Table& GlobalSymbolTable; /* note: do not invert declaration order */

 protected:
   DriveHandler DrvHdl;
   mutable OutputHandler OutHdl;
   
   /* Puntatore alla variabile Time nella GlobalSymbolTable */
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
#endif /* USE_STRUCT_NODES */

#if defined(HAVE_LOADABLE) && defined(HAVE_LTDL_H)
   bool loadableElemInitialized;
#endif /* HAVE_LOADABLE && HAVE_LTDL_H */

   flag fPrintDofStats;

   /* Parametri vari */
   char* sSimulationTitle;

 protected:  
   enum eRestart { NEVER, ATEND, ITERATIONS, TIME };
   eRestart RestartEvery;
   integer iRestartIterations;
   doublereal dRestartTime;
   
   integer iCurrRestartIter;
   doublereal dLastRestartTime;

   integer iOutputFrequency;
   mutable integer iOutputCount;

 public:
   enum ResType {
	   RES_NONE = 0x00,
	   RES_NATIVE = 0x01,
	   RES_ADAMS = 0x02,
	   RES_MOTIONVIEW = 0x04
   };

 protected:
   int ResMode;

#ifdef USE_ADAMS 
   char *sAdamsModelName;
   integer iAdamsOutputBlock;
   unsigned int iAdamsOutputNodes;
   unsigned int iAdamsOutputParts;
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW

#endif /* USE_MOTIONVIEW */
   
 private:
   /* chiamate dal costruttore per leggere i relativi articoli */
   void ReadControl(MBDynParser& HP, const char* sInputFileName, const char* sOutputFileName);
   void ReadNodes(MBDynParser& HP);
   void ReadElems(MBDynParser& HP);
   void ReadDrivers(MBDynParser& HP);
         
   /* read functions */
   friend Node* ReadStructNode(DataManager* pDM, MBDynParser& HP, DofOwner* pDO, unsigned int uLabel);   
   friend Elem** ReadOneElem(DataManager* pDM, MBDynParser& HP, unsigned int uLabel, int CurrType);   
   friend Elem* ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadJoint(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
   friend Elem* ReadGenel(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
   friend Elem* ReadElectric(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
   friend Elem* ReadBulk(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadForce(DataManager* pDM, MBDynParser& HP, unsigned int uLabel, flag = 0);  
   friend Elem* ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadBeam2(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadHBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadRotor(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
   friend Elem* ReadAerodynamicBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadAerodynamicBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadAerodynamicBeam2(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadAerodynamicModal(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
#ifdef USE_EXTERNAL
   friend Elem* ReadAerodynamicExternal(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
   friend Elem* ReadAerodynamicExternalModal(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
#endif /*  USE_EXTERNAL */
   friend Elem* ReadHydraulicElem(DataManager* pDM, MBDynParser& HP, const DofOwner* pDO, unsigned int uLabel);
   friend Drive* ReadFileDriver(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);   
   
 public:
   flag fReadOutput(MBDynParser& HP, enum Elem::Type t);
   flag fReadOutput(MBDynParser& HP, enum Node::Type t);

   doublereal dReadScale(MBDynParser& HP, enum DofOwner::Type t);
   
 private:
   /* legge i legami costitutivi */
   ConstitutiveLaw1D* ReadConstLaw1D(MBDynParser& HP, DefHingeType::Type& T);
   ConstitutiveLaw3D* ReadConstLaw3D(MBDynParser& HP, DefHingeType::Type& T);
   ConstitutiveLaw6D* ReadConstLaw6D(MBDynParser& HP, DefHingeType::Type& T);
   
   enum ConstLawDim {
      SCALAR = 0,
	DIM3,
	DIM6,
	LASTDIM
   };   
   
   /* chiamate a funzioni di inizializzazione */
#if defined(USE_STRUCT_NODES)   
   void InitialJointAssembly(void);
#endif /* USE_STRUCT_NODES */
   
   void DofOwnerSet(void);
   void DofOwnerInit(void);
   
   
 public:
   
   /* costruttore - legge i dati e costruisce le relative strutture */
   DataManager(MBDynParser& HP, 
	       doublereal dInitialTime,
	       const char* sInputFileName,
	       const char* sOutputFileName, 
	       flag fAbortAfterInput);
   
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

   /* Setta il valore della variabile tempo nella tabella dei simboli del
    * DataManager e nel DriveHandler */
   void SetTime(doublereal dTime);

   NamedValue *InsertSym(const char* const s, const Real& v, int redefine = 0);
   NamedValue *InsertSym(const char* const s, const Int& v, int redefine = 0);
   
   /* Collega il DataManager ed il DriveHandler ai vettori soluzione */
   void LinkToSolution(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr);

   /* Restituisce il numero di dof per la costruzione delle matrici ecc. */
   integer iGetNumDofs(void) const { 
      return iTotDofs; 
   };

   /* Restituisce il puntatore alla struttura dei dof */
   VecIter<Dof>& GetDofIterator(void) /* const */ { 
      return DofIter; 
   };

   /* Restituisce l'ostream al file di output, 
    * usato dai vari metodi per scrivere il log del calcolo */
   std::ostream& GetOutFile(void) { return OutHdl.Output(); };
   
   /* Restituisce il DriveHandler */
   const DriveHandler* pGetDrvHdl(void) const { return &DrvHdl; };
   
   /* Assembla lo jacobiano */
   virtual void AssJac(MatrixHandler& JacHdl, doublereal dCoef);
   
   /* Assembla le matrici per gli autovalori */
   virtual void AssEig(MatrixHandler& A_Hdl, MatrixHandler& B_Hdl);
   
   /* Assembla il residuo */
   virtual void AssRes(VectorHandler &ResHdl, doublereal dCoef);
   
   /* stampa i risultati */
   virtual void Output(bool force = false) const;
   virtual void Output(const VectorHandler& X, const VectorHandler& XP) const;
   virtual void Output_pch(std::ostream& pch) const;
   virtual void Output_f06(std::ostream& f06, const VectorHandler& X) const;
   virtual void Output_f06(std::ostream& f06, const VectorHandler& Xr, const VectorHandler& Xi) const;
   virtual void Output_OpenDX(std::ostream& dx, const VectorHandler& Xr, const VectorHandler& Xi) const;

#ifdef USE_ADAMS 
   /* MSC's ADAMS/View .res output */
   bool fAdamsOutput(void) const;
   void AdamsResOutputInit(void);
   void AdamsResOutput(integer iBlock, const char *type, const char *id) const;
#endif /* USE_ADAMS */

#ifdef USE_MOTIONVIEW 
   /* Altair's Motion View output */
   bool fMotionViewOutput(void) const;
   void MotionViewResOutputInit(const char *sOutputFileName);
   void MotionViewResOutput(integer iBlock, const char *type, const char *id) const;
   void MotionViewResOutputFini(void) const;
#endif /* USE_MOTIONVIEW */

   /* Prepara la soluzione con i valori iniziali */
   void SetValue(VectorHandler& X, VectorHandler& XP);
   
   /* Funzioni di aggiornamento dati durante la simulazione */
   virtual void MakeRestart(void);
   virtual void DerivativesUpdate(void) const;
   virtual void BeforePredict(VectorHandler& X, VectorHandler& XP, VectorHandler& XPrev, VectorHandler& XPPrev) const;        
   virtual void AfterPredict(void) const;        
   virtual void Update(void) const;
   virtual void AfterConvergence(void) const;        

   
   
   /* da ElemManager */
   friend class ElemManIterator;
   friend class InitialAssemblyIterator;
     
 public:
   enum DerivationTable { 
      ELEM = 0,  // pleonastico
	DOFOWNER = 1,
	GRAVITYOWNER = 2,
	AIRPROPOWNER = 4,
	INITIALASSEMBLY = 8
   };
   
 protected:
   
   /* struttura dei dati fondamentali degli elementi */
   struct ElemDataStructure {
      Elem** ppFirstElem; /* puntatore al puntatore al primo el. del tipo */
      unsigned int iNum;          /* numero di elementi del tipo */
      DofOwner::Type DofOwnerType; /* Tipo di DofOwner */
      unsigned int iDerivation;   /* Tabella delle derivazioni */
      flag fIsUnique;             /* Flag di elemento unico */
      flag fToBeUsedInAssembly;   /* Se deve essere usato nell'assemblaggio */
      flag fGeneratesInertialForces; /* Se genera forze d'inerzia */
      flag fUsesAirProperties;    /* Usa le proprieta' dell'aria */
      flag fDefaultOut;           /* Flag di default output */
      
      OutputHandler::OutFiles OutFile; /* Tipo di file in output */
      
   } ElemData[Elem::LASTELEMTYPE];
   
   mutable VecIter<Elem*> ElemIter;
   
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
   integer iMaxWorkNumRows; /* dimensioni massime della matrice di lavoro */
   integer iMaxWorkNumCols; /*    ''         ''   */
   integer iWorkIntSize;    /* dimensioni degli spazi di lavoro */
   integer iWorkDoubleSize;
   integer* piWorkIndex;    /* puntatore ad interi - array di lavoro */
   doublereal* pdWorkMat;   /* puntatore a double - matrice di lavoro */
  
   VariableSubMatrixHandler* pWorkMatA;  /* SubMatrix di lavoro */
   VariableSubMatrixHandler* pWorkMatB;
   


   /* ricerca elementi*/
   void* pFindElem(Elem::Type Typ, unsigned int uL) const;
   void* pFindElem(Elem::Type Typ, unsigned int uL, unsigned int iDeriv) const;
   void* pChooseElem(Elem* p, unsigned int iDeriv) const;
   
   Elem** ppFindElem(Elem::Type Typ, unsigned int uL) const;
   
   /* ricerca drives */
   void* pFindDrive(Drive::Type Typ, unsigned int uL) const;
   
   
   flag fGetDefaultOutputFlag(const Elem::Type& t) const;
   
 public:
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
   void ElemOutput_pch(std::ostream& pch) const;
   void ElemOutput_f06(std::ostream& f06, const VectorHandler& X) const;
   void ElemOutput_f06(std::ostream& f06, const VectorHandler& Xr, const VectorHandler& Xi) const;
   
   
   /* da NodeManager */
 protected:
   
   /* struttura dei dati dei nodi. Per ogni tipo: 
    * puntatore al puntatore al primo dato, numero degli item per tipo */
   struct {
      Node** ppFirstNode;
      unsigned int iNum;
      flag fDefaultOut;
      
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
   void NodeOutput_pch(std::ostream& pch) const;
   void NodeOutput_f06(std::ostream& f06, const VectorHandler& X) const;
   void NodeOutput_f06(std::ostream& f06, const VectorHandler& Xr, const VectorHandler& Xi) const;
   
   
   /* da DofManager */
 protected:

   /* struttura dei dati generali dei dof: numero totale per tipo, 
    * dimensione caratteristica (se esiste), puntatore al primo del tipo */
   struct {
      DofOwner* pFirstDofOwner;     /* punt. al primo DofOwner di ogni tipo */
      integer iNum;                 /* numero di DofOwners per ogni tipo */
      integer iSize;                /* numero di Dof (se fisso, es. nodi) */
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
   Elem::Type CurrType;
   Elem** ppCurr;
 
 public:
   InitialAssemblyIterator(const DataManager::ElemDataStructure 
			   (*pED)[Elem::LASTELEMTYPE]);
   InitialAssemblyElem* GetFirst(void) const;
   InitialAssemblyElem* GetNext(void) const;
};

/* InitialAssemblyIterator - end */



extern DriveCaller* ReadDriveData(const DataManager* pDM, MBDynParser& HP, const DriveHandler* pDH);
extern ScalarDof ReadScalarDof(const DataManager* pDM, MBDynParser& HP, flag fOrder);

#if (defined(USE_STRUCT_NODES) && defined(USE_AERODYNAMIC_ELEMS))
extern Shape* ReadShape(MBDynParser& HP);
#endif /* STRUCT && AERODYNAMIC */

#endif /* DATAMAN_H */

