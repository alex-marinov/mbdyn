/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>

extern "C" {
#include <strings.h>
#include <time.h>
}

#if defined(HAVE_LOADABLE) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif /* HAVE_LOADABLE && HAVE_LTDL_H */

#include <constltp.h>
#include <dataman_.h>
/* plugins per math parser */
#include <dofpgin.h>
#include <dummypgin.h>
#ifdef USE_TCL
#include <tclpgin.h>
#endif /* USE_TCL */

/* temporaneo? */
#include <beam.h>
		

/* DataManager - begin */

/* linka i singoli DriveCaller al DriveHandler posseduto dal DataManager */
extern void SetDrvHdl(DriveHandler*);

const flag fDefaultInitialJointAssemblyToBeMade = flag(0);
const flag fDefaultSkipInitialJointAssembly = flag(0);
const doublereal dDefaultInitialStiffness = 1.;
const flag fDefaultOmegaRotates = flag(0);
const doublereal dDefaultInitialAssemblyToll = 1.e-6;
const integer iDefaultMaxInitialIterations = 1;

const char sDefaultOutputFileName[] = "MBDyn";



/*
 * costruttore: inizializza l'oggetto, legge i dati e crea le strutture di
 * gestione di Dof, nodi, elementi e drivers.
 */

DataManager::DataManager(MBDynParser& HP, 
			 doublereal dInitialTime,
			 const char* sInputFileName, 
			 const char* sOutputFileName,
			 flag fAbortAfterInput)
:
MathPar(HP.GetMathParser()),
GlobalSymbolTable((HP.GetMathParser()).GetSymbolTable()),
DrvHdl((HP.GetMathParser()).GetSymbolTable()),
OutHdl(),
pTime(NULL),
pXCurr(NULL), pXPrimeCurr(NULL),
#if defined(USE_STRUCT_NODES)
fInitialJointAssemblyToBeDone(fDefaultInitialJointAssemblyToBeMade),
fSkipInitialJointAssembly(fDefaultSkipInitialJointAssembly),
fOutputFrames(0),
dInitialPositionStiffness(dDefaultInitialStiffness), 
dInitialVelocityStiffness(dDefaultInitialStiffness),
fOmegaRotates(fDefaultOmegaRotates),
dInitialAssemblyToll(dDefaultInitialAssemblyToll),
iMaxInitialIterations(iDefaultMaxInitialIterations),
#endif /* USE_STRUCT_NODES */
#if defined(HAVE_LOADABLE) && defined(HAVE_LTDL_H)
loadableElemInitialized(false),
#endif /* HAVE_LOADABLE && HAVE_LTDL_H */
fPrintDofStats(0),
sSimulationTitle(NULL),
RestartEvery(NEVER),
iRestartIterations(0),
dRestartTime(0.),
iCurrRestartIter(0),
dLastRestartTime(dInitialTime),
iOutputFrequency(1),
iOutputCount(0),
fAdamsResOutput(0),
sAdamsModelName(NULL),
iAdamsOutputBlock(1),
iAdamsOutputNodes(0),
iAdamsOutputParts(0),

/* ElemManager */
ElemIter(), 
ppElems(NULL), 
iTotElem(0), 
ppDrive(NULL), 
iTotDrive(0),
iMaxWorkNumRows(0),
iMaxWorkNumCols(0),
iWorkIntSize(0),
iWorkDoubleSize(0),
piWorkIndex(NULL),
pdWorkMat(NULL),
pWorkMatA(NULL),
pWorkMatB(NULL),

/* NodeManager */
NodeIter(), 
iTotNodes(0), 
ppNodes(NULL),

/* DofManager */
iTotDofOwners(0), 
pDofOwners(NULL), 
iTotDofs(0), 
pDofs(NULL), 
DofIter()
{
   DEBUGCOUTFNAME("DataManager::DataManager");

   /* pseudocostruttori */
   ElemManager(); 
   NodeManager();
   DofManager();

   /* registra il plugin per i dofs */
   HP.GetMathParser().RegisterPlugIn("dof", dof_plugin, this);

#ifdef USE_TCL
   /* registra il plugin per il tcl */
   HP.GetMathParser().RegisterPlugIn("tcl", tcl_plugin, NULL);
#else /* !USE_TCL */
   HP.GetMathParser().RegisterPlugIn("tcl", dummy_plugin,
		   (void *)"configure with --with-tcl to use tcl plugin");
#endif /* USE_TCL */
   
   /* Setta il tempo al valore iniziale */
   pTime = (Var *)GlobalSymbolTable.Get("Time");
   ASSERT(pTime != NULL);
   if(pTime == NULL) {
      DEBUGCERR("");
      std::cerr << "error in setting Time symbol" << std::endl;

      THROW(DataManager::ErrGeneric());
   }   
   SetTime(dInitialTime);
   
   DEBUGLCOUT(MYDEBUG_INIT, "GlobalSymbolTable:" 
	      << GlobalSymbolTable << std::endl);
   
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
	"controldata",
	"nodes",
	"elements",
	"drivers",
	"output"
   };
   
   /* enum delle parole chiave */
   enum KeyWords { 
      BEGIN = 0,
	CONTROLDATA,
	NODES,
	ELEMENTS,
	DRIVERS,
	OUTPUT,
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
      
   /* aggiorna la tabella del parser */
   HP.PutKeyTable(K);
      
   /* legge i dati di controllo */
   if(KeyWords(HP.GetDescription()) != BEGIN) {
      DEBUGCERR("");
      std::cerr << "<begin> expected at line "
	<< HP.GetLineData() << std::endl;

      THROW(DataManager::ErrGeneric());
   }
   if(KeyWords(HP.GetWord()) != CONTROLDATA) {
      DEBUGCERR("");
      std::cerr << "<begin: control data;> expected at line " 
	<< HP.GetLineData() << std::endl;

      THROW(DataManager::ErrGeneric());
   }

   ReadControl(HP, sInputFileName, sOutputFileName);
   HP.PutKeyTable(K);
      
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
   if(iTotNodes > 0) {	
      if(KeyWords(HP.GetDescription()) != BEGIN) {
	 DEBUGCERR("");
	 std::cerr << "<begin> expected at line " 
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != NODES) {
	 DEBUGCERR("");
	 std::cerr << "<begin: nodes;> expected at line " 
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
      
      ReadNodes(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("");
      std::cerr << "warning, no nodes are defined" << std::endl;
   }   
   /* fine lettura nodi */
   
   /*
    * Costruzione struttura ElemData, DriveData
    * e creazione array Elem* e Drive*
    */
   ElemDataInit();
   
   /* legge i drivers, crea la struttura ppDrive */
   if(iTotDrive > 0) {	
      if(KeyWords(HP.GetDescription()) != BEGIN) {
	 DEBUGCERR("");
	 std::cerr << "<begin> expected at line "
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != DRIVERS) {
	 DEBUGCERR("");
	 std::cerr << "<begin: drivers;> expected at line " 
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
      
      ReadDrivers(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("warning, no drivers are defined" << std::endl);
   }   
	
   /* fine lettura drivers */
   
   
   /*
    * legge gli elementi, crea la struttura ppElems
    * e contemporaneamente aggiorna i dof
    */
   if(iTotElem > 0) {	
      if(KeyWords(HP.GetDescription()) != BEGIN) {
	 DEBUGCERR("");
	 std::cerr << "<begin> expected at line " 
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != ELEMENTS) {
	 DEBUGCERR("");
	 std::cerr << "<begin: elements;> expected at line " 
	   << HP.GetLineData() << std::endl;

	 THROW(DataManager::ErrGeneric());
      }
	
      ReadElems(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("");
      std::cerr << "warning, no elements are defined" << std::endl;
   }

#ifdef USE_STRUCT_NODES
   if (fOutputFrames) {
      OutHdl.Open(OutputHandler::REFERENCEFRAMES);
      HP.OutputFrames(OutHdl.ReferenceFrames());
   }

#ifdef USE_AERODYNAMIC_ELEMS
   if (ElemData[Elem::AIRPROPERTIES].iNum > 0) {
      OutHdl.Open(OutputHandler::AIRPROPS);
   }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

   /* fine lettura elementi */     
   
   
   /*
    * Qui non deve leggere piu' niente, 
    * quindi chiude e apre i files di output
    */
   HP.Close();
   
#ifdef __HACK_BEAM__
   for (unsigned int i = 0; i < ElemData[Elem::BEAM].iNum; i++) {
      Beam* p = (Beam*)ElemData[Elem::BEAM].ppFirstElem[i]->pGet();
      OutHdl.Output()
	<< "beam " << p->GetLabel()
	<< " " << p->pGetNode(1)->GetLabel()
	<< " " << p->pGetNode(2)->GetLabel()
	<< " " << p->pGetNode(3)->GetLabel()
	<< std::endl;
	  
   }
#endif /* __HACK_BEAM__ */
   
   for (int i = 0; i < Node::LASTNODETYPE; i++) {
      if(NodeData[i].iNum > 0 
	 && NodeData[i].OutFile != OutputHandler::UNKNOWN) {
	 OutHdl.Open(NodeData[i].OutFile);
      }
   }

   for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
      if(ElemData[i].iNum > 0 && ElemData[i].OutFile != OutputHandler::UNKNOWN) {
	 OutHdl.Open(ElemData[i].OutFile);
      }
   }
   
   /* Inizializza il drive handler */
   DrvHdl.iRandInit(0);

   
   /* Verifica dei dati di controllo */
#ifdef DEBUG
   if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
      /* mostra in modo succinto il numero di DofOwner per tipo */
      for(int i = 0; i < DofOwner::LASTDOFTYPE; i++)
	std::cout << "DofType " << i << " (" << psDofOwnerNames[i] 
	<< "), n. of owners: " << DofData[i].iNum << std::endl;
      
      /* mostra in modo succinto il numero di nodi per tipo */
      for(int i = 0; i < Node::LASTNODETYPE; i++)
	std::cout << "NodeType " << i << " (" << psNodeNames[i] 
	<< "), n. of nodes: " << NodeData[i].iNum << std::endl;
      
      /* mostra in modo succinto il numero di elementi per tipo */
      for(int i = 0; i < Elem::LASTELEMTYPE; i++)
	std::cout << "Element Type " << i << " (" << psElemNames[i] 
	<< "), n. of elems: " << ElemData[i].iNum << std::endl;
      
      /* mostra in modo succinto il numero di drivers per tipo */
      for(int i = 0; i < Drive::LASTDRIVETYPE; i++)
	std::cout << "DriveType " << i << " (" << psDriveNames[i] 
	<< "), n. of drivers: " << DriveData[i].iNum << std::endl;
   }
#endif /* DEBUG */

   if(fAbortAfterInput) {
      silent_cout("Only input is required" << std::endl);
      return;
   }
   
   /* Se richiesto, inizializza il file di output AdamsRes */
   if (fAdamsOutput()) {
      AdamsResOutputInit();
      iAdamsOutputBlock = 1;
      AdamsResOutput(iAdamsOutputBlock, "INPUT", "MBDyn");
   }   
   
   /* Qui intercetto la struttura dei Dof prima che venga costruita e modifico
    * in modo che sia adatta all'assemblaggio dei vincoli; quindi la resetto 
    * e la lascio generare correttamente.
    * L'assemblaggio iniziale viene gestito dal DataManager perche' 
    * concettualmente la consistenza dei dati deve essere assicurata da lui,
    * che conosce gli aspetti fisici del problema */

#if defined(USE_STRUCT_NODES)   
   if (fInitialJointAssemblyToBeDone) {
      if (!fSkipInitialJointAssembly) {
	 InitialJointAssembly();
      }	else {
	 silent_cout("Skipping initial joints assembly" << std::endl);
      }
   } else {
      silent_cout("No initial assembly is required since no joints are defined"
		  << std::endl);
   }
#endif /* USE_STRUCT_NODES */
   
   
   /* Costruzione dei dati dei Dof definitivi da usare nella simulazione */
   
   /* Aggiornamento DofOwner degli elementi (e dei nodi) */
   DofOwnerSet();
   
   
   /* Verifica dei dati di controllo */
#ifdef DEBUG
   if (DEBUG_LEVEL_MATCH(MYDEBUG_INIT)) {
      /* mostra in modo succinto i DofOwners */
      int k = 0;
      for(int i = 0; i < DofOwner::LASTDOFTYPE; i++) {
	 std::cout << "DofType " << i << ':' << std::endl;
	 for(int j = 0; j < DofData[i].iNum; j++) {
	    std::cout << "DofOwner " << j << ", n. of dofs: "
	      << pDofOwners[k++].iNumDofs << std::endl;
	 }
      }
   }
#endif /* DEBUG */
   
   /* a questo punto i dof sono completamente scritti in quanto a numero. 
    * Rimane da sistemare il vettore con gli indici "interni" ed il tipo */   
   
   
   /* Costruzione array DofOwner e Dof */
   DofInit();
   
   /* Aggiornamento Dof di proprieta' degli elementi e dei nodi */
   DofOwnerInit();
   
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
   
   /* Se richiesto, esegue l'output delle condizioni iniziali*/
   if (fAdamsOutput()) {
      iAdamsOutputBlock = 2;
      AdamsResOutput(iAdamsOutputBlock, "INITIAL CONDITIONS", "MBDyn");
   }
   
} /* End of DataManager::DataManager() */




/* Distruttore: se richiesto, crea il file di restart; quindi
 * chiama i distruttori degli oggetti propri e libera la memoria di lavoro */

DataManager::~DataManager(void)
{
   /* Se e' richiesto il file di restart, il distruttore del DataManager
    * crea il file e forza gli oggetti a scrivere il loro contributo nel modo
    * opportuno */
   if(RestartEvery == ATEND) {
      MakeRestart();
   }
   
   if(sSimulationTitle != NULL) {
      SAFEDELETEARR(sSimulationTitle);
   }
   
   ElemManagerDestructor();
   NodeManagerDestructor();
   DofManagerDestructor();

#if defined(HAVE_LOADABLE) && defined(HAVE_LTDL_H)
   if (loadableElemInitialized) {
      if (lt_dlexit()) {
	 std::cerr << "lt_dlexit failed" << std::endl;
	 THROW(ErrGeneric());
      }
   }
#endif /* HAVE_LOADABLE && HAVE_LTDL_H */
} /* End of DataManager::DataManager() */


void DataManager::MakeRestart(void) 
{
  if (RestartEvery == ATEND) {	
   silent_cout("Making restart file ..." << std::endl);
   
   /* Inizializzazione del file di restart */
   OutHdl.RestartOpen();
  
   time_t tCurrTime(time(NULL));
   OutHdl.Restart() << "# Restart file prepared by Mbdyn, " 
     << ctime(&tCurrTime) << std::endl << std::endl;
   
   /* Dati iniziali */
   OutHdl.Restart() << "begin: data;" << std::endl
     << "# uncomment this line to use the default integrator" << std::endl
     << "  # integrator: multistep;" << std::endl
     << "end: data;" << std::endl << std::endl
     << "# the following block contains data for the multistep integrator" 
     << std::endl
     << "begin: multistep;" << std::endl
     << "# add data for the multistep integrator:" << std::endl
     << "  # initial time: " << pTime->GetVal().GetReal() << ';' << std::endl
     << "  # final time: " << pTime->GetVal().GetReal() << ';' << std::endl
     << "  # time step: 1.;" << std::endl
     << "  # method: ms, .0, .0;" << std::endl
     << "  # max iterations: 1;" << std::endl
     << "  # tolerance: 1.e-6;" << std::endl
     << "  # derivatives max iterations: 1;" << std::endl
     << "  # derivatives tolerance: 1.e-6;" << std::endl
     << "  # derivatives coefficient: 1.e-6;" << std::endl
     << "  # fictitious steps max iterations: 1;" << std::endl
     << "  # fictitious steps tolerance: 1.e-6;" << std::endl
     << "  # fictitious steps number: 2;" << std::endl
     << "  # fictitious steps ratio: 1.e-3;" << std::endl
     << "  # Newton Raphson: true;" << std::endl
     << "end: multistep;" << std::endl << std::endl;
   
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
      if (ElemData[iCnt].iNum > 0) {
	 if (ElemData[iCnt].fIsUnique == 1) {
	    OutHdl.Restart() << "  " << psReadControlElems[iCnt] 
	      << ';' << std::endl;
	 } else {		     		     
	    OutHdl.Restart() << "  " << psReadControlElems[iCnt] << ": "
	      << ElemData[iCnt].iNum << ';' << std::endl;
	 }		  
      }
   }	
   
   if (sSimulationTitle != NULL) {
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
     << "  # initial tolerance: " << dInitialAssemblyToll << ';' << std::endl
     << "  # max initial iterations: " << iMaxInitialIterations 
     << ';' << std::endl;
#endif // USE_STRUCT_NODES     
   OutHdl.Restart() << "# uncomment this line if restart file is to be made again" 
     << std::endl
     << "  # make restart file;" << std::endl
     << "# remember: it will replace the present file if you don't change its name"
     << std::endl;
   
   OutHdl.Restart() << "end: control data;" << std::endl << std::endl;
   
   /* Dati dei nodi */
   OutHdl.Restart() << "begin: nodes;" << std::endl;
   for (Node** ppTmpNode = ppNodes;
	ppTmpNode < ppNodes+iTotNodes;
	ppTmpNode++) {
      (*ppTmpNode)->Restart(OutHdl.Restart());
   }	
   OutHdl.Restart() << "end: nodes;" << std::endl << std::endl;
   
   /* Dati dei driver */
   if (iTotDrive > 0) {	     
      OutHdl.Restart() << "begin: drivers;" << std::endl;
      for (Drive** ppTmpDrv = ppDrive;
	   ppTmpDrv < ppDrive+iTotDrive;
	   ppTmpDrv++) {
	 OutHdl.Restart() 
	   << "  # file driver " << (*ppTmpDrv)->GetLabel()
	     << " is required" << std::endl;
      }	
      OutHdl.Restart() << "end: drivers;" << std::endl << std::endl;
   }	
   
   /* Dati degli elementi */
   OutHdl.Restart() << "begin: elements;" << std::endl;
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppLastElem; ppTmpEl++) {
      ASSERT(*ppTmpEl != NULL);
      (*ppTmpEl)->Restart(OutHdl.Restart());
   }
   OutHdl.Restart() << "end: elements;" << std::endl << std::endl << std::endl << std::endl;
  }
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

