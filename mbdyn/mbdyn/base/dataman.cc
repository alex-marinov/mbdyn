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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

const char sDMClassName[] = "DataManager";

/* linka i singoli DriveCaller al DriveHandler posseduto dal DataManager */
extern void SetDrvHdl(DriveHandler*);

#ifdef DEBUG_MEMMANAGER
clMemMan DMmm("DataManager");
#endif

const flag fDefaultInitialJointAssemblyToBeMade = flag(0);
const flag fDefaultSkipInitialJointAssembly = flag(0);
const doublereal dDefaultInitialStiffness = 1.;
const flag fDefaultOmegaRotates = flag(0);
const doublereal dDefaultInitialAssemblyToll = 1.e-6;
const integer iDefaultMaxInitialIterations = 1;

const char sDefaultOutputFileName[] = "Mbdyn";



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
fInitialJointAssemblyToBeMade(fDefaultInitialJointAssemblyToBeMade),
fSkipInitialJointAssembly(fDefaultSkipInitialJointAssembly),
dInitialPositionStiffness(dDefaultInitialStiffness), 
dInitialVelocityStiffness(dDefaultInitialStiffness),
fOmegaRotates(fDefaultOmegaRotates),
dInitialAssemblyToll(dDefaultInitialAssemblyToll),
iMaxInitialIterations(iDefaultMaxInitialIterations),
#endif /* USE_STRUCT_NODES */
fPrintDofStats(0),
sSimulationTitle(NULL),
RestartEvery(NEVER),
iRestartIterations(0),
dRestartTime(0.),
iCurrRestartIter(0),
dLastRestartTime(dInitialTime),
fAdamsResOutput(0),
sAdamsModelName(NULL),
iAdamsOutputBlock(1),
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
      cerr << "error in setting Time symbol" << endl;

      THROW(DataManager::ErrGeneric());
   }   
   SetTime(dInitialTime);
   
   DEBUGLCOUT(MYDEBUG_INIT, "GlobalSymbolTable:" 
	      << endl << GlobalSymbolTable << endl);
   
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
      cerr << "<begin> expected at line "
	<< HP.GetLineData() << endl;

      THROW(DataManager::ErrGeneric());
   }
   if(KeyWords(HP.GetWord()) != CONTROLDATA) {
      DEBUGCERR("");
      cerr << "<begin: control data;> expected at line " 
	<< HP.GetLineData() << endl;

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
	 cerr << "<begin> expected at line " 
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != NODES) {
	 DEBUGCERR("");
	 cerr << "<begin: nodes;> expected at line " 
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
      
      ReadNodes(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("");
      cerr << "warning, no nodes are defined" << endl;
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
	 cerr << "<begin> expected at line "
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != DRIVERS) {
	 DEBUGCERR("");
	 cerr << "<begin: drivers;> expected at line " 
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
      
      ReadDrivers(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("warning, no drivers are defined" << endl);
   }   
	
   /* fine lettura drivers */
   
   
   /*
    * legge gli elementi, crea la struttura ppElems
    * e contemporaneamente aggiorna i dof
    */
   if(iTotElem > 0) {	
      if(KeyWords(HP.GetDescription()) != BEGIN) {
	 DEBUGCERR("");
	 cerr << "<begin> expected at line " 
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
      if(KeyWords(HP.GetWord()) != ELEMENTS) {
	 DEBUGCERR("");
	 cerr << "<begin: elements;> expected at line " 
	   << HP.GetLineData() << endl;

	 THROW(DataManager::ErrGeneric());
      }
	
      ReadElems(HP);
      HP.PutKeyTable(K);
   } else {
      DEBUGCERR("");
      cerr << "warning, no elements are defined" << endl;
   }   
	
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
	<< endl;
	  
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
	cout << "DofType " << i << " (" << psDofOwnerNames[i] 
	<< "), n. of owners: " << DofData[i].iNum << endl;
      
      /* mostra in modo succinto il numero di nodi per tipo */
      for(int i = 0; i < Node::LASTNODETYPE; i++)
	cout << "NodeType " << i << " (" << psNodeNames[i] 
	<< "), n. of nodes: " << NodeData[i].iNum << endl;
      
      /* mostra in modo succinto il numero di elementi per tipo */
      for(int i = 0; i < Elem::LASTELEMTYPE; i++)
	cout << "Element Type " << i << " (" << psElemNames[i] 
	<< "), n. of elems: " << ElemData[i].iNum << endl;
      
      /* mostra in modo succinto il numero di drivers per tipo */
      for(int i = 0; i < Drive::LASTDRIVETYPE; i++)
	cout << "DriveType " << i << " (" << psDriveNames[i] 
	<< "), n. of drivers: " << DriveData[i].iNum << endl;
   }
#endif /* DEBUG */

   if(fAbortAfterInput) {
      silent_cout("Only input is required" << endl);
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
   if (fInitialJointAssemblyToBeMade) {
      if (!fSkipInitialJointAssembly) {
	 InitialJointAssembly();
      }	else {
	 silent_cout("Skipping initial joints assembly" << endl);
      }
   } else {
      silent_cout("No initial assembly is required since no joints are defined"
		  << endl);
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
	 cout << "DofType " << i << ':' << endl;
	 for(int j = 0; j < DofData[i].iNum; j++) {
	    cout << "DofOwner " << j << ", n. of dofs: "
	      << pDofOwners[k++].iNumDofs << endl;
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
	 cout << "Dof " << setw(4) << iCnt+1 << ": order "
	   << (pDofs+iCnt)->Order << endl;
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
} /* End of DataManager::DataManager() */


void DataManager::MakeRestart(void) 
{
  if (RestartEvery == ATEND) {	
   silent_cout("Making restart file ..." << endl);
   
   /* Inizializzazione del file di restart */
   OutHdl.RestartOpen();
  
   time_t tCurrTime(time(NULL));
   OutHdl.Restart() << "# Restart file prepared by Mbdyn, " 
     << ctime(&tCurrTime) << endl << endl;
   
   /* Dati iniziali */
   OutHdl.Restart() << "begin: data;" << endl
     << "# uncomment this line to use the default integrator" << endl
     << "  # integrator: multistep;" << endl
     << "end: data;" << endl << endl
     << "# the following block contains data for the multistep integrator" 
     << endl
     << "begin: multistep;" << endl
     << "# add data for the multistep integrator:" << endl
     << "  # initial time: " << pTime->GetVal().GetReal() << ';' << endl
     << "  # final time: " << pTime->GetVal().GetReal() << ';' << endl
     << "  # time step: 1.;" << endl
     << "  # method: ms, .0, .0;" << endl
     << "  # max iterations: 1;" << endl
     << "  # tolerance: 1.e-6;" << endl
     << "  # derivatives max iterations: 1;" << endl
     << "  # derivatives tolerance: 1.e-6;" << endl
     << "  # derivatives coefficient: 1.e-6;" << endl
     << "  # fictitious steps max iterations: 1;" << endl
     << "  # fictitious steps tolerance: 1.e-6;" << endl
     << "  # fictitious steps number: 2;" << endl
     << "  # fictitious steps ratio: 1.e-3;" << endl
     << "  # Newton Raphson: true;" << endl
     << "end: multistep;" << endl << endl;
   
   /* Dati di controllo */	
   OutHdl.Restart() << "begin: control data;" << endl;
   
   /* Nodi */
   for (int iCnt = 0; iCnt < Node::LASTNODETYPE; iCnt++) {	     
      if (NodeData[iCnt].iNum > 0) {		  
	 OutHdl.Restart() << "  " << psReadControlNodes[iCnt] << ": "
	   << NodeData[iCnt].iNum << ';' << endl;
      }
   }
   
   /* Drivers */
   for (int iCnt = 0; iCnt < Drive::LASTDRIVETYPE; iCnt++) {	     
      if (DriveData[iCnt].iNum > 0) {		  
	 OutHdl.Restart() << "  " 
	   << psReadControlDrivers[iCnt] << ": "
	   << DriveData[iCnt].iNum << ';' << endl;
      }
   }
   
   /* Elementi */
   for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {	     
      if (ElemData[iCnt].iNum > 0) {
	 if (ElemData[iCnt].fIsUnique == 1) {
	    OutHdl.Restart() << "  " << psReadControlElems[iCnt] 
	      << ';' << endl;
	 } else {		     		     
	    OutHdl.Restart() << "  " << psReadControlElems[iCnt] << ": "
	      << ElemData[iCnt].iNum << ';' << endl;
	 }		  
      }
   }	
   
   if (sSimulationTitle != NULL) {
      OutHdl.Restart() << "  title: \"" 
	<< sSimulationTitle << "\";" << endl;
   }

#if defined(USE_STRUCT_NODES)   
   OutHdl.Restart() << endl
     << "# comment this line if the model is to be modified!" << endl
     << "  skip initial joint assembly;" << endl
     << "# uncomment the following lines to improve the satisfaction of constraints"
     << endl
     << "  # initial stiffness: " << dInitialPositionStiffness << ", "
     << dInitialVelocityStiffness << ';' << endl
     << "  # initial tolerance: " << dInitialAssemblyToll << ';' << endl
     << "  # max initial iterations: " << iMaxInitialIterations 
     << ';' << endl;
#endif // USE_STRUCT_NODES     
   OutHdl.Restart() << "# uncomment this line if restart file is to be made again" 
     << endl
     << "  # make restart file;" << endl
     << "# remember: it will replace the present file if you don't change its name"
     << endl;
   
   OutHdl.Restart() << "end: control data;" << endl << endl;
   
   /* Dati dei nodi */
   OutHdl.Restart() << "begin: nodes;" << endl;
   for (Node** ppTmpNode = ppNodes;
	ppTmpNode < ppNodes+iTotNodes;
	ppTmpNode++) {
      (*ppTmpNode)->Restart(OutHdl.Restart());
   }	
   OutHdl.Restart() << "end: nodes;" << endl << endl;
   
   /* Dati dei driver */
   if (iTotDrive > 0) {	     
      OutHdl.Restart() << "begin: drivers;" << endl;
      for (Drive** ppTmpDrv = ppDrive;
	   ppTmpDrv < ppDrive+iTotDrive;
	   ppTmpDrv++) {
	 OutHdl.Restart() 
	   << "  # file driver " << (*ppTmpDrv)->GetLabel()
	     << " is required" << endl;
      }	
      OutHdl.Restart() << "end: drivers;" << endl << endl;
   }	
   
   /* Dati degli elementi */
   OutHdl.Restart() << "begin: elements;" << endl;
   Elem** ppLastElem = ppElems+iTotElem;
   for (Elem** ppTmpEl = ppElems; ppTmpEl < ppLastElem; ppTmpEl++) {
      ASSERT(*ppTmpEl != NULL);
      (*ppTmpEl)->Restart(OutHdl.Restart());
   }
   OutHdl.Restart() << "end: elements;" << endl << endl << endl << endl;
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

