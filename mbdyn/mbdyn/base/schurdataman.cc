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

/* Schur Data Manager */

/*
 * Copyright 1999-2014 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

/* libraries for mesh partitioning computation */
#include "metiswrap.h"
#include "chacowrap.h"

#include "schurdataman.h"
#include "mbcomm.h"
#include "mysleep.h"
#include "except.h"
#include "solver.h"
#undef min
#undef max
#include <vector>
#include <algorithm>

#include "rotor.h"

/* struttura contenente le distribuzioni di dofs */
struct Adjacency{
	int *pXadj;	/* x ogni nodo: pos. su pAdjncy lista connessioni */
	int *pAdjncy;	/* lista liste connessioni */
};

/* !!!!!!! Costanti da verificare !!!!!*/
const integer iDefaultMaxConnectionsPerVertex = 100;
const integer iDefaultMaxNodesPerElem = 50;
const integer iDefaultMaxInterfNodes = 50;
const int ADJ_UNDEFINED = -1;

#ifndef USE_MPI
SchurDataManager::SchurDataManager(MBDynParser& HP,
		unsigned OF,
		Solver* pS,
		doublereal dInitialTime,
		const char* sOutputFileName,
		const char* sInputFileName,
		bool bAbortAfterInput)
: DataManager(HP, OF, pS, dInitialTime, sOutputFileName, sInputFileName, bAbortAfterInput)
{
	silent_cerr("fatal error: you are building SchurDataManager, "
		"but mbdyn was compiled without MPI. "
		"Something weird is happening. "
		"Anyway, please compile with -DUSE_MPI "
		"to enable parallel solution" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

SchurDataManager::~SchurDataManager()
{
	NO_OP;
}

void
SchurDataManager::AssRes(VectorHandler&, doublereal) throw(ChangedEquationStructure)
{
	NO_OP;
}

void
SchurDataManager::AssJac(MatrixHandler&, doublereal)
{
	NO_OP;
}

void
SchurDataManager::DerivativesUpdate(void) const
{
	NO_OP;
}

void
SchurDataManager::BeforePredict(VectorHandler&, VectorHandler&,
		VectorHandler&, VectorHandler&) const
{
	NO_OP;
}

void
SchurDataManager::AfterPredict(void) const
{
	NO_OP;
}

void
SchurDataManager::Update(void) const
{
	NO_OP;
}

void
SchurDataManager::AfterConvergence(void) const
{
	NO_OP;
}

integer
SchurDataManager::HowManyDofs(DofType) const
{
	return 0;
}

integer *
SchurDataManager::GetDofsList(DofType) const
{
	return 0;
}

Dof *
SchurDataManager::pGetDofsList(void) const
{
	return 0;
}

void
SchurDataManager::Output(bool force) const
{
	NO_OP;
}
#else /* USE_MPI */

/* NOTE: define to use Wait/Waitall instead of Test/Testall
 * Apparently, this results in far better performances,
 * so we might want to extend it to all other communications */
#define USE_MPI_WAIT

/* Costruttore - begin */
SchurDataManager::SchurDataManager(MBDynParser& HP,
		unsigned OF,
		Solver *pS,
		doublereal dInitialTime,
		const char* sOutputFileName,
		const char* sInputFileName,
		bool bAbortAfterInput)
:DataManager(HP, OF, pS, dInitialTime, sOutputFileName, sInputFileName, bAbortAfterInput),
iTotVertices(0),
ppMyElems(NULL),
iNumLocElems(0),
ppMyIntElems(NULL),
iNumIntElems(0),
ppMyNodes(NULL),
iNumLocNodes(0),
pLocalDofs(NULL),
iNumLocDofs(0),
pLocalIntDofs(NULL),
iNumIntDofs(0),
ppIntNodes(NULL),
iNumIntNodes(0),
iNumMyInt(0),
pMyIntDofs(NULL),
pLabelsList(NULL),
wgtflag(WEIGHT_VERTICES),
pParAmgProcs(NULL),
Partitioner(PARTITIONER_DEFAULT),
pIndVelComm(NULL),
ppExpCntNodes(NULL),
ppExpCntElems(NULL),
iTotalExpConnections(0)
{
	DEBUGCOUT("Entering SchurDataManager" << std::endl);

	/* Inizializza il communicator */
	DataComm = MBDynComm.Dup();
	DataCommSize = DataComm.Get_size();
	MyRank = DataComm.Get_rank();

	DEBUGCOUT("Communicator Size: " << DataCommSize << std::endl);

	iTotVertices = iTotNodes + Elems.size();
	DEBUGCOUT("iTotVertices: " << iTotVertices << std::endl);

	/* parole chiave del blocco parallelizzazione */
	const char* sKeyWords[] = {
		"begin",
		"parallel",
		"connection",
		"number" "of" "connections",
		"element",
		"node",
		"force",
		"rigid" "body",
		"joint",
		"beam",
		"plate",
		"rotor",
		"aeromodal",
		"aerodynamic" "element",
		"electric" "bulk",
		"electric",
		"thermal",
		"genel",
		"hydraulic",
		"bulk",
		"loadable",
		"driven",
		"abstract",
		"structural",
		"parameter",
		"weights",
		"partition",
		"partitioner",
		"end",
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		BEGIN = 0,
		PARALLEL,
		CONNECTION,
		NUMBEROFCONNECTIONS,
		ELEMENT,
		NODE,
		FORCE,
		RIGIDBODY,
		JOINT,
		BEAM,
		PLATE,
		ROTOR,
		AEROMODAL,
		AERODYNAMICELEMENT,
		ELECTRICBULK,
		ELECTRIC,
		THERMAL,
		GENEL,
		HYDRAULIC,
		BULK,
		LOADABLE,
		DRIVEN,
		ABSTRACT,
		STRUCTURAL,
		PARAMETER,
		WEIGHTS,
		PARTITION,
		PARTITIONER,
		END,
		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* legge la distribuzione degli elementi sulle diverse CPU */

	try {
		if (KeyWords(HP.GetDescription()) != BEGIN) {
			pedantic_cerr("no explicit connections declared "
				"for this input file" << std::endl);
			return;
		}

	} catch (EndOfFile) {
		pedantic_cerr("no explicit connections declared "
			"for this input file" << std::endl);
		return;
	}

	int iNumElems = 0;
	int iNumNodes = 0;
	if (KeyWords(HP.GetWord()) != PARALLEL) {
		silent_cerr("Error: \"begin: parallel;\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	while (true) {
		switch (KeyWords(HP.GetDescription())) {
		case WEIGHTS:
			if (!HP.IsArg()) {
				silent_cerr("Error: Weight flag expected "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			while (HP.IsArg()) {
				if (HP.IsKeyWord("none")) {
					wgtflag = WEIGHT_NONE;

				} else if (HP.IsKeyWord("vertices")
						|| HP.IsKeyWord("computation"))
				{
					wgtflag |= WEIGHT_VERTICES;

				} else if (HP.IsKeyWord("communication")) {
					wgtflag |= WEIGHT_COMM;

				} else if (HP.IsKeyWord("edges")) {
					wgtflag |= WEIGHT_EDGES;

				} else {
					silent_cerr("invalid weight "
						" at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
			break;

		case PARTITION:
			SAFENEWARR(pParAmgProcs, int, iTotVertices);
			for (int i = 0; i < iTotVertices; i++) {
				if (!HP.IsArg()) {
					silent_cerr("the partition "
						"assignment is not complete; "
						"only " << i << " vertices "
						"input so far "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				pParAmgProcs[i] = HP.GetInt();
				if (pParAmgProcs[i] < 0 ||
						pParAmgProcs[i] >= DataCommSize)
				{
					silent_cerr("illegal value "
						<< pParAmgProcs[i]
						<< " for partition "
						"assignment[" << i << "] "
						"at line " << HP.GetLineData()
						<< "; must be between 0 and "
						<< DataCommSize - 1 << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
			break;

		case PARTITIONER:
			if (HP.IsKeyWord("metis")) {
#ifdef USE_METIS
				Partitioner = PARTITIONER_METIS;
#else /* ! USE_METIS */
				silent_cerr("METIS partitioner not available; aborting..." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_METIS */

			} else if (HP.IsKeyWord("chaco")) {
#ifdef USE_CHACO
				Partitioner = PARTITIONER_CHACO;
#else /* ! USE_CHACO */
				silent_cerr("CHACO partitioner not available; aborting..." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! USE_CHACO */
			} else if (HP.IsKeyWord("manual")) {
				Partitioner = PARTITIONER_MANUAL;

			} else {
				silent_cerr("unknown partitioner "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case END:
			if (KeyWords(HP.GetWord()) != PARALLEL) {
				silent_cerr("Error: \"end: parallel;\" expected "
					"at line " << HP.GetLineData()
					<< "; aborting..." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			goto endcycle;

		default:
			silent_cerr("Unknown input at line "
				<< HP.GetLineData() << "; aborting..." 
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		case NUMBEROFCONNECTIONS:
			iTotalExpConnections = HP.GetInt();
			SAFENEWARR(ppExpCntNodes, Node*, iTotalExpConnections);
			SAFENEWARR(ppExpCntElems, Elem*, iTotalExpConnections);

			Elem::Type CurrElType;
			Node::Type CurrNdType;
			unsigned int j = 0;

			for (int i = 0; i < iTotalExpConnections; i++) {
				if (KeyWords(HP.GetDescription()) != CONNECTION) {
					silent_cerr("Error: <Connection> "
						"expected at line "
						<< HP.GetLineData() 
						<< "; aborting..." 
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				for (int k = 0; k < 2; k++) {
					switch (KeyWords(HP.GetWord())) {
					case ELEMENT:
						switch (KeyWords(HP.GetWord())) {
						case FORCE:
							CurrElType = Elem::FORCE;
							break;

						case RIGIDBODY:
							CurrElType = Elem::BODY;
							break;

						case JOINT:
							CurrElType = Elem::JOINT;
							break;

						case BEAM:
							CurrElType = Elem::BEAM;
							break;

						case PLATE:
							CurrElType = Elem::PLATE;
							break;

						case ROTOR:
							CurrElType = Elem::INDUCEDVELOCITY;
							break;

						case AEROMODAL:
							CurrElType = Elem::AEROMODAL;
							break;

						case AERODYNAMICELEMENT:
							CurrElType = Elem::AERODYNAMIC;
							break;

						case ELECTRICBULK:
							CurrElType = Elem::ELECTRICBULK;
							break;

						case ELECTRIC:
							CurrElType = Elem::ELECTRIC;
							break;

						case THERMAL:
							CurrElType = Elem::THERMAL;
							break;

						case GENEL:
							CurrElType = Elem::GENEL;
							break;

						case HYDRAULIC:
							CurrElType = Elem::HYDRAULIC;
							break;

						case BULK:
							CurrElType = Elem::BULK;
							break;

						case LOADABLE:
							CurrElType = Elem::LOADABLE;
							break;

						case DRIVEN:
							CurrElType = Elem::DRIVEN;
							break;

						default:
							silent_cerr("Error: invalid element type "
								"at line " << HP.GetLineData()
								<< "; aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (HP.IsArg()) {
							j = HP.GetInt();

						} else {
							silent_cerr("Error: label expected "
								"for " << psElemNames[CurrElType]
								<< " element type at line "
								<< HP.GetLineData()
								<< "; aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						ppExpCntElems[iNumElems] =
							*(ppFindElem(CurrElType, j));
						if (ppExpCntElems[iNumElems] == NULL) {
							silent_cerr("Error at line "
								<< HP.GetLineData() << ": "
								<< psElemNames[CurrElType]
								<< "(" << j << ") undefined; "
								"aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						iNumElems++;
						break;

					case NODE:
						switch (KeyWords(HP.GetWord())) {
						case ABSTRACT:
							CurrNdType = Node::ABSTRACT;
							break;

						case STRUCTURAL:
							CurrNdType = Node::STRUCTURAL;
							break;

						case ELECTRIC:
							CurrNdType = Node::ELECTRIC;
							break;

						case THERMAL:
							CurrNdType = Node::THERMAL;
							break;

						case PARAMETER:
							CurrNdType = Node::PARAMETER;
							break;

						case HYDRAULIC:
							CurrNdType = Node::HYDRAULIC;
							break;

						default:
							silent_cerr("Error: invalid node type "
								"at line " << HP.GetLineData()
								<< "; aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						if (HP.IsArg()) {
							j = HP.GetInt();

						} else {
							silent_cerr("Error: label expected "
								"at line " << HP.GetLineData()
								<< "; aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}

						ppExpCntNodes[iNumNodes] =
							pFindNode(CurrNdType, j);
						if (ppExpCntNodes[iNumNodes] == NULL ) {
							silent_cerr("Error: at line "
								<< HP.GetLineData() << ":"
								<< psNodeNames[CurrNdType]
								<< "(" << j << ") undefined; "
								<< "aborting..." << std::endl);
							throw ErrGeneric(MBDYN_EXCEPT_ARGS);
						}
						iNumNodes++;
						break;

					default:
						silent_cerr("Unknown input at line "
							<< HP.GetLineData()
							<< "; aborting..." << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}
			}

			ASSERT(iNumNodes == iNumElems);
			if (iNumNodes + iNumElems !=  2*iTotalExpConnections) {
				silent_cerr("Error: total number of nodes and elements"
					" in the parallel section at line "
					<< HP.GetLineData()
					<< " is not consistent; aborting..." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;
		}
	}
endcycle:;
}
/* Costruttore - End */

/* Distruttore - begin */
SchurDataManager::~SchurDataManager(void)
{
	DEBUGCOUTFNAME("SchurDataManager::~SchurDataManager");

	if (ppMyNodes != NULL) {
		SAFEDELETEARR(ppMyNodes);
	}

	if (ppMyElems != NULL) {
		SAFEDELETEARR(ppMyElems);
	}

	if (ppMyIntElems != NULL) {
		SAFEDELETEARR(ppMyIntElems);
	}

	if (ppIntNodes != NULL) {
		SAFEDELETEARR(ppIntNodes);
	}

	if (pLocalDofs != NULL) {
		SAFEDELETEARR(pLocalDofs);
	}

	if (pLocalIntDofs != NULL) {
		SAFEDELETEARR(pLocalIntDofs);
	}

	if (ppExpCntNodes != NULL) {
		SAFEDELETEARR(ppExpCntNodes);
	}

	if (ppExpCntElems != NULL) {
		SAFEDELETEARR(ppExpCntElems);
	}
}

integer
SchurDataManager::HowManyDofs(DofType who) const
{
	switch(who) {
	case TOTAL:
		return iTotDofs;

	case LOCAL:
		return iNumLocDofs;

	case INTERNAL:
		return iNumIntDofs;

	case MYINTERNAL:
		return iNumMyInt;

	default:
		silent_cerr("SchurDataManager::HowManyDofs: "
			"illegal request (" << unsigned(who) << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

integer *
SchurDataManager::GetDofsList(DofType who) const
{
	switch(who) {
	case LOCAL:
		return pLocalDofs;

	case INTERNAL:
		return pLocalIntDofs;

	case MYINTERNAL:
		return pMyIntDofs;

	default:
		silent_cerr("SchurDataManager::GetDofsList: "
			"illegal request (" << unsigned(who) << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Dof*
SchurDataManager::pGetDofsList(void) const
{
	return pDofs;
}

/* Ripartisce i nodi fra i processi e all'interno di ogni singolo processo */
void
SchurDataManager::CreatePartition(void)
{
	DEBUGCOUT("Entering SchurDataManager::CreatePartition()" << std::endl);

	std::vector<std::vector<int> > adj(iTotVertices);	/* adjacency */
	std::vector<int> CommWgts(iTotVertices);		/* communication weights */

	Adjacency Vertices;	/* Struttura contenente le connessioni fra i vertici */
	std::vector<int> VertexWgts(iTotVertices);	/* Pesi dei vertici = dofs x ogni v. utile per METIS */
	Vertices.pXadj = 0;
	Vertices.pAdjncy = 0;
	integer iMax = 0;
	integer iRMax = 0;
	int iCount = 0;
	int iNumRt = 0;

	ASSERT(iTotVertices > 0);
	ASSERT(DataCommSize > 0);
	
	/* Costruisco e inizializzo le due strutture */
	SAFENEWARR(Vertices.pXadj, int, iTotVertices + 1);
	memset(Vertices.pXadj, 0, (iTotVertices + 1)*sizeof(int));

	/* Ciclo per la scrittura degli array delle connessioni.
	 * Il ciclo viene ripetuto se per un vertice si ha un numero
	 * di connessioni superiore al max consentito per default
	 * iDefaultMaxConnectionsPerVertex */
	int iMaxConnectionsPerVertex =
		(iTotVertices < iDefaultMaxConnectionsPerVertex) ? iTotVertices
		: iDefaultMaxConnectionsPerVertex;
	int GravityPos = 0, AirPropPos = 0;
	int* pRotPos = NULL;
	unsigned int* pRotLab = NULL;
	if (!ElemData[Elem::INDUCEDVELOCITY].ElemContainer.empty()) {
		unsigned iNum = ElemData[Elem::INDUCEDVELOCITY].ElemContainer.size();
		SAFENEWARR(pRotPos, int, iNum);
		SAFENEWARR(pRotLab, unsigned int, iNum);
		memset(pRotPos, 0, iNum*sizeof(int));
		memset(pRotLab, 0, iNum*sizeof(unsigned int));
	}
	SAFENEWARR(pLabelsList, unsigned int, iTotNodes);
	memset(pLabelsList, 0, iTotNodes*sizeof(unsigned int));

	std::vector<const Node *> connectedNodes;

	while (true) {
		InitList(Vertices.pXadj, iTotVertices + 1, 0);

		SAFENEWARR(Vertices.pAdjncy, int, iTotVertices*iMaxConnectionsPerVertex);
		InitList(Vertices.pAdjncy, iTotVertices*iMaxConnectionsPerVertex, ADJ_UNDEFINED);
		ASSERT(Elems.begin() != Elems.end());

		for (unsigned int i = 0; i < iTotNodes; i++) {
			pLabelsList[i] = ppNodes[i]->GetLabel();
		}

		/* ciclo sugli elementi per assemblare la struttura delle connessioni */
		Node** ppCurrNode = NULL;
		/* per numerare i nodi prendo la posizione del puntatore
		 * al nodo nell'array ppNodes */
		int position;
		iCount = iTotNodes;
		iNumRt = 0;

		for (ElemVecType::const_iterator pTmpEl = Elems.begin();
			pTmpEl != Elems.end();
			pTmpEl++, iCount++)
		{
			if ((*pTmpEl)->GetElemType() == Elem::GRAVITY) {
				GravityPos = iCount - iTotNodes;

			} else if ((*pTmpEl)->GetElemType() == Elem::AIRPROPERTIES) {
				AirPropPos = iCount - iTotNodes;

			} else if ((*pTmpEl)->GetElemType() == Elem::INDUCEDVELOCITY) {
				pRotPos[iNumRt] = iCount - iTotNodes;
				pRotLab[iNumRt] = (*pTmpEl)->GetLabel();
				iNumRt++;
			}

			(*pTmpEl)->GetConnectedNodes(connectedNodes);

			/* peso dell'elemento */
			integer dimA, dimB;
			(*pTmpEl)->WorkSpaceDim(&dimA, &dimB);
			VertexWgts[iCount] = dimA * dimB;

			CommWgts[iCount] = (*pTmpEl)->iGetNumDof()*(*pTmpEl)->iGetNumDof();

			for (std::vector<const Node *>::const_iterator i = connectedNodes.begin();
				i != connectedNodes.end();
				i++)
			{
				Vertices.pXadj[iCount + 1]++;

				/* trovo la pos. del nodo nella lista dei puntatori ai nodi */
				Node::Type type = (*i)->GetNodeType();
				unsigned label = (*i)->GetLabel();
				ppCurrNode = SearchNode(NodeData[type].ppFirstNode,
						NodeData[type].iNum, label);
				position = ppCurrNode - ppNodes;
				
				/* Aggiungo al peso dell'elemento il numero di dofs
				 * di ciascun nodo connesso */

#if 0   /* FIXME: i nodi hanno peso comp. nullo */
				VertexWgts[iCount] += (*ppCurrNode)->iGetNumDof();
#endif

				/* aggiungo fra le connessioni dell'elemento il nodo attuale */
				if ((iCount*iMaxConnectionsPerVertex) + Vertices.pXadj[iCount + 1] - 1 < iTotVertices*iMaxConnectionsPerVertex) {
					Vertices.pAdjncy[(iCount*iMaxConnectionsPerVertex)
						+ Vertices.pXadj[iCount + 1] - 1] = position;
				}

				/* aggiungo alle connessioni del nodo l'elemento attuale */
				Vertices.pXadj[position + 1]++;
				if ((position*iMaxConnectionsPerVertex) + Vertices.pXadj[position + 1] - 1 < iTotVertices*iMaxConnectionsPerVertex) {
					Vertices.pAdjncy[(position*iMaxConnectionsPerVertex)
						+ Vertices.pXadj[position + 1] - 1] = iCount;
				}

				/* peso (di comunicazione) del nodo */
				CommWgts[position] = (*ppCurrNode)->iGetNumDof();
			}
		}

		if (iTotalExpConnections != 0) {
			for (int i = 0; i < iTotalExpConnections; i++) {
				int iNdPos, iElPos;
				int j;

				for (j = 0; ppExpCntNodes[i] != ppNodes[j]; j++) {
					NO_OP;
				}
				iNdPos = j;

				for (j = 0; ppExpCntElems[i] != Elems[j]; j++) {
					NO_OP;
				}
				iElPos = j;

				Vertices.pXadj[iNdPos + 1]++;
				Vertices.pXadj[iTotNodes + iElPos + 1]++;
				Vertices.pAdjncy[(iNdPos*iMaxConnectionsPerVertex)
					+ Vertices.pXadj[iNdPos + 1] - 1] = iElPos + iTotNodes;
				Vertices.pAdjncy[((iTotNodes + iElPos)*iMaxConnectionsPerVertex)
					+ Vertices.pXadj[iTotNodes + iElPos + 1] - 1] = iNdPos;
			}
		}

		iMax = 0;
		for (int i = 0; i < iTotVertices; i++) {
			if (Vertices.pXadj[i] > iMaxConnectionsPerVertex) {
				iMax = Vertices.pXadj[i];
			}
		}

		if (iMax <= iMaxConnectionsPerVertex) {
			break;
		}

		iMaxConnectionsPerVertex = iMax;
		SAFEDELETEARR(Vertices.pAdjncy);
		Vertices.pAdjncy = 0;
	}

	for (int i = 1; i <= iTotVertices; i++) {
		Vertices.pXadj[i] += Vertices.pXadj[i - 1];
	}

	/* Compatta il vettore delle adiacenze */
	Pack(Vertices.pAdjncy, iTotVertices*iMaxConnectionsPerVertex);

	/* Chiamo la routine per ottere la partizione fra i diversi processi METIS.
	 * Se ne usano due diverse a seconda della dimensione della partizione */
	if (pParAmgProcs == NULL) {
		SAFENEWARR(pParAmgProcs, int, iTotVertices);

#ifdef DEBUG
		if (MyRank == 0) {
			std::ofstream ofPartition;

			ofPartition.open("partition.debug");
			ofPartition << "# METIS-like Input File" << std::endl
				<< "# Column 1: computational weight" << std::endl
				<< "# Column 2: communication weight" << std::endl
				<< "# Total Vertices: " << iTotVertices << std::endl
				<< "# Nodes" << std::endl;
			for (unsigned int i = 0; i < iTotNodes; i++) {
				ofPartition << "# " << i << "  Node Type: "
					<< "(" << psNodeNames[ppNodes[i]->GetNodeType()] << ")"
					<< " Label: " << ppNodes[i]->GetLabel() << std::endl
					<< VertexWgts[i] << " " << CommWgts[i];
				for (int j = Vertices.pXadj[i]; j < Vertices.pXadj[i + 1]; j++) {
					ofPartition << " " << Vertices.pAdjncy[j];
				}
				ofPartition << std::endl;
			}
			ofPartition << "# Elements" << std::endl;
			for (unsigned int i = 0; i < iTotElem; i++) {
				ofPartition << "# " << i + iTotNodes << "  Element Type: "
					<< "("  << psElemNames[Elems[i]->GetElemType()] << ")"
					<< " Label: " << Elems[i]->GetLabel() << std::endl
					<< VertexWgts[i + iTotNodes]
					<< " " << CommWgts[i + iTotNodes];
				for (int j = Vertices.pXadj[i + iTotNodes];
						j < Vertices.pXadj[i + iTotNodes + 1];
						j++) {
					ofPartition << " " << Vertices.pAdjncy[j];
				}
				ofPartition << std::endl;
			}
			ofPartition.close();
		}
#endif /* DEBUG */

		switch (Partitioner) {
		case PARTITIONER_MANUAL:
			ASSERT(0);
			break;
				
		case PARTITIONER_CHACO: {
			int	*vwgt = 0;
			int	*cwgt = 0;
			int	*ewgt = 0;

			if (wgtflag & WEIGHT_VERTICES) {
				vwgt = &VertexWgts[0];
			}

			if (wgtflag & WEIGHT_COMM) {
				/* unsupported */
				silent_cout("communication weights currently unsupported by CHACO; trying to emulate..." << std::endl);
				cwgt = &CommWgts[0];
			}

			if (wgtflag & WEIGHT_EDGES) {
				/* unsupported ... */
				silent_cout("edges weights currently unsupported by MBDyn" << std::endl);
			}

			chaco_interface(iTotVertices,
					Vertices.pXadj,
					Vertices.pAdjncy,
					vwgt,
					cwgt,
					ewgt,
					DataCommSize,
					pParAmgProcs);
			break;
		}

		case PARTITIONER_METIS: {
			int	*vwgt = 0;
			int	*cwgt = 0;

			if (wgtflag & WEIGHT_VERTICES) {
				vwgt = &VertexWgts[0];
			}

			if (wgtflag & WEIGHT_EDGES) {
				/* unsupported ... */
				silent_cout("edges weights currently unsupported by MBDyn" << std::endl);
			}

			if (wgtflag & WEIGHT_COMM) {
				cwgt = &CommWgts[0];
			}

			mbdyn_METIS_PartGraph(iTotVertices,
					Vertices.pXadj,
					Vertices.pAdjncy,
					vwgt,
					cwgt,
					NULL,	/* ewgt */
					DataCommSize,
					pParAmgProcs);
			break;
		}

		default:
			silent_cerr("no partition library is available; "
				"partition assignments must be provided "
				"manually" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* Stima del # vertici per blocco */
	int MyDim = iTotVertices/DataCommSize + DataCommSize;

	/* Lista dei nodi appartenenti a questo processo */
	Adjacency InterfNodes; /* nodi di interfaccia */

	/* Anche qui si usa un ciclo while per verificare
	 * se il numero di nodi di interfaccia per processo
	 * ipotizzato, iMaxInterfNodes, è sufficiente */
	integer iMaxInterfNodes = iDefaultMaxInterfNodes;
	InterfNodes.pXadj = NULL;
	InterfNodes.pAdjncy = NULL;
	SAFENEWARR(InterfNodes.pXadj, int, DataCommSize);
	SAFENEWARR(ppMyNodes, Node*, 2*MyDim);
	memset(InterfNodes.pXadj, 0, DataCommSize*sizeof(int));
	memset(ppMyNodes, 0, 2*MyDim*sizeof(Node*));

	while (true) {
		InitList(InterfNodes.pXadj, DataCommSize, 0);
		SAFENEWARR(InterfNodes.pAdjncy, int, DataCommSize*iMaxInterfNodes*8);
		memset(InterfNodes.pAdjncy, 0, DataCommSize*iMaxInterfNodes*8*sizeof(int));
		InitList(InterfNodes.pAdjncy, DataCommSize*iMaxInterfNodes*8, ADJ_UNDEFINED);

		iNumLocNodes = 0;
		iNumLocDofs = 0;
		iCount = 0;
		int TmpPrc = 0;
		bool bIsNInterf;

		for (unsigned int i = 0; i < iTotNodes; i++) {
			bIsNInterf = true;

			if (pParAmgProcs[i] == MyRank) {
				/* se uno dei nodi e' connesso ad un elemento non appartenente
				 * a questo processo e' un nodo di interfaccia */
				for (int j = Vertices.pXadj[i]; j < Vertices.pXadj[i + 1]; j++) {
					TmpPrc = pParAmgProcs[Vertices.pAdjncy[j]];
					if (TmpPrc != MyRank) {
						InterfNodes.pXadj[TmpPrc] += 1;
						int k = TmpPrc * iMaxInterfNodes * 2 + InterfNodes.pXadj[TmpPrc] - 1;
						if (k < DataCommSize * iMaxInterfNodes * 2) {
							InterfNodes.pAdjncy[k] = i;
							bIsNInterf = false;
						}
					}
				}
				
				if (bIsNInterf) {
					ppMyNodes[iCount] = ppNodes[i];
					iNumLocNodes++;
					iNumLocDofs += ppMyNodes[iCount]->iGetNumDof();
					iCount++;
				}
			}
		}

		iMax = 0;
		for (int i = 0; i < DataCommSize; i++) {
			if (InterfNodes.pXadj[i] > iMaxInterfNodes) {
				iMax = InterfNodes.pXadj[i];
			}
		}

		DataComm.Allreduce(&iMax, &iRMax, 1, MPI::INT, MPI::MAX);
		iMax = iRMax;
		if (iMax <= iMaxInterfNodes) {
			break;
		}

		iMaxInterfNodes = iMax;
		SAFEDELETEARR(InterfNodes.pAdjncy);
		InterfNodes.pAdjncy = 0;
	}

	/* Scambio i dati riguardo ai nodi di interfaccia */
	MPI::Request* pRReq = NULL;
	MPI::Request* pSReq = NULL;
	
	SAFENEWARRNOFILL(pRReq, MPI::Request, DataCommSize);
	SAFENEWARRNOFILL(pSReq, MPI::Request, DataCommSize);
	for (int i = 0; i < DataCommSize; i++) {
		pRReq[i] = MPI::REQUEST_NULL;
		pSReq[i] = MPI::REQUEST_NULL;
	}

	const int DIM_TAG = 10;

	for (int i = 0; i < DataCommSize; i++) {
#if 0
		if (i != MyRank)
#endif
		{
			pRReq[i] = DataComm.Irecv(InterfNodes.pAdjncy + iMaxInterfNodes + i*iMaxInterfNodes*2,
					iMaxInterfNodes, MPI::INT, i, DIM_TAG);
			pSReq[i] = DataComm.Isend(InterfNodes.pAdjncy + i*iMaxInterfNodes*2,
					iMaxInterfNodes, MPI::INT, i, DIM_TAG);
		}
	}

	/* lista degli elementi appartenenti a questo processo */
	SAFENEWARR(ppMyElems, Elem*, 2*MyDim);
	memset(ppMyElems, 0, 2*MyDim*sizeof(Elem*));

	/* Trattamento elementi particolari */
	int move = 0;

	/* Gravity */
	if (!ElemData[Elem::GRAVITY].ElemContainer.empty()) {
		/* FIXME: there's a better way to find GravityPos and so on... */
		ppMyElems[iNumLocElems] = Elems[GravityPos];
		iNumLocElems++;
		pParAmgProcs[GravityPos + iTotNodes] = -1;
		move++;
	}

	/* Air Properties */
	if (!ElemData[Elem::AIRPROPERTIES].ElemContainer.empty()) {
		ppMyElems[iNumLocElems] = Elems[AirPropPos];
		iNumLocElems++;
		pParAmgProcs[AirPropPos + iTotNodes] = -1;
		move++;
	}

	/* Induced Velocity elements */
	int iMyTotRot = 0;
	integer* pMyRot = NULL;
	integer iIVIsMine = 0;
	if (iNumRt  != 0) {
		SAFENEWARR(pMyRot, integer, iNumRt);
		for (ElemMapType::const_iterator i = ElemData[Elem::AERODYNAMIC].ElemContainer.begin();
			i != ElemData[Elem::AERODYNAMIC].ElemContainer.end();
			i++)
		{
			const AerodynamicElem *pAero = dynamic_cast<AerodynamicElem *>(i->second);
			ASSERT(pAero != NULL);
			const InducedVelocity *pIV = pAero->pGetInducedVelocity();

			if (pIV != NULL) {
				int pos = 0;
				unsigned int pTmpLab = pIV->GetLabel();
				
				for (int k = 0; k < iNumRt; k++) {
					if (pTmpLab == pRotLab[k]) {
						pos = pRotPos[k];
					}
				}

				for (int j = 0; j < iMyTotRot; j++) {
					if (pos != pMyRot[j]) {
						pMyRot[iMyTotRot] = pos;
						iMyTotRot++;
					}
				}

				if (!iMyTotRot) {
					pMyRot[iMyTotRot] = pos;
					iMyTotRot++;
				}
			}
		}

		/* Costruisco  i communicators per i rotori */
		SAFENEWARRNOFILL(pIndVelComm, MPI::Intracomm, iNumRt);
		
		int color, key;
		for (int i = 0; i < iNumRt; i++) {
			if (iMyTotRot == 0) {
				color = MPI::UNDEFINED;
				key = MyRank;
			} else {
				color = 1;
				if (pParAmgProcs[pMyRot[i] + iTotNodes] == MyRank) {
					silent_cout("InducedVelocity(" << Elems[pRotPos[i]]->GetLabel() << ")"
						<< " assigned to process "
						<< MyRank << std::endl);
					iIVIsMine = 1;
					key = 0;
				} else {
					key = MyRank + 1;
				}
			}
			pIndVelComm[i] = MBDynComm.Split(color, key);
#if 0
			IndVelComm[i] = MPI::COMM_WORLD.Split(color, key);
#endif
			InducedVelocity *r = dynamic_cast<InducedVelocity *>(Elems[pRotPos[i]]);
			ASSERT(r != 0);
			r->InitializeIndVelComm(pIndVelComm + i);
		}

		for (int i = 0; i < iMyTotRot; i++) {
			ppMyElems[iNumLocElems] = Elems[pMyRot[i]];
			if (pParAmgProcs[pMyRot[i] + iTotNodes] == MyRank) {
				iNumLocDofs += ppMyElems[iNumLocElems]->iGetNumDof();
			} else {
				move++;
			}
			iNumLocElems++;
			pParAmgProcs[pMyRot[i] + iTotNodes] = -1;
		}
	}

	for (unsigned int i = 0; i < iTotVertices - iTotNodes; i++) {
		if (pParAmgProcs[iTotNodes + i] == MyRank) {
			ppMyElems[iNumLocElems] = Elems[i];
			iNumLocDofs += ppMyElems[iNumLocElems]->iGetNumDof();
			iNumLocElems++;
		}
	}

	/* initialize local element iterator */
	MyElemIter.Init(ppMyElems, iNumLocElems);

#ifdef USE_MPI_WAIT
	MPI::Request::Waitall(DataCommSize, pRReq);
	MPI::Request::Waitall(DataCommSize, pSReq);
#else /* ! USE_MPI_WAIT */
	/* Verifico la ricezione dei nodi di interfaccia */
	bool bRecvFlag = false, bSentFlag = false;
	while (true) {
		if (!bRecvFlag) {
			bRecvFlag = MPI::Request::Testall(DataCommSize, pRReq);
		}

		if (!bSentFlag) {
			bSentFlag = MPI::Request::Testall(DataCommSize, pSReq);
		}
		
		if (bRecvFlag && bSentFlag) {
			break;
		}

		MYSLEEP(1000);
	}
#endif /* ! USE_MPI_WAIT */

	/* ordino i nodi interfaccia */
	std::sort(InterfNodes.pAdjncy,
			InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);
	
	/* elimino le ripetizioni */
	int* p = std::unique(InterfNodes.pAdjncy,
			InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);

	/* dimensione effettiva dell'interfaccia locale
	 * il -1 e' dovuto al primo valore che è sempre pari a -1 */
	/*
	 * InterfNodes.pAdjncy[0] == -1
	 * InterfNodes.pAdjncy[iNumIntNodes] == -1
	 */
	iNumIntNodes = p - &InterfNodes.pAdjncy[1];

	unsigned int* llabels = NULL;
	Node::Type* lTypes = NULL;
	SAFENEWARR(llabels, unsigned int, iNumIntNodes);
	SAFENEWARR(lTypes,  Node::Type, iNumIntNodes);
	SAFENEWARR(ppIntNodes, Node*, iNumIntNodes);
	for (int i = 0; i < iNumIntNodes; i++) {
		ppIntNodes[i] = ppNodes[InterfNodes.pAdjncy[i + 1]];
		llabels[i] = ppIntNodes[i]->GetLabel();
		lTypes[i] = ppIntNodes[i]->GetNodeType();
	}

	/* scrittura dei dofs degli stati interni di elementi che connettono
	 * un nodo interno ad uno di interfaccia, sulla lista dei dof
	 * di interfaccia */
	int* pPosIntElems = NULL;
	if (iNumLocElems != 0) {
		SAFENEWARR(pPosIntElems, int, iNumLocElems);
		SAFENEWARR(ppMyIntElems, Elem*, iNumLocElems);

		memset(pPosIntElems, 0, iNumLocElems*sizeof(int));
		memset(ppMyIntElems, 0, iNumLocElems*sizeof(Elem*));
	}

	for (int i = 0; i < iNumLocElems; i++) {
		if (ppMyElems[i]->iGetNumDof() != 0) {
			Elem::Type CType = ppMyElems[i]->GetElemType();
			switch (CType) {
			case Elem::INDUCEDVELOCITY:
				if (iIVIsMine == 1) {
					ppMyElems[i]->GetConnectedNodes(connectedNodes);
					for (std::vector<const Node *>::const_iterator j = connectedNodes.begin();
						j != connectedNodes.end();
						j++)
					{
						unsigned int* p = std::find(llabels, llabels + iNumIntNodes, (*j)->GetLabel());
						if (p != llabels + iNumIntNodes) {
							ppMyIntElems[iNumIntElems] =  ppMyElems[i];
							iNumIntDofs += ppMyElems[i]->iGetNumDof();
							iNumLocDofs -= ppMyElems[i]->iGetNumDof();
							pPosIntElems[iNumIntElems] = i;
							iNumIntElems++;
							break;
						}
					}
				}
				break;

			case Elem::HYDRAULIC:
			case Elem::GENEL:
			case Elem::ELECTRIC:
			case Elem::THERMAL:
				break;

			default:
				ppMyElems[i]->GetConnectedNodes(connectedNodes);
				for (std::vector<const Node *>::const_iterator j = connectedNodes.begin();
					j != connectedNodes.end();
					j++)
				{
					unsigned int* p = std::find(llabels, llabels + iNumIntNodes, (*j)->GetLabel());
					if (p != llabels + iNumIntNodes) {
						ppMyIntElems[iNumIntElems] =  ppMyElems[i];
						iNumIntDofs += ppMyElems[i]->iGetNumDof();
						iNumLocDofs -= ppMyElems[i]->iGetNumDof();
						pPosIntElems[iNumIntElems] = i;
						iNumIntElems++;
						break;
					}
				}
				break;
			}
		}
	}

	/* determina la liste dei dofs locali ed adiacenti suddivisa per processi,
	 * secondo la struttura Adjacency */
	SAFENEWARR(pLocalDofs, integer, iNumLocDofs);

	iCount = 0;
	for (int i = 0; i < iNumLocNodes; i++) {
		if (ppMyNodes[i]->iGetNumDof() != 0) {
			integer First = (ppMyNodes[i])->iGetFirstIndex();

			pLocalDofs[iCount] = First + 1;
			iCount++;
			for (unsigned int j = 1; j < ppMyNodes[i]->iGetNumDof(); j++) {
				pLocalDofs[iCount] = First + j + 1;
				iCount++;
			}
		}
	}

	/* Aggiungo i dofs degli elementi locali */
	integer TmpDofNum;
	int i2Count = 0;

	for (int i = move; i < iNumLocElems; i++) {
		TmpDofNum = ppMyElems[i]->iGetNumDof();
		if (TmpDofNum != 0) {
			if (i != pPosIntElems[i2Count]) {
				ElemWithDofs* pWithDofs = dynamic_cast<ElemWithDofs *>(ppMyElems[i]);
				integer First = (pWithDofs)->iGetFirstIndex();

				pLocalDofs[iCount] = First + 1;
				iCount++;
				for (int j = 1; j < TmpDofNum; j++) {
					pLocalDofs[iCount] = First + 1 + j;
					iCount++;
				}
				
			} else {
				i2Count++;
			}
		}
	}

	/* scrivo ora la lista dei dofs interfaccia */
	for (int i = 1; i < iNumIntNodes + 1; i++) {
		iNumIntDofs += ppNodes[InterfNodes.pAdjncy[i]]->iGetNumDof();
	}

	SAFENEWARR(pLocalIntDofs, integer, iNumIntDofs);
	SAFENEWARR(pMyIntDofs, integer, iNumIntDofs);

	iCount = 0;
	i2Count = 0;
	for (int i = 1; i <= iNumIntNodes; i++) {
		if (ppNodes[InterfNodes.pAdjncy[i]]->iGetNumDof() != 0) {
			integer First = ppNodes[InterfNodes.pAdjncy[i]]->iGetFirstIndex();

			pLocalIntDofs[iCount] = First + 1;
			iCount++;
			if (pParAmgProcs[InterfNodes.pAdjncy[i]] == MyRank) {
				pMyIntDofs[i2Count] = First + 1;
				i2Count++;
			}

			for (unsigned int j = 1;
					j < ppNodes[InterfNodes.pAdjncy[i]]->iGetNumDof();
					j++)
			{
				/* il - serve a distinguere questi dofs da quelli interni */
				pLocalIntDofs[iCount] = (First + 1 + j);
				iCount++;
				if (pParAmgProcs[InterfNodes.pAdjncy[i]] == MyRank) {
					pMyIntDofs[i2Count] = (First + 1 + j);
					i2Count++;
				}
			}
		}
	}

	/* Interfaccia degli elementi locali */
	for (int i = 0; i < iNumIntElems; i++) {
		TmpDofNum = ppMyIntElems[i]->iGetNumDof();
		ElemWithDofs* pWithDofs = dynamic_cast<ElemWithDofs *>(ppMyIntElems[i]);
		integer First = (pWithDofs)->iGetFirstIndex();
		pLocalIntDofs[iCount] = First + 1;
		iCount++;
		for (int j = 1; j < TmpDofNum; j++) {
			pLocalIntDofs[iCount] = First  + 1 + j;
			iCount++;
		}
	}

	iNumMyInt = i2Count;

	if ( Vertices.pXadj != NULL) {
		SAFEDELETEARR(Vertices.pXadj);
	}

	if ( Vertices.pAdjncy != NULL) {
		SAFEDELETEARR(Vertices.pAdjncy);
	}

	if ( InterfNodes.pXadj != NULL) {
		SAFEDELETEARR(InterfNodes.pXadj);
	}

	if ( InterfNodes.pAdjncy != NULL) {
		SAFEDELETEARR( InterfNodes.pAdjncy);
	}

	if ( pParAmgProcs != NULL) {
		SAFEDELETEARR(pParAmgProcs);
	}

	if ( pLabelsList != NULL) {
		SAFEDELETEARR(pLabelsList);
	}

	if ( pPosIntElems != NULL) {
		SAFEDELETEARR(pPosIntElems);
	}

	if ( lTypes != NULL) {
		SAFEDELETEARR(lTypes);
	}

	if ( llabels != NULL) {
		SAFEDELETEARR(llabels);
	}

	OutputPartition();
}

void
SchurDataManager::OutputPartition(void)
{
	silent_cout("Making partition file...");

	/* Inizializzazione */
	OutHdl.PartitionOpen();

	time_t tCurrTime(time(NULL));
	OutHdl.Partition()
		<< "# Partition file for MBDyn. Time: " << ctime(&tCurrTime)  << std::endl
		<< "# Partition produced ";
	switch (Partitioner) {
	case PARTITIONER_MANUAL:
		OutHdl.Partition()
			<< "manually";
		break;

	case PARTITIONER_METIS:
		OutHdl.Partition()
			<< "with METIS";
		break;

	case PARTITIONER_CHACO:
		OutHdl.Partition()
			<< "with CHACO";
		break;

	default:
		ASSERT(0);
	}
	OutHdl.Partition()
		<< std::endl;

	/* Dati */
	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "# Control data to verify the partition "
		<< std::endl
		<< std::endl
		<< "Total number of processes: " << DataCommSize << ";" << std::endl
		<< "Process #: " << MyRank << ";" << std::endl
		<< "Total number of Nodes: " << iTotNodes << ";"  << std::endl
		<< "Total number of Elements: " << Elems.size() << ";"  << std::endl
		<< "Total number of Dofs: " << iTotDofs << ";"  << std::endl;

	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Local Dofs number: " << iNumLocDofs << std::endl;

	OutHdl.Partition()
		<< std::endl
		<< "Local Nodes: " << iNumLocNodes << std::endl;
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		OutHdl.Partition()
			<< "Node Type: "
			<< "(" << psNodeNames[(ppMyNodes[i])->GetNodeType()] << ")"
			<< " Label: " << ppMyNodes[i]->GetLabel()
			<< " Dofs #: " << ppMyNodes[i]->iGetNumDof()
			<< std::endl;
	}

	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Local Elements: "<< iNumLocElems << std::endl;
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		OutHdl.Partition()
			<< "Element Type: "
			<< "("  << psElemNames[ppMyElems[i]->GetElemType()] << ")"
			<< " Label: " << ppMyElems[i]->GetLabel()
			<< " Dofs #: " << ppMyElems[i]->iGetNumDof()
			<< std::endl;
	}

	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Interface Dofs number: " << iNumIntDofs << std::endl
		<< std::endl
		<< "Interface Nodes: " << iNumIntNodes << std::endl;
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		OutHdl.Partition()
			<< "Node Type: "
			<< "(" << psNodeNames[(ppIntNodes[i])->GetNodeType()] << ")"
			<< " Label: " << ppIntNodes[i]->GetLabel()
			<< " Dofs #: " <<  ppIntNodes[i]->iGetNumDof()
			<< std::endl;
	}

	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Elements whose internal dofs are interface dofs: "
		<< iNumIntElems << std::endl;
	for (int i = 0; i < iNumIntElems; i++) {
		OutHdl.Partition()
			<< "Element Type: "
			<< "("  << psElemNames[(ppMyIntElems[i])->GetElemType()] << ")"
			<< " Label: " << ppMyIntElems[i]->GetLabel()
			<< " Dofs #: " << ppMyIntElems[i]->iGetNumDof()
			<< std::endl;
	}

#ifdef DEBUG
	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Local Dofs List:" << std::endl;

	for (int i = 0; i < iNumLocDofs; i++) {
		OutHdl.Partition() << pLocalDofs[i] << " ";
		if (!(iNumLocDofs%10)) {
			OutHdl.Partition() << std::endl;
		}
	}
	if (iNumLocDofs%10) {
		OutHdl.Partition()
			<< std::endl;
	}

	OutHdl.Partition()
		<< std::endl
		<< std::endl
		<< "Local Interface Dofs List:" <<std::endl;
	for (int i = 0; i < iNumIntDofs; i++) {
		OutHdl.Partition() << pLocalIntDofs[i] << " ";
		if (!(i%10)) {
			OutHdl.Partition() << std::endl;
		}
	}
	OutHdl.Partition() << std::endl;
#endif /* DEBUG */

	silent_cout("... partition file done" << std::endl);
}

/* compatta il vettore delle adiacenze */
void
SchurDataManager::Pack(int* pList, int dim)
{
	int* pOld = pList;
	int* pNew = pList;

	for (; pOld < pList + dim; pOld++) {
		if (*pOld != ADJ_UNDEFINED) {
			*pNew = *pOld;
			pNew++;
		}
	}
}

/* Inizializza le varie liste di interi usate in Create Partition */
void
SchurDataManager::InitList(int* list, int dim, int value)
{
	for (int i = 0; i < dim; i++) {
		list[i] = value;
	}
}

void
SchurDataManager::InitList(float* list, int dim, int value)
{
	for (int i = 0; i < dim; i++) {
		list[i] = value;
	}
}


/* Trova la posizione del nodo nell'array ppNodes */
Node**
SchurDataManager::SearchNode(Node** ppFirst, int dim, unsigned int& label)
{
	unsigned int* pbegin = pLabelsList + (ppFirst - ppNodes);
	unsigned int* p = std::find(pbegin, pbegin + dim, label);

	return ppNodes + (p - pLabelsList);
}


void
SchurDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef) throw(ChangedEquationStructure)
{
	DEBUGCOUT("Entering SchurDataManager::AssRes()" << std::endl);
	ASSERT(pWorkVec != NULL);

	try {
		DataManager::AssRes(ResHdl, dCoef, MyElemIter, *pWorkVec);

	} catch (ChangedEquationStructure) {
		Elem *pEl = NULL;

		if (MyElemIter.bGetCurr(pEl) == true) {
			silent_cerr("Jacobian reassembly requested by "
				<< psElemNames[pEl->GetElemType()]
				<< "(" << pEl->GetLabel() << "); "
				"currently unsupported by Schur data manager."
				<< std::endl);

		} else {
			silent_cerr("Jacobian reassembly requested "
				"by an element; currently unsupported "
				"by Schur data manager." << std::endl);
		}

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
SchurDataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
	DEBUGCOUT("Entering SchurDataManager::AssJac()" << std::endl);
	ASSERT(pWorkMat != NULL);

	DataManager::AssJac(JacHdl, dCoef, MyElemIter, *pWorkMat);
}
/* End of AssJac */

void
SchurDataManager::Update(void) const
{
	DEBUGCOUT("Entering SchurDataManager::Update()" << std::endl);

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		ppMyNodes[i]->Update(*pXCurr, *pXPrimeCurr);
	}

	/* Nodi di interfaccia */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		ppIntNodes[i]->Update(*pXCurr, *pXPrimeCurr);
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->Update(*pXCurr, *pXPrimeCurr);
	}
}
/* End of Update */

void
SchurDataManager::DerivativesUpdate(void) const
{
	DEBUGCOUT("Entering SchurDataManager::DerivativesUpdate()" << std::endl);

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		if (ppMyNodes[i]->GetNodeType() == Node::STRUCTURAL) {
			((StructNode*)ppMyNodes[i])->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
		} else {
			ppMyNodes[i]->Update(*pXCurr, *pXPrimeCurr);
		}
	}

	/* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		if (ppIntNodes[i]->GetNodeType() == Node::STRUCTURAL) {
			((StructNode*)ppIntNodes[i])->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
		} else {
			ppIntNodes[i]->Update(*pXCurr, *pXPrimeCurr);
		}
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->Update(*pXCurr, *pXPrimeCurr);
	}
}
/* End of DerivativeUpdate */


void
SchurDataManager::BeforePredict(VectorHandler& X, VectorHandler& XP,
		VectorHandler& XPrev, VectorHandler& XPPrev) const
{
	DEBUGCOUT("Entering SchurDataManager::BeforePredict()" << std::endl);

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		ppMyNodes[i]->BeforePredict(X, XP, XPrev, XPPrev);
	}

	/* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		ppIntNodes[i]->BeforePredict(X, XP, XPrev, XPPrev);
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->BeforePredict(X, XP, XPrev, XPPrev);
	}
}
/* End of BeforePredict */

void
SchurDataManager::AfterPredict(void) const
{
	DEBUGCOUT("Entering SchurDataManager::AfterPredict()" << std::endl);

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		ppMyNodes[i]->AfterPredict(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}

	/* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		ppIntNodes[i]->AfterPredict(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->AfterPredict(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}
}
/* End of AfterPredict */


void
SchurDataManager::AfterConvergence(void) const
{
	DEBUGCOUT("Entering SchurDataManager::AfterConvergence()" << std::endl);

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		ppMyNodes[i]->AfterConvergence(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}

	/* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		ppIntNodes[i]->AfterConvergence(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->AfterConvergence(*(VectorHandler*)pXCurr,
				*(VectorHandler*)pXPrimeCurr);
	}

	/* Restart condizionato */
	switch (RestartEvery) {
	case ITERATIONS:
		if (++iCurrRestartIter == iRestartIterations) {
			iCurrRestartIter = 0;
			((SchurDataManager*)this)->MakeRestart();
		}
		break;

	case TIME: {
		doublereal dT = dGetTime();
		if (dT - dLastRestartTime >= dRestartTime) {
			dLastRestartTime = dT;
			const_cast<SchurDataManager *>(this)->MakeRestart();
		}
		break;
	}

	default:
		break;
	}
}
/* End of AfterConvergence */


/* stampa i risultati */
void
SchurDataManager::Output(bool force) const
{
	DEBUGCOUT("Entering SchurDataManager::Output()" << std::endl);

	/* output only at multiples of iOutputFrequency */
	if ((!force) && !pOutputMeter->dGet()) {
		return;
	}

	/* Nodi */
	for (int i = 0; i < iNumLocNodes; i++) {
		ASSERT(ppMyNodes[i] != NULL);
		ppMyNodes[i]->Output(OutHdl);
	}

	/* Nodi di interfaccia */
	for (int i = 0; i < iNumIntNodes; i++) {
		ASSERT(ppIntNodes[i] != NULL);
		ppIntNodes[i]->Output(OutHdl);
	}

	/* Elementi */
	for (int i = 0; i < iNumLocElems; i++) {
		ASSERT(ppMyElems[i] != NULL);
		ppMyElems[i]->Output(OutHdl);
	}
}

/* End Output */

/* End SchurDataManager */

#endif /* USE_MPI */

