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

/* Schur Data Manager */

/* 
 * Copyright 1999-2003 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

//#ifdef USE_MPI
/* libreria per il calcolo delle partizioni */
#ifdef USE_METIS
extern "C" {
#include <metis.h>
}
#undef ASSERT /* kill annoying redefiniton message */
#endif /* USE_METIS */

#include <schurdataman.h>
#include <mbcomm.h>
#include <mysleep.h>
#include <except.h>

#undef min
#undef max
#include <algorithm>

#include <rotor.h>

/* struttura contenente le distribuzioni di dofs */
struct Adjacency{
    int* pXadj;            /* x ogni nodo: pos. su pAdjncy lista connessioni */
    int* pAdjncy;          /* lista liste connessioni */ 
};

/* !!!!!!! Costanti da verificare !!!!!*/
const integer iDefaultMaxConnectionsPerVertex = 100;
const integer iDefaultMaxNodesPerElem = 50;
const integer iDefaultMaxInterfNodes = 50;
const int ADJ_UNDEFINED = -1;

#ifndef USE_MPI
SchurDataManager::SchurDataManager(MBDynParser& HP, 
		unsigned OF,
		doublereal dInitialTime,
		const char* sInputFileName,
		const char* sOutputFileName, 
		bool bAbortAfterInput)
: DataManager(HP, OF, dInitialTime, sInputFileName, 
              sOutputFileName, bAbortAfterInput)
{
	std::cerr << "fatal error: you are building SchurDataManager,\n" <<
		"but mbdyn was compiled without MPI.\n" <<
		"Something weird is happening\n" <<
		"Anyway, please compile with -DUSE_MPI\n" <<
		"to enable paralle solution\n";
	THROW(ErrGeneric());
};

SchurDataManager::~SchurDataManager() {
	NO_OP;
};

void SchurDataManager::AssRes(VectorHandler&, double) { NO_OP; };
void SchurDataManager::AssJac(MatrixHandler&, double) { NO_OP; };
void SchurDataManager::DerivativesUpdate(void) const { NO_OP; };
void SchurDataManager::BeforePredict(VectorHandler&, VectorHandler&,
		VectorHandler&, VectorHandler&) const { NO_OP; };
void SchurDataManager::AfterPredict(void) const { NO_OP; };
void SchurDataManager::Update(void) const { NO_OP; };
void SchurDataManager::AfterConvergence(void) const { NO_OP; };
integer SchurDataManager::HowManyDofs(DofType) const { NO_OP; };
integer* SchurDataManager::GetDofsList(DofType) const { NO_OP; };
Dof* SchurDataManager::pGetDofsList(void) const { NO_OP; };
#else /* USE_MPI */

/* Costruttore - begin */
SchurDataManager::SchurDataManager(MBDynParser& HP, 
		unsigned OF,
		doublereal dInitialTime,
		const char* sInputFileName,
		const char* sOutputFileName, 
		bool bAbortAfterInput)
: DataManager(HP, OF, dInitialTime, sInputFileName, 
              sOutputFileName, bAbortAfterInput),
iTotVertices(0),
ppMyElems(NULL),
iNumLocElems(0),
ppMyIntElems(NULL),
iNumIntElems(0),
ppMyNodes(NULL),
iNumLocNodes(0),
LocalDofs(NULL),
iNumLocDofs(0),
LocalIntDofs(NULL),
iNumIntDofs(0),
ppIntNodes(NULL),
iNumIntNodes(0),
iNumMyInt(0),
pMyIntDofs(NULL),
pLabelsList(NULL),
wgtflag(1),
pParAmgProcs(NULL),
pRotorComm(NULL),
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

    iTotVertices = iTotNodes + iTotElem;
    DEBUGCOUT("iTotVertexes: " << iTotVertices << std::endl);

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
        END,
        LASTKEYWORD
    };

    /* tabella delle parole chiave */
    KeyTable K((int)LASTKEYWORD, sKeyWords);
  
    /* parser del blocco parallelizzazione */
    HP.PutKeyTable(K);
  
    /* legge la distribuzione degli elementi sulle diverse CPU */

    if (KeyWords(HP.GetDescription()) != BEGIN) {
    std::cerr << std::endl 
	 << "Warning: no explicit connection declared"
	 << " for this input file " << std::endl;
        return;
    } else {
        int iNumElems = 0;
        int iNumNodes = 0;
        if (KeyWords(HP.GetWord()) != PARALLEL) {
            std::cerr << std::endl 
	        << "Error: <begin: parallel;> expected at line "
	        << HP.GetLineData() << "; aborting ..." << std::endl;
            THROW(ErrGeneric());
        }
        while (true) {
            KeyWords CDesc = KeyWords(HP.GetDescription());
            switch (CDesc) {
            case WEIGHTS:
	        if (HP.fIsArg()) {
	            wgtflag = HP.GetInt();
	        } else {
	            std::cerr << std::endl << "Error: Number expected at line "  
	                << HP.GetLineData() << "; aborting ..." << std::endl;
	            THROW(ErrGeneric());
	        }
	        break;
	
            case PARTITION:
      	        SAFENEWARR(pParAmgProcs, int, iTotVertices);	
      	        for (int i = 0; i < iTotVertices; i++) {
	            if (HP.fIsArg()) {
	                pParAmgProcs[i] = HP.GetInt();
	            } else {
	                std::cerr << std::endl 
			    << "Error: the partition assignment"
			    " is not complete at line "  
	                    << HP.GetLineData() << "; aborting ..." << std::endl;
	                THROW(ErrGeneric());
	            }
	        }  
	        break;				
	
            case END:   
	        if (KeyWords(HP.GetWord()) != PARALLEL) {
	            std::cerr << std::endl 
		        << "Error: <end: parallel;> expected at line "
	                << HP.GetLineData() << "; aborting ..." << std::endl;
	            THROW(ErrGeneric());
	        }
	        goto endcycle;
    
            default:
	        std::cerr << "Unknown input at line "
	            << HP.GetLineData() << "; aborting ..." << std::endl;
	        THROW(ErrGeneric());
    
            case NUMBEROFCONNECTIONS: 
	        iTotalExpConnections = HP.GetInt();
	        SAFENEWARR(ppExpCntNodes, Node*, iTotalExpConnections);
	        SAFENEWARR(ppExpCntElems, Elem*, iTotalExpConnections);
      
	        Elem::Type ActualElType;
	        Node::Type ActualNdType;
	        unsigned int j = 0;
      
	        for (int i = 0; i < iTotalExpConnections; i++) {
	            if (KeyWords(HP.GetDescription()) != CONNECTION) {
	                std::cerr << std::endl 
			    << "Error: <Connection> expected at line "
		            << HP.GetLineData() << "; aborting ..." << std::endl;
	                THROW(ErrGeneric());
	            }
	  
	            for (int k = 0; k < 2; k++) {
	                KeyWords CurrDesc = KeyWords(HP.GetWord());
	                switch (CurrDesc) {
	                case ELEMENT: {
	                    KeyWords ElDesc = KeyWords(HP.GetWord());
	                    switch (ElDesc) {
	                    case FORCE:
		                ActualElType = Elem::FORCE;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		                    THROW(ErrGeneric());
		                }
		                ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
		                if (ppExpCntElems[iNumElems] == NULL ) {
		                    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		                    THROW(ErrGeneric());
		                }
		                iNumElems++;
		                break;
				
	                    case RIGIDBODY:
		                ActualElType = Elem::BODY;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		                    THROW(ErrGeneric());
		                }
		                ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		                    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		                    THROW(ErrGeneric());
		                }
		                iNumElems++;
		                break;
		
	                    case JOINT:
		                ActualElType = Elem::JOINT;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    std::cerr << std::endl 
				        << "Error: Label expected at line "  
				        << HP.GetLineData() 
				        << "; aborting ..." << std::endl;
				    THROW(ErrGeneric());
			       	}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
			  	    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case BEAM:
				ActualElType = Elem::BEAM;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				      	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;

	      		    case PLATE:
				ActualElType = Elem::PLATE;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ROTOR:
				ActualElType = Elem::ROTOR;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case AEROMODAL:
				ActualElType = Elem::AEROMODAL;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;

	      		    case AERODYNAMICELEMENT:
				ActualElType = Elem::AERODYNAMIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
			  	    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ELECTRICBULK:
				ActualElType = Elem::ELECTRICBULK;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ELECTRIC:
				ActualElType = Elem::ELECTRIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case GENEL:
				ActualElType = Elem::GENEL;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case HYDRAULIC:
				ActualElType = Elem::HYDRAULIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case BULK:
				ActualElType = Elem::BULK;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
			  	    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case LOADABLE:
				ActualElType = Elem::LOADABLE;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
				    std::cerr << "Error: at line " 
				     	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
				    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case DRIVEN:
				ActualElType = Elem::DRIVEN;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				      	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    std::cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    default:
				std::cerr << "Error: Element Type on line " 
				   << HP.GetLineData() 
				   << " not valid; aborting ..." << std::endl;
				THROW(ErrGeneric());
	      		    }
	      		    break;
	    		}
	    
	    		case NODE: {
	      		    KeyWords NdDesc = KeyWords(HP.GetWord());
	      		    switch (NdDesc) {		
	      		    case ABSTRACT:
				ActualNdType = Node::ABSTRACT;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				     	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    case STRUCTURAL:
				ActualNdType = Node::STRUCTURAL;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
				    std::cerr << std::endl 
				      	<< "Error: Label expected at line "  
		     			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
				    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    case ELECTRIC:
				ActualNdType = Node::ELECTRIC;
				
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				      	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	    		    case PARAMETER:
	      			ActualNdType = Node::PARAMETER;
	      			if (HP.fIsArg()) {
				    j = HP.GetInt();
	      			} else {
				    std::cerr << std::endl 
				    	<< "Error: Label expected at line "  
		     			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
				    THROW(ErrGeneric());
	      			}
	      			ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
	      			if (ppExpCntNodes[iNumNodes] == NULL ) {
				    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< std::endl;
				    THROW(ErrGeneric());
	      			}
	      			iNumNodes++;
	      			break;
	      
	      		    case HYDRAULIC:
				ActualNdType = Node::HYDRAULIC;
				
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    std::cerr << std::endl 
				    	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << std::endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    std::cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< std::endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    default:
				std::cerr << "Error: Node Type on line " 
				    << HP.GetLineData() << " not valid "
		     		    << "; aborting ..." << std::endl;
				THROW(ErrGeneric());
	      		    }  
	      		    break;
	    		}
	    
	    		default:
	      		    std::cerr << " Unknown input at line "
		   		<< HP.GetLineData() 
				<< "; aborting ..." << std::endl;
	      		    THROW(ErrGeneric());
	    		}
	  	    }
		}
   
		ASSERT(iNumNodes == iNumElems);
		if (iNumNodes + iNumElems !=  2*iTotalExpConnections) {
	  	    std::cerr << std::endl 
		     	<< "Error: Total number of Nodes and elements"
			" in the parallel section at line "
	       		<< HP.GetLineData() 
			<< " is not consistent; aborting ..." << std::endl;
	  	    THROW(ErrGeneric());
		}  
		break;
      	    }
    	}
endcycle:
    	return; 
    }
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
    if (LocalDofs != NULL) {
    	SAFEDELETEARR(LocalDofs);
    }
    if (LocalIntDofs != NULL) {
    	SAFEDELETEARR(LocalIntDofs);
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
		std::cerr << "SchurDataManager::HowManyDofs: "
			"illegal request (" << unsigned(who) << ")"
			<< std::endl;
		THROW(ErrGeneric());
	}
}

integer*
SchurDataManager::GetDofsList(DofType who) const
{
	switch(who) {
	case LOCAL:
		return LocalDofs;

	case INTERNAL:
		return LocalIntDofs;

	case MYINTERNAL: 
		return pMyIntDofs;

	default:
		std::cerr << "SchurDataManager::GetDofsList: "
			"illegal request (" << unsigned(who) << ")"
			<< std::endl;
		THROW(ErrGeneric());
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

    Adjacency Vertices;  /* Struttura contenente le connessioni fra i vertici */
    int* pVertexWgts;   /* Pesi dei vertici = dofs x ogni v. utile per METIS */
    int* pCommWgts;
    Vertices.pXadj = NULL;
    Vertices.pAdjncy = NULL;
    pVertexWgts = NULL;
    pCommWgts = NULL;
    integer iMax = 0;
    integer iRMax = 0;
    int iCount = 0;
    int  iNumRt = 0;
    
    ASSERT(iTotVertices > 0);
    ASSERT(DataCommSize > 0);
    /* Costruisco e inizializzo le due strutture */
    SAFENEWARR(Vertices.pXadj, int, iTotVertices+1);
    SAFENEWARR(pVertexWgts, int, iTotVertices*2);
    SAFENEWARR(pCommWgts, int, iTotVertices);
    
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
    int iNumberOfNodes;
    Node::Type* pMyTypes = NULL;
    unsigned int* pMyLabels = NULL;
    if (ElemData[Elem::ROTOR].iNum != 0) {
    	SAFENEWARR(pRotPos, int, ElemData[Elem::ROTOR].iNum);
    	SAFENEWARR(pRotLab, unsigned int, ElemData[Elem::ROTOR].iNum);
    }
    SAFENEWARR(pMyTypes, Node::Type, iDefaultMaxNodesPerElem);
    SAFENEWARR(pMyLabels, unsigned int, iDefaultMaxNodesPerElem);
    SAFENEWARR(pLabelsList, unsigned int, iTotNodes);

    while (true) {
        InitList(Vertices.pXadj, iTotVertices+1, 0);
        InitList(pVertexWgts, iTotVertices*2, 0);
       	InitList(pCommWgts, iTotVertices, 0);
       	SAFENEWARR(Vertices.pAdjncy, int, iTotVertices*iMaxConnectionsPerVertex);
       	InitList(Vertices.pAdjncy, iTotVertices*iMaxConnectionsPerVertex, ADJ_UNDEFINED);
       	ASSERT(ppElems != NULL);
       
       	/* ciclo sugli elementi per assemblare la struttura delle connessioni */
       	Node** ppActualNode = NULL;
       	/* per numerare i nodi prendo la posizione del puntatore 
         * al nodo nell'array ppNodes */
       	int position;
       	iCount = iTotNodes;
       	for (unsigned int i = 0; i < iTotNodes; i++) {
            pLabelsList[i] = ppNodes[i]->GetLabel();
       	}

       	iNumRt = 0;

       	for (Elem** ppTmpEl = ppElems; 
	     ppTmpEl < ppElems+iTotElem; ppTmpEl++, 
	     iCount++) {
            if ((*ppTmpEl)->GetElemType() == Elem::GRAVITY) {
             	GravityPos = ppTmpEl - ppElems;
            }
            if ((*ppTmpEl)->GetElemType() == Elem::AIRPROPERTIES) {
             	AirPropPos = ppTmpEl - ppElems;
            }
            if ((*ppTmpEl)->GetElemType() == Elem::ROTOR) {
             	pRotPos[iNumRt] = ppTmpEl - ppElems;
             	pRotLab[iNumRt] = (*ppTmpEl)->GetLabel();
             	iNumRt++;
            }
            (*ppTmpEl)->GetConnectedNodes(iNumberOfNodes, pMyTypes, pMyLabels);
            /* peso dell'elemento */
            integer dimA, dimB;
            (*ppTmpEl)->WorkSpaceDim(&dimA, &dimB);
            pVertexWgts[iCount] = dimA * dimB;

            pCommWgts[iCount] = 
	       	(*ppTmpEl)->iGetNumDof()*(*ppTmpEl)->iGetNumDof();

	    for (int i = 0; i <= iNumberOfNodes-1; i++) {
	     	Vertices.pXadj[iCount+1] += 1;
             	/* trovo la pos. del nodo nella lista dei puntatori ai nodi */
	     	ppActualNode = SearchNode(NodeData[pMyTypes[i]].ppFirstNode, 
	                                  NodeData[pMyTypes[i]].iNum, 
				          pMyLabels[i]);
	     	position = ppActualNode - ppNodes;
	     	/* Aggiungo al peso dell'elemento il numero di dofs 
	      	 * di ciascun nodo connesso */

#if 0   /* FIXME: i nodi hanno peso comp. nullo */
	     	pVertexWgts[iCount] += (*ppActualNode)->iGetNumDof();
#endif

	     	/* aggiungo fra le connessioni dell'elemento il nodo attuale */
	     	if ((iCount*iMaxConnectionsPerVertex) + Vertices.pXadj[iCount+1] - 1 < iTotVertices*iMaxConnectionsPerVertex) {
	            Vertices.pAdjncy[(iCount*iMaxConnectionsPerVertex) 
		  	+ Vertices.pXadj[iCount+1] - 1] = position;
	     	}
	     	/* aggiungo alle connessioni del nodo l'elemento attuale */
	     	Vertices.pXadj[position+1] += 1;
	     	if ((position*iMaxConnectionsPerVertex) + Vertices.pXadj[position+1] - 1 < iTotVertices*iMaxConnectionsPerVertex) {
	            Vertices.pAdjncy[(position*iMaxConnectionsPerVertex) 
		  	+ Vertices.pXadj[position+1] - 1] = iCount;
	     	}
	     	/* peso (di comunicazione) del nodo */
	     	pCommWgts[position] = (*ppActualNode)->iGetNumDof();
	    }
       	}

       	if (iTotalExpConnections != 0) {
            for (int i = 0; i < iTotalExpConnections; i++) {
             	int iNdPos, iElPos;
             	int j = 0;
             	while (ppExpCntNodes[i] != ppNodes[j]) {
                    j++;
             	}
             	iNdPos = j;
             	j = 0;
             	while (ppExpCntElems[i] != ppElems[j]) {
                    j++;
             	}
             	iElPos = j;
             	Vertices.pXadj[iNdPos+1] += 1;
             	Vertices.pXadj[iTotNodes + iElPos+1] += 1;
             	Vertices.pAdjncy[(iNdPos*iMaxConnectionsPerVertex) 
	       	    + Vertices.pXadj[iNdPos+1] - 1] = iElPos+iTotNodes;
             	Vertices.pAdjncy[((iTotNodes+iElPos)*iMaxConnectionsPerVertex) 
	       	    + Vertices.pXadj[iTotNodes+iElPos+1] - 1] = iNdPos;
            }
     	}
     
        iMax = 0;
        for (int i = 0; i < iTotVertices; i++) {
            if (Vertices.pXadj[i] > iMaxConnectionsPerVertex) {
             	iMax = Vertices.pXadj[i];
	    }
      	}
      	if (iMax > iMaxConnectionsPerVertex) {
            iMaxConnectionsPerVertex = iMax;
            SAFEDELETEARR(Vertices.pAdjncy);
      	} else {
            break;
      	}
    }
    
    for (int i = 1; i <= iTotVertices; i++) {
       	Vertices.pXadj[i] += Vertices.pXadj[i-1];
    }
    /* Compatta il vettore delle adiacenze */
    Pack(Vertices.pAdjncy, iTotVertices*iMaxConnectionsPerVertex);
    
    /* Chiamo la routine per ottere la partizione fra i diversi processi METIS.
     * Se ne usano due diverse a seconda della dimensione della partizione */
    int edgecut; /* e' un dato fornito in output da METIS non usato */
    
    if (pParAmgProcs == NULL) {   
      	SAFENEWARR(pParAmgProcs, int, iTotVertices);

#ifdef DEBUG
      	ofstream ofMetis;
      	if (MyRank == 0) {
      	    ofMetis.open("metis_conn.debug");
      	    ofMetis << "# METIS Input File" << std::endl
                << "Column 1 is for Computational weights " << std::endl
            	<< "Column 2 is for Comunicational weight " << std::endl
            	<< "Total Vertexes: " << iTotVertices << std::endl
            	<< "# Nodes" << std::endl;
      	    for (int i = 0; i < iTotNodes; i++) {
            	ofMetis << "# " << i << "  Node Type: "
          	    << "(" << psNodeNames[ppNodes[i]->GetNodeType()] << ")"
          	    << " Label: " << ppNodes[i]->GetLabel() << std::endl
          	    << pVertexWgts[i] << " " << pCommWgts[i];
            	for (int j = Vertices.pXadj[i]; j < Vertices.pXadj[i+1]; j++) {
          	    ofMetis << " " << Vertices.pAdjncy[j];
            	}
            	ofMetis << std::endl;
      	    }
            ofMetis << "# Elements" << std::endl;
      	    for (int i = 0; i < iTotElem; i++) {
            	ofMetis << "# " << i+iTotNodes << "  Element Type: "
          	    << "("  << psElemNames[ppElems[i]->GetElemType()] << ")"
          	    << " Label: " << ppElems[i]->GetLabel() << std::endl
          	    << pVertexWgts[i+iTotNodes] 
		    << " " << pCommWgts[i+iTotNodes];
            	for (int j = Vertices.pXadj[i+iTotNodes]; 
                     j < Vertices.pXadj[i+iTotNodes+1]; 
	             j++) {
          	    ofMetis << " " << Vertices.pAdjncy[j];
	    	}
            	ofMetis << std::endl;
      	    }
            ofMetis.close();
    	}
#endif /* DEBUG */
    
#ifdef USE_METIS
    	int numflag = 0;
    	int options = 0;      

    	METIS_PartGraphVKway(&iTotVertices,
                             Vertices.pXadj,
			     Vertices.pAdjncy,
			     pVertexWgts,
			     pCommWgts,
			     &wgtflag,
			     &numflag,
			     &DataCommSize,
			     &options,
			     &edgecut,
			     pParAmgProcs);    
#else /* !USE_METIS */
    	std::cerr <<"Sorry. You need to compile with -DUSE_METIS." << std::endl 
    	    << "No other partition library is implemented yet."
	    " Aborting ..." << std::endl;
#endif /* !USE_METIS */    
    }
 
    int MyDim = iTotVertices/DataCommSize; /* Stima del # vertici per blocco */
    /* Lista dei nodi appartenenti a questo processo */
    Adjacency InterfNodes; /* nodi di interfaccia */

    /* Anche qui si usa un ciclo while per verificare 
     * se il numero di nodi di interfaccia per processo
     * ipotizzato, iMaxInterfNodes, è sufficiente */
    integer iMaxInterfNodes = iDefaultMaxInterfNodes;
    InterfNodes.pXadj = NULL;
    InterfNodes.pAdjncy = NULL;
    SAFENEWARR(InterfNodes.pXadj, int, DataCommSize+1);
    SAFENEWARR(ppMyNodes, Node*, 2*MyDim);
    
    while (true) {
        InitList(InterfNodes.pXadj, DataCommSize+1, 0);
      	SAFENEWARR(InterfNodes.pAdjncy,int, DataCommSize*iMaxInterfNodes*8);
      	InitList(InterfNodes.pAdjncy, DataCommSize*iMaxInterfNodes*8, ADJ_UNDEFINED);

      	iNumLocNodes = 0;
      	iNumLocDofs = 0;
      	iCount = 0;
      	int TmpPrc = 0;
      	flag fIsNInterf;
      	for (unsigned int i = 0; i < iTotNodes; i++) {
            fIsNInterf = flag(1);

	    if (pParAmgProcs[i] == MyRank) {
	  	/* se uno dei nodi è connesso ad un elemento non appartenente 
	   	 * a questo processo e' un nodo di interfaccia */
	  	for (int j = Vertices.pXadj[i]; j < Vertices.pXadj[i+1]; j++) {
	    	    if ((TmpPrc = pParAmgProcs[Vertices.pAdjncy[j]]) != MyRank) {
	      		InterfNodes.pXadj[TmpPrc] += 1;
	      		if (TmpPrc * iMaxInterfNodes * 2 + InterfNodes.pXadj[TmpPrc] - 1
	          	    < DataCommSize * iMaxInterfNodes * 2) {
	        	    InterfNodes.pAdjncy[TmpPrc * iMaxInterfNodes * 2 
		  		+ InterfNodes.pXadj[TmpPrc] - 1] = i;
	        	    fIsNInterf = flag(0);
	      		}
	    	    }
	  	}
	  	if (fIsNInterf) {
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
      	if (iMax > iMaxInterfNodes) {
            iMaxInterfNodes = iMax;
            SAFEDELETEARR(InterfNodes.pAdjncy);
      	} else {
            break;
      	}
    }
    
    /* Scambio i dati riguardo ai nodi di interfaccia */
    
    MPI::Request* pRReq = NULL;
    MPI::Request* pSReq = NULL;
    SAFENEWARR(pRReq, MPI::Request, DataCommSize);
    SAFENEWARR(pSReq, MPI::Request, DataCommSize);
    for (int i=0; i < DataCommSize; i++) { 
    	pRReq[i] = MPI::REQUEST_NULL;
    	pSReq[i] = MPI::REQUEST_NULL;
    }
    const int DIM_TAG = 10;
    
    for (int i = 0; i <= DataCommSize-1; i++) {
      	if(i != MyRank) {
            pRReq[i] = DataComm.Irecv(InterfNodes.pAdjncy + iMaxInterfNodes + i*iMaxInterfNodes*2, iMaxInterfNodes, MPI::INT, i, DIM_TAG);
	    pSReq[i] = DataComm.Isend(InterfNodes.pAdjncy + i*iMaxInterfNodes*2, iMaxInterfNodes, MPI::INT, i, DIM_TAG);
      	}
    }  
 
    /* lista degli elementi appartenenti a questo processo */
    SAFENEWARR(ppMyElems, Elem*, 2*MyDim);
    
    /* Trattamento elementi particolari */
    int move = 0;
    /* Gravity */
#if 0
    Elem* pTmpElem;
#endif /* 0 */
    if (ElemData[Elem::GRAVITY].iNum != 0) {
        ppMyElems[iNumLocElems] = ppElems[GravityPos];
      	iNumLocElems += 1;
      	pParAmgProcs[GravityPos+iTotNodes] = -1;
      	move++;
    }
    
    /* Air Properties */
    if (ElemData[Elem::AIRPROPERTIES].iNum != 0) {
      	ppMyElems[iNumLocElems] = ppElems[AirPropPos];
      	iNumLocElems += 1;
      	pParAmgProcs[AirPropPos+iTotNodes] = -1;
      	move++;
    }
    
    /* Rotors */    
    int iMyTotRot = 0;
    integer* pMyRot = NULL;
    integer iRotorIsMine = 0;
    if (iNumRt  != 0) {
      	SAFENEWARR(pMyRot, integer, iNumRt);
      	for (unsigned int i = 0; i < ElemData[Elem::AERODYNAMIC].iNum; i++) {
	    const AerodynamicElem *pAero = 
	    	(ElemData[Elem::AERODYNAMIC].ppFirstElem[i])->pGetAerodynamicElem();
	    ASSERT(pAero != NULL);
	    const Rotor *pRotor = pAero->pGetRotor();

	    if (pRotor != NULL) {
	        int pos = 0;
	        unsigned int pTmpLab = pRotor->GetLabel();
	  	for (int k = 0; k < iNumRt; k++) {
	    	    if (pTmpLab == pRotLab[k]) {
	      		pos  = pRotPos[k];
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
      	SAFENEWARR(pRotorComm, MPI::Intracomm, iNumRt);
	int color, key;
      	for (int i = 0; i < iNumRt; i++) {
	    if (iMyTotRot == 0) {
	  	color = MPI::UNDEFINED;
	  	key = MyRank;
	    } else {
	  	color = 1;
	  	if (pParAmgProcs[pMyRot[i]+iTotNodes] == MyRank) {
	    	    silent_cout("Rotor " << ppElems[pRotPos[i]]->GetLabel() 
	      	                << " assigned to process " << MyRank << std::endl);
	    	    iRotorIsMine = 1;
	    	    key = 0;
	  	} else {
	    	    key = MyRank+1;
	  	}
	    }
	    pRotorComm[i] = MBDynComm.Split(color,key);
	    /* RotorComm[i] = MPI::COMM_WORLD.Split(color,key); */
	    Rotor *r = (Rotor *)ppElems[pRotPos[i]]->pGet();
	    r->InitializeRotorComm(pRotorComm + i);
	
      	}
        for (int i = 0; i < iMyTotRot; i++) {
	    ppMyElems[iNumLocElems] = ppElems[pMyRot[i]];
	    if (pParAmgProcs[pMyRot[i]+iTotNodes] == MyRank) { 
	  	iNumLocDofs += (ppMyElems[iNumLocElems])->iGetNumDof();
	    } else { 
		move++;
	    }
	    iNumLocElems += 1;
	    pParAmgProcs[pMyRot[i]+iTotNodes] = -1;
      	}
    }

    for (int i = iTotNodes; i < iTotVertices; i++) {
      	if (pParAmgProcs[i] == MyRank) {
	    ppMyElems[iNumLocElems] = ppElems[i-iTotNodes];
	    iNumLocDofs += (ppMyElems[iNumLocElems])->iGetNumDof();
	    iNumLocElems += 1;
      	}
    }

    /* Verifico la ricezione dei nodi di interfaccia */
    flag RecvFlag, SentFlag;
    while (true) {
      	RecvFlag = MPI::Request::Testall(DataCommSize, pRReq);
      	SentFlag = MPI::Request::Testall(DataCommSize, pSReq);
      	if (RecvFlag && SentFlag) {
	    break;
      	} else {
	    MYSLEEP(1000);
      	}
    }
    
    /* ordino  i nodi interfaccia */
    std::sort(InterfNodes.pAdjncy, 
         InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);
    /* elimino le ripetizioni */
    int* p = std::unique(InterfNodes.pAdjncy, 
                    InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);
    /* dimensione effettiva dell'interfaccia locale 
     * il -1 e' dovuto al primo valore che è sempre pari a -1 */
    iNumIntNodes = p - InterfNodes.pAdjncy;

    unsigned int* llabels = NULL;
    Node::Type* lTypes = NULL;
    SAFENEWARR(llabels, unsigned int, iNumIntNodes);
    SAFENEWARR(lTypes,  Node::Type, iNumIntNodes);
    SAFENEWARR(ppIntNodes, Node*, iNumIntNodes); 
    for (int i = 0; i < iNumIntNodes-1; i++) {
      	ppIntNodes[i] = ppNodes[InterfNodes.pAdjncy[i+1]];
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
    }
    
    for (int i = 0; i < iNumLocElems; i++) {
      	if (ppMyElems[i]->iGetNumDof() != 0) {
            Elem::Type CType = (ppMyElems[i])->GetElemType();
            switch (CType) {
            case Elem::ROTOR:
          	if (iRotorIsMine == 1) {
            	    (ppMyElems[i])->GetConnectedNodes(iNumberOfNodes, pMyTypes, pMyLabels);
                    for (int j = 0; j < iNumberOfNodes; j++) {
              	    	unsigned int* p = 
			    std::find(llabels, llabels+iNumIntNodes, pMyLabels[j]);
              		if (p != llabels+iNumIntNodes) {
                	    ppMyIntElems[iNumIntElems] =  ppMyElems[i];
                	    iNumIntDofs += ppMyElems[i]->iGetNumDof();
                	    iNumLocDofs -= ppMyElems[i]->iGetNumDof();
                	    pPosIntElems[iNumIntElems] = i;
                	    iNumIntElems++;
                	    j = iNumberOfNodes;
              		}
            	    }
          	}
                break;

            case Elem::HYDRAULIC:
            case Elem::GENEL:
            case Elem::ELECTRIC:
          	break;

            default:
          	(ppMyElems[i])->GetConnectedNodes(iNumberOfNodes, 
						  pMyTypes, 
						  pMyLabels);
                for (int j = 0; j < iNumberOfNodes; j++) {
            	    unsigned int* p = 
		    	std::find(llabels, llabels+iNumIntNodes, pMyLabels[j]);
            	    if (p != llabels+iNumIntNodes) {
              		ppMyIntElems[iNumIntElems] =  ppMyElems[i];
              		iNumIntDofs += ppMyElems[i]->iGetNumDof();
              		iNumLocDofs -= ppMyElems[i]->iGetNumDof();
              		pPosIntElems[iNumIntElems] = i;
              		iNumIntElems++;
              		j = iNumberOfNodes;
            	    }
	  	}
	  	break;
	    }
      	}
    }
    
    /* determina la liste dei dofs locali ed adiacenti suddivisa per processi,
     * secondo la struttura Adjacency */
    SAFENEWARR(LocalDofs, integer ,iNumLocDofs);
    
    iCount = 0;
    for (int i = 0; i < iNumLocNodes; i++) {
      	if ((ppMyNodes[i])->iGetNumDof() != 0) {
            integer First = (ppMyNodes[i])->iGetFirstIndex();
            LocalDofs[iCount] = First+1;
            iCount++;
            for (unsigned int j = 1; j < (ppMyNodes[i])->iGetNumDof(); j++) {
          	LocalDofs[iCount] = First + j + 1;
          	iCount++;
            }
      	}
    }
    
    /* Aggiungo i dofs degli elementi locali */
    integer TmpDofNum;
    int i2Count = 0;
    for (int i = move; i < iNumLocElems; i++) {
      	if ((TmpDofNum = (ppMyElems[i])->iGetNumDof()) != 0) {
	    if (i != pPosIntElems[i2Count]) { 
	  	ElemWithDofs* pWithDofs = (ppMyElems[i])->pGetElemWithDofs();
	  	integer First = (pWithDofs)->iGetFirstIndex();
	  	LocalDofs[iCount] = First+1;
	  	iCount++;
	  	for (int j = 1; j < TmpDofNum; j++) {
	    	    LocalDofs[iCount] = First + 1 + j;
	    	    iCount++;
	  	}
	    } else {
	  	i2Count++;
	    }
      	}
    }

    /* scrivo ora la lista dei dofs interfaccia */    
    for (int i = 1; i < iNumIntNodes; i++) {
      	iNumIntDofs += (ppNodes[InterfNodes.pAdjncy[i]])->iGetNumDof();
    }

    SAFENEWARR(LocalIntDofs, integer, iNumIntDofs);
    SAFENEWARR(pMyIntDofs, integer, iNumIntDofs);
    
    iCount = 0;
    i2Count = 0;
    for (int i = 1; i < iNumIntNodes; i++) {
      	if ((ppNodes[InterfNodes.pAdjncy[i]])->iGetNumDof() != 0) {
            integer First = (ppNodes[InterfNodes.pAdjncy[i]])->iGetFirstIndex();
            LocalIntDofs[iCount] = (First + 1);
            iCount++;
            if (pParAmgProcs[InterfNodes.pAdjncy[i]] == MyRank) {
          	pMyIntDofs[i2Count] = (First + 1);
	  	i2Count++;
	    }
	    for (unsigned int j = 1; 
	         j < (ppNodes[InterfNodes.pAdjncy[i]])->iGetNumDof(); 
	         j++) {
	  	/* il - serve a distinguere questi dofs da quelli interni */
	  	LocalIntDofs[iCount] = (First  + 1 + j);
	  	iCount++;
	  	if (pParAmgProcs[InterfNodes.pAdjncy[i]] == MyRank) {
	    	    pMyIntDofs[i2Count] = (First  + 1 + j);
	    	    i2Count++;
	  	}
	    }
      	}
    }

    /* Interfaccia degli elementi locali */
    for (int i = 0; i < iNumIntElems; i++) {
      	TmpDofNum = (ppMyIntElems[i])->iGetNumDof();
      	ElemWithDofs* pWithDofs = (ppMyIntElems[i])->pGetElemWithDofs();
      	integer First = (pWithDofs)->iGetFirstIndex();
      	LocalIntDofs[iCount] = (First + 1);
      	iCount++;
      	for (int j = 1; j < TmpDofNum; j++) {
	    LocalIntDofs[iCount] = (First  + 1 + j);
	    iCount++;
      	}
    }

    iNumMyInt = i2Count;
    iNumIntNodes = iNumIntNodes-1;

    if (pMyTypes != NULL) {  
      	SAFEDELETE(pMyTypes);
    }
    if (pMyLabels != NULL) {  
      	SAFEDELETE(pMyLabels);
    }
    if ( Vertices.pXadj != NULL) {
      	SAFEDELETEARR(Vertices.pXadj);
    }
    if ( Vertices.pAdjncy != NULL) {
     	SAFEDELETEARR(Vertices.pAdjncy);
    }
    if ( InterfNodes.pXadj != NULL) {
     	SAFEDELETEARR( InterfNodes.pXadj);
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
    if ( pVertexWgts != NULL) {
      	SAFEDELETEARR(pVertexWgts);
    }
    if ( pCommWgts != NULL) {
      	SAFEDELETEARR(pCommWgts);
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
    silent_cout("Making partition file ...");
  
    /* Inizializzazione */
    OutHdl.PartitionOpen();
  
    time_t tCurrTime(time(NULL));
    OutHdl.Partition() 
    	<< "# Partition file for Mbdyn. Time: " << ctime(&tCurrTime)  << std::endl
	<< "# Partition produced with METIS" << std::endl << std::endl << std::endl;
  
    /* Dati */
  
    OutHdl.Partition()
    	<< "# Control data useful to verify the partition "
    	<< std::endl << std::endl 
    	<< "Total number of processes: " << DataCommSize << ";" << std::endl
    	<< " Process #: " << MyRank << ";" << std::endl
    	<< "Total number of Nodes: " << iTotNodes << ";"  << std::endl
    	<< "Total number of Elements: " << iTotElem << ";"  << std::endl
    	<< "Total number of Dofs: " << iTotDofs << ";"  << std::endl
    	<< std::endl << std::endl;

    OutHdl.Partition() << " Local Dofs number: " << iNumLocDofs << std::endl;
    
    OutHdl.Partition() << " Local Nodes: " << iNumLocNodes << std::endl;
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	OutHdl.Partition() 
            << "Node Type: "  
            << "(" << psNodeNames[(ppMyNodes[i])->GetNodeType()] << ")"
            << " Label: " << (ppMyNodes[i])->GetLabel()
            << " Dofs #: " << (ppMyNodes[i])->iGetNumDof()
            << std::endl;
    }

    OutHdl.Partition() 
    	<< std::endl << std::endl << "Local Elements: "<< iNumLocElems << std::endl << std::endl;
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	OutHdl.Partition() 
            << "Element Type: " 
      	    << "("  << psElemNames[(ppMyElems[i])->GetElemType()] << ")"
      	    << " Label: " << (ppMyElems[i])->GetLabel()
      	    << " Dofs #: " << (ppMyElems[i])->iGetNumDof()
      	    << std::endl;
    }

    OutHdl.Partition() 
    	<< std::endl << std::endl << " Interface Dofs number: " << iNumIntDofs << std::endl;
    OutHdl.Partition() 
    	<< std::endl << "Interface Nodes: " << iNumIntNodes << std::endl;
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	OutHdl.Partition() 
      	    << "Node Type: " 
      	    << "(" << psNodeNames[(ppIntNodes[i])->GetNodeType()] << ")" 
      	    << " Label: " << (ppIntNodes[i])->GetLabel()
      	    << " Dofs #: " <<  (ppIntNodes[i])->iGetNumDof()
      	    << std::endl;
    }
  
    OutHdl.Partition() 
    	<< std::endl << std::endl << "Elements whose internal dofs are interface dofs: "
	<< iNumIntElems << std::endl << std::endl;
    for (int i = 0; i < iNumIntElems; i++) {
    	OutHdl.Partition() 
      	    << "Element Type: " 
      	    << "("  << psElemNames[(ppMyIntElems[i])->GetElemType()] << ")"
      	    << " Label: " << (ppMyIntElems[i])->GetLabel()
      	    << " Dofs #: " << (ppMyIntElems[i])->iGetNumDof()
      	    << std::endl;
    }

#ifdef DEBUG
    OutHdl.Partition() << std::endl << "Local Dofs List:" << std::endl;
    int j = 0;
    for (int i = 0; i < iNumLocDofs; i++) {
    	OutHdl.Partition() << LocalDofs[i] << " ";
    	j++;
    	if (j > 10) {
      	    OutHdl.Partition() << std::endl;
      	    j = 0;
    	}
    }
    OutHdl.Partition() << std::endl;

    OutHdl.Partition() << std::endl << "Local Interface Dofs List:" <<std::endl;
    j = 0;
    for (int i = 0; i < iNumIntDofs; i++) {
	OutHdl.Partition() << LocalIntDofs[i] << " ";
    	j++;
    	if (j > 10) {
      	    OutHdl.Partition() << std::endl;
      	    j = 0;
    	}
    }
    OutHdl.Partition() << std::endl;
  
#endif /* DEBUG */   

    silent_cout("done" << std::endl);
}

/* compatta il vettore delle adiacenze */
void
SchurDataManager::Pack(int* pList, int dim)
{
    int iCount = 0;
    int* pOld = pList;
    int* pNew = pList;
    for (; pOld < pList + dim; pOld++) {
    	if (*pOld != ADJ_UNDEFINED) {
      	    *pNew++ = *pOld;
      	    iCount++;
    	}
    }
}

/* Inizializza le varie liste di interi usate in Create Partition */
void
SchurDataManager::InitList(int* list, int dim, int value)
{ 
    for (int i = 0; i <= dim-1; i++) {
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
SchurDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
    DEBUGCOUT("Entering SchurDataManager::AssRes()" << std::endl);
   
    /* Vedi quanto scritto per lo jacobiano */   
    ASSERT(iWorkIntSize >= iWorkDoubleSize);
    MySubVectorHandler WorkVec(iWorkDoubleSize, piWorkIndex, pdWorkMat);
    for (Elem** ppTmpEl = ppMyElems; 
         ppTmpEl < ppMyElems+iNumLocElems; 
	 ppTmpEl++) {
        ResHdl += (*ppTmpEl)->AssRes(WorkVec, dCoef, *pXCurr, *pXPrimeCurr);
    }
}

void 
SchurDataManager::AssJac(MatrixHandler& JacHdl, doublereal dCoef)
{
    DEBUGCOUT("Entering SchurDataManager::AssJac()" << std::endl);
    ASSERT(pWorkMatA != NULL);   
    ASSERT(ppElems != NULL);

    for (Elem** ppTmpEl = ppMyElems; 
         ppTmpEl < ppMyElems+iNumLocElems; 
	 ppTmpEl++) {
    	JacHdl += (*ppTmpEl)->AssJac(*pWorkMatA, dCoef, *pXCurr, *pXPrimeCurr);
    }
}	   
/* End of AssJac */

void
SchurDataManager::Update(void) const
{  
    DEBUGCOUT("Entering SchurDataManager::Update()" << std::endl);

    /* Nodi */   
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	(ppMyNodes[i])->Update(*pXCurr, *pXPrimeCurr);
    }

    /* Nodi di interfaccia */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	(ppIntNodes[i])->Update(*pXCurr, *pXPrimeCurr);
    }
 
    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->Update(*pXCurr, *pXPrimeCurr);
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
    	if ((ppMyNodes[i])->GetNodeType() == Node::STRUCTURAL) {
      	    ((StructNode*)ppMyNodes[i])->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
    	} else {	 
      	    (ppMyNodes[i])->Update(*pXCurr, *pXPrimeCurr);
    	}
    }

    /* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	if ((ppIntNodes[i])->GetNodeType() == Node::STRUCTURAL) {
      	    ((StructNode*)ppIntNodes[i])->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
     	} else {	 
      	    (ppIntNodes[i])->Update(*pXCurr, *pXPrimeCurr);
      	}
    }

    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->Update(*pXCurr, *pXPrimeCurr);
    }        
}
/* End of DerivativeUpdate */


void 
SchurDataManager::BeforePredict(VectorHandler& X,
                                VectorHandler& XP,
				VectorHandler& XPrev, 
				VectorHandler& XPPrev) const
{
    DEBUGCOUT("Entering SchurDataManager::BeforePredict()" << std::endl);

    /* Nodi */   
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	(ppMyNodes[i])->BeforePredict(X, XP, XPrev, XPPrev);
    }
  
    /* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	(ppIntNodes[i])->BeforePredict(X, XP, XPrev, XPPrev);
    }
  
    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->BeforePredict(X, XP, XPrev, XPPrev);
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
    	(ppMyNodes[i])->AfterPredict(*(VectorHandler*)pXCurr,
				     *(VectorHandler*)pXPrimeCurr);
    }
  
    /* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	(ppIntNodes[i])->AfterPredict(*(VectorHandler*)pXCurr,
				      *(VectorHandler*)pXPrimeCurr);
    }
  
    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->AfterPredict(*(VectorHandler*)pXCurr,
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
    	(ppMyNodes[i])->AfterConvergence(*(VectorHandler*)pXCurr,
				     *(VectorHandler*)pXPrimeCurr);
    }
  
    /* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	(ppIntNodes[i])->AfterConvergence(*(VectorHandler*)pXCurr,
				      *(VectorHandler*)pXPrimeCurr);
    }
  
    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->AfterConvergence(*(VectorHandler*)pXCurr,
				     *(VectorHandler*)pXPrimeCurr);
    }        
}
/* End of AfterConvergence */


 /* stampa i risultati */
void
SchurDataManager::Output(void) const
{
    DEBUGCOUT("Entering SchurDataManager::Output()" << std::endl);

    /* Dati intestazione */
    OutHdl.Output() 
     	<< "Time: " 
	<< std::setw(16) << std::setprecision(8) << DrvHdl.dGetTime() << std::endl;

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
  
    /* Restart condizionato */
    switch (RestartEvery) {
    case ITERATIONS:
    	if (++((integer&)iCurrRestartIter) == iRestartIterations) {
      	    (integer&)iCurrRestartIter = 0;
      	    ((SchurDataManager*)this)->MakeRestart();
    	}
    	break;
  
    case TIME:
    	ASSERT(pTime != NULL);
    	if (pTime->GetVal().GetReal()-dLastRestartTime >= dRestartTime) {
      	    (doublereal&)dLastRestartTime = pTime->GetVal().GetReal();
      	    ((SchurDataManager*)this)->MakeRestart();
    	}	   
    	break;
	
    default:
    	break;
    }
}

/* End Output */

/* End SchurDataManager */

#endif /* USE_MPI */

