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

/* Schur Data Manager */

/* 
 * Copyright 1999-2000 Giuseppe Quaranta <giuquaranta@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_MPI

#include <schurdataman.h>
#include <except.h>

/* libreria per il calcolo delle partizioni */
#ifdef USE_METIS
extern "C" {
#include <metis.h>
}
#endif /* USE_METIS */

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

 /* Costruttore - begin */
SchurDataManager::SchurDataManager(MBDynParser& HP, 
				   doublereal dInitialTime,
				   const char* sInputFileName,
				   const char* sOutputFileName, 
				   flag fAbortAfterInput)
: DataManager(HP, dInitialTime, sInputFileName, 
              sOutputFileName,fAbortAfterInput),
iTotVertexs(0),
pLabelsList(NULL),
ppMyNodes(NULL),
iNumLocNodes(0),
ppMyElems(NULL),
iNumLocElems(0),
ppMyIntElems(NULL),
iNumIntElems(0),
LocalDofs(NULL),
iNumLocDofs(0),
LocalIntDofs(NULL),
iNumIntDofs(0),
ppIntNodes(NULL),
iNumIntNodes(0),
iNumMyInt(0),
pMyIntDofs(NULL),
iTotalExpConnections(0),
ppExpCntNodes(NULL),
ppExpCntElems(NULL),
pParAmgProcs(NULL),
wgtflag(1),
pRotorComm(NULL)
{
    DEBUGCOUT("Entering SchurDataManager" << endl);

    /* Inizializza il communicator */ 
    DataComm = MPI::COMM_WORLD.Dup();
    DataCommSize = DataComm.Get_size();
    MyRank = DataComm.Get_rank();

    DEBUGCOUT("Communicator Size: " << DataCommSize << endl);

    iTotVertexs = iTotNodes + iTotElem;
    DEBUGCOUT("iTotVertexes: " << iTotVertexs << endl);

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
        cerr << endl 
	     << "Warning: no explicit connection declared"
	     " for this input file " << endl;
        return;
    } else {
        int iNumElems = 0;
        int iNumNodes = 0;
        if (KeyWords(HP.GetWord()) != PARALLEL) {
            cerr << endl 
	        << "Error: <begin: parallel;> expected at line "
	        << HP.GetLineData() << "; aborting ..." << endl;
            THROW(ErrGeneric());
        }
        while (1) {
            KeyWords CDesc = KeyWords(HP.GetDescription());
            switch (CDesc) {
            case WEIGHTS:
	        if (HP.fIsArg()) {
	            wgtflag = HP.GetInt();
	        } else {
	            cerr << endl << "Error: Number expected at line "  
	                << HP.GetLineData() << "; aborting ..." << endl;
	            THROW(ErrGeneric());
	        }
	        break;
	
            case PARTITION:
      	        SAFENEWARR(pParAmgProcs, int, iTotVertexs, DMmm);	
      	        for (int i = 0; i < iTotVertexs; i++) {
	            if (HP.fIsArg()) {
	                pParAmgProcs[i] = HP.GetInt();
	            } else {
	                cerr << endl 
			    << "Error: the partition assignment"
			    " is not complete at line "  
	                    << HP.GetLineData() << "; aborting ..." << endl;
	                THROW(ErrGeneric());
	            }
	        }  
	        break;				
	
            case END:   
	        if (KeyWords(HP.GetWord()) != PARALLEL) {
	            cerr << endl 
		        << "Error: <end: parallel;> expected at line "
	                << HP.GetLineData() << "; aborting ..." << endl;
	            THROW(ErrGeneric());
	        }
	        goto endcycle;
    
            default:
	        cerr << "Unknown input at line "
	            << HP.GetLineData() << "; aborting ..." << endl;
	        THROW(ErrGeneric());
    
            case NUMBEROFCONNECTIONS: 
	        iTotalExpConnections = HP.GetInt();
	        SAFENEWARR(ppExpCntNodes, Node*, iTotalExpConnections, DMmm);
	        SAFENEWARR(ppExpCntElems, Elem*, iTotalExpConnections, DMmm);
      
	        ElemType::Type ActualElType;
	        NodeType::Type ActualNdType;
	        unsigned int j = 0;
      
	        for (int i = 0; i < iTotalExpConnections; i++) {
	            if (KeyWords(HP.GetDescription()) != CONNECTION) {
	                cerr << endl 
			    << "Error: <Connection> expected at line "
		            << HP.GetLineData() << "; aborting ..." << endl;
	                THROW(ErrGeneric());
	            }
	  
	            for (int k = 0; k < 2; k++) {
	                KeyWords CurrDesc = KeyWords(HP.GetWord());
	                switch (CurrDesc) {
	                case ELEMENT: {
	                    KeyWords ElDesc = KeyWords(HP.GetWord());
	                    switch (ElDesc) {
	                    case FORCE:
		                ActualElType = ElemType::FORCE;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		                    THROW(ErrGeneric());
		                }
		                ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
		                if (ppExpCntElems[iNumElems] == NULL ) {
		                    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		                    THROW(ErrGeneric());
		                }
		                iNumElems++;
		                break;
				
	                    case RIGIDBODY:
		                ActualElType = ElemType::BODY;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		                    THROW(ErrGeneric());
		                }
		                ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		                    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		                    THROW(ErrGeneric());
		                }
		                iNumElems++;
		                break;
		
	                    case JOINT:
		                ActualElType = ElemType::JOINT;
		
		                if (HP.fIsArg()) {
		                    j = HP.GetInt();
		                } else {
		                    cerr << endl 
				        << "Error: Label expected at line "  
				        << HP.GetLineData() 
				        << "; aborting ..." << endl;
				    THROW(ErrGeneric());
			       	}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
			  	    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case BEAM:
				ActualElType = ElemType::BEAM;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				      	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;

	      		    case PLATE:
				ActualElType = ElemType::PLATE;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ROTOR:
				ActualElType = ElemType::ROTOR;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case AERODYNAMICELEMENT:
				ActualElType = ElemType::AERODYNAMIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
			  	    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ELECTRICBULK:
				ActualElType = ElemType::ELECTRICBULK;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case ELECTRIC:
				ActualElType = ElemType::ELECTRIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case GENEL:
				ActualElType = ElemType::GENEL;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case HYDRAULIC:
				ActualElType = ElemType::HYDRAULIC;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case BULK:
				ActualElType = ElemType::BULK;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
					<< HP.GetLineData() 
					<< "; aborting ..." << endl;
			  	    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case LOADABLE:
				ActualElType = ElemType::LOADABLE;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				        << "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
				    cerr << "Error: at line " 
				     	<< HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
				    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    case DRIVEN:
				ActualElType = ElemType::DRIVEN;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				      	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntElems[iNumElems] = 
				    *(ppFindElem(ActualElType, j));
				if (ppExpCntElems[iNumElems] == NULL ) {
		  		    cerr << "Error: at line " 
				        << HP.GetLineData() 
					<< " undefined element; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumElems++;
				break;
		
	      		    default:
				cerr << "Error: Element Type on line " 
				   << HP.GetLineData() 
				   << " not valid; aborting ..." << endl;
				THROW(ErrGeneric());
	      		    }
	      		    break;
	    		}
	    
	    		case NODE: {
	      		    KeyWords NdDesc = KeyWords(HP.GetWord());
	      		    switch (NdDesc) {		
	      		    case ABSTRACT:
				ActualNdType = NodeType::ABSTRACT;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				     	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    case STRUCTURAL:
				ActualNdType = NodeType::STRUCTURAL;
		
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
				    cerr << endl 
				      	<< "Error: Label expected at line "  
		     			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
				    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    case ELECTRIC:
				ActualNdType = NodeType::ELECTRIC;
				
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				      	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	    		    case PARAMETER:
	      			ActualNdType = NodeType::PARAMETER;
	      			if (HP.fIsArg()) {
				    j = HP.GetInt();
	      			} else {
				    cerr << endl 
				    	<< "Error: Label expected at line "  
		     			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
				    THROW(ErrGeneric());
	      			}
	      			ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
	      			if (ppExpCntNodes[iNumNodes] == NULL ) {
				    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< endl;
				    THROW(ErrGeneric());
	      			}
	      			iNumNodes++;
	      			break;
	      
	      		    case HYDRAULIC:
				ActualNdType = NodeType::HYDRAULIC;
				
				if (HP.fIsArg()) {
		  		    j = HP.GetInt();
				} else {
		  		    cerr << endl 
				    	<< "Error: Label expected at line "  
		       			<< HP.GetLineData() 
					<< "; aborting ..." << endl;
		  		    THROW(ErrGeneric());
				}
				ppExpCntNodes[iNumNodes] = 
				    pFindNode(ActualNdType, j);
				if (ppExpCntNodes[iNumNodes] == NULL ) {
		  		    cerr << "Error: at line " 
				    	<< HP.GetLineData() 
					<< " undefined node; aborting ..." 
					<< endl;
		  		    THROW(ErrGeneric());
				}
				iNumNodes++;
				break;
		
	      		    default:
				cerr << "Error: Node Type on line " 
				    << HP.GetLineData() << " not valid "
		     		    << "; aborting ..." << endl;
				THROW(ErrGeneric());
	      		    }  
	      		    break;
	    		}
	    
	    		default:
	      		    cerr << " Unknown input at line "
		   		<< HP.GetLineData() 
				<< "; aborting ..." << endl;
	      		    THROW(ErrGeneric());
	    		}
	  	    }
		}
   
		ASSERT(iNumNodes == iNumElems);
		if (iNumNodes + iNumElems !=  2*iTotalExpConnections) {
	  	    cerr << endl 
		     	<< "Error: Total number of Nodes and elements"
			" in the parallel section at line "
	       		<< HP.GetLineData() 
			<< " is not consistent; aborting ..." << endl;
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
    	SAFEDELETEARR(ppMyNodes, DMmm);
    }
    if (ppMyElems != NULL) {
    	SAFEDELETEARR(ppMyElems, DMmm);
    } 
    if (ppMyIntElems != NULL) {
    	SAFEDELETEARR(ppMyIntElems, DMmm);
    } 
    if (ppIntNodes != NULL) {
    	SAFEDELETEARR(ppIntNodes, DMmm);
    }
    if (LocalDofs != NULL) {
    	SAFEDELETEARR(LocalDofs, DMmm);
    }
    if (LocalIntDofs != NULL) {
    	SAFEDELETEARR(LocalIntDofs, DMmm);
    }
    if (ppExpCntNodes != NULL) {
    	SAFEDELETEARR(ppExpCntNodes, DMmm);
    }
    if (ppExpCntElems != NULL) {
    	SAFEDELETEARR(ppExpCntElems, DMmm);
    }
}

/* Ripartisce i nodi fra i processi e all'interno di ogni singolo processo */
void 
SchurDataManager::CreatePartition(void)
{ 
    const char sFuncName[] = "CreatePartition()";
    DEBUGCOUT("Entering " << sFuncName << endl);

    Adjacency Vertexs;  /* Struttura contenente le connessioni fra i vertici */
    int* pVertexWgts;   /* Pesi dei vertici = dofs x ogni v. utile per METIS */
    int* pCommWgts;
    Vertexs.pXadj = NULL;
    Vertexs.pAdjncy = NULL;
    pVertexWgts = NULL;
    pCommWgts = NULL;
    integer iMax = 0;
    integer iRMax = 0;
    int iCount = 0;
    int  iNumRt = 0;
    
    ASSERT(iTotVertexs > 0);
    ASSERT(DataCommSize > 0);
    /* Costruisco e inizializzo le due strutture */
    SAFENEWARR(Vertexs.pXadj, int, iTotVertexs+1, DMmm);
    SAFENEWARR(pVertexWgts, int, iTotVertexs*2, DMmm);
    SAFENEWARR(pCommWgts, int, iTotVertexs, DMmm);
    
    /* Ciclo per la scrittura degli array delle connessioni. 
     * Il ciclo viene ripetuto se per un vertice si ha un numero
     * di connessioni superiore al max consentito per default
     * iDefaultMaxConnectionsPerVertex */
    int iMaxConnectionsPerVertex = 
      (iTotVertexs < iDefaultMaxConnectionsPerVertex) ? iTotVertexs 
      : iDefaultMaxConnectionsPerVertex;
    int GravityPos, AirPropPos;
    int* pRotPos = NULL;
    integer* pRotLab = NULL;
    int iNumberOfNodes;
    NodeType::Type* pMyTypes = NULL;
    unsigned int* pMyLabels = NULL;
    SAFENEWARR(pRotPos, int, ElemData[ElemType::ROTOR].iNum, DMmm);
    SAFENEWARR(pRotLab, integer, ElemData[ElemType::ROTOR].iNum, DMmm);

    SAFENEWARR(pMyTypes, NodeType::Type, iDefaultMaxNodesPerElem, DMmm);
    SAFENEWARR(pMyLabels, unsigned int, iDefaultMaxNodesPerElem, DMmm);
    SAFENEWARR(pLabelsList, unsigned int, iTotNodes, DMmm);

    while (1) {
        InitList(Vertexs.pXadj, iTotVertexs+1, 0);
        InitList(pVertexWgts, iTotVertexs*2, 0);
       	InitList(pCommWgts, iTotVertexs, 0);
       	SAFENEWARR(Vertexs.pAdjncy, int, iTotVertexs*iMaxConnectionsPerVertex, DMmm);
       	InitList(Vertexs.pAdjncy, iTotVertexs*iMaxConnectionsPerVertex, ADJ_UNDEFINED);
       	ASSERT(ppElems != NULL);
       
       	/* ciclo sugli elementi per assemblare la struttura delle connessioni */
       	Node** ppActualNode = NULL;
       	/* per numerare i nodi prendo la posizione del puntatore 
         * al nodo nell'array ppNodes */
       	int position;
       	iCount = iTotNodes;
       	for (int i = 0; i < iTotNodes; i++) {
            pLabelsList[i] = (ppNodes[i])->GetLabel();
       	}

       	iNumRt = 0;

       	for (Elem** ppTmpEl = ppElems; 
	     ppTmpEl < ppElems+iTotElem; ppTmpEl++, 
	     iCount++) {
            if ((*ppTmpEl)->GetElemType() == ElemType::GRAVITY) {
             	GravityPos = ppTmpEl - ppElems;
            }
            if ((*ppTmpEl)->GetElemType() == ElemType::AIRPROPERTIES) {
             	AirPropPos = ppTmpEl - ppElems;
            }
            if ((*ppTmpEl)->GetElemType() == ElemType::ROTOR) {
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
	     	Vertexs.pXadj[iCount+1] += 1;
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
	     	if ((iCount*iMaxConnectionsPerVertex) + Vertexs.pXadj[iCount+1] - 1 < iTotVertexs*iMaxConnectionsPerVertex) {
	            Vertexs.pAdjncy[(iCount*iMaxConnectionsPerVertex) 
		  	+ Vertexs.pXadj[iCount+1] - 1] = position;
	     	}
	     	/* aggiungo alle connessioni del nodo l'elemento attuale */
	     	Vertexs.pXadj[position+1] += 1;
	     	if ((position*iMaxConnectionsPerVertex) + Vertexs.pXadj[position+1] - 1 < iTotVertexs*iMaxConnectionsPerVertex) {
	            Vertexs.pAdjncy[(position*iMaxConnectionsPerVertex) 
		  	+ Vertexs.pXadj[position+1] - 1] = iCount;
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
             	Vertexs.pXadj[iNdPos+1] += 1;
             	Vertexs.pXadj[iTotNodes + iElPos+1] += 1;
             	Vertexs.pAdjncy[(iNdPos*iMaxConnectionsPerVertex) 
	       	    + Vertexs.pXadj[iNdPos+1] - 1] = iElPos+iTotNodes;
             	Vertexs.pAdjncy[((iTotNodes+iElPos)*iMaxConnectionsPerVertex) 
	       	    + Vertexs.pXadj[iTotNodes+iElPos+1] - 1] = iNdPos;
            }
     	}
     
        iMax = 0;
        for (int i = 0; i < iTotVertexs; i++) {
            if (Vertexs.pXadj[i] > iMaxConnectionsPerVertex) {
             	iMax = Vertexs.pXadj[i];
	    }
      	}
      	if (iMax > iMaxConnectionsPerVertex) {
            iMaxConnectionsPerVertex = iMax;
            SAFEDELETEARR(Vertexs.pAdjncy, DMmm);
      	} else {
            break;
      	}
    }
    
    for (int i = 1; i <= iTotVertexs; i++) {
       	Vertexs.pXadj[i] += Vertexs.pXadj[i-1];
    }
    /* Compatta il vettore delle adiacenze */
    Pack(Vertexs.pAdjncy, iTotVertexs*iMaxConnectionsPerVertex);
    
    /* Chiamo la routine per ottere la partizione fra i diversi processi METIS.
     * Se ne usano due diverse a seconda della dimensione della partizione */
    int edgecut; /* e' un dato fornito in output da METIS non usato */
    
    if (pParAmgProcs == NULL) {   
      	SAFENEWARR(pParAmgProcs, int, iTotVertexs, DMmm);

#ifdef DEBUG
      	ofstream ofMetis;
      	if (MyRank == 0) {
      	    ofMetis.open("metis_conn.debug");
      	    ofMetis << "# METIS Input File" << endl
                << "Column 1 is for Computational weights " << endl
            	<< "Column 2 is for Comunicational weight " << endl
            	<< "Total Vertexes: " << iTotVertexs << endl
            	<< "# Nodes" << endl;
      	    for (int i = 0; i < iTotNodes; i++) {
            	ofMetis << "# " << i << "  Node Type: "
          	    << "(" << psNodeNames[(ppNodes[i])->GetNodeType()] << ")"
          	    << " Label: " << (ppNodes[i])->GetLabel() << endl
          	    << pVertexWgts[i] << " " << pCommWgts[i];
            	for (int j = Vertexs.pXadj[i]; j < Vertexs.pXadj[i+1]; j++) {
          	    ofMetis << " " << Vertexs.pAdjncy[j];
            	}
            	ofMetis << endl;
      	    }
            ofMetis << "# Elements" << endl;
      	    for (int i = 0; i < iTotElem; i++) {
            	ofMetis << "# " << i+iTotNodes << "  Element Type: "
          	    << "("  << psElemNames[(ppElems[i])->GetElemType()] << ")"
          	    << " Label: " << (ppElems[i])->GetLabel() << endl
          	    << pVertexWgts[i+iTotNodes] 
		    << " " << pCommWgts[i+iTotNodes];
            	for (int j = Vertexs.pXadj[i+iTotNodes]; 
                     j < Vertexs.pXadj[i+iTotNodes+1]; 
	             j++) {
          	    ofMetis << " " << Vertexs.pAdjncy[j];
	    	}
            	ofMetis << endl;
      	    }
            ofMetis.close();
    	}
#endif /* DEBUG */
    
#ifdef USE_METIS
    	int numflag = 0;
    	int options = 0;      

    	METIS_PartGraphVKway(&iTotVertexs,
                             Vertexs.pXadj,
			     Vertexs.pAdjncy,
			     pVertexWgts,
			     pCommWgts,
			     &wgtflag,
			     &numflag,
			     &DataCommSize,
			     &options,
			     &edgecut,
			     pParAmgProcs);    
#else /* !USE_METIS */
    	cerr <<"Sorry. You need to compile with -DUSE_METIS." << endl 
    	    << "No other partition library is implemented yet."
	    " Aborting ..." << endl;
#endif /* !USE_METIS */    
    }
 
    int MyDim = iTotVertexs/DataCommSize; /* Stima del # vertici per blocco */
    /* Lista dei nodi appartenenti a questo processo */
    Adjacency InterfNodes; /* nodi di interfaccia */

    /* Anche qui si usa un ciclo while per verificare 
     * se il numero di nodi di interfaccia per processo
     * ipotizzato, iMaxInterfNodes, � sufficiente */
    integer iMaxInterfNodes = iDefaultMaxInterfNodes;
    InterfNodes.pXadj = NULL;
    InterfNodes.pAdjncy = NULL;
    SAFENEWARR(InterfNodes.pXadj, int, DataCommSize+1, DMmm);
    SAFENEWARR(ppMyNodes, Node*, 2*MyDim, DMmm);
    
    while (1) {
        InitList(InterfNodes.pXadj, DataCommSize+1, 0);
      	SAFENEWARR(InterfNodes.pAdjncy,int, DataCommSize*iMaxInterfNodes*8, DMmm);
      	InitList(InterfNodes.pAdjncy, DataCommSize*iMaxInterfNodes*8, ADJ_UNDEFINED);

      	iNumLocNodes = 0;
      	iNumLocDofs = 0;
      	iCount = 0;
      	int TmpPrc = 0;
      	flag fIsNInterf;
      	for (int i = 0; i < iTotNodes; i++) {
            fIsNInterf = flag(1);

	    if (pParAmgProcs[i] == MyRank) {
	  	/* se uno dei nodi � connesso ad un elemento non appartenente 
	   	 * a questo processo e' un nodo di interfaccia */
	  	for (int j = Vertexs.pXadj[i]; j < Vertexs.pXadj[i+1]; j++) {
	    	    if ((TmpPrc = pParAmgProcs[Vertexs.pAdjncy[j]]) != MyRank) {
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
	    	    iNumLocDofs += (ppMyNodes[iCount])->iGetNumDof();
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
            SAFEDELETEARR(InterfNodes.pAdjncy, DMmm);
      	} else {
            break;
      	}
    }
    
    /* Scambio i dati riguardo ai nodi di interfaccia */
    
    MPI::Request* pRReq = NULL;
    MPI::Request* pSReq = NULL;
    SAFENEWARR(pRReq, MPI::Request, DataCommSize, DMmm);
    SAFENEWARR(pSReq, MPI::Request, DataCommSize, DMmm);
    pRReq[MyRank] = MPI::REQUEST_NULL;
    pSReq[MyRank] = MPI::REQUEST_NULL;
    const int DIM_TAG = 10;
    
    for (int i = 0; i <= DataCommSize-1; i++) {
      	if(i != MyRank) {
            pRReq[i] = DataComm.Irecv(InterfNodes.pAdjncy + iMaxInterfNodes + i*iMaxInterfNodes*2, iMaxInterfNodes, MPI::INT, i, DIM_TAG);
	    pSReq[i] = DataComm.Isend(InterfNodes.pAdjncy + i*iMaxInterfNodes*2, iMaxInterfNodes, MPI::INT, i, DIM_TAG);
      	}
    }  
 
    /* lista degli elementi appartenenti a questo processo */
    SAFENEWARR(ppMyElems, Elem*, 2*MyDim, DMmm);
    
    /* Trattamento elementi particolari */
    int move = 0;
    /* Gravity */
    Elem* pTmpElem;
    if (ElemData[ElemType::GRAVITY].iNum != 0) {
        ppMyElems[iNumLocElems] = ppElems[GravityPos];
      	iNumLocElems += 1;
      	pParAmgProcs[GravityPos+iTotNodes] = -1;
      	move++;
    }
    
    /* Air Properties */
    if (ElemData[ElemType::AIRPROPERTIES].iNum != 0) {
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
      	SAFENEWARR(pMyRot, integer, iNumRt, DMmm);
      	for (int i = 0; i < ElemData[ElemType::AERODYNAMIC].iNum; i++) {
	    integer pTmpLab = 
	    	((ElemData[ElemType::AERODYNAMIC].ppFirstElem)[i])->GetRotor();
	    int pos;
	    if (pTmpLab != -1) {
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
      	SAFENEWARR(pRotorComm, MPI::Intracomm, iNumRt, DMmm);
	int color, key;
      	for (int i = 0; i < iNumRt; i++) {
	    if (iMyTotRot == 0) {
	  	color = MPI::UNDEFINED;
	  	key = MyRank;
	    } else {
	  	color = 1;
	  	if (pParAmgProcs[pMyRot[i]+iTotNodes] == MyRank) {
	    	    silent_cout("Rotor " << ppElems[pRotPos[i]]->GetLabel() 
	      	                << " assigned to process " << MyRank << endl);
	    	    iRotorIsMine = 1;
	    	    key = 0;
	  	} else {
	    	    key = MyRank+1;
	  	}
	    }
	    pRotorComm[i] = MPI::COMM_WORLD.Split(color,key);
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

    for (int i = iTotNodes; i < iTotVertexs; i++) {
      	if (pParAmgProcs[i] == MyRank) {
	    ppMyElems[iNumLocElems] = ppElems[i-iTotNodes];
	    iNumLocDofs += (ppMyElems[iNumLocElems])->iGetNumDof();
	    iNumLocElems += 1;
      	}
    }

    /* Verifico la ricezione dei nodi di interfaccia */
    flag RecvFlag, SentFlag;
    while (1) {
      	RecvFlag = MPI::Request::Testall(DataCommSize, pRReq);
      	SentFlag = MPI::Request::Testall(DataCommSize, pSReq);
      	if (RecvFlag && SentFlag) {
	    break;
#ifdef USE_MYSLEEP
      	} else {
	    mysleep(1000);
#endif /* USE_MYSLEEP */
      	}
    }
    
    /* ordino  i nodi interfaccia */
    sort(InterfNodes.pAdjncy, 
         InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);
    /* elimino le ripetizioni */
    int* p = unique(InterfNodes.pAdjncy, 
                    InterfNodes.pAdjncy + iMaxInterfNodes * 2 * DataCommSize);
    /* dimensione effettiva dell'interfaccia locale 
     * il -1 e' dovuto al primo valore che � sempre pari a -1 */
    iNumIntNodes = p - InterfNodes.pAdjncy;

    unsigned int* llabels = NULL;
    NodeType::Type* lTypes = NULL;
    SAFENEWARR(llabels, unsigned int, iNumIntNodes, DMmm);
    SAFENEWARR(lTypes,  NodeType::Type, iNumIntNodes, DMmm);
    SAFENEWARR(ppIntNodes, Node*, iNumIntNodes, DMmm); 
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
      	SAFENEWARR(pPosIntElems, int, iNumLocElems, DMmm);
      	SAFENEWARR(ppMyIntElems, Elem*, iNumLocElems, DMmm);
    }
    
    for (int i = 0; i < iNumLocElems; i++) {
      	if (ppMyElems[i]->iGetNumDof() != 0) {
            ElemType::Type CType = (ppMyElems[i])->GetElemType();
            switch (CType) {
            case ElemType::ROTOR:
          	if (iRotorIsMine == 1) {
            	    (ppMyElems[i])->GetConnectedNodes(iNumberOfNodes, pMyTypes, pMyLabels);
                    for (int j = 0; j < iNumberOfNodes; j++) {
              	    	unsigned int* p = 
			    find(llabels, llabels+iNumIntNodes, pMyLabels[j]);
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

            case ElemType::HYDRAULIC:
          	break;

            case ElemType::GENEL:
	  	break;

            case ElemType::ELECTRIC:
          	break;

            default:
          	(ppMyElems[i])->GetConnectedNodes(iNumberOfNodes, 
						  pMyTypes, 
						  pMyLabels);
                for (int j = 0; j < iNumberOfNodes; j++) {
            	    unsigned int* p = 
		    	find(llabels, llabels+iNumIntNodes, pMyLabels[j]);
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
    SAFENEWARR(LocalDofs, integer ,iNumLocDofs, DMmm);
    
    iCount = 0;
    for (int i = 0; i < iNumLocNodes; i++) {
      	if ((ppMyNodes[i])->iGetNumDof() != 0) {
            integer First = (ppMyNodes[i])->iGetFirstIndex();
            LocalDofs[iCount] = First+1;
            iCount++;
            for (int j = 1; j < (ppMyNodes[i])->iGetNumDof(); j++) {
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

    SAFENEWARR(LocalIntDofs, integer, iNumIntDofs, DMmm);
    SAFENEWARR(pMyIntDofs, integer, iNumIntDofs, DMmm);
    
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
	    for (int j = 1; 
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
      	SAFEDELETE(pMyTypes, DMmm);
    }
    if (pMyLabels != NULL) {  
      	SAFEDELETE(pMyLabels, DMmm);
    }
    if (pRReq != NULL) {  
      	SAFEDELETEARR(pRReq, DMmm);
    }
    if (pSReq != NULL) {  
      	SAFEDELETEARR(pSReq, DMmm);
    }
    if ( Vertexs.pXadj != NULL) {
      	SAFEDELETEARR(Vertexs.pXadj, DMmm);
    }
    if ( Vertexs.pAdjncy != NULL) {
     	SAFEDELETEARR(Vertexs.pAdjncy, DMmm);
    }
    if ( InterfNodes.pXadj != NULL) {
     	SAFEDELETEARR( InterfNodes.pXadj, DMmm);
    }
    if ( InterfNodes.pAdjncy != NULL) {
     	SAFEDELETEARR( InterfNodes.pAdjncy, DMmm);
    }
    if ( pParAmgProcs != NULL) {
      	SAFEDELETEARR(pParAmgProcs, DMmm);
    }
    if ( pLabelsList != NULL) {
      	SAFEDELETEARR(pLabelsList, DMmm);
    }
    if ( pVertexWgts != NULL) {
      	SAFEDELETEARR(pVertexWgts, DMmm);
    }
    if ( pCommWgts != NULL) {
      	SAFEDELETEARR(pCommWgts, DMmm);
    } 
    if ( pPosIntElems != NULL) {
      	SAFEDELETEARR(pPosIntElems, DMmm);
    }
    if ( lTypes != NULL) {
      	SAFEDELETEARR(lTypes, DMmm);
    }
    if ( llabels != NULL) {
      	SAFEDELETEARR(llabels, DMmm);
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
    	<< "# Partition file for Mbdyn. Time: " << ctime(&tCurrTime)  << endl
	<< "# Partition produced with METIS" << endl << endl << endl;
  
    /* Dati */
  
    OutHdl.Partition()
    	<< "# Control data useful to verify the partition "
    	<< endl << endl 
    	<< "Total number of processes: " << DataCommSize << ";" << endl
    	<< " Process #: " << MyRank << ";" << endl
    	<< "Total number of Nodes: " << iTotNodes << ";"  << endl
    	<< "Total number of Elements: " << iTotElem << ";"  << endl
    	<< "Total number of Dofs: " << iTotDofs << ";"  << endl
    	<< endl << endl;

    OutHdl.Partition() << " Local Dofs number: " << iNumLocDofs << endl;
    
    OutHdl.Partition() << " Local Nodes: " << iNumLocNodes << endl;
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	OutHdl.Partition() 
            << "Node Type: "  
            << "(" << psNodeNames[(ppMyNodes[i])->GetNodeType()] << ")"
            << " Label: " << (ppMyNodes[i])->GetLabel()
            << " Dofs #: " << (ppMyNodes[i])->iGetNumDof()
            << endl;
    }

    OutHdl.Partition() 
    	<< endl << endl << "Local Elements: "<< iNumLocElems << endl << endl;
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	OutHdl.Partition() 
            << "Element Type: " 
      	    << "("  << psElemNames[(ppMyElems[i])->GetElemType()] << ")"
      	    << " Label: " << (ppMyElems[i])->GetLabel()
      	    << " Dofs #: " << (ppMyElems[i])->iGetNumDof()
      	    << endl;
    }

    OutHdl.Partition() 
    	<< endl << endl << " Interface Dofs number: " << iNumIntDofs << endl;
    OutHdl.Partition() 
    	<< endl << "Interface Nodes: " << iNumIntNodes << endl;
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	OutHdl.Partition() 
      	    << "Node Type: " 
      	    << "(" << psNodeNames[(ppIntNodes[i])->GetNodeType()] << ")" 
      	    << " Label: " << (ppIntNodes[i])->GetLabel()
      	    << " Dofs #: " <<  (ppIntNodes[i])->iGetNumDof()
      	    << endl;
    }
  
    OutHdl.Partition() 
    	<< endl << endl << "Elements whose internal dofs are interface dofs: "
	<< iNumIntElems << endl << endl;
    for (int i = 0; i < iNumIntElems; i++) {
    	OutHdl.Partition() 
      	    << "Element Type: " 
      	    << "("  << psElemNames[(ppMyIntElems[i])->GetElemType()] << ")"
      	    << " Label: " << (ppMyIntElems[i])->GetLabel()
      	    << " Dofs #: " << (ppMyIntElems[i])->iGetNumDof()
      	    << endl;
    }

#ifdef DEBUG
    OutHdl.Partition() << endl << "Local Dofs List:" << endl;
    int j = 0;
    for (int i = 0; i < iNumLocDofs; i++) {
    	OutHdl.Partition() << LocalDofs[i] << " ";
    	j++;
    	if (j > 10) {
      	    OutHdl.Partition() << endl;
      	    j = 0;
    	}
    }
    OutHdl.Partition() << endl;

    OutHdl.Partition() << endl << "Local Interface Dofs List:" <<endl;
    j = 0;
    for (int i = 0; i < iNumIntDofs; i++) {
	OutHdl.Partition() << LocalIntDofs[i] << " ";
    	j++;
    	if (j > 10) {
      	    OutHdl.Partition() << endl;
      	    j = 0;
    	}
    }
    OutHdl.Partition() << endl;
  
#endif /* DEBUG */   

    silent_cout("done" << endl);
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
    unsigned int* p = find(pbegin, pbegin + dim, label);
    return ppNodes + (p - pLabelsList);
}


void
SchurDataManager::AssRes(VectorHandler& ResHdl, doublereal dCoef)
{
    DEBUGCOUT("Entering SchurDataManager::AssRes()" << endl);
   
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
    const char sFuncName[] = "SchurDataManager::AssJac()";
  
    DEBUGCOUT("Entering " << sFuncName << endl);
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
    DEBUGCOUT("Entering SchurDataManager::Update()" << endl);

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
    DEBUGCOUT("Entering SchurDataManager::DerivativesUpdate()" << endl);

    /* Nodi */   
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	if ((ppMyNodes[i])->GetNodeType() == NodeType::STRUCTURAL) {
      	    ((StructNode*)ppMyNodes[i])->DerivativesUpdate(*pXCurr, *pXPrimeCurr);
    	} else {	 
      	    (ppMyNodes[i])->Update(*pXCurr, *pXPrimeCurr);
    	}
    }

    /* Nodi adiacenti i cui valori influenzano gli assemblaggi degli elementi */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	if ((ppIntNodes[i])->GetNodeType() == NodeType::STRUCTURAL) {
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
    DEBUGCOUT("Entering SchurDataManager::BeforePredict()" << endl);

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
    DEBUGCOUT("Entering SchurDataManager::AfterPredict()" << endl);

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


 /* stampa i risultati */
void
SchurDataManager::Output(void) const
{
    DEBUGCOUT("Entering SchurDataManager::Output()" << endl);

    /* Dati intestazione */
    ((OutputHandler&)OutHdl).Output() 
     	<< "Time: " 
	<< setw(16) << setprecision(8) << DrvHdl.dGetTime() << endl;

    OutputHandler& OH = (OutputHandler&)OutHdl;
    /* Nodi */
    for (int i = 0; i < iNumLocNodes; i++) {
    	ASSERT(ppMyNodes[i] != NULL);
    	(ppMyNodes[i])->Output(OH);
    }

    /* Nodi di interfaccia */
    for (int i = 0; i < iNumIntNodes; i++) {
    	ASSERT(ppIntNodes[i] != NULL);
    	(ppIntNodes[i])->Output(OH);
    }

    /* Elementi */
    for (int i = 0; i < iNumLocElems; i++) {
    	ASSERT(ppMyElems[i] != NULL);
    	(ppMyElems[i])->Output(OH);
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

        
/* End MakeRestart */

/* End SchurDataManager */

#endif /* USE_MPI */

