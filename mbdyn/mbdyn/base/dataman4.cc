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

/* Lettura elementi */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <dataman_.h>

/* Elementi */
#ifdef USE_STRUCT_NODES
#include <autostr.h>   /* Elementi automatici associati ai nodi dinamici */
#include <gravity.h>   /* Elemento accelerazione di gravita' */
#ifdef USE_AERODYNAMIC_ELEMS
#include <aerodyn.h>   /* Classe di base degli elementi aerodinamici */
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

#include <driven.h>    /* driven elements */
#include <j2p.h>       /* bind di elementi a nodi parametrici */

#include <drive.h>
#include <tpldrive.h>
#include <tpldrive_.h>

#ifdef HAVE_LOADABLE
#include <loadable.h>
#endif /* HAVE_LOADABLE */

static int iNumTypes[ElemType::LASTELEMTYPE];

/* enum delle parole chiave */
enum KeyWords {
   UNKNOWNKEYWORD = -1,
   END = 0,
   ELEMENTS,
   
   GRAVITY,
   BODY,
   AUTOMATICSTRUCTURAL,
   JOINT,
   COUPLE,
   BEAM,
   
   AIRPROPERTIES,
   ROTOR,
   AERODYNAMICBODY,
   AERODYNAMICBEAM,
   
   FORCE,
   
   GENEL,
   ELECTRIC,
   
   HYDRAULIC,
   
   BULK,
   LOADABLE,
   DRIVEN,
   
   EXISTING,
   OUTPUT,
   BIND,
   
   LASTKEYWORD
};

void DataManager::ReadElems(MBDynParser& HP)
{
   DEBUGCOUTFNAME("DataManager::ReadElems");
   
   /* parole chiave del blocco degli elementi */
   const char* sKeyWords[] = { 
      "end",
      
      "elements",
      
      "gravity",
      "body",
      "automatic" "structural",
      "joint",
      "couple",
      "beam",
      
      "airproperties",
      "rotor",
      "aerodynamicbody",
      "aerodynamicbeam",
      
      "force",
      
      "genel",
      "electric",
      
      "hydraulic",
      
      "bulk",
      "loadable",
      "driven",
      
      "existing",
      "output",
      "bind"
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
      
   /* strutture di conteggio degli elementi letti */
   for (int i = 0; i < ElemType::LASTELEMTYPE; 
	iNumTypes[i] = ElemData[i++].iNum) { 
      NO_OP; 
   }   
         
   int iMissingElems = iTotElem;
   DEBUGLCOUT(MYDEBUG_INPUT, "Expected elements: " << iMissingElems << endl);
   
#ifdef USE_STRUCT_NODES
   /* Aggiunta degli elementi strutturali automatici legati ai nodi dinamici */
   if (ElemData[ElemType::AUTOMATICSTRUCTURAL].iNum > 0) {	
      StructNode** ppTmpNod = 
	(StructNode**)NodeData[NodeType::STRUCTURAL].ppFirstNode;
      int iTotNod = NodeData[NodeType::STRUCTURAL].iNum;
      
      Elem** ppTmpEl = 
	ElemData[ElemType::AUTOMATICSTRUCTURAL].ppFirstElem;
      for (StructNode** ppTmp = ppTmpNod; ppTmp < ppTmpNod+iTotNod; ppTmp++) {
	 if ((*ppTmp)->GetStructNodeType() == StructNodeType::DYNAMIC) {
	    
	    SAFENEWWITHCONSTRUCTOR(*ppTmpEl, AutomaticStructElem,
				   AutomaticStructElem((DynamicStructNode*)(*ppTmp)), DMmm);
	    
	    ppTmpEl++;
	    iMissingElems--;
	    iNumTypes[ElemType::AUTOMATICSTRUCTURAL]--;
	    DEBUGLCOUT(MYDEBUG_INPUT, 
		       "Initialising automatic structural element linked to node " 
		       << (*ppTmp)->GetLabel() << endl);
	 }
      }
   }
#endif /* USE_STRUCT_NODES */
   
   KeyWords CurrDesc;   
   while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
      
      if (CurrDesc == OUTPUT) {
	 DEBUGLCOUT(MYDEBUG_INPUT, "Elements to be output: ");
	 ElemType::Type Typ;
	 switch (KeyWords(HP.GetWord())) {
	    
#ifdef USE_STRUCT_NODES
	  case BODY: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "bodies" << endl);
	     Typ = ElemType::BODY;
	     break;
	  }
	  case AUTOMATICSTRUCTURAL: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "automatic structural" << endl);
	     Typ = ElemType::AUTOMATICSTRUCTURAL;
	     break;
	  }
	  case JOINT: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "joints" << endl);
	     Typ = ElemType::JOINT;
	     break;
	  }
	  case BEAM: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "beams" << endl);
	     Typ = ElemType::BEAM;
	     break;
	  }
#ifdef USE_AERODYNAMIC_ELEMS
	  case ROTOR: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "rotors" << endl);
	     Typ = ElemType::ROTOR;
	     break;
	  }
	  case AERODYNAMICBODY:
	  case AERODYNAMICBEAM: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "aerodynamic" << endl);
	     Typ = ElemType::AERODYNAMIC;
	     break;
	  }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

	  case FORCE:
#ifdef USE_STRUCT_NODES
	  case COUPLE:
#endif /* USE_STRUCT_NODES */
	  {		 
	     DEBUGLCOUT(MYDEBUG_INPUT, "forces" << endl);
	     Typ = ElemType::FORCE;
	     break;
	  }	
#ifdef USE_ELECTRIC_NODES
	  case GENEL: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "genels" << endl);
	     Typ = ElemType::GENEL;
	     break;
	  }
	  case ELECTRIC: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "electric" << endl);
	     Typ = ElemType::ELECTRIC;
	     break;
	  }
#endif /* USE_ELECTRIC_NODES */

#ifdef USE_HYDRAULIC_NODES
	  case HYDRAULIC: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic elements" << endl);
	     Typ = ElemType::HYDRAULIC;
	     break;
	  }
#endif /* USE_HYDRAULIC_NODES */

	  case BULK: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "bulk" << endl);
	     Typ = ElemType::BULK;
	     break;
	  }

#ifdef HAVE_LOADABLE
	  case LOADABLE: {
	     DEBUGLCOUT(MYDEBUG_INPUT, "loadable" << endl);
	     Typ = ElemType::LOADABLE;
	     break;
	  }
#endif /* HAVE_LOADABLE */
	    
	  case UNKNOWNKEYWORD: {
	     cerr << "Error: unknown element type, cannot modify output" << endl;	     
	     THROW(DataManager::ErrGeneric());
	  }
	    
	  default: {
	     cerr << "Error: element type " << sKeyWords[CurrDesc] 
	       << " at line " << HP.GetLineData() << " is not allowed"
	       << endl;
	     THROW(DataManager::ErrGeneric());	    
	  }
	 }
	 
	 /* Elements list */
	 while (HP.fIsArg()) {
	    unsigned int uL = (unsigned int)HP.GetInt();	    	
	    Elem* pE = (Elem*)pFindElem(Typ, uL);
	    if (pE == NULL) {
	       cerr << "Error: " << psElemNames[Typ] << "(" << uL
		 << ") is not defined; output cannot be modified" << endl;
	    } else {
	       DEBUGLCOUT(MYDEBUG_INPUT, "element " << uL << endl);
	       pE->SetOutputFlag(flag(1));
	    }
	 }

      } else if (CurrDesc == BIND) {
	 /* Label dell'elemento */
	 unsigned int uL = HP.GetInt();
	 
	 /* Tipo dell'elemento */
	 ElemType::Type t = ElemType::UNKNOWN;	 
	 switch (KeyWords(HP.GetWord())) {
	  case BODY:
	    t = ElemType::BODY;
	    break;
	  case AUTOMATICSTRUCTURAL:
	    t = ElemType::AUTOMATICSTRUCTURAL;
	    break;
	  case JOINT:
	    t = ElemType::JOINT;
	    break;
	  case FORCE:
	  case COUPLE:
	    t = ElemType::FORCE;
	    break;
	  case BEAM:
	    t = ElemType::BEAM;
	    break;
	  case ROTOR:
	    t = ElemType::ROTOR;
	    break;
	  case AERODYNAMICBODY:
	  case AERODYNAMICBEAM:	
	    t = ElemType::AERODYNAMIC;
	    break;
	  case GENEL:
	    t = ElemType::GENEL;
	    break;
	  case ELECTRIC:
	    t = ElemType::ELECTRIC;
	    break;
	  case HYDRAULIC:     
	    t = ElemType::HYDRAULIC;
	    break;
	  case BULK:
	    t = ElemType::BULK;
	    break;	 
	  case LOADABLE:
	    t = ElemType::LOADABLE;
	    break;
	  default:
	    THROW(ErrGeneric());
	 }
	 
	 Elem* pEl = ((Elem*)pFindElem(t, uL));
	 if (pEl == NULL) {
	    cerr << "can't find " << psElemNames[t] << " (" << uL 
	      << ") at line " << HP.GetLineData() << endl;
	    THROW(ErrGeneric());
	 }
      
	 /* Label del nodo parameter */
	 uL = HP.GetInt();
	 
	 Elem2Param* pNd = ((Elem2Param*)pFindNode(NodeType::PARAMETER, uL));
	 if (pNd == NULL) {
	    cerr << "can't find parameter node (" << uL
	      << ") at line " << HP.GetLineData() << endl;
	    THROW(ErrGeneric());
	 }

	 /* Numero d'ordine del dato privato a cui fare il binding */
	 unsigned int i = HP.GetInt();
	 
	 /* indice del dato a cui il parametro e' bound */
	 if (i <= 0 || i > pEl->iGetNumPrivData()) {
	    cerr << "error in private data number " << i << " for element "
	      << psElemNames[t] << " (" << pEl->GetLabel() 
	      << ") at line " << HP.GetLineData() << endl;
	    THROW(ErrGeneric());
	 }
	 
	 /* fa il binding del ParameterNode all'elemento */
	 DEBUGLCOUT(MYDEBUG_INPUT, "Binding " << psElemNames[t] 
		    << " (" << pEl->GetLabel() 
		    << ") to Parameter " << pNd->GetLabel() << endl);
	 pNd->Bind(pEl, i);

	 
      /* gestisco a parte gli elementi automatici strutturali, perche'
       * sono gia' stati costruiti altrove e li devo solo inizializzare;
       * eventualmente si puo' fare altrimenti */
      } else if (CurrDesc == AUTOMATICSTRUCTURAL) {
	 unsigned int uLabel = HP.GetInt();
	 Elem* pEl = (Elem*)pFindElem(ElemType::AUTOMATICSTRUCTURAL, uLabel);
	 if (pEl == NULL) {
	    cerr << "line " << HP.GetLineData() 
	      << ": unable to find automatic structural element " 
	      << uLabel << endl;
	    THROW(ErrGeneric());
	 }
	 
	 DEBUGCOUT("reading automatic structural element " << uLabel << endl);

	 /* forse e' il caso di usare il riferimento del nodo? */
	 
	 /* nota: i primi due sono gestiti direttamente 
	  * dagli elementi con inerzia, e quindi non sono usati */
	 Vec3 q(HP.GetVecAbs(AbsRefFrame));
	 Vec3 g(HP.GetVecAbs(AbsRefFrame));
	 Vec3 qp(HP.GetVecAbs(AbsRefFrame));
	 Vec3 gp(HP.GetVecAbs(AbsRefFrame));
	 
	 DEBUGCOUT("Q  = " << q << endl
		   << "G  = " << g << endl
		   << "Qp = " << qp << endl
		   << "Gp = " << gp << endl);
	   
	 AutomaticStructElem* pAuto = (AutomaticStructElem*)pEl->pGet();
	 pAuto->Init(q, g, qp, gp);

      /* default: leggo un elemento e lo creo */
      } else {
      
      
	 /* puntatore al puntatore all'elemento */
	 Elem** ppE = NULL;
	 
	 unsigned int uLabel;
	 switch (CurrDesc) {
	    /* Qui vengono elencati gli elementi unici, che non richiedono label
	     * (per ora: accelerazione di gravita' e proprieta' dell'aria */
	    
#ifdef USE_STRUCT_NODES
	    /* Accelerazione di gravita' */
	  case GRAVITY: {
	     silent_cout("Reading gravity acceleration" << endl);
	     
	     if (iNumTypes[ElemType::GRAVITY]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": gravity acceleration element is not defined" << endl;
		
		THROW(DataManager::ErrGeneric());
	     }
	     
	     ppE = ElemData[ElemType::GRAVITY].ppFirstElem;
	     uLabel = 1;
	     
	      
	     TplDriveCaller<Vec3>* pDC 
	       = ReadTplDrive(this, HP, &DrvHdl, Vec3(0.));
	     HP.PutKeyTable(K);
	     
	     flag fOut = fReadOutput(HP, ElemType::GRAVITY);
	     
	     SAFENEWWITHCONSTRUCTOR(*ppE,
				    Gravity,
				    Gravity(pDC, fOut),
				    DMmm);
	     
	     break;
	  }
	    
#ifdef USE_AERODYNAMIC_ELEMS
	    /* Elementi aerodinamici: proprieta' dell'aria */
	  case AIRPROPERTIES: {
	     silent_cout("Reading air properties" << endl);
	     
	     if(iNumTypes[ElemType::AIRPROPERTIES]-- <= 0) {
		DEBUGCERR("");
		cerr << "line " << HP.GetLineData() 
		  << ": air properties element is not defined" << endl;
		
		THROW(DataManager::ErrGeneric());
	     }
	     
	     ppE = ElemData[ElemType::AIRPROPERTIES].ppFirstElem;
	     uLabel = 1;
	     
	     doublereal dRho = HP.GetReal();
	     DEBUGLCOUT(MYDEBUG_INPUT, "Air density: " << dRho << endl);
	     if (dRho <= 0.) {
		cerr 
		  << "illegal null or negative air density at line "
		  << HP.GetLineData() << endl;
		
		THROW(DataManager::ErrGeneric());
	     }
	     
	     doublereal dSS = HP.GetReal();
	     DEBUGLCOUT(MYDEBUG_INPUT, "Sound speed: " << dSS << endl);
	     if (dSS <= 0.) {
		cerr
		  << "illegal null or negative sound speed at line "
		  << HP.GetLineData() << endl;
		
		THROW(DataManager::ErrGeneric());
	     }	      
	     
	     /* Driver multiplo */	   
	     TplDriveCaller<Vec3>* pDC 
	       = ReadTplDrive(this, HP, &DrvHdl, Vec3(0.));	      	      
	     HP.PutKeyTable(K);	      	      
	     
	     flag fOut = fReadOutput(HP, ElemType::AIRPROPERTIES);
	     
	     SAFENEWWITHCONSTRUCTOR(*ppE, 
				    AirProperties,
				    AirProperties(pDC, dRho, dSS, fOut),
				    DMmm);
	     
	     break;
	  }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
	    
	    /* Elemento generico: legge la label e fa un nuovo switch */
	  default: {	      
	     /* legge la label */	
	     uLabel = unsigned(HP.GetInt());

	     /* in base al tipo, avviene l'allocazione */
	     switch (CurrDesc) {
		/* corpo rigido */
		
		
	      case DRIVEN: {
		 /* Reads the driver */
		 DriveCaller* pDC = ReadDriveData(this, HP, &DrvHdl);
		 HP.PutKeyTable(K);
		 
		 HP.ExpectDescription();
		 KeyWords CurrDriven = KeyWords(HP.GetDescription());
		 
#ifdef DEBUG
		 switch (CurrDriven) {
		  case FORCE:
#ifdef USE_STRUCT_NODES
		  case BODY:
		  case JOINT:
		  case COUPLE:
		  case BEAM:
#ifdef USE_AERODYNAMIC_ELEMS
		  case ROTOR:
		  case AERODYNAMICBODY:
		  case AERODYNAMICBEAM:
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
#ifdef USE_ELECTRIC_NODES
		  case GENEL:
		  case ELECTRIC:
#endif /* USE_ELECTRIC_NODES */
#ifdef USE_HYDRAULIC_NODES
		  case HYDRAULIC:
#endif /* USE_HYDRAULIC_NODES */
		  case BULK:
#ifdef HAVE_LOADABLE
		  case LOADABLE:
#endif /* HAVE_LOADABLE */
		  case EXISTING: {
		     DEBUGLCOUT(MYDEBUG_INPUT, "OK, this element can be driven" << endl);
		     break;
		  }
		    
		  default: {
		     DEBUGCERR("warning, this element can't be driven" << endl);
		     break;
		  }
		 }		     		  
#endif /* DEBUG */
		 
		 if (CurrDriven == EXISTING) {
		    iMissingElems++;
		    CurrDriven = KeyWords(HP.GetWord());
		    unsigned int uL = (unsigned int)HP.GetInt();
		    if (uL != uLabel) {
		       cerr << "Error: the driving element must have the same label of the driven" << endl;
		       
		       THROW(DataManager::ErrGeneric());
		    }
		    
		    switch (CurrDriven) {
		     case FORCE: {
			ppE = ppFindElem(ElemType::FORCE, uLabel);
			break;
		     }
#ifdef USE_STRUCT_NODES
		     case BODY: {			 
			ppE = ppFindElem(ElemType::BODY, uLabel);
			break;
		     }
		     case JOINT: {		       
			ppE = ppFindElem(ElemType::JOINT, uLabel);
			break;
		     }		    
		     case COUPLE: {
			ppE = ppFindElem(ElemType::FORCE, uLabel);
			break;
		     }
		     case BEAM: {
			ppE = ppFindElem(ElemType::BEAM, uLabel);
			break;
		     }
#ifdef USE_AERODYNAMIC_ELEMS
		     case ROTOR: {
			ppE = ppFindElem(ElemType::ROTOR, uLabel);
			break;
		     }
		     case AERODYNAMICBODY:
		     case AERODYNAMICBEAM: {
			ppE = ppFindElem(ElemType::AERODYNAMIC, uLabel);
			break;
		     }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
#ifdef USE_ELECTRIC_NODES
		     case GENEL: {		       
			ppE = ppFindElem(ElemType::GENEL, uLabel);
			break;
		     }
		     case ELECTRIC: {		       
			ppE = ppFindElem(ElemType::ELECTRIC, uLabel);
			break;
		     }
#endif /* USE_ELECTRIC_NODES */
#ifdef USE_HYDRAULIC_NODES
		     case HYDRAULIC: {		       
			ppE = ppFindElem(ElemType::HYDRAULIC, uLabel);
			break;
		     }
#endif /* USE_HYDRAULIC_NODES */
		     case BULK: {		       
			ppE = ppFindElem(ElemType::BULK, uLabel);
			break;
		     }
#ifdef HAVE_LOADABLE
		     case LOADABLE: {		       
			ppE = ppFindElem(ElemType::LOADABLE, uLabel);
			break;
		     }
#endif /* HAVE_LOADABLE */
		       
		     default: {
			DEBUGCERR("warning, this element can't be driven" << endl);
			break;
		     }			
		    }
		    
		    if (ppE == NULL) {
		       cerr << "Error: element " << uLabel 
			 << "cannot be driven since it doesn't exist" << endl;
		       
		       THROW(DataManager::ErrGeneric());
		    }
		    
		    
		    flag fOut = fReadOutput(HP, (*ppE)->GetElemType());
		    (*ppE)->SetOutputFlag(fOut);
		    
		    
		 } else {
		    unsigned int uDummy = (unsigned int)HP.GetInt();
		    if (uDummy != uLabel) {
		       cerr << "Error: the element label must be the same of the driving element" << endl;
		       
		       THROW(DataManager::ErrGeneric());
		    }
		    
		    /* Reads the true element */
		    ppE = ReadOneElem(this, HP, uLabel, CurrDriven);
		    HP.PutKeyTable(K);
		    
		    if (*ppE == NULL) {
		       DEBUGCERR("");
		       cerr << "error in allocation of element "
			 << uLabel << endl;
		       
		       THROW(ErrMemory());
		    }		  
		 }
		 
		 
		 /* Creates the driver for the element */
		 Elem* pEl = NULL;
		 SAFENEWWITHCONSTRUCTOR(pEl,
					DrivenElem,
					DrivenElem(pDC, *ppE),
					DMmm);
		 		 
		 /* Substitutes the element with the driver */
		 *ppE = pEl;
		 
		 break;
	      }
		
		
		/* Normal element */
	      case FORCE:
#ifdef USE_STRUCT_NODES
	      case BODY:
	      case JOINT:
	      case COUPLE:
	      case BEAM:
#ifdef USE_AERODYNAMIC_ELEMS
	      case ROTOR:
	      case AERODYNAMICBODY:
	      case AERODYNAMICBEAM:
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
#ifdef USE_ELECTRIC_NODES
	      case GENEL:
	      case ELECTRIC:
#endif /* USE_ELECTRIC_NODES */
#ifdef USE_HYDRAULIC_NODES
	      case HYDRAULIC:
#endif /* USE_HYDRAULIC_NODES */
	      case BULK:
#ifdef HAVE_LOADABLE
	      case LOADABLE:
#endif /* HAVE_LOADABLE */
		  {		 
		 /* Nome dell'elemento */
             	 const char *sName = NULL;
             	 if (HP.IsKeyWord("name")) {
	            const char *sTmp = HP.GetStringWithDelims();
	            SAFESTRDUP(sName, sTmp, DMmm);
	         }
		 
		 ppE = ReadOneElem(this, HP, uLabel, CurrDesc);
		 HP.PutKeyTable(K);
		 
		 if (sName != NULL) {
		    (*ppE)->PutName(sName);
   		    SAFEDELETEARR(sName, DMmm);
		 }

		 break;
	      }
		
		
		/* in caso di tipo sconosciuto */
	      case UNKNOWNKEYWORD: {
		 DEBUGCERR("");
		 cerr << "error - unknown element type at line " 
		   << HP.GetLineData() << endl;
		 
		 THROW(DataManager::ErrGeneric());
	      }
		
	      default: {
		 DEBUGCERR("");
		 cerr << "error - element type " << sKeyWords[CurrDesc]
		   << " at line " << HP.GetLineData() 
		   << " is not allowed " << endl;
		 
		 THROW(DataManager::ErrGeneric());
	      }        
	     }	      	      	      
	  }        
	 }
	 
	 /* verifica dell'allocazione */
	 ASSERT(*ppE != NULL);
	 
	 /* Aggiorna le dimensioni massime degli spazi di lavoro 
	  * (qui va bene perche' il puntatore e' gia' stato verificato) */
	 integer iNumRows = 0;
	 integer iNumCols = 0;
	 (*ppE)->WorkSpaceDim(&iNumRows, &iNumCols);
	 if (iNumRows > iMaxWorkNumRows) {
	    iMaxWorkNumRows = iNumRows;
	    DEBUGLCOUT(MYDEBUG_INIT, "Current max work rows number: " 
		       << iMaxWorkNumRows << endl);
	 }
	 if (iNumCols > iMaxWorkNumCols) {
	    iMaxWorkNumCols = iNumCols;
	    DEBUGLCOUT(MYDEBUG_INIT, "Current max work cols number: "
		       << iMaxWorkNumCols << endl);
	 }
	 
	 /* decrementa il totale degli elementi mancanti */
	 iMissingElems--;	  
      }
   }
   
   if (KeyWords(HP.GetWord()) != ELEMENTS) {
      DEBUGCERR("");
      cerr << "<end: elements;> expected at line" 
	<< HP.GetLineData() << endl;
           
      THROW(DataManager::ErrGeneric());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      DEBUGCERR("");
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;
     
      THROW(DataManager::ErrGeneric());
   }   
   
   if (iMissingElems > 0) {
      DEBUGCERR("");
      cerr << "warning: " << iMissingElems
	<< " elements are missing;" << endl;
      for (int iCnt = 0; iCnt < ElemType::LASTELEMTYPE; iCnt++) {
	 if (iNumTypes[iCnt] > 0) {
	    cerr << "  " << iNumTypes[iCnt] 
	      << ' ' << psElemNames[iCnt] << endl;
	 }	 
      }      
      
      THROW(DataManager::ErrGeneric());
   }
   
#ifdef USE_STRUCT_NODES
   /* Linka gli elementi che generano forze d'inerzia all'elemento 
    * accelerazione di gravita' */
   if (ElemData[ElemType::GRAVITY].iNum > 0) {
      Gravity* pGrav = 
	(Gravity*)(*(ElemData[ElemType::GRAVITY].ppFirstElem))->pGet();
      for (int iCnt = 0; iCnt < ElemType::LASTELEMTYPE; iCnt++) {
	 if ((ElemData[iCnt].fGeneratesInertialForces == 1)
	     && (ElemData[iCnt].iNum > 0)) {
	    
	    Elem** ppTmp = ElemData[iCnt].ppFirstElem;
	    Elem** ppLastEl = ppTmp+ElemData[iCnt].iNum;
	    while (ppTmp < ppLastEl) {	       
	       ASSERT((*ppTmp)->pGetElemGravityOwner() != NULL);
	       (*ppTmp)->pGetElemGravityOwner()->PutGravity(pGrav);
	       ppTmp++;
	    }	 
	 }      
      }   
   }

#ifdef USE_AERODYNAMIC_ELEMS
   /* Linka gli elementi che usano le proprieta' dell'aria all'elemento
    * proprieta' dell'aria */
   if (ElemData[ElemType::AIRPROPERTIES].iNum > 0) {      
      AirProperties* pProp =
	(AirProperties*)(*(ElemData[ElemType::AIRPROPERTIES].ppFirstElem))->pGet();
      
      for (int iCnt = 0; iCnt < ElemType::LASTELEMTYPE; iCnt++) {
	 if ((ElemData[iCnt].fUsesAirProperties == 1)
	     && (ElemData[iCnt].iNum > 0)) {	    
	    
	    Elem** ppTmp = ElemData[iCnt].ppFirstElem;	    
	    Elem** ppLastEl = ppTmp+ElemData[iCnt].iNum;	    
	    while (ppTmp < ppLastEl) {	      
	       ASSERT((*ppTmp)->pGetAerodynamicElem() != NULL);
	       (*ppTmp)->pGetAerodynamicElem()->PutAirProperties(pProp);
	       ppTmp++;
	    }	 	     
	 }      
      }
      
   } else {
      /* Esegue un controllo per vedere se esistono elementi aerodinamici
       * ma non sono definite le proprieta' dell'aria, nel qual caso
       * il calcolo deve essere arrestato */
      flag fStop(0);
      
      for (int iCnt = 0; iCnt < ElemType::LASTELEMTYPE; iCnt++) {
	 if ((ElemData[iCnt].fUsesAirProperties == 1)
	     && (ElemData[iCnt].iNum > 0)) {
	    if (fStop == 0) {
	       cerr << "warning, the following aerodynamic elements are defined: " << endl;	       
	       fStop = 1;
	    }
	    cerr << ElemData[iCnt].iNum << " " 
	      << psElemNames[iCnt] << endl;
	 }
      }
      
      if (fStop) {
	 cerr << "while no air properties are defined; aborting ..." << endl;
	 
	 THROW(DataManager::ErrGeneric());
      }
   }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
   
   DEBUGLCOUT(MYDEBUG_INPUT, "End of elements data" << endl);
} /* End of DataManager::ReadElems() */



Elem** ReadOneElem(DataManager* pDM, 
		   MBDynParser& HP,
		   unsigned int uLabel,
		   int CurrType)
{

   Elem** ppE = NULL;
   
   switch(KeyWords(CurrType)) {

      /* forza */
    case FORCE:
#ifdef USE_STRUCT_NODES
    case COUPLE:
#endif /* USE_STRUCT_NODES */
	{
       int iForceType;
       if (KeyWords(CurrType) == FORCE) {
	  iForceType = 0;
	  silent_cout("Reading force " << uLabel << endl);
       } else /* if(KeyWords(CurrType) == COUPLE) */ {
	  iForceType = 1;
	  silent_cout("Reading couple " << uLabel << endl);
       }
       
       if (iNumTypes[ElemType::FORCE]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": force " << uLabel
	    << " exceedes force elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::FORCE, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": force " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */		      
       int i = pDM->ElemData[ElemType::FORCE].iNum
	 -iNumTypes[ElemType::FORCE]-1;
       ppE = pDM->ElemData[ElemType::FORCE].ppFirstElem+i;
       
       *ppE = ReadForce(pDM, HP, uLabel, iForceType);
       
       break;
    }
      
#ifdef USE_STRUCT_NODES
    case BODY: {
       silent_cout("Reading rigid body " << uLabel << endl);
       
       if(iNumTypes[ElemType::BODY]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": rigid body " << uLabel
	    << " exceedes rigid body elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* verifica che non sia gia' definito */		    
       if(pDM->pFindElem(ElemType::BODY, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": rigid body " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::BODY].iNum-iNumTypes[ElemType::BODY]-1;
       ppE = pDM->ElemData[ElemType::BODY].ppFirstElem+i;
       
       *ppE = ReadBody(pDM, HP, uLabel);     
       
       break;
    }
      
      /* vincoli */
    case JOINT: {    
       silent_cout("Reading joint " << uLabel << endl);
       
       if (iNumTypes[ElemType::JOINT]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": joint " << uLabel
	    << " exceedes joint elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */		      
       if (pDM->pFindElem(ElemType::JOINT, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": joint " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::JOINT].iNum
	 -iNumTypes[ElemType::JOINT]-1;
       ppE = pDM->ElemData[ElemType::JOINT].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::JOINT].pFirstDofOwner+i;
       
       *ppE = ReadJoint(pDM, HP, pDO, uLabel);     
       
       /* attenzione: i Joint aggiungono DofOwner e quindi devono
	* completare la relativa struttura */
       break;
    }
      
      /* trave */
    case BEAM: {      
       silent_cout("Reading beam " << uLabel << endl);
       
       if (iNumTypes[ElemType::BEAM]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": beam " << uLabel
	    << " exceedes beam elements number" << endl;
	 
	  THROW(DataManager::ErrGeneric());
       }
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::BEAM, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": beam " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::BEAM].iNum
	 -iNumTypes[ElemType::BEAM]-1;
       ppE = pDM->ElemData[ElemType::BEAM].ppFirstElem+i;
       
       *ppE = ReadBeam(pDM, HP, uLabel);
       
       break;
    }
      
#ifdef USE_AERODYNAMIC_ELEMS
      /* Elementi aerodinamici: rotori */
    case ROTOR: {
       silent_cout("Reading rotor " << uLabel << endl);
       
       if (iNumTypes[ElemType::ROTOR]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": rotor " << uLabel
	    << " exceedes rotor elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::ROTOR, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": rotor " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::ROTOR].iNum
	 -iNumTypes[ElemType::ROTOR]-1;
       ppE = pDM->ElemData[ElemType::ROTOR].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::ROTOR].pFirstDofOwner+i;
       
       *ppE = ReadRotor(pDM, HP, pDO, uLabel);     
       
       break;
    }	 
      
      /* Elementi aerodinamici: rotori */
    case AERODYNAMICBODY:
    case AERODYNAMICBEAM: {
       silent_cout("Reading aerodynamic element " << uLabel << endl);
       
       if (iNumTypes[ElemType::AERODYNAMIC]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": aerodynamic element " << uLabel
	    << " exceedes aerodynamic elements number" << endl;
	 
	  THROW(DataManager::ErrGeneric());
       }
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::AERODYNAMIC, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": aerodynamic element " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::AERODYNAMIC].iNum
	 -iNumTypes[ElemType::AERODYNAMIC]-1;
       ppE = pDM->ElemData[ElemType::AERODYNAMIC].ppFirstElem+i;
       
       switch(KeyWords(CurrType)) {
	case AERODYNAMICBODY: {
	   *ppE = ReadAerodynamicBody(pDM, HP, uLabel);
	   break;
	}
	  
	case AERODYNAMICBEAM: {
	   *ppE = ReadAerodynamicBeam(pDM, HP, uLabel);
	   break;
	}
	  
	default: {
	   ASSERTMSG(0, "You shouldn't have reached this point");
	   break;
	}
       }     
       
       break;
    }
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */

      
#ifdef USE_ELECTRIC_NODES
      /* genel */
    case GENEL: {
       silent_cout("Reading genel " << uLabel << endl);

       if(iNumTypes[ElemType::GENEL]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": genel " << uLabel
	    << " exceedes genel elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */		      
       if(pDM->pFindElem(ElemType::GENEL, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": genel " << uLabel
	    << " already defined" << endl;
	 
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::GENEL].iNum
	 -iNumTypes[ElemType::GENEL]-1;
       ppE = pDM->ElemData[ElemType::GENEL].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::GENEL].pFirstDofOwner+i;
       
       *ppE = ReadGenel(pDM, HP, pDO, uLabel);
       
       break;
    }
      
#ifdef USE_HYDRAULIC_NODES
      /* elementi idraulici */
    case HYDRAULIC: {
       silent_cout("Reading hydraulic element " << uLabel << endl);
       
       if(iNumTypes[ElemType::HYDRAULIC]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": hydraulic element " << uLabel
	    << " exceedes hydraulic elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::HYDRAULIC, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": hydraulic element " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::HYDRAULIC].iNum
	 -iNumTypes[ElemType::HYDRAULIC]-1;
       ppE = pDM->ElemData[ElemType::HYDRAULIC].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::HYDRAULIC].pFirstDofOwner+i;
       
       *ppE = ReadHydraulicElem(pDM, HP, pDO, uLabel);
       
       /* attenzione: gli elementi elettrici aggiungono DofOwner 
	* e quindi devono completare la relativa struttura */
       
       break;
    }
#endif /* USE_HYDRAULIC_NODES */

      
      /* elementi elettrici */
    case ELECTRIC: {
       silent_cout("Reading electric element " << uLabel << endl);
       
       if(iNumTypes[ElemType::ELECTRIC]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": electric element " << uLabel
	    << " exceedes electric elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::ELECTRIC, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": electric element " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */
       int i = pDM->ElemData[ElemType::ELECTRIC].iNum
	 -iNumTypes[ElemType::ELECTRIC]-1;
       ppE = pDM->ElemData[ElemType::ELECTRIC].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::ELECTRIC].pFirstDofOwner+i;
       
       *ppE = ReadElectric(pDM, HP, pDO, uLabel);
       
       /* attenzione: gli elementi elettrici aggiungono DofOwner 
	* e quindi devono completare la relativa struttura */
       
       break;
    }
#endif /* USE_ELECTRIC_NODES */
      
      /* elementi bulk */
    case BULK: {
       silent_cout("Reading bulk element " << uLabel << endl);
       
       if (iNumTypes[ElemType::BULK]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": bulk element " << uLabel
	    << " exceedes bulk elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::BULK, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": bulk element " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       /* allocazione e creazione */		     
       int i = pDM->ElemData[ElemType::BULK].iNum
	 -iNumTypes[ElemType::BULK]-1;
       ppE = pDM->ElemData[ElemType::BULK].ppFirstElem+i;
       
       *ppE = ReadBulk(pDM, HP, uLabel);
       /* HP.PutKeyTable(K); */
       
       break;
    }		 		                     
    
#ifdef HAVE_LOADABLE
      /* elementi loadable */
    case LOADABLE: {
       silent_cout("Reading loadable element " << uLabel << endl);
       
       if (iNumTypes[ElemType::LOADABLE]-- <= 0) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": loadable element " << uLabel
	    << " exceedes loadable elements number" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }	  
       
       /* verifica che non sia gia' definito */
       if (pDM->pFindElem(ElemType::LOADABLE, uLabel) != NULL) {
	  DEBUGCERR("");
	  cerr << "line " << HP.GetLineData() 
	    << ": loadable element " << uLabel
	    << " already defined" << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* allocazione e creazione */		     
       int i = pDM->ElemData[ElemType::LOADABLE].iNum
	 -iNumTypes[ElemType::LOADABLE]-1;
       ppE = pDM->ElemData[ElemType::LOADABLE].ppFirstElem+i;
       DofOwner* pDO = pDM->DofData[DofType::LOADABLE].pFirstDofOwner+i;
              
       *ppE = ReadLoadable(pDM, HP, pDO, uLabel);
       
       break;
    }		 		                     
#endif /* defined(HAVE_LOADABLE) */
      
      
     
      /* In case the element type is not correct */
    default: {
       cerr << "You shouldn't be here" << endl;
       
       THROW(DataManager::ErrGeneric());
    }
   }
   
   /* Ritorna il puntatore al puntatore all'elemento appena costruito */
   return ppE;
}
