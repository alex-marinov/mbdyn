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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include <mbconfig.h>

#include <dataman.h>
#include <preselem.h>

#include <actuator.h>
#include <hfluid.h>
#include <hminor.h>
#include <hutils.h>
#include <pipe.h>
#include <valve.h>

/* HydraulicElem - begin */

HydraulicElem::HydraulicElem(unsigned int uL, const DofOwner* pDO,
			     HydraulicFluid* hf, flag fOut)
: Elem(uL, ElemType::HYDRAULIC, fOut), 
ElemWithDofs(uL, ElemType::HYDRAULIC, pDO, fOut), 
HF(hf)
{
   if (HF == NULL) {
      cerr << "HydraulicElem::HydraulicElem(" << GetLabel() 
	<< "): NULL hydraulic fluid pointer (FIXME)" << endl; 
   }
   ASSERT(HF != NULL);
} 


HydraulicElem::~HydraulicElem(void) 
{
   if (HF != NULL) {
      SAFEDELETE(HF, DMmm);
   }
}


/* Tipo dell'elemento (usato per debug ecc.) */
ElemType::Type HydraulicElem::GetElemType(void) const
{
   return ElemType::HYDRAULIC;
}


/* Contributo al file di restart 
 * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
 * relativa alla classe derivata */
ostream& HydraulicElem::Restart(ostream& out) const 
{
   return out << "  hydraulic: " << GetLabel();
}

/* Output */
void HydraulicElem::Output(OutputHandler& OH) const 
{
   NO_OP;
}

void HydraulicElem::SetInitialValue(VectorHandler& /* X */ ) const {
   NO_OP; 
}

/* HydraulicElem - end */





Elem* ReadHydraulicElem(DataManager* pDM,
			MBDynParser& HP, 
			const DofOwner* pDO, 
			unsigned int uLabel)
{
   DEBUGCOUT("ReadHydraulicElem()");
   
   const char* sKeyWords[] = {
      "minor" "loss",
      "control" "valve",
      "dynamic" "control" "valve",
      "pressure" "flow" "control",
      "pressure" "valve",
      "flow" "valve",
      "orifice",
      "accumulator",
      "tank",
      "pipe",
      "dynamic" "pipe",
      "actuator"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
      MINOR_LOSS = 0, 
      CONTROL_VALVE,
      DYNAMIC_CONTROL_VALVE,
      PRESSURE_FLOW_CONTROL_VALVE,
      PRESSURE_VALVE,
      FLOW_VALVE,
      ORIFICE,
      ACCUMULATOR,
      TANK,
      PIPE,
      DYNAMIC_PIPE,
      ACTUATOR,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   /* lettura del tipo di elemento elettrico */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      cout << "hydraulic element type: "
	<< sKeyWords[CurrKeyWord] << endl;
   }   
#endif   
   
   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {
      
    case ACTUATOR: {
#if defined(USE_STRUCT_NODES)
       /* lettura dei dati specifici */
       /* due nodi idraulici e due nodi strutturali */
       
       /* nodo idraulico 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNodeHyd1;
       if ((pNodeHyd1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo idraulico 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNodeHyd2;
       if ((pNodeHyd2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo strutturale 1 */
       unsigned int uNode1 = (unsigned int)HP.GetInt();
       
       DEBUGCOUT("Linked to Node " << uNode1 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNodeStr1;
       if ((pNodeStr1 = pDM->pFindStructNode(uNode1)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode1
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 f1(HP.GetPosRel(ReferenceFrame(pNodeStr1)));
       DEBUGCOUT("Offset 1: " << f1 << endl);
       
       /* nodo strutturale 2 */
       unsigned int uNode2 = (unsigned int)HP.GetInt();
       DEBUGCOUT("Linked to Node " << uNode2 << endl);
       
       /* verifica di esistenza del nodo */
       StructNode* pNodeStr2;
       if ((pNodeStr2 = pDM->pFindStructNode(uNode2)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": structural node " << uNode2
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       Vec3 f2(HP.GetPosRel(ReferenceFrame(pNodeStr2)));
       DEBUGCOUT("Offset 2: " << f2 << endl);  
       
       ReferenceFrame RF(pNodeStr1);
       Vec3 axis(0., 0., 1.); 
       if (HP.IsKeyWord("direction")) {
	  axis = HP.GetVecRel(RF);
	  doublereal d = axis.Norm();
	  if (d < DBL_EPSILON) {
	     cerr << "need a definite direction, not " << axis << "!" << endl;
	     THROW(ErrGeneric());
	  }
	  axis /= d;
       } 
       
       /* Area nodo1 */
       doublereal area1 = HP.GetReal();
       if (area1 <= DBL_EPSILON) {		  
	  cerr << "null or negative area1 in actuator"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Area1: " << area1 << endl);
       
       /* Area nodo2 */
       doublereal area2 = HP.GetReal();
       if (area2 <= DBL_EPSILON) {		  
	  cerr << "null or negative area2 in actuator"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Area2: " << area1 << endl);
       
       /* lunghezza cilindro (a meno dello spessore */
       doublereal dl = HP.GetReal();
       if (dl <= DBL_EPSILON) {		  
	  cerr << "null or negative dl in actuator"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("dl: " << area1 << endl);
       
       
       HydraulicFluid* hf1 = HP.GetHydraulicFluid();
       ASSERT(hf1 != NULL);
       
       HydraulicFluid* hf2 = NULL; 
       if (HP.IsKeyWord("same")) {
	  hf2 = hf1->pCopy();
       } else {
	  hf2 = HP.GetHydraulicFluid();
       }
       ASSERT(hf2 != NULL);

       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Actuator,
			      Actuator(uLabel, pDO, 
				       pNodeHyd1, pNodeHyd2, 
				       pNodeStr1, pNodeStr2,
				       f1, f2, axis, hf1, hf2, 
				       area1, area2, dl,
				       fOut),
			      DMmm);
       
       break;
    }	
      
#else /* defined(USE_STRUCT_NODES) */
      cerr << "you are not allowed to use actuators" << endl;
      THROW(ErrGeneric());
#endif /* defined(USE_STRUCT_NODES) */
      
    case MINOR_LOSS: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Kappa1 diretto */
       doublereal dKappa1 = HP.GetReal();
       if (dKappa1 < 0.) {		  
	  cerr << "negative Kappa1 in minor loss"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Kappa1: " << dKappa1 << endl);
       
       /* Kappa2 inverso */
       doublereal dKappa2 = HP.GetReal();
       if (dKappa2 < 0.) {		  
	  cerr << "negative Kappa2 in minor loss"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Kappa2: " << dKappa2 << endl);
       
       /* Area */
       doublereal area = HP.GetReal();
       if (area <= DBL_EPSILON) {		  
	  cerr << "null or negative area in minor loss"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Area: " << area << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);

       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Minor_loss,
                              Minor_loss(uLabel, pDO, hf, pNode1, pNode2, 
					 dKappa1, dKappa2, area, fOut), 
			      DMmm);
       
       break;
    }

    case CONTROL_VALVE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 3 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode3;
       if ((pNode3 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 4 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode4;
       if ((pNode4 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Area massima della valvola */
       doublereal area_max = HP.GetReal();
       if (area_max <= 0.) {		  
	  cerr << "null or negative area_max in control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_max: " << area_max << endl);
       
       /* Area di trafilamento in % sull'area massima:valore di default = 1.e-6 */
       doublereal loss_area = 1.e-6; 
       if (HP.IsKeyWord("loss")) {
	  loss_area = HP.GetReal();
	  if (loss_area  < 0.) {		  
	     cerr << "negative loss_area in control valve "
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	  
	  DEBUGCOUT("Loss_area in %= " << loss_area << endl); 
       }
       
       /* Stato */
       DriveCaller* pDC =ReadDriveData(pDM, HP,pDM->pGetDrvHdl());
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Control_valve,
			      Control_valve(uLabel, pDO, hf, pNode1, pNode2, 
					    pNode3, pNode4, area_max, loss_area, pDC,
					    fOut), 
			      DMmm);
       
       break;
    }
      
    case DYNAMIC_CONTROL_VALVE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 3 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode3;
       if ((pNode3 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 4 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode4;
       if ((pNode4 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Forza */
       DriveCaller* pDC = ReadDriveData(pDM, HP,pDM->pGetDrvHdl());
       
       /* spostamento iniziale */
       doublereal start = HP.GetReal();
       DEBUGCOUT("Start: " << start << endl);
       
       /* Spostamento massimo della valvola */
       doublereal s_max = HP.GetReal();
       if (s_max < 0.) {		  
	  cerr << "negative s_max in dynamic control valve"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("S_max: " << s_max << endl);
       
       /* Larghezza del condotto */
       doublereal width = HP.GetReal();
       if (width <= 0.) {		  
	  cerr << "null or negative width in dynamic control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Width: " << width << endl);
       
       /* Area di trafilamento in % sull'area massima(==width*s_max):valore di default = 1.e-6 */
       doublereal loss_area =1.e-6;
       if (HP.IsKeyWord("loss")) {
	  loss_area = HP.GetReal();
	  if (loss_area < 0.) {		  
	     cerr << "negative loss_area in dynamic control valve "
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  DEBUGCOUT("Loss_area in %= " << loss_area << endl); 
       }
       
       /* Diametro della valvola */
       doublereal valve_diameter = HP.GetReal();
       if (valve_diameter <= 0.) {		  
	  cerr << "null or negative valve diameter in dynamic control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve diameter: " << valve_diameter << endl);
       
       /* Densita' del corpo della valvola */
       doublereal valve_density = HP.GetReal();
       if (valve_density <= 0.) {		  
	  cerr << "null or negative valve density in dynamic control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve density: " << valve_density << endl);
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       /* c della velocita' */
       doublereal c_vel = HP.GetReal();
       DEBUGCOUT("c_vel: " << c_vel << endl);
       
       /* c della accelerazione */
       doublereal c_acc = HP.GetReal();
       DEBUGCOUT("c_acc: " << c_acc << endl);
       	
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Dynamic_control_valve,
			      Dynamic_control_valve(uLabel, pDO, hf, 
						    pNode1, pNode2, 
						    pNode3, pNode4, 
						    pDC, start,
						    s_max, width, 
						    loss_area, 
						    valve_diameter, 
						    valve_density,
						    c_spost, c_vel, c_acc,
						    fOut), 
			      DMmm);
       break;
    }

    case PRESSURE_FLOW_CONTROL_VALVE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 3 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode3;
       if ((pNode3 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 4 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode4;
       if ((pNode4 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
         /* nodo 5 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode5;
       if ((pNode5 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
         /* nodo 6 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode6;
       if ((pNode6 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       
       
       /* Forza */
       DriveCaller* pDC =ReadDriveData(pDM, HP,pDM->pGetDrvHdl());
       
       /* spostamento iniziale */
       doublereal start = HP.GetReal();
       if (start < 0.) {		  
	  cerr << "negative start in pressure flow control valve"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       } 
       DEBUGCOUT("Start: " << start << endl);
       
       /* Spostamento massimo della valvola */
       doublereal s_max = HP.GetReal();
       if (s_max < 0.) {		  
	  cerr << "negative s_max in pressure flow control valve"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("S_max: " << s_max << endl);
       
       /* Larghezza del condotto */
       doublereal width = HP.GetReal();
       if (width <= 0.) {		  
	  cerr << "null or negative width in pressure flow control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Width: " << width << endl);
       
       /* Area di trafilamento in % sull'area massima(==width*s_max):valore di default = 1.e-6 */
       doublereal loss_area =1.e-6;
       if (HP.IsKeyWord("loss")) {
	  loss_area = HP.GetReal();
	  if (loss_area < 0.) {		  
	     cerr << "negative loss_area in pressure flow control valve "
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  DEBUGCOUT("Loss_area in %= " << loss_area << endl); 
       }
       
       /* Diametro della valvola */
       doublereal valve_diameter = HP.GetReal();
       if (valve_diameter <= 0.) {		  
	  cerr << "null or negative valve diameter in pressure flow control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve diameter: " << valve_diameter << endl);
       
       /* Densita' del corpo della valvola */
       doublereal valve_density = HP.GetReal();
       if (valve_density <= 0.) {		  
	  cerr << "null or negative valve density in pressure flow control valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve density: " << valve_density << endl);
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       /* c della velocita' */
       doublereal c_vel = HP.GetReal();
       DEBUGCOUT("c_vel: " << c_vel << endl);
       
       /* c della accelerazione */
       doublereal c_acc = HP.GetReal();
       DEBUGCOUT("c_acc: " << c_acc << endl);
       	
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Pressure_flow_control_valve,
			      Pressure_flow_control_valve(uLabel, pDO, hf, 
						    pNode1, pNode2, 
						    pNode3, pNode4, 
						    pNode5, pNode6, 
						    pDC, start,
						    s_max, width, 
						    loss_area, 
						    valve_diameter, 
						    valve_density,
						    c_spost, c_vel, c_acc,
						    fOut), 
			      DMmm);
       break;
    }
      
      
    case PRESSURE_VALVE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Area diaframma */
       doublereal area_diaf = HP.GetReal();
       if (area_diaf <= 0.) {		  
	  cerr << "null or negative area_diaf in pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_diaf: " << area_diaf << endl);
       
       /* Massa valvola */
       doublereal mass = HP.GetReal();
       if (mass <= 0.) {		  
	  cerr << "null or negative valve mass in pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve mass: " << mass << endl);
       
       /* Area massima della valvola */
       doublereal area_max = HP.GetReal();
       if (area_max <= 0.) {		  
	  cerr << "null or negative area_max in pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_max: " << area_max << endl);
       
       /* Spostamento massimo della valvola */
       doublereal s_max = HP.GetReal();
       if (s_max < 0.) {		  
	  cerr << "negative s_max in pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("S_max: " << s_max << endl);
       
       /* Kappa : costante della molla */
       doublereal Kappa = HP.GetReal();
       if (Kappa < 0.) {		  
	  cerr << "negative Kappa in  pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Kappa: " << Kappa << endl);
       
       /* Forza0: precarico della molla */
       doublereal force0 = HP.GetReal();
       if (force0 < 0.) {		  
	  cerr << "negative force0 in  pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Force0: " << force0 << endl);
       
       /* Larghezza luce di passaggio */
       doublereal width = HP.GetReal();
       if (width <= 0.) {		  
	  cerr << "null or negative width in pressure valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Width: " << width << endl);
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       /* c della velocita' */
       doublereal c_vel = HP.GetReal();
       DEBUGCOUT("c_vel: " << c_vel << endl);
       
       /* c della accelerazione */
       doublereal c_acc = HP.GetReal();
       DEBUGCOUT("c_acc: " << c_acc << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
	       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Pressure_valve,
                              Pressure_valve(uLabel, pDO, hf, pNode1, pNode2, 
					     area_diaf, mass, area_max, 
					     s_max, Kappa, force0, width,
					     c_spost, c_vel, c_acc,
					     fOut), 
			      DMmm);
       
       break;
    }
      
    case FLOW_VALVE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* nodo 3 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode3;
       if ((pNode3 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Area diaframma */
       doublereal area_diaf = HP.GetReal();
       if (area_diaf <= 0.) {		  
	  cerr << "null or negative area_diaf in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_diaf: " << area_diaf << endl);
       
       /* Massa valvola */
       doublereal mass = HP.GetReal();
       if (mass <= 0.) {
	  cerr << "null or negative valve mass in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Valve mass: " << mass << endl);
       
       /* Area tubo */
       doublereal area_pipe = HP.GetReal();
       if (area_pipe <= 0.) {
	  cerr << "null or negative area_pipe in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_pipe: " << area_pipe << endl);
            
       /* Area massima della valvola */
       doublereal area_max = HP.GetReal();
       if (area_max <= 0.) {
	  cerr << "null or negative area_max in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_max: " << area_max << endl);
       
       /* Kappa : costante della molla */
       doublereal Kappa = HP.GetReal();
       if (Kappa <= 0.) {
	  cerr << "null or negative Kappa in  flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Kappa: " << Kappa << endl);
       
       /* Forza0: precarico della molla */
       doublereal force0 = HP.GetReal();
       if (force0 < 0.) {		  
	  cerr << "negative force0 in  flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Force0: " << force0 << endl);
       
       /* Larghezza luce di passaggio */
       doublereal width = HP.GetReal();
       if (width <= 0.) {		  
	  cerr << "null or negative width in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Width: " << width << endl);
       
       /* Corsa massima della valvola */
       doublereal s_max = HP.GetReal();
       if (s_max < 0.) {		  
	  cerr << "negative s_max in flow valve "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("s_max: " << s_max << endl);
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       /* c della velocita' */
       doublereal c_vel = HP.GetReal();
       DEBUGCOUT("c_vel: " << c_vel << endl);
       
       /* c della accelerazione */
       doublereal c_acc = HP.GetReal();
       DEBUGCOUT("c_acc: " << c_acc << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Flow_valve,
                              Flow_valve(uLabel, pDO, hf, 
					 pNode1, pNode2, pNode3,
					 area_diaf, mass,area_pipe, area_max,
					 Kappa, force0, width, s_max,
					 c_spost, c_vel, c_acc,
					 fOut),
			      DMmm);
       
       break;
    }
      
    case ORIFICE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Diametro */
       doublereal diameter = HP.GetReal();
       if (diameter <= 0.) {		  
	  cerr << "null or negative diameter in orifice"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Diameter: " << diameter << endl);
       
       /* Area diaframma */
       doublereal area_diaf = HP.GetReal();
       if (area_diaf <= 0.) {		  
	  cerr << "null or negative area_diaf in orifice"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Area_diaf: " << area_diaf << endl);
  
       /* Area del tubo */
       doublereal area_pipe = diameter*diameter*0.785;
       if (HP.IsKeyWord("area")) 
	 {
	    area_pipe = HP.GetReal();
	    if (area_pipe <= 0.) 
	      {		  
		 cerr << "null or negative area_pipe in orifice"
		   << uLabel << " at " << HP.GetLineData() << endl;
		 THROW(DataManager::ErrGeneric());
	      }	 
	 }
       DEBUGCOUT("Area_pipe: " << area_pipe << endl);
       
       doublereal ReCr = 10;
       if (HP.IsKeyWord("ReCr")) 
	 {
	    ReCr = HP.GetReal();
	    if (ReCr <= 0.) 
	      {		  
		 cerr << "null or negative Reynold's number in orifice"
		   << uLabel << " at " << HP.GetLineData() << endl;
		 THROW(DataManager::ErrGeneric());
	      }	 
	 }
       DEBUGCOUT("Reynold critico: " << ReCr << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Orifice,
			      Orifice(uLabel, pDO, hf, 
				      pNode1, pNode2, 
				      diameter, 
				      area_diaf, area_pipe, ReCr, fOut),
			      DMmm);
       break;
    }
      
    case ACCUMULATOR: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* Corsa pistone */
       doublereal stroke = HP.GetReal();
       if (stroke <= 0.) {		  
	  cerr << "null or negative stroke in accumulator "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Stroke: " << stroke << endl);
          
       doublereal start = 0.;
       if (HP.IsKeyWord("start")) {	       
	  // Corsa iniziale del setto    	   
	  start = HP.GetReal();
	  if (start > stroke) 
	    {		  
	       cerr << "Accumulator: stroke minor then inizial position"
		 << uLabel << " at " << HP.GetLineData() << endl;
	       THROW(DataManager::ErrGeneric());
	    }
       }	    	  
       DEBUGCOUT("start: " << start << endl);
       
       /* Area stantuffo */
       doublereal area = HP.GetReal();
       if (area <= 0.) {		  
	  cerr << "null or negative area in accumulator "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area: " << area << endl);
       
       /* Area pipe */
       doublereal area_pipe = HP.GetReal();
       if (area_pipe <= 0.) {		  
	  cerr << "null or negative area_pipe in accumulator "
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("area_pipe: " << area_pipe << endl); 
       
       /* Massa stantuffo */
       doublereal mass = HP.GetReal();
       if (mass <= 0.) {		  
	  cerr << "null or negative mass in accumulator"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	 
       DEBUGCOUT("Mass: " << mass << endl);
         
       doublereal h_in = 1;
       if (HP.IsKeyWord("lossin")) {	       
	  // Perdita di carico entrata
	  h_in = HP.GetReal();
	  if (h_in < 0.) 
	    {		  
	       cerr << "Negative loss_in in accumulator"
		      << uLabel << " at " << HP.GetLineData() << endl;
	       THROW(DataManager::ErrGeneric());
	    }
       }	    	  
       DEBUGCOUT("Loss_in: " << h_in << endl);
       
       doublereal h_out = 0.5;
	    if (HP.IsKeyWord("lossout")) {	       
	       // Perdita di carico uscita    	   
	       h_out = HP.GetReal();
 	       if (h_out < 0.) 
		 {		  
		    cerr << "Negative loss_out in accumulator"
		      << uLabel << " at " << HP.GetLineData() << endl;
		    THROW(DataManager::ErrGeneric());
		 }
	    }	    	  
	    DEBUGCOUT("loss_out: " << h_out << endl);
       
       doublereal press0   = 0.;
       doublereal press_max= 0.;
       doublereal Kappa    = 0.;
       
       if (HP.IsKeyWord("gas")) {
	  
	  /* Pressione gas accumulatore scarico */
	  press0 = HP.GetReal();
	  if (press0 <= 0.) {		  
	     cerr << "null or negative pressure0 in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	    	  
	  DEBUGCOUT("press0: " << press0 << endl);
	  
	  /* Pressione massima del gas */
	  press_max = HP.GetReal();
	  if (press_max <= 0.) {		  
	     cerr << "null or negative pressure max in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	     
	  DEBUGCOUT("Pressure max: " << press_max << endl);
	  
	  Kappa = HP.GetReal();
	  if (Kappa < 0.) {		  
	     cerr << "negative Kappa in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	 
	  DEBUGCOUT("Kappa: " << Kappa << endl);
       }
       
       doublereal weight = 0.;
       if (HP.IsKeyWord("weight")) {
	  weight = HP.GetReal();
	  if (weight <= 0.) {		  
	     cerr << "null or negative weight in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	  
	  DEBUGCOUT("weight: " << weight << endl);
       }
       
       doublereal spring = 0.;
       doublereal force0 = 0.;
       if (HP.IsKeyWord("spring")) {
	  spring = HP.GetReal();
	  if (spring < 0.) {		  
	     cerr << "negative spring in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  
	  force0 = HP.GetReal();
	  if (force0 < 0.) {		  
	     cerr << "negative force0 in accumulator"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }
	  DEBUGCOUT("spring: " << spring << endl);
	  DEBUGCOUT("force0: " << force0 << endl);
       }
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       /* c della velocita' */
       doublereal c_vel = HP.GetReal();
       DEBUGCOUT("c_vel: " << c_vel << endl);
       
       /* c della accelerazione */
       doublereal c_acc = HP.GetReal();
       DEBUGCOUT("c_acc: " << c_acc << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Accumulator,
			      Accumulator(uLabel, pDO, hf, pNode1, 
					  stroke, start, area, area_pipe, 
					  mass,h_in, h_out,
					  press0, press_max,
					  Kappa, weight, spring, force0, 
					  c_spost, c_vel, c_acc, fOut),
			      DMmm);
       break;
    }
      
    case TANK: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Pressione serbatoio */
       doublereal press = HP.GetReal();
       if (press <= 0.) {		  
	  cerr << "null or negative pressure in tank"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Pressure: " << press << endl);
       
       /* Area pipe */
       doublereal area_pipe = HP.GetReal();
       if (area_pipe <= 0.) {		  
	  cerr << "null or negative area_pipe in tank"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area_pipe: " << area_pipe << endl); 
       
       /* Area serbatoio */
       doublereal area_serb = HP.GetReal();
       if (area_serb <= 0.) {		  
	  cerr << "null or negative area_serb in tank"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Area serbatoio: " << area_serb << endl);
       
       /* Livello massimo dell'olio */
       doublereal s_max = HP.GetReal();
       if (s_max < 0.) {		  
	  cerr << "negative s_max in tank"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Livello massimo dell'olio: " << s_max << endl);
       
       /* Livello iniziale */
       doublereal level= .5*s_max; /* valore di default 50% del massimo */
       
       if (HP.IsKeyWord("startlevel")) {
	  level = HP.GetReal();
	  if (level < 0.) {		  
	     cerr << "negative level in tank"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	     
	  DEBUGCOUT("Livello iniziale: " << level << endl);
       }
       
       /* Soglia di allarme */
       doublereal s_min = .1*s_max; /* valore di default 10% del massimo */
       if (HP.IsKeyWord("alarmlevel")) {
	  doublereal s_min = HP.GetReal();
	  if (s_min < 0.) {
	     cerr << "negative s_min in tank"
	       << uLabel << " at " << HP.GetLineData() << endl;
	     THROW(DataManager::ErrGeneric());
	  }	     
	  DEBUGCOUT("Soglia di allarme: " << s_min << endl);
       }
       
       /* c dello spostamento */
       doublereal c_spost = HP.GetReal();
       DEBUGCOUT("c_spost: " << c_spost << endl);
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Tank,
			      Tank (uLabel, pDO, hf, pNode1,pNode2, press,
				    area_pipe, area_serb,
				    level, s_max, s_min, c_spost, fOut),
			      DMmm);
       break;
    }
      
    case PIPE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Diametro */
       doublereal diameter = HP.GetReal();
       if (diameter <= 0.) {		  
	  cerr << "null or negative diameter in pipe"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Diameter: " << diameter << endl);
       
       // Area      	   
       doublereal area = diameter*diameter*0.785;
       if (HP.IsKeyWord("area")) 
	 {
	    area = HP.GetReal();
	    if (area <= 0.) 
	      {		  
		 cerr << "null or negative area in pipe"
		   << uLabel << " at " << HP.GetLineData() << endl;
		 THROW(DataManager::ErrGeneric());
	      }	
	      }
       DEBUGCOUT("Area: " << area << endl);
       
       /* Lunghezza */
       doublereal lenght = HP.GetReal();
       if (lenght <= 0.) {		  
	  cerr << "null or negative lenght in pipe"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Lenght: " << lenght << endl); 
       
       /* Transizione se e' 0 parto da laminare se e' 1 parto da turbolento */
       flag turbulent = 0;
       if (HP.IsKeyWord("turbulent")) {
	  turbulent = 1;
	  DEBUGCOUT("Turbulent" << endl); 
       }
       doublereal q0 = 0.;
       if (HP.IsKeyWord("initialvalue")) {
	  q0 = HP.GetReal();
	  DEBUGCOUT("Initial q = " << q0 << endl); 
       }
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Pipe,
			      Pipe(uLabel, pDO, hf, pNode1, pNode2, 
				   diameter, 
				   area, lenght, turbulent, q0, fOut),
			      DMmm);
       break;
    }
      
    case DYNAMIC_PIPE: {
       
       /* nodo 1 */
       unsigned int uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode1;
       if ((pNode1 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }		  
       
       /* nodo 2 */
       uNode = (unsigned int)HP.GetInt();	     
       DEBUGCOUT("Linked to Hydraulic Node " << uNode << endl);
       
       /* verifica di esistenza del nodo idraulico */
       PressureNode* pNode2;
       if ((pNode2 = (PressureNode*)pDM->pFindNode(NodeType::HYDRAULIC, uNode)) == NULL) {
	  cerr << "line " << HP.GetLineData() 
	      << ": hydraulic node " << uNode
	    << " not defined" << endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       doublereal diameter = HP.GetReal();
       if (diameter <= 0.) {		  
	  cerr << "null or negative diameter in dynamic pipe"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Diameter: " << diameter << endl);
       
       // Area      	   
       doublereal area = diameter*diameter*0.785;
       if (HP.IsKeyWord("area")) 
	      {
		 area = HP.GetReal();
		 if (area <= 0.) 
		   {		  
		      cerr << "null or negative area in pipe"
			<< uLabel << " at " << HP.GetLineData() << endl;
		      THROW(DataManager::ErrGeneric());
		   }	
	      }
       DEBUGCOUT("Area: " << area << endl);
       
       /* Lunghezza */
       doublereal lenght = HP.GetReal();
       if (lenght <= 0.) {		  
	  cerr << "null or negative lenght in dynamic pipe"
	    << uLabel << " at " << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }	     
       DEBUGCOUT("Lenght: " << lenght << endl); 
       
       /* Transizione se e' 0 parto da laminare se e' 1 parto da turbolento */
       flag turbulent = 0;
       if (HP.IsKeyWord("turbulent")) {
	  turbulent = 1;
	  DEBUGCOUT("Turbulent" << endl); 
       }
       doublereal q0 = 0.;
       if (HP.IsKeyWord("initialvalue")) {
	  q0 = HP.GetReal();
	  DEBUGCOUT("Initial q = " << q0 << endl); 
       }
       
       HydraulicFluid* hf = HP.GetHydraulicFluid();
       ASSERT(hf != NULL);
       
       flag fOut = pDM->fReadOutput(HP, ElemType::HYDRAULIC);

#if 0
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      Dynamic_pipe,
			      Dynamic_pipe(uLabel, pDO, hf,
					   pNode1, pNode2, diameter, 
					   area, lenght, turbulent, q0, fOut),
			      DMmm);
#else
       SAFENEWWITHCONSTRUCTOR(pEl,
			      DynamicPipe,
			      DynamicPipe(uLabel, pDO, hf,
					   pNode1, pNode2, diameter, 
					   area, lenght, turbulent, q0, fOut),
			      DMmm);
#endif
       break;
    }	   
      
      /* Aggiungere altri elementi idraulici */
      
    default: {
       cerr << "unknown hydraulic element type in hydraulic element " << uLabel
	 << " at line " << HP.GetLineData() << endl;       
       THROW(DataManager::ErrGeneric());
    }	
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << "semicolon expected at line " << HP.GetLineData() << endl;     
      THROW(DataManager::ErrGeneric());
   }
   
   return pEl;
}
