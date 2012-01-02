/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef HUTILS_H
#define HUTILS_H

#include "preselem.h"

/* Accumulator - begin */

class  Accumulator : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   doublereal stroke;     /* corsa pistone */
   doublereal start;      /* posizione iniziale*/
   doublereal area;       /* area stantufo */
   doublereal area_pipe;  /* area del tubo */
   doublereal mass;       /* massa valvola */
   doublereal press0;     /* pressione del gas con l'accumulatore scarico */
   doublereal press_max;  /* pressione massima del gas nell'accumulatore */
   doublereal Kappa;      /* costante della trasformazione */
   doublereal weight;     /* peso */
   doublereal spring;     /* costante della molla */
   doublereal force0;     /* precarico della molla */
   doublereal h_in;       /* coef per la perdita di carico concentrata */
   doublereal h_out;      /*  in entrata & uscita accumulatore */
   doublereal h;          /* quella giusta delle due*/
   doublereal density;    /* densita' del fluido */    
    
 
   doublereal s_min_gas;  /* spazio minimo del gas */
   doublereal s_max;      /* corsa massima della valvola */
   doublereal ratio2;      /* rapporto (area accumulatore/area tubo)^2 */
   doublereal c_spost;    /* coef dello spostamento */
   doublereal c_vel;      /* coef della velocita' */
   doublereal c_acc;      /* coef dell'accelerazione */
   doublereal c1;         /* coef dello spostamento */
   doublereal c2;         /* coef della velocita' */
   doublereal c3;         /* coef dell'accelerazione */
   doublereal c4;         /* coef. per il termine costante */
   doublereal cf1;        /* coef dello spostamento  finale */
   doublereal cf2;        /* coef della velocita'  finale */
   doublereal cf3;        /* coef dell'accelerazione finale */
   doublereal cf4;        /* coef. per il termine costante finale */
   
   doublereal s;          /* spostamento diaframma    (per l'output) */
   doublereal v;          /* velocita' diaframma      (per l'output) */
   doublereal sp;         /* velocita' diaframma       */
   doublereal vp;         /* accelerazione diaframma  (per l'output) */
   doublereal pgas;       /* pressione gas            (per l'output) */
   doublereal flow;       /* portata nodo 1           (per l'output) */
   
 public:
   Accumulator(unsigned int uL, const DofOwner* pD,
	       HydraulicFluid* hf, const PressureNode* p1, 
	       doublereal St, doublereal start, doublereal As, doublereal A_pipe, 
	       doublereal ms, doublereal h_in, doublereal h_out,
	       doublereal P0, doublereal Pmax, doublereal k, 
	       doublereal Wg, doublereal Kspr, doublereal force0, 
	       doublereal cs, doublereal cv, doublereal ca, flag fOut);

   ~Accumulator(void);
   
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicElem::Type GetHydraulicType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void Output(OutputHandler& OH) const;
   
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode1;
   };
  /* ************************************************ */ 	 
};

/* Accumulator - end */


/* Tank - begin */

class Tank : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   doublereal press;
   doublereal area_pipe;
   doublereal area_serb;
   doublereal level;
   doublereal s_max;   /* livello massimo dell'olio */
   doublereal s_min;   /* soglia d'allarme */
   doublereal Kappa1;  /* coef.di perdita concentrata fluido entrante nel serbatoio */
   doublereal Kappa2;  /* coef.di perdita concentrata fluido uscente dal serbatoio */
   doublereal c_spost; 
   doublereal s;       /* livello */
   doublereal sp;      /* tasso di crescita del livello */
   doublereal flow1;   /* portata al nodo di carico  (per l'output) */
   doublereal flow2;   /* portata al nodo di scarico (per l'output) */
   
 public:
   Tank(unsigned int uL, const DofOwner* pD,
	HydraulicFluid* hf,
	const PressureNode* p1, const PressureNode* p2, 
	doublereal press, doublereal A_pipe, doublereal A_serb,
	doublereal lev, doublereal s_mx, doublereal s_mn,
	doublereal c_s, flag fOut);
   
   ~Tank(void);
   
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicElem::Type GetHydraulicType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void Output(OutputHandler& OH) const;
   
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

  /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */
};

/* Tank - end */

#endif
