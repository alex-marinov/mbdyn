/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef PIPE_H
#define PIPE_H

#include "preselem.h"
// #include "hfluid.h"

/* Pipe - begin */

class Pipe : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   doublereal diameter;
   doublereal viscosity; /* e' globale perche' non dipende dalla pressione */
   doublereal area;
   doublereal length;
   flag turbulent;
   doublereal q0;

   doublereal vel;
   doublereal flow;
   doublereal Re;
   doublereal ktrb;   /* coefficiente per il moto turbolento */
   doublereal klam;   /* coefficiente per il moto laminare */
   doublereal ktra;   /* coefficiente per il moto di transizione */
   
   
 public:
   Pipe(unsigned int uL, const DofOwner* pD,
	HydraulicFluid* hf,const PressureNode* p1, const PressureNode* p2, 
	doublereal Dh, doublereal A, 
	doublereal L, flag transition, doublereal q0, flag fOut);
   
   ~Pipe(void);
   
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
   
   virtual void AfterConvergence(const VectorHandler& X, 
		   const VectorHandler& XP);
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

/* Pipe - end */


/* Dynamic_pipe Tubo ad elementi finiti - begin */

class Dynamic_pipe : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   doublereal diameter;
   doublereal area;
   doublereal length;
   flag turbulent;
   doublereal q0;              /* portata iniziale */
   
   doublereal flow1;           /* portata nodo 1  utilizzata per l'output */
   doublereal flow2;           /* portata nodo 2  utilizzata per l'output */
   doublereal Re;              /* numero di Reynolds medio */
   doublereal ktrb;            /* coefficiente per il moto turbolento */
   doublereal klam;            /* coefficiente per il moto laminare */
   doublereal ktra;            /* coefficiente per il moto di transizione */
   doublereal densitySt;       /* densita' dipendente dalla pressione inizio per l'output */
   doublereal densityMe;       /* densita' dipendente dalla pressione meta' per l'output*/
   doublereal densityEn;       /* densita' dipendente dalla pressione fine per l'output */
   doublereal densityDPres;    /* densita' diviso beta */
   
   doublereal viscosity;       /* viscosita' */

   doublereal VelS;            /* velocita' all'inizio del tubo */
   doublereal VelM;            /* velocita' a meta' del tubo */
   doublereal VelE;            /* velocita' alla fine del tubo */
   doublereal fa;              /* coefficiente di attrito */
   
   doublereal pp;              /* output */
   
 public:
   Dynamic_pipe(unsigned int uL, const DofOwner* pD, HydraulicFluid* hf,
	const PressureNode* p1, const PressureNode* p2, 
	doublereal Dh,
	doublereal A, doublereal L, 
	flag transition, doublereal q0, flag fOut);
   
   ~Dynamic_pipe(void);
   
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
   
   virtual void AfterConvergence(const VectorHandler& X, 
		   const VectorHandler& XP);
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

/* Dynamic_pipe - end */


/* DynamicPipe (tubo ad elementi finiti) - begin */

class DynamicPipe : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   
   doublereal diameter;
   doublereal area;
   doublereal length;
   
   flag turbulent;
   doublereal q0;              /* portata iniziale */
   
   doublereal Re;              /* numero di Reynolds medio */
   
   doublereal p1;
   doublereal p2;
   doublereal p1p;
   doublereal p2p;
   
   doublereal q1;
   doublereal q2;
   doublereal q1p;
   doublereal q2p;

   doublereal density1;       /* densita' dipendente dalla pressione inizio per l'output */
   doublereal density0;       /* densita' dipendente dalla pressione meta' per l'output*/
   doublereal density2;       /* densita' dipendente dalla pressione fine per l'output */
   doublereal densityDPres1;    /* densita' diviso beta */
   doublereal densityDPres2;    /* densita' diviso beta */
   
   doublereal viscosity;       /* viscosita' */
   
   doublereal dKlam;
   doublereal dKtrb;
   doublereal dKtra;
   
 public:
   DynamicPipe(unsigned int uL, const DofOwner* pD, HydraulicFluid* hf,
	       const PressureNode* p1, const PressureNode* p2, 
	       doublereal Dh,
	       doublereal A, doublereal L,
	       flag transition, doublereal q0, flag fOut);
   
   ~DynamicPipe(void);
   
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicElem::Type GetHydraulicType(void) const;
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   virtual DofOrder::Order GetEqType(unsigned int i) const;
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void AfterConvergence(const VectorHandler& X, 
		   const VectorHandler& XP);
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

/* Dynamic_pipe - end */

#endif
