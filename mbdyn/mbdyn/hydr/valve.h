/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef VALVE_H
#define VALVE_H

#include <preselem.h>

/* Control_valve - begin */

class Control_valve : virtual public Elem, public HydraulicElem, public DriveOwner {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   const PressureNode* pNode3;
   const PressureNode* pNode4;
   doublereal Stato;
   doublereal Cd;         /* coefficiente di perdita */
   doublereal area_max;   /* larghezza del condotto: A=x*area_max */
   doublereal loss_area;  /* area di trafilamento in % sull'area max*/
   doublereal s_max; 
   doublereal flow1;      /* portata nodo 1 (per l'output) */
   doublereal flow2;      /* portata nodo 2 (per l'output) */
   doublereal flow3;      /* portata nodo 3 (per l'output) */
   doublereal flow4;      /* portata nodo 4 (per l'output) */
   doublereal A1;
   doublereal A2;
   doublereal A3;
   doublereal A4;
   doublereal A1min;
   doublereal A2min;
   doublereal A3min;
   doublereal A4min;
   
 public:
   Control_valve(unsigned int uL, const DofOwner* pD, HydraulicFluid* hf,
		 const PressureNode* p1, const PressureNode* p2, 
		 const PressureNode* p3, const PressureNode* p4, 
		 doublereal A_max, doublereal Loss_a, const DriveCaller* pDC,flag fOut);
   
   ~Control_valve(void);
   
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
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(4);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
     connectedNodes[2] = pNode3;
     connectedNodes[3] = pNode4;
   };
   /* ************************************************ */ 	
};

/* Control_valve - end */


/* Control_valve2 - begin */

#define VALVE_6
// #undef VALVE_6

class Control_valve2 
: virtual public Elem, public HydraulicElem, public DriveOwner {
private:
	 /*
	  *
	  * Node 1    Node 4
	  *          __
	  *   |  \     | / \
	  *   |   \  /    |
	  *   |    \/     |
	  *   |    /\     |
	  *  \ /  / __|   |
	  *  
	  * Node 2     Node 3
	  * 
	  */
	enum {
		N1 = 0,
		N2 = 1,
		N3 = 2,
		N4 = 3,
		LAST_N = 4
	};
	enum {
		Q12 = 0,
		Q34 = 1,
		Q13 = 2,
		Q24 = 3,
#ifdef VALVE_6
		Q14 = 4,
		Q23 = 5,
		LAST_Q = 6
#else /* !VALVE_6 */
		LAST_Q = 4
#endif /* !VALVE_6 */
	};

	const PressureNode* pNode[LAST_N];
	doublereal q[LAST_Q];
	doublereal dp[LAST_Q];
	doublereal density;

	doublereal Stato;
	doublereal Cd;		/* coefficiente di perdita */
	doublereal area_max;	/* larghezza del condotto: A=x*area_max */
	doublereal loss_area;	/* area di trafilamento in % sull'area max */
	doublereal area_min;	/* area di trafilamento = area_max*loss_area */
	doublereal s_max; 
	doublereal f[LAST_N];
	doublereal A[LAST_Q];

	/* gathers data from nodes and computes intermediate member data */
	void Prepare(void);
   
public:
	Control_valve2(unsigned int uL, const DofOwner* pD, HydraulicFluid* hf,
			const PressureNode* p1, const PressureNode* p2, 
			const PressureNode* p3, const PressureNode* p4, 
			doublereal A_max, doublereal Loss_a, 
			const DriveCaller* pDC,flag fOut);

	~Control_valve2(void);
	
	/* Tipo di elemento idraulico (usato solo per debug ecc.) */
	virtual HydraulicElem::Type GetHydraulicType(void) const;
	
	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef, const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
	
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef, const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
	
	virtual void Output(OutputHandler& OH) const;
	
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);
	
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(4);
		connectedNodes[0] = pNode[N1];
		connectedNodes[1] = pNode[N2];
		connectedNodes[2] = pNode[N3];
		connectedNodes[3] = pNode[N4];
	};
	/* ************************************************ */ 	
};

/* Control_valve2 - end */


/* Control_valve_dinamica - begin */

class Dynamic_control_valve : virtual public Elem, public HydraulicElem, public DriveOwner {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   const PressureNode* pNode3;
   const PressureNode* pNode4;
   doublereal start;        /* posizione iniziale della valvola */
   doublereal c_spost;      /* coef dello spostamento */
   doublereal c_vel;        /* coef della velocita' */
   doublereal c_acc;        /* coef dell'accelerazione */
   doublereal Force;        /* forza applicata alla valvola */
    
   doublereal Cd;             /* coefficiente di perdita */
   doublereal width;          /* larghezza del condotto: A=x*width */
   doublereal loss_area;      /* area di trafilamento in % sull'area max*/
   doublereal valve_diameter; /* diametro della valvola */
   doublereal valve_density;  /* densita' del corpo della valvola */
  
   doublereal s_max;        /* corsa massima della valvola */

   doublereal A1;           /* area di collegamento nodo 1 & nodo 2 */
   doublereal A2;           /* area di collegamento nodo 1 & nodo 3 */
   doublereal A3;           /* area di collegamento nodo 3 & nodo 4 */
   doublereal A4;           /* area di collegamento nodo 2 & nodo 4 */
  
   doublereal Mass;         /* massa del corpo della valvola */
   doublereal flow1;        /* portata nodo 1 (per l'output) */
   doublereal flow2;        /* portata nodo 2 (per l'output) */
   doublereal flow3;        /* portata nodo 3 (per l'output) */
   doublereal flow4;        /* portata nodo 4 (per l'output) */
   doublereal s;            /* spostamento valvola    (per l'output) */
   doublereal v;            /* velocita' valvola      (per l'output) */
   doublereal sp;           /* velocita' valvola */
   doublereal vp;           /* accelerazione valvola  (per l'output) */
   doublereal c1;           /* coef dello spostamento */
   doublereal c2;           /* coef della velocita' */
   doublereal c3;           /* coef dell'accelerazione */
   doublereal cf1;          /* coef dello spostamento finale */
   doublereal cf2;          /* coef della velocita' finale */
   doublereal cf3;          /* coef dell'accelerazione finale */
   doublereal deltaP;       /* salto di pressione p1-p2+p3 */
   
 public:
   Dynamic_control_valve(unsigned int uL, const DofOwner* pD,
			 HydraulicFluid* hf,
			 const PressureNode* p1, const PressureNode* p2, 
			 const PressureNode* p3, const PressureNode* p4, 
			 const DriveCaller* pDC, doublereal s0, 
			 doublereal s_mx, doublereal W, doublereal Loss_A, 
			 doublereal Valve_d, doublereal Valve_rho, 
			 doublereal cs, doublereal cv, doublereal ca, flag fOut);
   
   ~Dynamic_control_valve(void);
   
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
     connectedNodes.resize(4);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
     connectedNodes[2] = pNode3;
     connectedNodes[3] = pNode4;
   };
   /* ************************************************ */
};

/* Dynamic_control_valve - end */

/* Pressure flow control_valve_dinamica - begin */

class Pressure_flow_control_valve : virtual public Elem, public HydraulicElem, public DriveOwner {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   const PressureNode* pNode3;
   const PressureNode* pNode4;
   const PressureNode* pNode5;
   const PressureNode* pNode6;
   doublereal start;        /* posizione iniziale della valvola */
   doublereal c_spost;      /* coef dello spostamento */
   doublereal c_vel;        /* coef della velocita' */
   doublereal c_acc;        /* coef dell'accelerazione */
   doublereal Force;        /* forza applicata alla valvola */
    
   doublereal Cd;             /* coefficiente di perdita */
   doublereal width;          /* larghezza del condotto: A=x*width */
   doublereal loss_area;      /* area di trafilamento in % sull'area max*/
   doublereal valve_diameter; /* diametro della valvola */
   doublereal valve_density;  /* densita' del corpo della valvola */
   doublereal valve_area;     /* area della valvola */  
   doublereal s_max;        /* corsa massima della valvola */

   doublereal A1;           /* area di collegamento nodo 1 & nodo 2 */
   doublereal A2;           /* area di collegamento nodo 1 & nodo 3 */
   doublereal A3;           /* area di collegamento nodo 3 & nodo 4 */
   doublereal A4;           /* area di collegamento nodo 2 & nodo 4 */
  
   doublereal Mass;         /* massa del corpo della valvola */
   doublereal flow1;        /* portata nodo 1 (per l'output) */
   doublereal flow2;        /* portata nodo 2 (per l'output) */
   doublereal flow3;        /* portata nodo 3 (per l'output) */
   doublereal flow4;        /* portata nodo 4 (per l'output) */
   doublereal flow5;        /* portata nodo 5 (per l'output) */
   doublereal flow6;        /* portata nodo 6 (per l'output) */
   
   doublereal s;            /* spostamento valvola    (per l'output) */
   doublereal v;            /* velocita' valvola      (per l'output) */
   doublereal sp;           /* velocita' valvola */
   doublereal vp;           /* accelerazione valvola  (per l'output) */
   doublereal c1;           /* coef dello spostamento */
   doublereal c2;           /* coef della velocita' */
   doublereal c3;           /* coef dell'accelerazione */
   doublereal cf1;          /* coef dello spostamento finale */
   doublereal cf2;          /* coef della velocita' finale */
   doublereal cf3;          /* coef dell'accelerazione finale */
   doublereal deltaP;       /* salto di pressione p1-p2+p3 */
   
 public:
   Pressure_flow_control_valve(unsigned int uL, const DofOwner* pD,
			 HydraulicFluid* hf,
			 const PressureNode* p1, const PressureNode* p2, 
			 const PressureNode* p3, const PressureNode* p4, 
			 const PressureNode* p5, const PressureNode* p6, 
		         const DriveCaller* pDC, doublereal s0, 
			 doublereal s_mx, doublereal W, doublereal Loss_A, 
			 doublereal Valve_d, doublereal Valve_rho, 
			 doublereal cs, doublereal cv, doublereal ca, flag fOut);
   
   ~Pressure_flow_control_valve(void);
   
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
     connectedNodes.resize(6);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
     connectedNodes[2] = pNode3;
     connectedNodes[3] = pNode4;
     connectedNodes[4] = pNode5;
     connectedNodes[5] = pNode6;
  };
   /* ************************************************ */
};

/* Pressure_flow_control_valve - end */


/* Pressure_valve - begin */

class Pressure_valve : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   doublereal area_diaf;  /* area diaframma */
   doublereal mass;       /* massa valvola */
   doublereal area_max;   /* area maggiore della valvola */
   doublereal Kappa;      /* costante della molla */
   doublereal force0;     /* precarico della molla */
   doublereal width;      /* larghezza luce di passaggio */
   doublereal s_max;      /* corsa massima della valvola */
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
   doublereal sp;         /* velocita' diaframma  */
   doublereal vp;         /* accelerazione diaframma  (per l'output) */
   doublereal flow1;      /* portata nodo 1           (per l'output) */
   doublereal flow2;      /* portata nodo 2           (per l'output) */
  
 public:
   Pressure_valve(unsigned int uL, const DofOwner* pD,
		  HydraulicFluid* hf,
		  const PressureNode* p1, const PressureNode* p2, 
		  doublereal A_dia, doublereal mv, 
		  doublereal A_max, doublereal s_mx, doublereal K, 
		  doublereal F0, doublereal w,
		  doublereal cs, doublereal cv, doublereal ca, flag fOut);
   
   ~Pressure_valve(void);
  
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

/* Pressure_valve - end */


/* Flow_valve - begin */

class  Flow_valve : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   const PressureNode* pNode3;
   doublereal area_diaf;  /* area diaframma */
   doublereal mass;       /* massa valvola */
   doublereal area_pipe;  /* area del tubo */
   doublereal area_max;   /* area maggiore della valvola */
   doublereal Kappa;      /* costante della molla */
   doublereal force0;     /* precarico della molla */
   doublereal width;      /* larghezza luce di passaggio */
   doublereal s_max;      /* corsa massima della valvola */
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

   doublereal h;          /* perdita di carico concentrata tra i nodi 1 e 2 (smorza il moto della valvola) */
   doublereal s;          /* spostamento valvola  (per l'output) */
   doublereal sp;         /* velocita' valvola  (per l'output) */
   doublereal v;          /* velocita' valvola  (per l'output) */
   doublereal vp;         /* accelerazione valvola  (per l'output) */
   doublereal flow1;      /* portata nodo 1  (per l'output) */
   doublereal flow2;      /* portata nodo 2  (per l'output) */
   doublereal flow3;      /* portata nodo 3  (per l'output) */
   
 public:
   Flow_valve(unsigned int uL, const DofOwner* pD, HydraulicFluid* hf,
              const PressureNode* p1, const PressureNode* p2, 
	      const PressureNode* p3, doublereal A_dia, 
	      doublereal mv, doublereal A_pipe,doublereal A_max, 
	      doublereal K, doublereal F0, doublereal w,doublereal s_mx,
	      doublereal cs, doublereal cv, doublereal ca,  flag fOut);
   
   ~Flow_valve(void);
   
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
     connectedNodes.resize(3);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
     connectedNodes[2] = pNode3;
   };
  /* ************************************************ */ 	 
};

/* Flow_valve - end */

#endif /* VALVE_H */

