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

#ifndef HMINOR_H
#define HMINOR_H

#include "preselem.h"

/* MinorLoss - begin */

class MinorLoss : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* m_pNode1;
   const PressureNode* m_pNode2;
   
   doublereal m_dKappa12;
   doublereal m_dKappa21;
   doublereal m_Area;

   doublereal flow;  /* utilizzato per l'output */
   doublereal vel;   /* utilizzato per l'output */
   doublereal m_dKappa;
 
 public:
   MinorLoss(unsigned int uL, const DofOwner* pD,
		HydraulicFluid* hf,
		const PressureNode* p1, const PressureNode* p2,
		doublereal dK12, doublereal dK21, doublereal A, flag fOut);
   
   ~MinorLoss(void);
   
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
		connectedNodes.resize(2);
		connectedNodes[0] = m_pNode1;
		connectedNodes[1] = m_pNode2;
	};
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* MinorLoss - end */


/* ThreeWayMinorLoss - begin */

class ThreeWayMinorLoss : virtual public Elem, public HydraulicElem {
private:
	const PressureNode* m_pNode0;
	const PressureNode* m_pNode1;
	const PressureNode* m_pNode2;
	const PressureNode* m_pNodeN;

	doublereal m_dKappa12;
	doublereal m_dKappa21;
	doublereal m_Area1;
	doublereal m_Area2;
	doublereal m_Area;
	
	doublereal flow;  /* utilizzato per l'output */
	doublereal vel;   /* utilizzato per l'output */
	doublereal m_dKappa;

public:
	ThreeWayMinorLoss(unsigned int uL, const DofOwner* pD,
		HydraulicFluid* hf, const PressureNode* p0,
		const PressureNode* p1, const PressureNode* p2,
		doublereal dK12, doublereal dK21, 
		doublereal A1, doublereal A2, flag fOut);

	~ThreeWayMinorLoss(void);

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
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(3);
		connectedNodes[0] = m_pNode0;
		connectedNodes[1] = m_pNode1;
		connectedNodes[2] = m_pNode2;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* ThreeWayMinorLoss - end */


/* Orifice - begin */

class Orifice : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* m_pNode1;
   const PressureNode* m_pNode2;
   doublereal diameter;
   doublereal viscosity;
   doublereal m_Area_diaf;
   doublereal m_Area_pipe;
   doublereal ReCr;
  
   doublereal CriticJump;
   doublereal delta;
   doublereal Cd;
   doublereal flow;  /* utilizzato per l'output */
   doublereal vel;   /* utilizzato per l'output */
   doublereal Re;    /* utilizzato per l'output */
   bool turbulent;
      
 public:
   Orifice(unsigned int uL, const DofOwner* pD,
	   HydraulicFluid* hf,
	   const PressureNode* p1, const PressureNode* p2, 
	   doublereal Dh,
	   doublereal A_diaf, doublereal A_pipe, doublereal ReCR, flag fOut);
   
   ~Orifice(void);
   
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
     connectedNodes.resize(2);
     connectedNodes[0] = m_pNode1;
     connectedNodes[1] = m_pNode2;
   };
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* Orifice - end */

#endif
