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

#ifndef HMINOR_H
#define HMINOR_H

#include "preselem.h"

/* Minor_loss - begin */

class Minor_loss : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   
   doublereal dKappa1;
   doublereal dKappa2;
   doublereal area;

   doublereal flow;  /* utilizzato per l'output */
   doublereal vel;   /* utilizzato per l'output */
   doublereal dKappa;
 
 public:
   Minor_loss(unsigned int uL, const DofOwner* pD,
		HydraulicFluid* hf,
		const PressureNode* p1, const PressureNode* p2,
		doublereal dK1,	doublereal dK2,  doublereal A, flag fOut);
   
   ~Minor_loss(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicType::Type GetHydraulicType(void) const;

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order SetDof(unsigned int i) const;
   
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
   
   virtual void SetValue(VectorHandler& X, VectorHandler& XP ) const;

  /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, NodeType::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */
};

/* Minor_loss - end */


/* Orifice - begin */

class Orifice : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNode1;
   const PressureNode* pNode2;
   doublereal diameter;
   doublereal viscosity;
   doublereal area_diaf;
   doublereal area_pipe;
   doublereal ReCr;
  
   doublereal CriticJump;
   doublereal delta;
   doublereal flow;  /* utilizzato per l'output */
   doublereal vel;   /* utilizzato per l'output */
   doublereal Re;    /* utilizzato per l'output */
      
 public:
   Orifice(unsigned int uL, const DofOwner* pD,
	   HydraulicFluid* hf,
	   const PressureNode* p1, const PressureNode* p2, 
	   doublereal Dh,
	   doublereal A_diaf, doublereal A_pipe, doublereal ReCR, flag fOut);
   
   ~Orifice(void);
   
   virtual inline void* pGet(void) const {
      return (void*)this;
   };
   
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicType::Type GetHydraulicType(void) const;

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order SetDof(unsigned int i) const;
   
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
   
   virtual void SetValue(VectorHandler& X, VectorHandler& XP ) const;

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, NodeType::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pNode1->GetNodeType();
     NdLabels[0] = pNode1->GetLabel();
     NdTyps[1] = pNode2->GetNodeType();
     NdLabels[1] = pNode2->GetLabel();
   };
   /* ************************************************ */
};

/* Orifice - end */

#endif
