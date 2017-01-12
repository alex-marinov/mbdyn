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

/* Vincolo di contatto in un piano */

#ifndef INPLCNTJ_H
#define INPLCNTJ_H

#include "joint.h"

/* InPlaneContact - begin */

/* Vincolo di contatto nel piano:
 * analogo al vincolo di giacenza nel piano per quanto riguarda la struttura,
 * tuttavia in questo caso i due corpi non sono obbligati ad essere 
 * in contatto. La loro distanza deve essere positiva nella direzione 
 * del vettore v o nulla. Il prodotto della distanza per la reazione vincolare
 * deve essere nullo sempre. Questo implica che la reazione puo' essere 
 * presente solo quando i due corpi sono in contatto. 
 * Valutare la possibilita' di introdurre anche una forza di attrito per 
 * strisciamento (anche nel caso del semplice contatto) */

class InPlaneContactJoint : public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 v;
   Vec3 p;
   doublereal dF;
   
 public:
   /* Costruttore banale */
   InPlaneContactJoint(unsigned int uL, const DofOwner* pDO);
   
   /* Costruttore non banale */
   InPlaneContactJoint(unsigned int uL, const DofOwner* pDO,
		       const StructNode* pN1, const StructNode* pN2, 
		       const Vec3& vTmp, const Vec3& pTmp);
   
   ~InPlaneContactJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const 
     { return Joint::INPLANECONTACT; };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
      
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 1);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 13;
      *piNumCols = 13; 
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void Output(OutputHandler& OH) const;

   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { return 2; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const 
     { *piNumRows = 26; *piNumCols = 26; };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X);

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

/* InPlaneContactJoint - end */

#endif
