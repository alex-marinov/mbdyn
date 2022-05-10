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

/* Giunti prismatici */

#ifndef PRISMJ_H
#define PRISMJ_H

#include "joint.h"
#include "drive.h"

/* PrismaticJoint - begin */

class PrismaticJoint : virtual public Elem, public Joint {
 private:
   /* Giunto prismatico - vincola due corpi alla traslazione 
    * l'uno rispetto all'altro. Il sistema di riferimento del prisma e' dato
    * rispetto a quelli dei nodi associati ai corpi.
    * In particolare rispetto al nodo 1 la trasformazione dal sistema
    * di riferimento del prisma al sistema globale e': R1*R1h, mentre per
    * il nodo 2 la medesima trasformazion e': R2*R2h.
    * Il vettore M esprie le reazioni vincolari di coppia. */
   const StructNode* pNode1;
   const StructNode* pNode2;
   Mat3x3 R1h;
   Mat3x3 R2h;
   Vec3 M;
   
 public:
   /* Costruttore non banale */
   PrismaticJoint(unsigned int uL, const DofOwner* pDO,
		  const StructNode* pN1, const StructNode* pN2,		 
		  const Mat3x3& R1hTmp, const Mat3x3& R2hTmp, flag fOut);
   
   /* Distruttore */
   ~PrismaticJoint(void);

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::PRISMATIC; 
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 3;
   };
   
   DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 5);
      return DofOrder::ALGEBRAIC; 
   };

   DofOrder::Order GetEqType(unsigned int i) const;

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 9;
      *piNumCols = 9; 
   };
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void OutputPrepare(OutputHandler& OH);
   void Output(OutputHandler& OH) const;
 
   void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 6;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows, 
				    integer* piNumCols) const {
      *piNumRows = 18; 
      *piNumCols = 18; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);   

#ifdef DEBUG
   virtual const char* sClassName(void) const { 
      return "PrismaticJoint";
   };
#endif 


   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;

   /* describes the dimension of components of equation */
   virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* PrismaticJoint - end */

#endif
