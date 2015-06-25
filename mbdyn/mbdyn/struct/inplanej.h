/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* Vincolo di giacitura in un piano */

#ifndef INPLANEJ_H
#define INPLANEJ_H

#include "joint.h"

/* InPlaneJoint - begin */

/* Vincolo di giacenza in un piano:
 * si impone che il corpo b giaccia in un piano solidale al corpo a, 
 * definito dalla posizione di un punto p
 * e da una direzione a partire da tale punto, v.
 * La reazione vincolare e' una forza di modulo F e direzione v 
 * applicata nel punto di contatto. Quindi per ipotesi la forza e' applicata
 * al sorpo b e non da' momento su di esso, mentre e' applicata nel punto 
 * di contatto al corpo a e da' un momento pari alla distanza tra il punto
 * di contatto ed il corpo a moltiplicato per la forza.
 */

class InPlaneJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   
   const Vec3 v;
   const Vec3 p;
   
   doublereal dF;
   
 public:
   /* Costruttore non banale */
   InPlaneJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2, 
		const Vec3& vTmp, const Vec3& pTmp, flag fOut);
   
   ~InPlaneJoint(void);

   /* Tipo di joint */
   virtual Joint::Type GetJointType(void) const {
      return INPLANE;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
   
#ifdef DEBUG
   virtual DofOrder::Order GetDofType(unsigned int i) const
#else
   virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i == 0);
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

#ifdef DEBUG
   virtual const char* sClassName(void) const { return "InPlaneJoint"; };
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

};

/* InPlaneJoint - end */


/* InPlaneWithOffsetJoint - begin */

/* Vincolo di giacenza in un piano, con offset:
 * si impone che il punto q del corpo b giaccia in un piano solidale 
 * al corpo a, definito dalla posizione di un punto p
 * e da una direzione a partire da tale punto, v.
 * La reazione vincolare e' una forza di modulo F e direzione v 
 * applicata nel punto di contatto. 
 */

class InPlaneWithOffsetJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   
   const Vec3 v; /* Direzione normale al piano */
   const Vec3 p; /* Offset del piano dal nodo 1 */
   const Vec3 q; /* Offset del punto di contatto dal nodo 2 */
   
   doublereal dF;
   
 public:
   /* Costruttore non banale */
   InPlaneWithOffsetJoint(unsigned int uL, const DofOwner* pDO,
                          const StructNode* pN1, const StructNode* pN2,
                          const Vec3& vT, const Vec3& pT, const Vec3& qT,
			  flag fOut);
   
   ~InPlaneWithOffsetJoint(void);

   /* Tipo di joint */
   virtual Joint::Type GetJointType(void) const {
      return INPLANE;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
   
#ifdef DEBUG
   virtual DofOrder::Order GetDofType(unsigned int i) const
#else
   virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const
#endif
   {
      ASSERT(i >= 0 && i < 1);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 13; *piNumCols = 13; };
   
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

#ifdef DEBUG
   virtual const char* sClassName(void) const { return "InPlaneWithOffsetJoint"; };
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
};

/* InPlaneWithOffsetJoint - end */

#endif
