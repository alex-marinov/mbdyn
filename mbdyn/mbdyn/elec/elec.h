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

/* elementi elettrici, tipo: Elem::ELECTRIC */

#ifndef ELEC_H
#define ELEC_H


/* include per derivazione della classe */

#include "elem.h"
#include "elecnode.h"
#include "strnode.h"
#include "drive.h"
      
extern const char* psElectricNames[];


/* Electric - begin */

class Electric : virtual public Elem, public ElemWithDofs {
 public:
   /* Tipi di elementi elettrici */
   enum Type {
      UNKNOWN = -1,
      
      ACCELEROMETER = 0,
      DISPLACEMENT,
      DISCRETECONTROL,
      
      LASTELECTRICTYPE
   };

 private:
   Electric::Type ElecT;
   
 public:
   Electric(unsigned int uL, Electric::Type T, 
	    const DofOwner* pDO, flag fOut);
   virtual ~Electric(void);
   
   /* Contributo al file di restart 
    * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
    * relativa alla classe derivata */
   virtual ostream& Restart(ostream& out) const;

   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual Elem::Type GetElemType(void) const;

   /* Tipo di elemento elettrico (usato solo per debug ecc.) */
   virtual Electric::Type GetElectric(void) const = 0;

   /* Output */
   virtual void Output(OutputHandler& OH) const;
};

/* Electric - end */


#if defined(USE_STRUCT_NODES)

/* Accelerometer - begin */

class Accelerometer : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   doublereal dOmega;
   doublereal dTau;
   doublereal dCsi;
   doublereal dKappa;
   
 public:
   Accelerometer(unsigned int uL, const DofOwner* pD, 
		 const StructNode* pS, const AbstractNode* pA,
		 const Vec3& TmpDir,
		 doublereal dO, doublereal dT, doublereal dC, doublereal dK,
		 flag fOut);
   ~Accelerometer(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectric(void) const {
      return Electric::ACCELEROMETER;
   };
   
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
   
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const;

 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pStrNode->GetNodeType();
     NdLabels[0] = pStrNode->GetLabel();
     NdTyps[1] = pAbsNode->GetNodeType();
     NdLabels[1] = pAbsNode->GetLabel();
   };
   /* ************************************************ */
};


inline void* Accelerometer::pGet(void) const
{
   return (void*)this;
}

/* Accelerometer - end */


/* TraslAccel - begin */

class TraslAccel : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   Vec3 f;
   
 public:
   TraslAccel(unsigned int uL, const DofOwner* pD, 
	      const StructNode* pS, const AbstractNode* pA,
	      const Vec3& TmpDir, const Vec3& Tmpf,
	      flag fOut);
   ~TraslAccel(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectric(void) const {
      return Electric::ACCELEROMETER;
   };
   
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
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const;

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pStrNode->GetNodeType();
     NdLabels[0] = pStrNode->GetLabel();
     NdTyps[1] = pAbsNode->GetNodeType();
     NdLabels[1] = pAbsNode->GetLabel();
   };
   /* ************************************************ */
};


inline void* TraslAccel::pGet(void) const 
{
   return (void*)this;
}
   
/* TraslAccel - end */


/* RotAccel - begin */

class RotAccel : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   
 public:
   RotAccel(unsigned int uL, const DofOwner* pD, 
	      const StructNode* pS, const AbstractNode* pA,
	      const Vec3& TmpDir,
	      flag fOut);
   ~RotAccel(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectric(void) const {
      return Electric::ACCELEROMETER;
   };
   
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
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const;

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pStrNode->GetNodeType();
     NdLabels[0] = pStrNode->GetLabel();
     NdTyps[1] = pAbsNode->GetNodeType();
     NdLabels[1] = pAbsNode->GetLabel();
   };
   /* ************************************************ */
};


inline void* RotAccel::pGet(void) const 
{
   return (void*)this;
}
   
/* RotAccel - end */

/* DispMeasure - begin */

class DispMeasure : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode1;
   const StructNode* pStrNode2;
   const AbstractNode* pAbsNode;
   Vec3 f1;
   Vec3 f2;
   
 public:
   DispMeasure(unsigned int uL, const DofOwner* pD,
	       const StructNode* pS1, const StructNode* pS2, 
	       const AbstractNode* pA,
	       const Vec3& Tmpf1, const Vec3& Tmpf2,
	       flag fOut);
   ~DispMeasure(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectric(void) const {
      return Electric::DISPLACEMENT;
   };
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
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
   
   /* Setta i valori iniziali delle variabili (e fa altre cose)
    * prima di iniziare l'integrazione */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 3;
     NdTyps[0] = pStrNode1->GetNodeType();
     NdLabels[0] = pStrNode1->GetLabel();
     NdTyps[1] = pStrNode2->GetNodeType();
     NdLabels[1] = pStrNode2->GetLabel();
     NdTyps[2] = pAbsNode->GetNodeType();
     NdLabels[2] = pAbsNode->GetLabel();
   };
   /* ************************************************ */
};


inline void* DispMeasure::pGet(void) const 
{
   return (void*)this;
}
   
/* DispMeasure - end */

#endif // USE_STRUCT_NODES

class DataManager;
class MBDynParser;

extern Elem* ReadElectric(DataManager* pDM,
			  MBDynParser& HP, 
			  const DofOwner* pDO, 
			  unsigned int uLabel);

#endif
