/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* GENEL filtro scalare analogico */

#ifndef GENFILT_H
#define GENFILT_H

#include "genel.h"

#if 0
/* GenelFilter - begin */

class GenelFilter : public Genel {
 protected:
   ScalarDof SD_y; /* uscita */
   ScalarDof SD_u; /* ingresso */
   
   unsigned int Na;
   unsigned int Nb;
   
   unsigned int iNumDofs;
   
   doublereal* pdP;
   doublereal* pdTau;
   
 public:
   GenelFilter(unsigned int uLabel, const DofOwner* pDO, 
	       const ScalarDof& y, const ScalarDof& u,
	       unsigned int na, unsigned int nb,
	       doublereal* p, doublereal* tau,	     
	       flag fOutput);
   virtual ~GenelFilter(void);
   virtual inline void* pGet(void) const;
   
   virtual unsigned int iGetNumDof(void) const;   
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ );

   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SD_y.pNode->GetNodeType();
     NdLabels[0] = SD_y.pNode->GetLabel();
     NdTyps[1] = SD_u.pNode->GetNodeType();
     NdLabels[1] = SD_u.pNode->GetLabel();
   };
   /* ************************************************ */
};


inline void* GenelFilter::pGet(void) const
{
   return (void*)this;
}

/* GenelFilter - end */
#endif 


/* GenelFilterEq - begin */

class GenelFilterEq : public Genel {
 protected:
   ScalarDof SD_y; /* uscita */
   ScalarDof SD_u; /* ingresso */
   
   unsigned int Na;
   unsigned int Nb;
   
   doublereal* pdA;
   doublereal* pdB;
   
   doublereal* pdAlpha;
   doublereal* pdBeta;
   
   flag fSteady;
   
 public:
   GenelFilterEq(unsigned int uLabel, const DofOwner* pDO, 
		 const ScalarDof& y, const ScalarDof& u,
		 unsigned int na, unsigned int nb,
		 doublereal* pa, doublereal* pb,
		 flag fSt, flag fOutput);
   virtual ~GenelFilterEq(void);
   virtual inline void* pGet(void) const;
   
   virtual unsigned int iGetNumDof(void) const;   
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::SCALARFILTER; 
   };

   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ );

   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);

   /* Setta i valori iniziali delle variabili (e fa altre cose) 
    * prima di iniziare l'integrazione */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SD_y.pNode->GetNodeType();
     NdLabels[0] = SD_y.pNode->GetLabel();
     NdTyps[1] = SD_u.pNode->GetNodeType();
     NdLabels[1] = SD_u.pNode->GetLabel();
   };
   /* ************************************************ */

};


inline void* GenelFilterEq::pGet(void) const
{
   return (void*)this;
}

/* GenelFilterEq - end */


/* GenelStateSpaceSISO - begin */

class GenelStateSpaceSISO : public Genel {
 protected:
   ScalarDof SD_y; /* uscita */
   ScalarDof SD_u; /* ingresso */
   
   unsigned int iNumDofs;
   
   doublereal* pdA;
   doublereal* pdB;
   doublereal* pdC;
   doublereal dD;
   
   doublereal* pdX;
   doublereal* pdXP;
   
 public:
   GenelStateSpaceSISO(unsigned int uLabel, const DofOwner* pDO, 
		       const ScalarDof& y, const ScalarDof& u,
		       unsigned int Order,
		       doublereal* pA, doublereal* pB,	     
		       doublereal* pC, doublereal D,
		       flag fOutput);
   
   virtual ~GenelStateSpaceSISO(void);
   
   virtual inline void* pGet(void) const {
      return (void*)this;
   };
   
   virtual unsigned int iGetNumDof(void) const;
   
   /* esegue operazioni sui dof di proprieta' dell'elemento */
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::STATESPACESISO; 
   };
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ );

   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = SD_y.pNode->GetNodeType();
     NdLabels[0] = SD_y.pNode->GetLabel();
     NdTyps[1] = SD_u.pNode->GetNodeType();
     NdLabels[1] = SD_u.pNode->GetLabel();
   };
   /* ************************************************ */
};

/* GenelStateSpaceSISO - end */


/* GenelStateSpaceMIMO - begin */

class GenelStateSpaceMIMO : public Genel {
 protected:
   unsigned int iNumOutputs;
   unsigned int iNumInputs;
   ScalarDof* pvSD_y; /* uscite */
   ScalarDof* pvSD_u; /* ingressi */
   
   unsigned int iNumDofs;
   
   doublereal* pdA;
   doublereal* pdB;
   doublereal* pdC;
   doublereal* pdD;
   
   doublereal* pdX;
   doublereal* pdXP;
   
 public:
   GenelStateSpaceMIMO(unsigned int uLabel, const DofOwner* pDO,
		       unsigned int iNumOut, const ScalarDof* y,
		       unsigned int iNumIn, const ScalarDof* u,
		       unsigned int Order,
		       doublereal* pA, doublereal* pB,	     
		       doublereal* pC, doublereal* pD,
		       flag fOutput);
   
   virtual ~GenelStateSpaceMIMO(void);
   
   virtual inline void* pGet(void) const {
      return (void*)this;
   };
   
   virtual unsigned int iGetNumDof(void) const;
   
   /* esegue operazioni sui dof di proprieta' dell'elemento */
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const { 
      return Genel::STATESPACEMIMO; 
   };

   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ );

   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);

 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = iNumInputs +iNumOutputs;
     for(unsigned int i= 0; i <= iNumOutputs-1; i++) { 
       NdTyps[i] = (pvSD_y[i].pNode)->GetNodeType();
       NdLabels[i] = (pvSD_y[i].pNode)->GetLabel();
     }
     for(unsigned int i= 0; i <= iNumInputs-1; i++) {
       NdTyps[iNumOutputs+i] = (pvSD_u[i].pNode)->GetNodeType();
       NdLabels[iNumOutputs+i] = (pvSD_u[i].pNode)->GetLabel();
     }
   };
   /* ************************************************ */

};

/* GenelStateSpaceMIMO - end */

#endif
