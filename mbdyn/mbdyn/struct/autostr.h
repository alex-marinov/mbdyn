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

/* Elemento strutturale automatico */

#ifndef AUTOSTR_H
#define AUTOSTR_H

#include "elem.h"
#include "strnode.h"
#include "matvec3.h"

/* AutomaticStructElem - begin */

class AutomaticStructElem : virtual public Elem {
   friend class DynamicStructNode;
   
 protected:
   const DynamicStructNode* pNode;
   Vec3 Q;
   Vec3 G;
   Vec3 QP;
   Vec3 GP;

   /* Accesso ai suoi dati */
   virtual inline const Vec3& GetQCurr(void) const { return Q; };
   virtual inline const Vec3& GetGCurr(void) const { return G; };     
   virtual inline const Vec3& GetQPCurr(void) const { return QP; };
   virtual inline const Vec3& GetGPCurr(void) const { return GP; };   
   
 public:
   AutomaticStructElem(const DynamicStructNode* pN);
   
   virtual ~AutomaticStructElem(void) { 
      NO_OP;
   };
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* inizializza i dati */
   void Init(const Vec3& q, const Vec3& g, const Vec3& qp, const Vec3& gp);

   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual ElemType::Type GetElemType(void) const { 
      return ElemType::AUTOMATICSTRUCTURAL; 
   };
   
   
   /* funzioni di servizio */

   /* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
    * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
    * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
    * Viene usato nella costruzione dei DofOwner e quindi deve essere 
    * indipendente da essi. In genere non comporta overhead in quanto il 
    * numero di dof aggiunti da un tipo e' una costante e non richede dati 
    * propri.
    * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
    * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
    * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
    * dof propri.
    * Il metodo SetDof() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */
   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 12; 
      *piNumCols = 6; 
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* assemblaggio eig */   
   void AssEig(VariableSubMatrixHandler& WorkMatA,
	       VariableSubMatrixHandler& WorkMatB,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);

   /* output; si assume che ogni tipo di elemento sappia, attraverso
    * l'OutputHandler, dove scrivere il proprio output */
   virtual void Output(OutputHandler& OH) const;
   
   /* Setta i valori iniziali delle variabili (e fa altre cose) 
    * prima di iniziare l'integrazione */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;  

    /* *******PER IL SOLUTORE PARALLELO******** */
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, NodeType::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();
   };
   /* ************************************************ */    
};

/* AutomaticStructElem - end */

#endif
