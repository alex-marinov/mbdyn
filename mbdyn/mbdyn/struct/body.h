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

/* elementi di massa, tipo: Elem::Type BODY */

#ifndef BODY_H
#define BODY_H

/* include per derivazione della classe */

#include "elem.h"
#include "strnode.h"
#include "gravity.h"

/* Body - begin */

class Body : 
virtual public Elem, public ElemGravityOwner, public InitialAssemblyElem {   
  private:
    const StructNode* pNode;
    doublereal dMass;
    Vec3 Xgc;
    Vec3 S0;
    Mat3x3 J0;
    
    mutable Vec3 S;
    mutable Mat3x3 J;
   
    /* Assembla le due matrici necessarie per il calcolo degli
     * autovalori e per lo jacobiano */  
    void AssMat_(FullSubMatrixHandler& WorkMatA,
		 FullSubMatrixHandler& WorkMatB,
		 doublereal dCoef,
		 const VectorHandler& XCurr,
		 const VectorHandler& XPrimeCurr);

    /* momento statico */
    Vec3 _GetS(void) const;

    /* momento d'inerzia */
    Mat3x3 _GetJ(void) const;
      
  public:
    /* Costruttore definitivo (da mettere a punto) */
    Body(unsigned int uL, const StructNode* pNodeTmp, 
	 doublereal dMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp, 
	 flag fOut);
   
    /* virtual ~Body(void); */
   
    virtual inline void* pGet(void) const { 
        return (void*)this;
    };

    /* massa totale */
    doublereal dGetM(void) const {
        return dMass;
    };
   
    /* Scrive il contributo dell'elemento al file di restart */
    virtual std::ostream& Restart(std::ostream& out) const;
   
    /* Tipo dell'elemento (usato solo per debug ecc.) */
    virtual Elem::Type GetElemType(void) const {
        return Elem::BODY; 
    };
   
    void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
        *piNumRows = 12; 
        *piNumCols = 6; 
    };
   
    VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				     doublereal dCoef,
				     const VectorHandler& XCurr, 
				     const VectorHandler& XPrimeCurr);
   
    void AssMats(VariableSubMatrixHandler& WorkMatA,
	        VariableSubMatrixHandler& WorkMatB,
	        const VectorHandler& XCurr,
	        const VectorHandler& XPrimeCurr);
   
    SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			     doublereal dCoef,
			     const VectorHandler& XCurr, 
			     const VectorHandler& XPrimeCurr);
   
    /* Numero gdl durante l'assemblaggio iniziale */
    virtual unsigned int iGetInitialNumDof(void) const { 
        return 0; 
    };
   
    /* Dimensione del workspace durante l'assemblaggio iniziale.
     * Occorre tener conto del numero di dof che l'elemento definisce
     * in questa fase e dei dof dei nodi che vengono utilizzati.
     * Sono considerati dof indipendenti la posizione e la velocita'
     * dei nodi */
    virtual void 
    InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
        *piNumRows = 12; 
        *piNumCols = 6; 
    };
   
    /* Contributo allo jacobiano durante l'assemblaggio iniziale */
    virtual VariableSubMatrixHandler& 
    InitialAssJac(VariableSubMatrixHandler& WorkMat,
                  const VectorHandler& XCurr);
   
    /* Contributo al residuo durante l'assemblaggio iniziale */   
    virtual SubVectorHandler& 
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

    /* Usata per inizializzare la quantita' di moto */
    virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;

    virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

    /******** PER IL SOLUTORE PARALLELO *********/        
    /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
     * utile per l'assemblaggio della matrice di connessione fra i dofs */
    virtual void 
    GetConnectedNodes(int& NumNodes, 
                      Node::Type* NdTyps, 
		      unsigned int* NdLabels) {
        NumNodes = 1;
        NdTyps[0] = pNode->GetNodeType();
        NdLabels[0] = pNode->GetLabel();
    };
    /**************************************************/
};

/* Body - end */

class DataManager;
class MBDynParser;

extern Elem* ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* BODY_H */

