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

/* Elemento modale */

/* 
 * Copyright 1999-2000 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef MODAL_H
#define MODAL_H

#include <ac/fstream>

#include <joint.h>

/* Modal - begin */

/* 
 * ATTENZIONE! 
 * per ora e' derivato da Joint; 
 * puo' darsi che venga creata una classe apposta
 */


class Modal : virtual public Elem, public Joint {
 protected:
   const ModalNode* pModalNode;
   unsigned int NModes;
   unsigned int NStrNodes;
   unsigned int NFemNodes;
   doublereal dMass;
   Vec3   Inv2;
   Mat3x3 Inv7;
   MatNxN *pModalMass;
   MatNxN *pModalStiff;
   MatNxN *pModalDamp;
   unsigned int* IdFemNodes;
   unsigned int* IntNodes;
   Mat3xN *pXYZFemNodes;
   Mat3xN *pOffsetNodes;   
   const StructNode** pInterfaceNodes;
 
   Mat3xN *pPHIt;
   Mat3xN *pPHIr;
   
   Mat3xN *pModeShapest;
   Mat3xN *pModeShapesr;
 
   Mat3xN *pInv3;
   Mat3xN *pInv4;
   Mat3xN *pInv5;
   Mat3xN *pInv8;
   Mat3xN *pInv9;

   Mat3xN *pInv10;
   Mat3xN *pInv11;
   
   Vec3   Inv3jaj;
   Vec3   Inv3jaPj;
   Mat3x3 Inv8jaj;
   Mat3x3 Inv8jaPj;
   Mat3xN Inv5jaj;
   Vec3 Inv4j, VInv5jaj, VInv5jaPj;
   Mat3x3 Inv8j, Inv9jkak, Inv9jkajaPk;

   VecN a;
   VecN b;
   VecN bPrime;

   Vec3* pd1tot;
   Vec3* pd2;
   Mat3x3* pR1tot;
   Mat3x3* pR2;
   Vec3* pF;
   Vec3* pM;
   mutable std::ofstream fOutFlex;
   
   integer iModalIndex;
      
 public:
   /* Costruttore non banale */
   Modal(unsigned int uL,
         const ModalNode* pModalNodeTmp, 
	 const DofOwner* pDO,
	 unsigned int N,
	 unsigned int NS,
         unsigned int NFN,
         doublereal dMass,
         const Vec3& STmp,
         const Mat3x3& JTmp,
         MatNxN *pGenMass,
         MatNxN *pGenStiff,
         MatNxN *pGenDamp,
         unsigned int* IdFemNodes,
         unsigned int* IntNodes,
	 Mat3xN *pN,
         Mat3xN *pXYZOffsetNodes,
         const StructNode** pN2,
         Mat3xN *pPHIt,
         Mat3xN *pPHIr,
	 Mat3xN *pModeShapest,
         Mat3xN *pModeShapesr,
         Mat3xN *pInv3,
         Mat3xN *pInv4,
         Mat3xN *pInv5,
         Mat3xN *pInv8,
         Mat3xN *pInv9,
         Mat3xN *pInv10,
         Mat3xN *pInv11,
	 VecN *a,
	 VecN *aP,
         const char *sFileMod,       
         DataManager* pDM,
         MBDynParser& HP,
         flag fOut);

   /* Distruttore */
   ~Modal(void);
   
   virtual inline void* pGet(void) const;
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const;
   DofOrder::Order SetDof(unsigned int i) const;
   
   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   void Output(OutputHandler& OH) const;

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const;
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const;
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   
   void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const;
  
   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual doublereal dGetPrivData(unsigned int i) const;   
   
   /* Funzioni che restituiscono dati che possono servire ad
    * altri elementi (ad es. agli elementi aerodinamici modali)
    */

   Mat3xN* pGetPHIt(void) {
      return pModeShapest;
   };

   Mat3xN* pGetPHIr(void) {
      return pModeShapesr;
   };
   
   Mat3xN* pGetFemNodesPosition(void) {
      return pXYZFemNodes;
   };
   
   integer uGetNModes(void) {
      return NModes;
   };
   
   integer uGetNFemNodes(void) {
      return NFemNodes;
   };
   
   integer iGetModalIndex(void) const {
      return iGetFirstIndex();
   };
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
    * utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
      NumNodes = 1+NStrNodes;
      NdTyps[0] = pModalNode[0].GetNodeType();
      NdLabels[0] = pModalNode[0].GetLabel();
      for(unsigned int j=0; j < NStrNodes; j++) {
	 NdTyps[1+j] = pInterfaceNodes[j]->GetNodeType();
	 NdLabels[1+j] = pInterfaceNodes[j]->GetLabel();
      }
   };
   /* ************************************************ */
};


inline void* Modal::pGet(void) const
{
   return (void*)this;
}

/* Modal - end */

class DataManager;
class MBDynParser;

extern Joint* ReadModal(DataManager* pDM, MBDynParser& HP, const DofOwner* pD0, unsigned int uLabel);

#endif /* MODAL_H */

