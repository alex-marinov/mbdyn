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

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */


#ifndef BEAMSLIDER_H
#define BEAMSLIDER_H

#include <myassert.h>
#include <except.h>

#include <strnode.h>
#include <elem.h>
#include <beam.h>
#include <joint.h>

/* BeamSliderJoint - begin */

class BeamSliderJoint: virtual public Elem, public Joint {
 private:
   unsigned int nRotConstr;

   const StructNode* pNode;
   const Beam** ppBeam;
   Vec3 f;
   Mat3x3 R;
   Vec3 F;
   Vec3 M;
   
 public:
   /* Costruttore non banale */
   BeamSliderJoint(unsigned int uL, const DofOwner* pDO,
		   const StructNode* pN, const Beam* pB,
		   const Vec3& fTmp, const Mat3x3& RTmp, flag fOut);
   
   /* Distruttore */
   ~BeamSliderJoint(void);

   virtual inline void* pGet(void) const { 
      return (void*)this;
   };

   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::BEAMSLIDER;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 3+nRotConstr;
   };
   
   DofOrder::Order SetDof(unsigned int i) const {
      ASSERT(i >= 0 && i < 3);
      return DofOrder::ALGEBRAIC;
   }
   
   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 24+iGetNumDof();
      *piNumCols = 24+iGetNumDof(); 
   };
   
      
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
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 6+2*nRotConstr;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const {
      *piNumRows = 48+iGetInitialNumDof(); 
      *piNumCols = 48+iGetInitialNumDof();
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const {
      return 6;
   };   
   
   virtual doublereal dGetPrivData(unsigned int i = 0) const;
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 1+1*Beam::NUMNODES;
     NdTyps[0] = pNode->GetNodeType();
     NdLabels[0] = pNode->GetLabel();

     /* for each beam */
     for (int i = 0; i < 1; i++) {
	/* for each node */
	for (int j = 1; j <= Beam::NUMNODES; j++) {
            const StructNode *pN = ppBeam[i]->pGetNode(j);
            NdTyps[Beam::NUMNODES*i+j] = pN->GetNodeType();
            NdLabels[Beam::NUMNODES*i+j] = pN->GetLabel();
	}
     }
   };
   /* ************************************************ */
};

#endif /* BEAMSLIDER_H */

