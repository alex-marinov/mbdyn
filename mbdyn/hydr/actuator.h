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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef ACTUATOR_H
#define ACTUATOR_H

#include "preselem.h"

/* Actuator - begin */

class Actuator : virtual public Elem, public HydraulicElem {
 private:
   const PressureNode* pNodeHyd1;
   const PressureNode* pNodeHyd2;
   HydraulicFluid* HF2;
   doublereal area1;              /* area su cui agisce p1 */
   doublereal area2;              /* area su cui agisce p2 */
   
   doublereal dl;                 /* lunghezza della camera 
				   * (a meno dello spessore del pistone) */
   
   const Vec3 axis;               /* asse dell'attuatore */
   
   doublereal dp1;                /* pressioni e derivate stati interni */
   doublereal dp2;
   doublereal dpP1;
   doublereal dpP2;
   
   doublereal flow1;
   doublereal flow2;
   doublereal Vol1;
   doublereal Vol2;
   doublereal density1;
   doublereal density2;
 
   const StructNode* pNodeStr1;   /* nodo strutturale 1 */
   const StructNode* pNodeStr2;   /* nodo strutturale 2 */
   const Vec3 f1;                 /* offset nodo1 */
   const Vec3 f2;                 /* offset nodo2 */
  
 public:
   Actuator(unsigned int uL, const DofOwner* pD,
	    const PressureNode* p1, const PressureNode* p2,  
	    const StructNode* pN1, const StructNode* pN2,
	    const Vec3& f1Tmp, const Vec3& f2Tmp,
	    const Vec3& axisTmp,
	    HydraulicFluid* hf1,
	    HydraulicFluid* hf2, 
	    doublereal A_1, doublereal A_2, 
	    doublereal l,
	    flag fOut);
   
   ~Actuator(void);
  
   /* Tipo di elemento idraulico (usato solo per debug ecc.) */
   virtual HydraulicElem::Type GetHydraulicType(void) const;

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;
   virtual DofOrder::Order GetEqType(unsigned int i) const;
   
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
   
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(4);
     connectedNodes[0] = pNodeHyd1;
     connectedNodes[1] = pNodeHyd2;
     connectedNodes[2] = pNodeStr1;
     connectedNodes[3] = pNodeStr2;
   };
   /* ************************************************ */
};

/* Actuator - end */
			      
#endif
