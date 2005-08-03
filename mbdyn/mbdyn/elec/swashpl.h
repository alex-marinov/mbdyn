/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

/* Genel piatto oscillante */

#ifndef SWASHPL_H
#define SWASHPL_H

#include "genel.h"
#include "drive.h"
#include "elecnode.h"

/* SwashPlate - begin */

class SwashPlate : virtual public Elem, public Genel {
protected:   
	const AbstractNode* pCollectiveIn;
	// DriveOwner Collective;
	const AbstractNode* pLongitudinalIn;
	// DriveOwner Longitudinal;
	const AbstractNode* pLateralIn;
	// DriveOwner Lateral;

	const AbstractNode* pNode1;
	const AbstractNode* pNode2;
	const AbstractNode* pNode3;

	doublereal dDynamicCoef;
	doublereal dCyclicFactor;
	doublereal dCollectiveFactor;

	flag fCollLimits;
	doublereal dCollMax;
	doublereal dCollMin;

	flag fForeAftLimits;
	doublereal dForeAftMax;
	doublereal dForeAftMin;

	flag fLatLimits;
	doublereal dLatMax;
	doublereal dLatMin;

public:
	SwashPlate(unsigned int uL, const DofOwner* pDO,
		const AbstractNode* pCollIn, // const DriveCaller* pColl, 
		const AbstractNode* pLongIn, // const DriveCaller* pLong, 
		const AbstractNode* pLatIn,  // const DriveCaller* pLat,
		const AbstractNode* pN1,
		const AbstractNode* pN2,
		const AbstractNode* pN3,
		doublereal dDynCoef,
		doublereal dCyclFact,
		doublereal dCollFact,
		flag fCL,
		doublereal dCMin,
		doublereal dCMax,
		flag fFL,
		doublereal dFMin,
		doublereal dFMax,
		flag fLL,
		doublereal dLMin,	      
		doublereal dLMax,
		flag fOut);

	virtual ~SwashPlate(void);
	virtual inline void* pGet(void) const { 
		return (void*)this;
	};

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Genel */
	virtual Genel::Type GetGenelType(void) const { 
		return Genel::SWASHPLATE; 
	};

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	virtual void SetInitialValue(VectorHandler& X) const;
	virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
		NumNodes = 6;
		NdTyps[0] = pCollectiveIn->GetNodeType();
		NdLabels[0] = pCollectiveIn->GetLabel();
		NdTyps[1] = pLongitudinalIn->GetNodeType();
		NdLabels[1] = pLongitudinalIn->GetLabel();
		NdTyps[2] = pLateralIn->GetNodeType();
		NdLabels[2] = pLateralIn->GetLabel();
		NdTyps[3] = pNode1->GetNodeType();
		NdLabels[3] = pNode1->GetLabel();
		NdTyps[4] = pNode2->GetNodeType();
		NdLabels[4] = pNode2->GetLabel();
		NdTyps[5] = pNode3->GetNodeType();
		NdLabels[5] = pNode3->GetLabel();
	};
	/* ************************************************ */
};

/* SwashPlate - end */

#endif /* SWASHPL_H */

