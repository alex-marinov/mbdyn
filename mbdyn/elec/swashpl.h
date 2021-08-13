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

/* Genel piatto oscillante */

#ifndef SWASHPL_H
#define SWASHPL_H

#include "genel.h"
#include "drive.h"
#include "elecnode.h"

/* SwashPlate - begin */

class SwashPlate : virtual public Elem, public Genel {
protected:
	const ScalarDifferentialNode* pCollectiveIn;
	// DriveOwner Collective;
	const ScalarDifferentialNode* pLongitudinalIn;
	// DriveOwner Longitudinal;
	const ScalarDifferentialNode* pLateralIn;
	// DriveOwner Lateral;

	const ScalarDifferentialNode* pNode1;
	const ScalarDifferentialNode* pNode2;
	const ScalarDifferentialNode* pNode3;

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
		const ScalarDifferentialNode* pCollIn, // const DriveCaller* pColl,
		const ScalarDifferentialNode* pLongIn, // const DriveCaller* pLong,
		const ScalarDifferentialNode* pLatIn,  // const DriveCaller* pLat,
		const ScalarDifferentialNode* pN1,
		const ScalarDifferentialNode* pN2,
		const ScalarDifferentialNode* pN3,
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

	virtual void SetInitialValue(VectorHandler& X);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(6);
		connectedNodes[0] = pCollectiveIn;
		connectedNodes[1] = pLongitudinalIn;
		connectedNodes[2] = pLateralIn;
		connectedNodes[3] = pNode1;
		connectedNodes[4] = pNode2;
		connectedNodes[5] = pNode3;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;

	/* describes the dimension of components of equation */
    virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* SwashPlate - end */

#endif /* SWASHPL_H */

