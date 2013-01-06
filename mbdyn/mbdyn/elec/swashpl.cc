/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Swash plate */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>

#include "swashpl.h"


/* SwashPlate - begin */

SwashPlate::SwashPlate(unsigned int uL, const DofOwner* pDO,
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
		flag fOut)
: Elem(uL, fOut),
Genel(uL, pDO, fOut),
pCollectiveIn(pCollIn),   // Collective(pColl),
pLongitudinalIn(pLongIn), // Longitudinal(pLong),
pLateralIn(pLatIn),       // Lateral(pLat),
pNode1(pN1), pNode2(pN2), pNode3(pN3),
dDynamicCoef(dDynCoef),
dCyclicFactor(dCyclFact),
dCollectiveFactor(dCollFact),
fCollLimits(fCL), dCollMax(dCMax), dCollMin(dCMin),
fForeAftLimits(fFL), dForeAftMax(dFMax), dForeAftMin(dFMin),
fLatLimits(fLL), dLatMax(dLMax), dLatMin(dLMin)
{
	ASSERT(pCollectiveIn != NULL);
	ASSERT(pCollectiveIn->GetNodeType() == Node::ABSTRACT);
	ASSERT(pLongitudinalIn != NULL);
	ASSERT(pLongitudinalIn->GetNodeType() == Node::ABSTRACT);
	ASSERT(pLateralIn != NULL);
	ASSERT(pLateralIn->GetNodeType() == Node::ABSTRACT);

	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::ABSTRACT);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::ABSTRACT);
	ASSERT(pNode3 != NULL);
	ASSERT(pNode3->GetNodeType() == Node::ABSTRACT);

	ASSERT(dCyclicFactor != 0.);
	ASSERT(dCollectiveFactor != 0.);
	ASSERT(dDynamicCoef >= 0.);
}


SwashPlate::~SwashPlate(void)
{
	NO_OP;
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
SwashPlate::Restart(std::ostream& out) const
{
	Genel::Restart(out) << ", swash plate, "
		<< pCollectiveIn->GetLabel() << ", ";
	// Collective.pGetDriveCaller()->Restart(out) << ", ";
	if (fCollLimits) {
		out << "limits, " << dCollMin << ", " << dCollMax << ", ";
	}
	out << pLongitudinalIn->GetLabel() << ", ";
	// Longitudinal.pGetDriveCaller()->Restart(out) << ", ";
	if (fForeAftLimits) {
		out << "limits, " << dForeAftMin << ", " << dForeAftMax << ", ";
	}
	out << pLateralIn->GetLabel() << ", ";
	// Lateral.pGetDriveCaller()->Restart(out) << ", ";
	if (fLatLimits) {
		out << "limits, " << dLatMin << ", " << dLatMax << ", ";
	}
	out << pNode1->GetLabel() << ", "
		<< pNode2->GetLabel() << ", "
		<< pNode3->GetLabel() << ", "
		<< dDynamicCoef << ", "
		<< dCyclicFactor << ", " << dCollectiveFactor << ';' << std::endl;

	return out;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
SwashPlate::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering SwashPlate::AssJac()" << std::endl);

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(14, 0);

	integer iCollFirstIndex = pCollectiveIn->iGetFirstIndex() + 1;
	integer iLongFirstIndex = pLongitudinalIn->iGetFirstIndex() + 1;
	integer iLatFirstIndex = pLateralIn->iGetFirstIndex() + 1;

	integer iNode1FirstIndex = pNode1->iGetFirstIndex() + 1;
	integer iNode2FirstIndex = pNode2->iGetFirstIndex() + 1;
	integer iNode3FirstIndex = pNode3->iGetFirstIndex() + 1;

	WM.PutItem(1, iCollFirstIndex, iCollFirstIndex, dCoef);
	WM.PutItem(2, iLongFirstIndex, iLongFirstIndex, dCoef);
	WM.PutItem(3, iLatFirstIndex, iLatFirstIndex, dCoef);

	doublereal d = dDynamicCoef + dCoef;

	WM.PutItem(4, iNode1FirstIndex, iNode1FirstIndex, d);
	WM.PutItem(5, iNode2FirstIndex, iNode2FirstIndex, d);
	WM.PutItem(6, iNode3FirstIndex, iNode3FirstIndex, d);

	d = dCollectiveFactor;

	WM.PutItem(7, iNode1FirstIndex, iCollFirstIndex, -d);
	WM.PutItem(8, iNode2FirstIndex, iCollFirstIndex, -d);
	WM.PutItem(9, iNode3FirstIndex, iCollFirstIndex, -d);

	d = dCollectiveFactor*dCyclicFactor;

	WM.PutItem(10, iNode1FirstIndex, iLongFirstIndex, d);

	d /= 2.;

	WM.PutItem(11, iNode2FirstIndex, iLongFirstIndex, -d);
	WM.PutItem(12, iNode3FirstIndex, iLongFirstIndex, -d);

	d *= sqrt(3.);

	WM.PutItem(13, iNode2FirstIndex, iLongFirstIndex, d);
	WM.PutItem(14, iNode3FirstIndex, iLongFirstIndex, -d);

	return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler&
SwashPlate::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering SwashPlate::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	WorkVec.ResizeReset(6);

	integer iCollFirstIndex = pCollectiveIn->iGetFirstIndex() + 1;
	integer iLongFirstIndex = pLongitudinalIn->iGetFirstIndex() + 1;
	integer iLatFirstIndex = pLateralIn->iGetFirstIndex() + 1;

	integer iNode1FirstIndex = pNode1->iGetFirstIndex() + 1;
	integer iNode2FirstIndex = pNode2->iGetFirstIndex() + 1;
	integer iNode3FirstIndex = pNode3->iGetFirstIndex() + 1;

	WorkVec.PutRowIndex(1, iCollFirstIndex);
	WorkVec.PutRowIndex(2, iLongFirstIndex);
	WorkVec.PutRowIndex(3, iLatFirstIndex);

	WorkVec.PutRowIndex(4, iNode1FirstIndex);
	WorkVec.PutRowIndex(5, iNode2FirstIndex);
	WorkVec.PutRowIndex(6, iNode3FirstIndex);

	doublereal dXColl = XCurr(iCollFirstIndex);
	doublereal dXLong = XCurr(iLongFirstIndex);
	doublereal dXLat = XCurr(iLatFirstIndex);

	WorkVec.PutCoef(1, -dXColl);
	WorkVec.PutCoef(2, -dXLong);
	WorkVec.PutCoef(3, -dXLat);

	/* Limits on pitch angles */
	if (fCollLimits) {
		if (dXColl > dCollMax) {
			dXColl = dCollMax;
		} else if (dXColl < dCollMin) {
			dXColl = dCollMin;
		}
	}

	if (fForeAftLimits) {
		if (dXLong > dForeAftMax) {
			dXLong = dForeAftMax;
		} else if (dXLong < dForeAftMin) {
			dXLong = dForeAftMin;
		}
	}

	if (fLatLimits) {
		if (dXLat > dLatMax) {
			dXLat = dLatMax;
		} else if (dXLat < dLatMin) {
			dXLat = dLatMin;
		}
	}

	WorkVec.PutCoef(4, dCollectiveFactor*(dXColl - dCyclicFactor*dXLong)
			- XCurr(iNode1FirstIndex) - dDynamicCoef*XPrimeCurr(iNode1FirstIndex));
	WorkVec.PutCoef(5, dCollectiveFactor*(dXColl + dCyclicFactor*(.5*dXLong - sqrt(3.)/2.*dXLat))
			- XCurr(iNode2FirstIndex) - dDynamicCoef*XPrimeCurr(iNode2FirstIndex));
	WorkVec.PutCoef(6, dCollectiveFactor*(dXColl + dCyclicFactor*(.5*dXLong + sqrt(3.)/2.*dXLat))
			- XCurr(iNode3FirstIndex) - dDynamicCoef*XPrimeCurr(iNode3FirstIndex));

	return WorkVec;
}

void
SwashPlate::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}

/* SwashPlate - end */

