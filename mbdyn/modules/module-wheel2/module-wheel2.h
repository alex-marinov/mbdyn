/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#ifndef MODULE_WHEEL2_H
#define MODULE_WHEEL2_H

#include "dataman.h"
#include "userelem.h"

class Wheel2
: virtual public Elem, public UserDefinedElem
{
private:

	// wheel node
	StructNode *pWheel;

	// wheel axle direction (wrt/ wheel node)
	Vec3 WheelAxle;

	// (flat) ground node
	StructNode *pGround;

	// ground position/orientation (wrt/ ground node)
	Vec3 GroundPosition;
	Vec3 GroundDirection;

	// wheel geometry
	doublereal dRadius;
	doublereal dInternalRadius;
	doublereal dVolCoef;

	doublereal dRefArea;
	doublereal dRNP;		/* R+nG'*pG */
	doublereal dV0;

	// tyre properties
	doublereal dP0;
	doublereal dGamma;
	doublereal dHystVRef;

	// friction data
	bool bSlip;
	DriveCaller *pMuX0;
	DriveCaller *pMuY0;
	DriveCaller *pMuY1;
	doublereal dvThreshold;

	// output data
	Vec3 F;
	Vec3 M;
	doublereal dInstRadius;
	doublereal dDeltaL;
	doublereal dVn;
	doublereal dSr;
	doublereal dAlpha;
	doublereal dAlphaThreshold;
	doublereal dMuX;
	doublereal dMuY;
	doublereal dVa;
	doublereal dVc;

public:
	Wheel2(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~Wheel2(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

#endif // ! MODULE_WHEEL2_H

