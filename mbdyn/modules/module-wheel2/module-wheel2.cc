/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "userelem.h"

#include <iostream>
#include <limits>
#include <cfloat>
#include <limits>

#include "module-wheel2.h"

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

Wheel2::Wheel2(unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	wheel2							\n"
"Author: 	Stefania Gualdi <gualdi@aero.polimi.it>			\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"Connects 2 structural nodes:						\n"
"     -	Wheel								\n"
"     -	Ground								\n"
"									\n"
"Note: 									\n"
"     -	The Axle and the Wheel structural nodes must be connected 	\n"
"	by a joint that allows relative rotations only about 		\n"
"	one axis (the axle)						\n"
"     -	The center of the wheel is assumed coincident with 		\n"
"	the position of the wheel structural node			\n"
"     -	The Ground structural node supports a plane defined		\n"
"	a point and a direction orthogonal to the plane (future 	\n"
"	versions might use an arbitrary, deformable surface)		\n"
"     -	The forces are applied at the \"contact point\", that 		\n"
"	is defined according to geometrical properties 			\n"
"	of the system and according to the relative position 		\n"
"	and orientation of the Wheel and Ground structural nodes	\n"
"									\n"
"     -	Input:								\n"
"		<wheel structural node label> ,				\n"
"		<wheel axle direction> ,				\n"
"		<ground structural node label> ,			\n"
"		<reference point position of the ground plane> ,	\n"
"		<direction orthogonal to the ground plane> ,		\n"
"		<wheel radius> ,					\n"
"		<torus radius> ,					\n"
"		<volume coefficient (black magic?)> ,			\n"
"		<tire pressure> ,					\n"
"		<tire polytropic exponent> ,				\n"
"		<reference velocity for tire hysteresis>		\n"
"		[ slip ,						\n"
"		<longitudinal friction coefficient drive>		\n"
"		<lateral friction coefficient drive for s.r.=0>		\n"
"		<lateral friction coefficient drive for s.r.=1>		\n"
"		[ , threshold , <slip ratio velocity threshold> , 	\n"
"			<slip angle velocity threshold> ] ]		\n"
"									\n"
"     -	Output:								\n"
"		1)	element label					\n"
"		2-4)	tire force in global reference frame		\n"
"		5-7)	tire couple in global reference frame		\n"
"		8)	effective radius				\n"
"		9)	tire radial deformation				\n"
"		10)	tire radial deformation velocity		\n"
"		11)	slip ratio					\n"
"		12)	slip angle					\n"
"		13)	longitudinal friction coefficient		\n"
"		14)	lateral friction coefficient			\n"
"		15)	axis relative tangential velocity		\n"
"		16)	point of contact relative tangential velocity	\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// read wheel node
	pWheel = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);

	// read wheel axle
	ReferenceFrame RF = ReferenceFrame(pWheel);
	WheelAxle = HP.GetVecRel(RF);

	// read ground node
	pGround = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);

	// read ground position/orientation
	RF = ReferenceFrame(pGround);
	GroundPosition = HP.GetPosRel(RF);
	GroundDirection = HP.GetVecRel(RF);

	// normalize ground orientation
	doublereal d = GroundDirection.Dot();
	if (d <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Wheel2(" << uLabel << "): "
			"null direction at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	GroundDirection /= sqrt(d);

	// wheel geometry
	dRadius = HP.GetReal();
	dInternalRadius = HP.GetReal();
	dVolCoef = HP.GetReal();

	// reference area
	dRefArea = 3.7*dInternalRadius*std::sqrt(dInternalRadius*(2.*dRadius - dInternalRadius));

	// parameter required to compute Delta L
	dRNP = dRadius + GroundPosition*GroundDirection;

	// reference volume
	dV0 = 2.*M_PI*(dRadius - dInternalRadius)
		*M_PI*dInternalRadius*dInternalRadius*dVolCoef;

	// tyre properties
	dP0 = HP.GetReal();
	dGamma = HP.GetReal();
	dHystVRef = HP.GetReal();

	// friction
	bSlip = false;
	if (HP.IsKeyWord("slip")) {
		bSlip = true;

		/*
		 * Parametri di attrito
		 */
		pMuX0 = HP.GetDriveCaller();
		pMuY0 = HP.GetDriveCaller();
		pMuY1 = HP.GetDriveCaller();
	
		dvThreshold = 0.;
		dAlphaThreshold = 0.;
		if (HP.IsKeyWord("threshold")) {
			dvThreshold = HP.GetReal();
			if (dvThreshold < 0.) {
				silent_cerr("Wheel2(" << uLabel << "): "
					"illegal velocity threshold " << dvThreshold
					<< " at line " << HP.GetLineData() << std::endl);
				dvThreshold = std::abs(dvThreshold);
			}

			dAlphaThreshold = HP.GetReal();
			if (dvThreshold < 0.) {
				silent_cerr("Wheel2(" << uLabel << "): "
					"illegal slip angle threshold " << dAlphaThreshold
					<< " at line " << HP.GetLineData() << std::endl);
				dAlphaThreshold = std::abs(dAlphaThreshold);
			}
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
	
	std::ostream& out = pDM->GetLogFile();
	out << "wheel2: " << uLabel
		<< " " << pWheel->GetLabel()	//node label
		<< " " << WheelAxle		//wheel axle
		<< " " << pGround->GetLabel()	//ground label
		<< " " << GroundDirection	//ground direction
		<< " " << dRadius		//wheel radius
		<< " " << dInternalRadius	//wheel internal radius
		<< std::endl;
}

Wheel2::~Wheel2(void)
{
	NO_OP;
}

void
Wheel2::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel()	// 1:	label
			<< " " << F			// 2-4:	force
			<< " " << M			// 5-7:	moment
			<< " " << dInstRadius		// 8:	inst. radius
			<< " " << dDeltaL		// 9:	radial deformation
			<< " " << dVn 			// 10:	radial deformation velocity
			<< " " << dSr			// 11:	slip ratio
			<< " " << 180./M_PI*dAlpha	// 12:	slip angle
			<< " " << dMuX			// 13:	longitudinal friction coefficient
			<< " " << dMuY			// 14:	lateral friction coefficient
			<< " " << dVa			// 15:	axis relative velocity
			<< " " << dVc			// 16:	POC relative velocity
			<< std::endl;
	}
}

void
Wheel2::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler& 
Wheel2::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	WorkMat.SetNullMatrix();
	
	return WorkMat;
}

SubVectorHandler& 
Wheel2::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// ground orientation in the absolute frame
	Vec3 n = pGround->GetRCurr()*GroundDirection;

	// wheel-ground distance in the absolute frame
	Vec3 x = pWheel->GetXCurr() - pGround->GetXCurr();

	// contact when dDeltaL > 0
	dDeltaL = dRNP - x*n;

	// reset output data
	dInstRadius = dRadius - dDeltaL;
	
	dSr = 0.;
	dAlpha = 0.;

	dMuX = 0.;
	dMuY = 0.;

	dVa = 0.;
	dVc = 0.;

	if (dDeltaL < 0.) {
		
		F = Zero3;
		M = Zero3;

		dInstRadius = dRadius;
		dDeltaL = 0.;

		// no need to assemble residual contribution
		WorkVec.Resize(0);

		return WorkVec;
	}

	// resize residual
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.ResizeReset(iNumRows);

	integer iGroundFirstMomIndex = pGround->iGetFirstMomentumIndex();
	integer iWheelFirstMomIndex = pWheel->iGetFirstMomentumIndex();

	// equations indexes
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iGroundFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iWheelFirstMomIndex + iCnt);
	}

	// relative speed between wheel axle and ground
	Vec3 va = pWheel->GetVCurr() - pGround->GetVCurr()
		- pGround->GetWCurr().Cross(pGround->GetRCurr()*GroundPosition);

	dVa = (va - n*(n*va)).Norm();
	
	// relative speed between wheel and ground at contact point
	Vec3 v = va - (pWheel->GetWCurr()).Cross(n*dInstRadius);
	
	// normal component of relative speed
	// (positive when wheel departs from ground)
	dVn = n*v;
	
	dVc = (v - n*dVn).Norm();
	
	// estimated contact area
	doublereal dA = dRefArea*(dDeltaL/dInternalRadius);
	
	// estimated intersection volume
	doublereal dDeltaV = .5*dA*dDeltaL;

	// estimated tyre pressure (polytropic)
	doublereal dP = dP0*pow(dV0/(dV0 - dDeltaV), dGamma);

	// we assume the contact point lies at the root of the normal
	// to the ground that passes through the center of the wheel
	// this means that the axle needs to remain parallel to the ground
	Vec3 pc = pWheel->GetXCurr() - (n*dInstRadius);

	// force
	doublereal dFn = (dA*dP*(1. - tanh(dVn/dHystVRef)));
	F = n*dFn;

	if (bSlip) {
		// "forward" direction: axle cross normal to ground
		Vec3 fwd = (pWheel->GetRCurr()*WheelAxle).Cross(n);
		doublereal d = fwd.Dot();
		if (d < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Wheel2(" << GetLabel() << "): "
				"wheel axle is (neraly) orthogonal "
				"to the ground" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		fwd /= sqrt(d);

		// slip ratio
		doublereal dvx = fwd.Dot(v);
		doublereal sgn = copysign(1., dvx);
		doublereal dfvx = std::abs(dvx);
		doublereal dvax = fwd.Dot(va);
		doublereal dfvax = std::abs(dvax);

		// FIXME: if vax ~ 0 (e.g. because the vehicle stopped)
		// the "sleep" ratio needs to be small, right?
		dSr = dfvx/(dfvax + dvThreshold);
		if (dSr > 1.) {
			dSr = 1.;
		}

		// "lateral" direction: normal to ground cross forward
		Vec3 lat = n.Cross(fwd);

		// lateral speed of wheel center
		doublereal dvay = lat.Dot(va);

		// wheel center drift angle
		// NOTE: the angle is restricted to the first quadran
		// the angle goes to zero when the velocity is below threshold
		dAlpha = atan2(dvay, std::abs(dvax));
		if (dAlphaThreshold > 0.) {
			doublereal dtmp = tanh(sqrt(dvax*dvax + dvay*dvay)/dAlphaThreshold);
			dAlpha *= dtmp*dtmp;
		}

		// longitudinal friction coefficient
		doublereal dMuX0 = pMuX0->dGet(dSr);

		// NOTE: -1 < alpha/(pi/2) < 1
		dMuX = dMuX0*sgn*(1. - std::abs(dAlpha)/M_PI_2);

		// force correction: the longitudinal friction coefficient
		// is used with the sign of the contact point relative speed
		F -= fwd*dFn*dMuX;

		if (dvay != 0.) {
			doublereal dMuY0 = pMuY0->dGet(dAlpha);
			doublereal dMuY1 = pMuY1->dGet(dAlpha);
			
			dMuY = dMuY0 + (dMuY1 - dMuY0)*dSr;

			// force correction
			F -= lat*dFn*dMuY;
		}
	}

	// moment
	M = (pc - pWheel->GetXCurr()).Cross(F);

	WorkVec.Sub(1, F);
	WorkVec.Sub(4, (pc - pGround->GetXCurr()).Cross(F));
	WorkVec.Add(7, F);
	WorkVec.Add(10, M);

	return WorkVec;
}

unsigned int
Wheel2::iGetNumPrivData(void) const
{
	return 0;
}

int
Wheel2::iGetNumConnectedNodes(void) const
{
	// wheel + ground
	return 2;
}

void
Wheel2::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	// wheel + ground
	connectedNodes.resize(2);

	connectedNodes[0] = pWheel;
	connectedNodes[1] = pGround;
}

void 
Wheel2::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
   	NO_OP;
}

std::ostream& 
Wheel2::Restart(std::ostream& out) const
{
   	return out << "# not implemented yet" << std::endl;
}

unsigned
Wheel2::iGetInitialNumDof(void) const
{
	return 0;
}

void 
Wheel2::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
Wheel2::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler& 
Wheel2::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	WorkVec.Resize(0);
	return WorkVec;
}

bool
wheel2_set(void)
{
	UserDefinedElemRead *rf = new UDERead<Wheel2>;

	if (!SetUDE("wheel2", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!wheel2_set()) {
		silent_cerr("Wheel2: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}
	return 0;
}
#endif // ! STATIC_MODULES
