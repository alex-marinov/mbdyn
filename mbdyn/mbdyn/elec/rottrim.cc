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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "genel.h"
#include "drive.h"
#include "dataman.h"
#include "constltp.h"
#include "strnode.h"
#include "rottrim.h"

RotorTrimBase::RotorTrimBase(unsigned int uL,
	const DofOwner* pDO,
	const ScalarDifferentialNode* pNode1,
	const ScalarDifferentialNode* pNode2,
	const ScalarDifferentialNode* pNode3,
	const DriveCaller* pDThrust,
	const DriveCaller* pDRollMoment,
	const DriveCaller* pDPitchMoment,
	const doublereal& dG,
	const doublereal& dp,
	const doublereal& dT0,
	const doublereal& dT1,
	const doublereal& dK0,
	const doublereal& dK1,
	const DriveCaller *pTrigger,
	flag fOut)
: Elem(uL, fOut),
Genel(uL, pDO, fOut),
dRadius(-1.),
DThrust(pDThrust),
DRollMoment(pDRollMoment),
DPitchMoment(pDPitchMoment),
Trigger(pTrigger),
dGamma(dG),
dP(dp),
dP2(dP*dP),
dC(8.*(dP*dP - 1.)/dGamma),
dC2(dC*dC),
dTau0(dT0),
dTau1(dT1),
dKappa0(dK0),
dKappa1(dK1)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode3 != NULL);
	ASSERT(dGamma > 0.);
	ASSERT(dP > 0.);

	pvNodes[0] = pNode1;
	pvNodes[1] = pNode2;
	pvNodes[2] = pNode3;
}

RotorTrimBase::~RotorTrimBase(void)
{
	NO_OP;
}

unsigned int
RotorTrimBase::iGetNumDof(void) const
{
 	return 0;
}

/* Dimensioni del workspace */
void
RotorTrimBase::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 3;
	*piNumCols = 3;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
RotorTrimBase::AssJac(VariableSubMatrixHandler& WorkMat,
	  doublereal dCoef,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering RotorTrimBase::AssJac()" << std::endl);

	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.Resize(3, 0);

	doublereal	d[3];
	if (Trigger.dGet() == 0.) {
		d[0] = dCoef;
		d[1] = dCoef;
		d[2] = dCoef;

	} else {
		d[0] = dTau0 + dCoef;
		d[1] = dTau1 + dCoef;
		d[2] = dTau1 + dCoef;
	}

	for (int iIdx = 0; iIdx < 3; iIdx++ ) {
		integer iRowIndex = pvNodes[iIdx]->iGetFirstRowIndex() + 1;
		integer iColIndex = pvNodes[iIdx]->iGetFirstColIndex() + 1;
       		WM.PutItem(iIdx + 1, iRowIndex, iColIndex, d[iIdx]);
	}

        return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
RotorTrimBase::AssRes(SubVectorHandler& WorkVec,
		  doublereal dCoef,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering RotorTrimBase::AssRes()" << std::endl);

	WorkVec.Resize(3);

	WorkVec.PutRowIndex(1, pvNodes[0]->iGetFirstRowIndex() + 1);
	WorkVec.PutRowIndex(2, pvNodes[1]->iGetFirstRowIndex() + 1);
	WorkVec.PutRowIndex(3, pvNodes[2]->iGetFirstRowIndex() + 1);

	doublereal dX1 = pvNodes[0]->dGetX();
	doublereal dX2 = pvNodes[1]->dGetX();
	doublereal dX3 = pvNodes[2]->dGetX();

	if (Trigger.dGet() == 0.) {
		WorkVec.PutCoef(1, -dX1);
		WorkVec.PutCoef(2, -dX2);
		WorkVec.PutCoef(3, -dX3);

		return WorkVec;
	}

	doublereal dX1Prime = pvNodes[0]->dGetXPrime();
	doublereal dX2Prime = pvNodes[1]->dGetXPrime();
	doublereal dX3Prime = pvNodes[2]->dGetXPrime();

	doublereal dThrust;
	doublereal dRollMoment;
	doublereal dPitchMoment;
	doublereal dRho;
	doublereal dOmega;
	doublereal dMu;

	GetData(dThrust, dRollMoment, dPitchMoment, dRho, dOmega, dMu);

	doublereal d = M_PI*pow(dRadius, 4)*dRho*dOmega*dOmega;
	dThrust /= d;
	d *= dRadius;
	dRollMoment /= d;
	dPitchMoment /= d;

	doublereal dMu2 = dMu*dMu;

	const doublereal f = dC/(1. + dC2);

	Mat3x3 m((1. + 3./2.*dMu2)/6.,
		-f*dMu/6.*(dC - dGamma/(16.*dP2)),
		f*dMu/6.*(1. + dC*dGamma/(16.*dP2)),
		2./9.*dMu,
		-f/16.*(dC*(1. + 3./2.*dMu2) - 2./9.*dMu2*dGamma/dP2),
		f/16.*((1. + 2.*dMu2) + 2./9.*dMu2*dC*dGamma/dP2),
		0.,
		-f/16.,
		-f*dC/16.*(1. + 1./2.*dMu2));
	Vec3 v(DThrust.dGet() - dThrust,
		DRollMoment.dGet() - dRollMoment,
		DPitchMoment.dGet() - dPitchMoment);
	v = m.Solve(v);

	WorkVec.PutCoef(1, v(1)*dKappa0 - dX1 - dTau0*dX1Prime);
	WorkVec.PutCoef(2, v(2)*dKappa1 - dX2 - dTau1*dX2Prime);
	WorkVec.PutCoef(3, v(3)*dKappa1 - dX3 - dTau1*dX3Prime);

	return WorkVec;
}

void
RotorTrimBase::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(3);
	for (int i = 0; i < 3; i++) {
		connectedNodes[i] = pvNodes[i];
	}
}

static unsigned iRotorTz = 0;
static unsigned iRotorMx = 0;
static unsigned iRotorMy = 0;

RotorTrim::RotorTrim(unsigned int uL,
	const DofOwner* pDO,
	const Rotor* pRot,
	const ScalarDifferentialNode* pNode1,
	const ScalarDifferentialNode* pNode2,
	const ScalarDifferentialNode* pNode3,
	const DriveCaller* pDThrust,
	const DriveCaller* pDRollMoment,
	const DriveCaller* pDPitchMoment,
	const doublereal& dG,
	const doublereal& dp,
	const doublereal& dT0,
	const doublereal& dT1,
	const doublereal& dK0,
	const doublereal& dK1,
	const DriveCaller *pTrigger,
	flag fOut)
: Elem(uL, fOut),
RotorTrimBase(uL, pDO, pNode1, pNode2, pNode3,
	pDThrust, pDRollMoment, pDPitchMoment,
	dG, dp, dT0, dT1, dK0, dK1, pTrigger, fOut),
pRotor(pRot)
{
	ASSERT(pRotor != NULL);

	dRadius = pRotor->dGetRadius();

	if (iRotorTz == 0) {
		/* first RotorTrim */
		iRotorTz = pRotor->iGetPrivDataIdx("Tz");
		iRotorMx = pRotor->iGetPrivDataIdx("Mx");
		iRotorMy = pRotor->iGetPrivDataIdx("My");
	}
}

RotorTrim::~RotorTrim(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
RotorTrim::Restart(std::ostream& out) const
{
	return out << " /* rotor trim not implemented yet */ ";
}

void
RotorTrim::GetData(doublereal &dThrust,
	doublereal &dRollMoment,
	doublereal &dPitchMoment,
	doublereal &dRho,
	doublereal &dOmega,
	doublereal &dMu) const
{
	dThrust = pRotor->dGetPrivData(iRotorTz);
	dRollMoment = pRotor->dGetPrivData(iRotorMx);
	dPitchMoment = pRotor->dGetPrivData(iRotorMy);
	dRho = pRotor->dGetAirDensity(pRotor->GetXCurr());
	dOmega = pRotor->dGetOmega();
	dMu = pRotor->dGetMu();
}

void
RotorTrim::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	pRotor->GetConnectedNodes(connectedNodes);
	int NumNodes = connectedNodes.size();
	connectedNodes.resize(NumNodes + 3);
	for (int i = 0; i < 3; i++) {
		connectedNodes[NumNodes + i] = pvNodes[i];
	}
}

RotorTrimGeneric::RotorTrimGeneric(unsigned int uL,
	const DofOwner *pDO,
	const StructNode *pStrNode,
	const DriveCaller *pThrust,
	const DriveCaller *pRollMoment,
	const DriveCaller *pPitchMoment,
	const AirProperties *pAP,
	doublereal dRadius,
	const DriveCaller *pOmega,
	const DriveCaller *pMu,
	const ScalarDifferentialNode* pNode1,
	const ScalarDifferentialNode* pNode2,
	const ScalarDifferentialNode* pNode3,
	const DriveCaller* pDThrust,
	const DriveCaller* pDRollMoment,
	const DriveCaller* pDPitchMoment,
	const doublereal& dG,
	const doublereal& dp,
	const doublereal& dT0,
	const doublereal& dT1,
	const doublereal& dK0,
	const doublereal& dK1,
	const DriveCaller *pTrigger,
	flag fOut)
: Elem(uL, fOut),
RotorTrimBase(uL, pDO, pNode1, pNode2, pNode3,
	pDThrust, pDRollMoment, pDPitchMoment,
	dG, dp, dT0, dT1, dK0, dK1, pTrigger, fOut),
pStrNode(pStrNode),
Thrust(pThrust),
RollMoment(pRollMoment),
PitchMoment(pPitchMoment),
pAP(pAP),
Omega(pOmega),
Mu(pMu)
{
	this->dRadius = dRadius;
}

RotorTrimGeneric::~RotorTrimGeneric(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
RotorTrimGeneric::Restart(std::ostream& out) const
{
	return out << " /* rotor trim not implemented yet */ ";
}

void
RotorTrimGeneric::GetData(doublereal &dThrust,
	doublereal &dRollMoment,
	doublereal &dPitchMoment,
	doublereal &dRho,
	doublereal &dOmega,
	doublereal &dMu) const
{
	dThrust = Thrust.dGet();
	dRollMoment = RollMoment.dGet();
	dPitchMoment = PitchMoment.dGet();
	dRho = pAP->dGetAirDensity(pStrNode ? pStrNode->GetXCurr() : Zero3);
	dOmega = Omega.dGet();
	dMu = Mu.dGet();
}

