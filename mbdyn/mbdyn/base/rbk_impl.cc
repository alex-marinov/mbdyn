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

/* Rigid body kinematics: structure, handling etc. */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "rbk_impl.h"
#include "Rot.hh"

ConstRigidBodyKinematics::ConstRigidBodyKinematics(void)
: X(Zero3), R(Eye3), V(Zero3), W(Zero3), XPP(Zero3), WP(Zero3)
{
	NO_OP;
}

ConstRigidBodyKinematics::ConstRigidBodyKinematics(const Vec3& X,
	const Mat3x3& R,
	const Vec3& V,
	const Vec3& W,
	const Vec3& XPP,
	const Vec3& WP)
: X(X), R(R), V(V), W(W), XPP(XPP), WP(WP)
{
	NO_OP;
}

ConstRigidBodyKinematics::~ConstRigidBodyKinematics(void)
{
	NO_OP;
}

const Vec3&
ConstRigidBodyKinematics::GetX(void) const
{
	return X;
}

const Mat3x3&
ConstRigidBodyKinematics::GetR(void) const
{
	return R;
}

const Vec3&
ConstRigidBodyKinematics::GetV(void) const
{
	return V;
}

const Vec3&
ConstRigidBodyKinematics::GetW(void) const
{
	return W;
}

const Vec3&
ConstRigidBodyKinematics::GetXPP(void) const
{
	return XPP;
}

const Vec3&
ConstRigidBodyKinematics::GetWP(void) const
{
	return WP;
}

DriveRigidBodyKinematics::DriveRigidBodyKinematics(
	const TplDriveCaller<Vec3> *pXDrv,
	const TplDriveCaller<Vec3> *pThetaDrv,
	const TplDriveCaller<Vec3> *pVDrv,
	const TplDriveCaller<Vec3> *pWDrv,
	const TplDriveCaller<Vec3> *pXPPDrv,
	const TplDriveCaller<Vec3> *pWPDrv)
: XDrv(pXDrv),
ThetaDrv(pThetaDrv),
VDrv(pVDrv),
WDrv(pWDrv),
XPPDrv(pXPPDrv),
WPDrv(pWPDrv)
{
	Update();
}

DriveRigidBodyKinematics::~DriveRigidBodyKinematics(void)
{
	NO_OP;
}

void
DriveRigidBodyKinematics::Update(void)
{
	bool bR(false);

	if (ThetaDrv.pGetDriveCaller()) {
		R = RotManip::Rot(ThetaDrv.Get());
		bR = true;
	}

	if (XDrv.pGetDriveCaller()) {
		if (bR) {
			X = R.MulTV(XDrv.Get());

		} else {
			X = XDrv.Get();
		}
	}

	if (VDrv.pGetDriveCaller()) {
		if (bR) {
			V = R.MulTV(VDrv.Get());

		} else {
			V = VDrv.Get();
		}
	}

	if (WDrv.pGetDriveCaller()) {
		if (bR) {
			W = R.MulTV(WDrv.Get());

		} else {
			W = WDrv.Get();
		}
	}

	if (XPPDrv.pGetDriveCaller()) {
		if (bR) {
			XPP = R.MulTV(XPPDrv.Get());

		} else {
			XPP = XPPDrv.Get();
		}
	}

	if (WPDrv.pGetDriveCaller()) {
		if (bR) {
			WP = R.MulTV(WPDrv.Get());

		} else {
			WP = WPDrv.Get();
		}
	}
}
