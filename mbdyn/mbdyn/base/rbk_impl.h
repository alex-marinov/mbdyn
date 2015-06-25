/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#ifndef RBK_IMPL_H
#define RBK_IMPL_H

#include "rbk.h"
#include "tpldrive.h"

class ConstRigidBodyKinematics : public RigidBodyKinematics {
protected:
	Vec3 X;
	Mat3x3 R;
	Vec3 V;
	Vec3 W;
	Vec3 XPP;
	Vec3 WP;

public:
	ConstRigidBodyKinematics(void);
	ConstRigidBodyKinematics(const Vec3& X,
		const Mat3x3& R,
		const Vec3& V,
		const Vec3& W,
		const Vec3& XPP,
		const Vec3& WP);
	virtual ~ConstRigidBodyKinematics(void);
	virtual const Vec3& GetX(void) const;
	virtual const Mat3x3& GetR(void) const;
	virtual const Vec3& GetV(void) const;
	virtual const Vec3& GetW(void) const;
	virtual const Vec3& GetXPP(void) const;
	virtual const Vec3& GetWP(void) const;
};

class DriveRigidBodyKinematics : public ConstRigidBodyKinematics {
private:
	TplDriveOwner<Vec3> XDrv;
	TplDriveOwner<Vec3> ThetaDrv;
	TplDriveOwner<Vec3> VDrv;
	TplDriveOwner<Vec3> WDrv;
	TplDriveOwner<Vec3> XPPDrv;
	TplDriveOwner<Vec3> WPDrv;

public:
	DriveRigidBodyKinematics(
		const TplDriveCaller<Vec3> *pXDrv,
		const TplDriveCaller<Vec3> *pThetaDrv,
		const TplDriveCaller<Vec3> *pVDrv,
		const TplDriveCaller<Vec3> *pWDrv,
		const TplDriveCaller<Vec3> *pXPPDrv,
		const TplDriveCaller<Vec3> *pWPDrv);
	virtual ~DriveRigidBodyKinematics(void);
	virtual void Update(void);
};

#endif // RBK_IMPL_H

