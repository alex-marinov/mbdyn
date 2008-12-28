/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "rbk.h"

RigidBodyKinematics::~RigidBodyKinematics(void)
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

