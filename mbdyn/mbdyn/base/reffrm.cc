/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Reference frame: structure, handling etc. */

#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <iomanip>

#include "reffrm.h"
#include "Rot.hh"

ReferenceFrame::ReferenceFrame(void)
: WithLabel(0),
x(::mb_zero<Vec3>()),
R(::mb_deye<Mat3x3>(1.)),
v(::mb_zero<Vec3>()),
w(::mb_zero<Vec3>()),
od(UNKNOWN_ORIENTATION_DESCRIPTION)
{
	NO_OP;
}

ReferenceFrame::ReferenceFrame(unsigned int uLabel, 
		const Vec3& xIn, const Mat3x3& RIn,
		const Vec3& vIn, const Vec3& wIn,
		const OrientationDescription& ood)
: WithLabel(uLabel), x(xIn), R(RIn), v(vIn), w(wIn), od(ood)
{
	ASSERT(Eye3.IsSame(R.MulTM(R), 1e-12));
}

ReferenceFrame::ReferenceFrame(const RigidBodyKinematics* pRBK)
: WithLabel(0),
x(pRBK->GetX()), R(pRBK->GetR()),
v(pRBK->GetV()), w(pRBK->GetW()),
od(EULER_123)
{
	ASSERT(Eye3.IsSame(R.MulTM(R), 1e-12));
}

ReferenceFrame::~ReferenceFrame(void) { 
	NO_OP;
}
	
const Vec3&
ReferenceFrame::GetX(void) const
{
	return x; 
}
	
const Mat3x3&
ReferenceFrame::GetR(void) const
{
	return R;
}
	
const Vec3&
ReferenceFrame::GetV(void) const
{
	return v;
}
	
const Vec3&
ReferenceFrame::GetW(void) const
{
	return w;
}

const Vec3&
ReferenceFrame::GetXPP(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

const Vec3&
ReferenceFrame::GetWP(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

ReferenceFrame&
ReferenceFrame::operator = (const ReferenceFrame& rf)
{
	PutLabel(rf.GetLabel());
	x = rf.x;
	R = rf.R;
	v = rf.v;
	w = rf.w;
	od = rf.od;
	return *this;
}

std::ostream&
ReferenceFrame::Output(std::ostream& out) const
{
	out 
		<< std::setw(8) << GetLabel() << " "
		<< x << " ";

	switch (od) {
	case EULER_123:
		out << MatR2EulerAngles123(R)*dRaDegr;
		break;

	case EULER_313:
		out << MatR2EulerAngles313(R)*dRaDegr;
		break;

	case EULER_321:
		out << MatR2EulerAngles321(R)*dRaDegr;
		break;

	case ORIENTATION_VECTOR:
		out << RotManip::VecRot(R);
		break;

	case ORIENTATION_MATRIX:
		out << R;
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return out << " "
		<< v << " " << w << " " << std::endl;
}

