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

#ifndef REFFRM_H
#define REFFRM_H

#include <iostream>

#include "myassert.h"
#include "mynewmem.h"

#include "withlab.h"
#include "matvec3.h"
#include "rbk.h"

class ReferenceFrame : public WithLabel, public RigidBodyKinematics {
private:
	Vec3 x;
	Mat3x3 R;
	Vec3 v;
	Vec3 w;

	OrientationDescription od;

public:
	ReferenceFrame(void);

	ReferenceFrame(unsigned int uLabel, 
			const Vec3& xIn, const Mat3x3& RIn,
			const Vec3& vIn, const Vec3& wIn,
			const OrientationDescription& ood);

	ReferenceFrame(const RigidBodyKinematics* pRBK);

	~ReferenceFrame(void);
	
	const Vec3& GetX(void) const;

	const Mat3x3& GetR(void) const;

	const Vec3& GetV(void) const;

	const Vec3& GetW(void) const;

	const Vec3& GetXPP(void) const;

	const Vec3& GetWP(void) const;

	ReferenceFrame& operator = (const ReferenceFrame& rf);

	std::ostream& Output(std::ostream& out) const;
};

#endif /* REFFRM_H */

