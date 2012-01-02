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

#ifndef GEOMDATA_H
#define GEOMDATA_H

#include "matvec3.h"

/* Geometry */
struct Geometry {
	unsigned uLabel;

	// kinematics
	Vec3 X;
	Mat3x3 R;
	Vec3 V;
	Vec3 W;

	// optional kinematics
	Vec3 XPP;
	Vec3 WP;

	// optional forces
	Vec3 F;
	Vec3 M;
};

struct GeometryData {
	enum Flags {
		X			= 0x0001U,

		R			= 0x0002U,
		RT			= 0x0004U,
		PHI			= 0x0008U,

		ORIENTATION_MASK	= (R | RT | PHI),

		V			= 0x0010U,
		W			= 0x0020U,

		XPP			= 0x0040U,
		WP			= 0x0080U,

		ACCELERATIONS_MASK	= (XPP | WP),

		F			= 0x0100U,
		M			= 0x0200U,

		FORCES_MASK		= (F | M)
	};
	unsigned uFlags;

	std::vector<Geometry> data;
};

#endif // GEOMDATA_H

