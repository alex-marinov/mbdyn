/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2008
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

#ifndef MODULE_AERODYN_H
#define MODULE_AERODYN_H

//#define USE_DOUBLE_PRECISION
#define USE_SINGLE_PRECISION
#include "AeroDyn.h"

struct AeroNode {
	StructNode	*pNode;
	Vec3		f;     // offset of the aero point wrt./ the node.
	Mat3x3		Ra;    // aerodynamic orientation of the aero point wrt./ the node
	doublereal	dBuiltinPitch;
	Vec3		F;     // Force on Node
	Vec3		M;     // Moment on Node.

	doublereal	PITNOW; // Node pitch angle.
};

/*
 * user-defined struct
 */
typedef struct module_aerodyn_t {
	/**
	 * Nacelle node; requirements:
	 * - axis 3 is the shaft axis
	 * - axis 3 in wind direction
	 */
	StructNode	*pNacelle;
	StructNode	*pHub;

	integer		nblades; // the number of blades.
	integer		nelems;  // the number of elements per blade.

	doublereal	Hub_Tower_xy_distance;

	/*
	 * node data
	 */
	std::vector<AeroNode>	nodes;		// nodes
	std::vector<Mat3x3>	bladeR;		// reference orientation of each blade in the hub reference frame

	std::string ofname;
	std::ofstream out;

	/* 
	 * internal states to access the Variables which is defined in the 
	 * common MODULEs of AeroDyn
	 */
	F_LOGICAL	FirstLoop;
	F_INTEGER       elem; // use to identify the current element in the interface module!
	F_INTEGER       c_elem; // use to identify the current element in AeroDyn!
	F_INTEGER	c_blade; // use to identify the current blade!

	bool bFirst;
	DriveOwner	Time;		// time drive
	doublereal	dOldTime;	// old time
	doublereal	dDT;		// time step
} module_aerodyn_t;

#endif /* MODULE_AERODYN_H */
