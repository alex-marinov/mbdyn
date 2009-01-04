/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2009
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

// uncomment to enable debug output
#define MODULE_AERODYN_DEBUG

// keep this consistent with AeroDyn build
// #define USE_DOUBLE_PRECISION
#define USE_SINGLE_PRECISION

#include "AeroDyn.h"

struct AeroNode {
	StructNode	*pNode;
	Vec3		f;	// offset of the aero point wrt./ the node,
				// constant in the node's reference frame
	Mat3x3		Ra;	// aerodynamic orientation of the aero point
				// wrt./ the node

	doublereal	dBuiltInTwist;

	Vec3		F;	// Force acting on Node
	Vec3		M;	// Moment acting on Node, with respect
				// to node's position

 	doublereal      FN;    // Normal Force on each blade element (Not being used now!).
 	doublereal      FT;    // Tangental force on each blade element.(Not being used now!)
 	doublereal      AM;    // Aerodynamic moment on each blade element.(Not being used now!)
 
	doublereal	PITNOW; // Node pitch angle
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
	std::vector<Mat3x3>	bladeR;		// orientation matrix of each blade root in the hub reference frame

	std::string ofname;
	std::ofstream out;

	/*
	 * Total aerodynamic data
	 */
	Vec3		TF;     // Total Force on the rotor in the absolute frame
	Vec3		TM;     // Total Moment on the rotor in the absolute frame.
	Vec3		TF_h;     // Total Force on the rotor in the hub frame
	Vec3		TM_h;     // Total Moment on the rotor in the hub frame.
	doublereal      Thrust;   // Rotor Thrust.
	doublereal      Torque;   // Rotor Torque.
	doublereal      Rotor_speed;   // Rotor angular velocity.
	/* 
	 * internal states to access the Variables which is defined in the 
	 * common MODULEs of AeroDyn
	 */
	F_LOGICAL	FirstLoop;
	F_INTEGER       elem;    // use to identify the current element in the interface module!
	F_INTEGER       c_elem;  // use to identify the current element in AeroDyn!
	F_INTEGER	c_blade; // use to identify the current blade in AeroDyn!
        F_REAL          rlocal;  // use to identify the current element position 
        F_REAL          r_hub;   // use to identify the hub radius.

	bool bFirst;
	DriveOwner	Time;		// time drive
	doublereal	dOldTime;	// old time
	doublereal      dCurTime;       // current time
	F_REAL          dDT;		// time step
} module_aerodyn_t;

#endif /* MODULE_AERODYN_H */
