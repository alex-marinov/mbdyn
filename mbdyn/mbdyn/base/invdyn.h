/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef INVDYN_H
#define INVDYN_H

class InverseDynamics {
public:
	enum Order {
		UNDEFINED = -2,
		INVERSE_DYNAMICS = -1,
		POSITION = 0,
		VELOCITY = 1,
		ACCELERATION = 2
	};

	enum Type {
		FULLY_ACTUATED = 0x01U,
		FULLY_DETERMINED = 0x02U,
		COLLOCATED = 0x04U,

		// fully actuated collocated
		FULLY_ACTUATED_COLLOCATED = (FULLY_ACTUATED|FULLY_DETERMINED|COLLOCATED),
		// fully actuated non collocated
		FULLY_ACTUATED_NON_COLLOCATED = (FULLY_ACTUATED|FULLY_DETERMINED),
		// underdetermined underactuated collocated
		UNDERDETERMINED_UNDERACTUATED_COLLOCATED = (COLLOCATED),
		// underdetermined fully actuated
		UNDERDETERMINED_FULLY_ACTUATED = (FULLY_ACTUATED),


		PRESCRIBED_MOTION = 0x10U,
		TORQUE = 0x20U,
		ERGONOMY = 0x40U,
		RIGHT_HAND_SIDE = 0x80U
	};
};

extern const char *invdyn2str(InverseDynamics::Order iOrder);

#endif // INVDYN_H
