/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307	 USA
 */

 /*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2022) all rights reserved.

	The copyright of this code is transferred
	to Pierangelo Masarati and Paolo Mantegazza
	for use in the software MBDyn as described
	in the GNU Public License version 2.1
  */

#ifndef __SOLVER_BASE_H__INCLUDED__
#define __SOLVER_BASE_H__INCLUDED__

struct SolverBase {
	enum StepIntegratorType {
		INT_CRANKNICOLSON,
		INT_MODCRANKNICOLSON,
		INT_MS2,
		INT_MS3,
		INT_MS4,
		INT_SS2,
		INT_SS3,
		INT_SS4,
		INT_HOPE,
		INT_BATHE,
		INT_MSSTC3,
		INT_MSSTH3,
		INT_MSSTC4,
		INT_MSSTH4,
		INT_MSSTC5,
		INT_MSSTH5,
		INT_DIRK33,
		INT_DIRK43,
		INT_DIRK54,
		INT_THIRDORDER,
		INT_IMPLICITEULER,
		INT_UNKNOWN,
		INT_DEFAULT = INT_UNKNOWN,
		INT_HYBRID
        };

        static constexpr auto INT_COUNT = INT_UNKNOWN + 1;
};

#endif
