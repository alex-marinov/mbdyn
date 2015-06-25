/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
/*
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Tingnan Zhang <tingnan1986@gatech.edu>
 */

#ifndef CUDATEST_H
#define CUDATEST_H

#include "matvec3.h"

class CUDATest {
public:
	virtual ~CUDATest(void) {};
	virtual void GetForce(Vec3& F, const Vec3& X, const Vec3& V) = 0;
};

extern "C" void *
mbdyn_CUDATest_init(unsigned n, double *pk, double *pr);

#endif // CUDATEST_H


