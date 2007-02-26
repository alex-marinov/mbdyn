/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#include <vector>
#include "ac/f2c.h"
#include "myassert.h"

/*
 * Compute spline coefficients
 */
void spline (	
	const std::vector<doublereal>& x, 
	const std::vector<doublereal>& y, 
	std::vector<doublereal>& b, 
	std::vector<doublereal>& c, 
	std::vector<doublereal>& d);
/*
 * Evaluate spline
*/
doublereal seval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const std::vector<doublereal>& b,
	const std::vector<doublereal>& c,
	const std::vector<doublereal>& d,
	const int diff = 0);
/*
 * Evaluate multilinear function
*/
doublereal leval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const int diff = 0);

