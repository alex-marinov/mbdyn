/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <fstream>
#include <cmath>

#include "myassert.h"
#include "solman.h"
#include "harwrap.h"

#include "dae-intg.h"

/*
 * Rationale: experiment with inverse dynamics for nonlinear control
 */

doublereal k = 1.;
doublereal m1 = 1.;
doublereal m2 = 1.;
doublereal g = 1.e-3;
doublereal gP = 1.;

static int
read(void** pp, const char* user_defined)
{
	return 0;
}

static int
size(void* p)
{
	return 7;
}

static int
init(void* p, VectorHandler& X, VectorHandler& XP)
{
	X.Reset();
	XP.Reset();

	return 0;
}

static int
grad(void* p, MatrixHandler& J, MatrixHandler& JP,
     const VectorHandler& X, const VectorHandler& XP, const doublereal& t)
{
	J(1, 2) = -1.;
	J(2, 1) = k;
	J(2, 3) = -k;
	J(3, 4) = -1.;
	J(4, 1) = -k;
	J(4, 3) = k;
	J(5, 1) = 1.;
	J(6, 3) = 1.;
	J(7, 4) = 1./m2;

	JP(1, 1) = m1;
	JP(2, 2) = 1.;
	JP(3, 3) = m2;
	JP(3, 7) = 1.;
	JP(4, 4) = 1.;
	JP(2, 5) = 1.;
	JP(5, 6) = -1.;
	JP(7, 7) = gP;

	return 0;
}

static int
func(void* p, VectorHandler& R, const VectorHandler& X, 
		const VectorHandler& XP, const doublereal& t)
{
	doublereal x1 = X(1);
	doublereal q1 = X(2);
	doublereal x2 = X(3);
	doublereal q2 = X(4);

	doublereal x1P = XP(1);
	doublereal q1P = XP(2);
	doublereal x2P = XP(3);
	doublereal q2P = XP(4);
	doublereal lambda = XP(5);
	doublereal xs = XP(6);
	doublereal mu = XP(7);

	R(1) = q1 - m1*x1P;
	R(2) = -q1P - k*(x1 - x2) - lambda;
	R(3) = q2 - m2*x2P - mu;
	R(4) = -q2P - k*(x2 - x1);
	R(5) = xs - x1;
	R(6) = -x2;
	R(7) = -q2/m2 - gP*mu;
	if (t > 1.) {
		R(6) += 1. - std::cos(2.*M_PI*(t - 1.));
		R(7) += 2.*M_PI*std::sin(2.*M_PI*(t - 1.));
	}

	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
   	return o << X(1) << " " << X(2)
		<< " " << X(3) << " " << X(4)
		<< " " << X(5) << " " << X(6)
		<< " " << X(7)
     		<< " " << XP(1) << " " << XP(2)
     		<< " " << XP(3) << " " << XP(4)
     		<< " " << XP(5) << " " << XP(6)
		<< " " << XP(7);
}

/* simboli da esportare */
static funcs _ff = {
	read,
	NULL,
	init,
	size,
	grad,
	func,
	out,
	NULL
};

void* ff = (void *)&_ff;

