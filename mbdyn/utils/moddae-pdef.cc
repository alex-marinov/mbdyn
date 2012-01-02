/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include <cmath>

#include <fstream>
#include <iostream>

#include "myassert.h"
#include "solman.h"

#include "dae-intg.h"

struct private_data {
	doublereal m;
	doublereal c;
	doublereal k;
	doublereal l;
	doublereal g;
	doublereal x[4];
	doublereal xP[4];
};

static int
read(void** pp, const char* user_defined)
{
	private_data* pd = new private_data;

	if (user_defined != NULL) {
		std::ifstream in(user_defined);
		if (!in) {
			 std::cerr << "unable to open file "
				 "\"" << user_defined << "\"" << std::endl;
			 exit(EXIT_FAILURE);
		}
		in >> pd->m >> pd->c >> pd->k >> pd->l >> pd->g
			>> pd->x[0] >> pd->x[1] >> pd->x[2] >> pd->x[3];

	} else {
		pd->m = 1.;
		pd->c = 0.;
		pd->k = 1.;
		pd->l = 1.;
		pd->g = 9.81;
		pd->x[0] = 0.;
		pd->x[1] = 0.;
		pd->x[2] = 1.;
		pd->x[3] = 0.;
	}

	pd->xP[0] = pd->x[2];
	pd->xP[1] = pd->x[3];
	pd->xP[2] = -(2.*pd->x[2]*pd->x[3] + pd->g*sin(pd->x[0]))/pd->l;
	pd->xP[3] = pd->l*pd->x[2]*pd->x[2] + pd->g*cos(pd->x[0])
		- (pd->k*pd->x[1] + pd->c*pd->x[3])/pd->m;

	std::cerr << "m=" << pd->m << ", c=" << pd->c << ", k=" << pd->k
			<< std::endl
		<< "l=" << pd->l << ", g=" << pd->g
			<< std::endl
		<< "x={" << pd->x[0] << "," << pd->x[1] << ","
			<< pd->x[2] << "," << pd->x[3] << "}" << std::endl
		<< "xP={" << pd->xP[0] << "," << pd->xP[1] << ","
			<< pd->xP[2] << "," << pd->xP[3] << "}" << std::endl;

	*pp = (void*)pd;

	return 0;
}

static std::ostream&
help(void *p, std::ostream& o)
{
	return o << "deformable pendulum" << std::endl;
}

static int
size(void* p)
{
	return 4;
}

static int
init(void* p, VectorHandler& X, VectorHandler& XP)
{
	private_data* pd = (private_data*)p;

	X.Reset();
	XP.Reset();
	for (int i = 1; i <= size(p); i++) {
		/* derivata iniziale */
		XP(i) = pd->xP[i - 1];
		/* posizione iniziale */
		X(i) = pd->x[i - 1];
	}

	return 0;
}

static int
grad(void* p, MatrixHandler& J, MatrixHandler& JP,
		const VectorHandler& X, const VectorHandler& XP,
		const doublereal& t)
{
	private_data* pd = (private_data*)p;

	doublereal theta = X(1);
	doublereal u = X(2);
	doublereal phi = X(3);
	doublereal w = X(4);

	doublereal ctheta = cos(theta);
	doublereal stheta = sin(theta);
	doublereal m = pd->m;
	doublereal c = pd->c;
	doublereal k = pd->k;
	doublereal l0 = pd->l;
	doublereal l = l0+u;
	doublereal g = pd->g;

	J(1, 3) = -1.;
	J(2, 4) = -1.;
	J(3, 1) = g*ctheta/l;
	J(3, 2) = -(2.*phi*w + g*stheta)/(l*l);
	J(3, 3) = 2.*w/l;
	J(3, 4) = 2.*phi/l;
	J(4, 1) = g*stheta;
	J(4, 2) = k/m - phi*phi;
	J(4, 3) = -2*l*phi;
	J(4, 4) = c/m;

	for (int i = 1; i <= 4; i++) {
		JP(i, i) = 1.;
	}

	return 0;
}

static int
func(void* p, VectorHandler& R,
		const VectorHandler& X, const VectorHandler& XP,
		const doublereal& t)
{
	private_data* pd = (private_data*)p;

	doublereal theta = X(1);
	doublereal u = X(2);
	doublereal phi = X(3);
	doublereal w = X(4);

	doublereal ctheta = cos(theta);
	doublereal stheta = sin(theta);
	doublereal m = pd->m;
	doublereal c = pd->c;
	doublereal k = pd->k;
	doublereal l0 = pd->l;
	doublereal l = l0+u;
	doublereal g = pd->g;

	R(1) = phi - XP(1);
	R(2) = w - XP(2);
	R(3) = - (2.*phi*w + g*stheta)/l - XP(3);
	R(4) = (phi*phi*l - (k*u + c*w)/m + g*ctheta) - XP(4);

	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
	private_data* pd = (private_data*)p;

	doublereal theta = X(1);
	doublereal ctheta = cos(theta);
	doublereal stheta = sin(theta);
	doublereal u = X(2);
	doublereal phi = X(3);
	doublereal w = X(4);
	doublereal m = pd->m;
	doublereal k = pd->k;
	doublereal l = pd->l+u;
	doublereal g = pd->g;
	doublereal x = l*stheta;
	doublereal y = -l*ctheta;
	doublereal xP = w*stheta+l*ctheta*phi;
	doublereal yP = -w*ctheta+l*stheta*phi;

	doublereal E = .5*m*(xP*xP+yP*yP)+m*g*y+.5*k*u*u;

	return o
		<< theta	/*  3 */
		<< " " << u	/*  4 */
		<< " " << phi	/*  5 */
		<< " " << w	/*  6 */
		<< " " << XP(1)	/*  7 */
		<< " " << XP(2)	/*  8 */
		<< " " << XP(3)	/*  9 */
		<< " " << XP(4)	/* 10 */
		<< " " << x 	/* 11 */
		<< " " << y 	/* 12 */
		<< " " << xP 	/* 11 */
		<< " " << yP 	/* 12 */
		<< " " << E;	/* 13 */
}

static int
destroy(void** p)
{
	private_data* pd = (private_data*)(*p);

	delete pd;
	*p = NULL;

	return 0;
}

static struct funcs funcs_handler = {
	read,
	help,
	init,
	size,
	grad,
	func,
	out,
	destroy
};

/* de-mangle name */
extern "C" {
void *ff = &funcs_handler;
}

