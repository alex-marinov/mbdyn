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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <fstream>
#include <iostream>
#include <cmath>

#include "myassert.h"
#include "solman.h"
#include "dae-intg.h"

struct private_data {
	doublereal l;
	doublereal g;
	doublereal thetaold;
	doublereal x[5];
	doublereal xP[5];
};

static int
read(void** pp, const char* user_defined)
{
	*pp = (void*)new private_data;
	private_data* pd = (private_data*)*pp;
   
	if (user_defined != NULL) {
		std::ifstream in(user_defined);
		if (!in) {
			std::cerr << "unable to open file "
				"\"" << user_defined << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
		in >> pd->l >> pd->g >> pd->x[0] >> pd->x[1]
			>> pd->x[2] >> pd->x[3] >> pd->xP[4];
	} else {
		pd->l = 1.;
		pd->g = 9.81;     
		pd->x[0] = 0.;
		pd->x[1] = 1.;
		pd->x[2] = 1.;
		pd->x[3] = 0.;
		pd->xP[4] = 0.;
	}

	pd->x[4] = 0.;
	pd->thetaold = atan2(pd->x[0], -pd->x[1]);

	pd->xP[0] = pd->x[2];
	pd->xP[1] = pd->x[3];

	doublereal phi = sqrt(pd->x[0]*pd->x[0] + pd->x[1]*pd->x[1]) - pd->l;
	if (fabs(phi) > 1.e-9) {
		std::cerr << "constraint violated" << std::endl;
		exit(EXIT_FAILURE);
	}

	doublereal psi = (pd->x[0]*pd->xP[0] + pd->x[1]*pd->xP[1])/pd->l;
	if (fabs(psi) > 1.e-9) {
		std::cerr << "constraint derivative violated" << std::endl;
		exit(EXIT_FAILURE);
	}

	doublereal lambda = pd->l*(pd->x[2]*pd->x[2] + pd->x[3]*pd->x[3]
			- pd->x[1]*pd->g)/pd->l;

	if (fabs(lambda - pd->xP[4]) > 1.e-9) {
		std::cerr << "constraint reaction incorrect" << std::endl;
		pd->xP[4] = lambda;
	}

	pd->xP[2] = - (pd->x[0]*pd->xP[4])/pd->l;
	pd->xP[3] = - ((pd->x[1]*pd->xP[4])/pd->l + pd->g);

	std::cerr 
		<< "l=" << pd->l << ", g=" << pd->g << std::endl
		<< "x={" << pd->x[0] << "," << pd->x[1] << "," << pd->x[2]
		<< "," << pd->x[3] << "," << pd->x[4] << "}" << std::endl
		<< "xP={" << pd->xP[0] << "," << pd->xP[1] << "," << pd->xP[2]
		<< "," << pd->xP[3] << "," << pd->xP[4] << "}" << std::endl;

	return 0;
}

static int
size(void* p)
{
	return 5;
}

static int
init(void* p, VectorHandler& X, VectorHandler& XP)
{
	private_data* pd = (private_data*)p;
	
	X.Reset();
	XP.Reset();
	
	for (int i = 1; i <= size(p); i++) {
		XP(i) = pd->xP[i - 1]; /* posiz. iniziale */
		X(i) = pd->x[i - 1]; /* posiz. iniziale */
	}

	return 0;
}

static int
grad(void* p, MatrixHandler& J, MatrixHandler& JP, 
		const VectorHandler& X, const VectorHandler& XP,
		const doublereal& t)
{
	private_data* pd = (private_data*)p;

	doublereal x = X(1);
	doublereal y = X(2);
	doublereal lambda = XP(5);

	doublereal l = pd->l;

	J(1, 3) = -1.;
	J(2, 4) = -1.;
	J(3, 1) = lambda/l;
	J(4, 2) = lambda/l;
	J(5, 1) = x/l;
	J(5, 2) = y/l;

	for (int i = 1; i <= 4; i++) {
		JP(i, i) = 1.;
	}

	JP(3, 5) = x/l;
	JP(4, 5) = y/l;

	return 0;
}

static int
func(void* p, VectorHandler& R,
		const VectorHandler& X, const VectorHandler& XP,
		const doublereal& t)
{
	private_data* pd = (private_data*)p;

	doublereal x = X(1);
	doublereal y = X(2);
	doublereal u = X(3);
	doublereal v = X(4);
	doublereal xP = XP(1);
	doublereal yP = XP(2);
	doublereal uP = XP(3);
	doublereal vP = XP(4);
	doublereal lambda = XP(5);

	doublereal l = pd->l;
	doublereal g = pd->g;

	R(1) = u - xP;
	R(2) = v - yP;
	R(3) = - (uP + lambda*x/l);
	R(4) = - (vP + lambda*y/l + g);
	R(5) = l - sqrt(x*x + y*y);

	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
	private_data* pd = (private_data*)p;

	doublereal x = X(1);
	doublereal y = X(2);
	doublereal xP = X(3);
	doublereal yP = X(4);
	doublereal lambda = XP(5);

	doublereal theta = atan2(x, -y);
	while (theta > pd->thetaold + M_PI_2) {
		theta -= M_PI;
	}

	while (theta < pd->thetaold - M_PI_2) {
		theta += M_PI;
	}

	pd->thetaold = theta;

	doublereal phi = 0.;
	if (fabs(x) > fabs(y)) {
		phi = yP/x;
	} else {
		phi = -xP/y;
	}
	doublereal g = pd->g;

	doublereal E = .5*(xP*xP+yP*yP)+g*y;

	return o
		<< theta		/*  3 */
		<< " " << phi		/*  4 */
		<< " " << XP(1)		/*  5 */
		<< " " << XP(2)		/*  6 */
		<< " " << x 		/*  7 */
		<< " " << y 		/*  8 */
		<< " " << xP		/*  9 */
		<< " " << yP		/* 10 */
		<< " " << E		/* 11 */
		<< " " << lambda;	/* 12 */
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
	0,
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

