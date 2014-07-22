/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
	doublereal x[2];
	doublereal xP[2];
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

		in >> pd->l >> pd->g >> pd->x[0] >> pd->x[1];

	} else {
		pd->l = 1.;
		pd->g = 9.81;
		pd->x[0] = 0.;
		pd->x[1] = 1.;
	}

	pd->xP[0] = pd->x[1];
	pd->xP[1] = -pd->g*sin(pd->x[0])/pd->l;

	std::cerr
		<< "l=" << pd->l << ", g=" << pd->g << std::endl
		<< "x={" << pd->x[0] << "," << pd->x[1] << "}" << std::endl
		<< "xP={" << pd->xP[0] << "," << pd->xP[1] << "}" << std::endl;

	return 0;
}

static int
size(void* p)
{
	return 2;
}

static int
init(void* p, VectorHandler& X, VectorHandler& XP)
{
	private_data* pd = (private_data*)p;

	X.Reset();
	XP.Reset();

	for (int i = 1; i <= size(p); i++) {
		/* velocita' iniziale */
		XP(i) = pd->xP[i-1];
		/* posizione iniziale */
		X(i) = pd->x[i-1];
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

	doublereal l = pd->l;
	doublereal g = pd->g;

	J(1, 2) = -1.;
	J(2, 1) = g*cos(theta)/l;

	for (int i = 1; i <= 2; i++) {
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
	doublereal phi = X(2);

	doublereal l = pd->l;
	doublereal g = pd->g;

	R(1) = phi - XP(1);
	R(2) = -g*sin(theta)/l - XP(2);

	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
	private_data* pd = (private_data*)p;

	doublereal theta = X(1);
	doublereal ctheta = cos(theta);
	doublereal stheta = sin(theta);
	doublereal phi = X(2);
	doublereal l = pd->l;
	doublereal g = pd->g;
	doublereal x = l*stheta;
	doublereal y = -l*ctheta;
	doublereal xP = l*ctheta*phi;
	doublereal yP = l*stheta*phi;

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
		<< " " << E;		/* 11 */
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

