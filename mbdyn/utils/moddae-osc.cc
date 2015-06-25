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

#include <iostream>
#include <fstream>
#include <cmath>

#include "myassert.h"
#include "solman.h"
#include "harwrap.h"

#include "dae-intg.h"

struct private_data {
   	doublereal m;
   	doublereal c;
   	doublereal k;
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
			std::cerr << "unable to open file \"" << user_defined
				<< "\"" << std::endl;
			exit(EXIT_FAILURE);
		}
		in >> pd->m >> pd->c >> pd->k >> pd->x[0] >> pd->x[1];
	} else {
		pd->m = 1.;
		pd->c = 1.e-2;
		pd->k = 1.;     
		pd->x[0] = 0.;
		pd->x[1] = 0.;
	}

	pd->xP[0] = pd->x[1];
	pd->xP[1] = - (pd->k*pd->x[0] + pd->c*pd->x[1])/pd->m;
	
	return 0;
}

static std::ostream&
help(void *p, std::ostream& o)
{
	private_data* pd = (private_data*)p;

	return o
		<< "harmonic oscillator:" << std::endl
		<< std::endl
		<< "\tm \\ddot{x} + c \\dot{x} + k x = 0" << std::endl
		<< std::endl
		<< "\tm  = " << pd->m << std::endl
		<< "\tc  = " << pd->c << std::endl
		<< "\tk  = " << pd->k << std::endl
		<< "\tx0 = " << pd->x[0] << std::endl
		<< "\tv0 = " << pd->x[1] << std::endl
		<< std::endl;
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
	for (int i = 1; i <= size(p); i++) {      
		/* posiz. iniziale */
		X(i) = pd->x[i-1];
	}
	return 0;
}

static int
grad(void* p, MatrixHandler& J, MatrixHandler& JP,
     const VectorHandler& X, const VectorHandler& XP, const doublereal& t)
{
	private_data* pd = (private_data*)p;

	J(1, 2) = -1.;
	J(2, 1) = pd->k/pd->m;
	J(2, 2) = pd->c/pd->m;

	JP(1, 1) = 1.;
	JP(2, 2) = 1.;

	return 0;
}

static int
func(void* p, VectorHandler& R, const VectorHandler& X, 
		const VectorHandler& XP, const doublereal& t)
{
	private_data* pd = (private_data*)p;

	doublereal x = X(1);
	doublereal v = X(2);
	doublereal xP = XP(1);
	doublereal vP = XP(2);

	R(1) = v - xP;
	R(2) = -vP - (pd->k*x + pd->c*v)/pd->m;
	
	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
#if 0
   	private_data* pd = (private_data*)p;
#endif

   	return o << X(1) << " " << X(2)
     		<< " " << XP(1) << " " << XP(2);
}

static int
destroy(void** p)
{
   	private_data* pd = (private_data*)(*p);

   	delete pd;
   	*p = NULL;

   	return 0;
}

/* simboli da esportare */
static funcs _ff = {
	read,
	help,
	init,
	size,
	grad,
	func,
	out,
	destroy
};

void* ff = (void *)&_ff;

