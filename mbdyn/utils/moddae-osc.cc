/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/iostream>
#include <ac/fstream>
#include <math.h>

#include <myassert.h>
#include <solman.h>
#include <harwrap.h>

#include <dae-intg.h>

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
		// std::cerr << "opening file \"" << user_defined 
		// 	<< "\"" << std::endl;
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
	
	// std::cerr << "m=" << pd->m << ", c=" << pd->c 
	// 	<< ", k=" << pd->k << std::endl
	//      << "x={" << pd->x[0] << "," << pd->x[1] << "}" << std::endl;
	
	return 0;
}

static int
size(void* p)
{
	// private_data* pd = (private_data*)p;
	return 2;
}

static int
init(void* p, VectorHandler& X, VectorHandler& XP)
{
	private_data* pd = (private_data*)p;
	X.Reset(0.);
	for (int i = 1; i <= size(p); i++) {      
		X.fPutCoef(i, pd->x[i-1]); /* posiz. iniziale */
	}
	return 0;
}

static int
grad(void* p, MatrixHandler& J, MatrixHandler& JP,
     const VectorHandler& X, const VectorHandler& XP, const doublereal& t)
{
	private_data* pd = (private_data*)p;
	J.fPutCoef(1, 2, -1.);
	J.fPutCoef(2, 1, pd->k/pd->m);
	J.fPutCoef(2, 2, pd->c/pd->m);

	JP.fPutCoef(1, 1, 1.);
	JP.fPutCoef(2, 2, 1.);
	return 0;
}

static int
func(void* p, VectorHandler& R, const VectorHandler& X, 
		const VectorHandler& XP, const doublereal& t)
{
	private_data* pd = (private_data*)p;
	doublereal x = X.dGetCoef(1);
	doublereal v = X.dGetCoef(2);
	doublereal xP = XP.dGetCoef(1);
	doublereal vP = XP.dGetCoef(2);
	R.fPutCoef(1, v - xP);
	R.fPutCoef(2, -vP - (pd->k*x + pd->c*v)/pd->m);
	return 0;
}

static std::ostream&
out(void* p, std::ostream& o, const VectorHandler& X, const VectorHandler& XP)
{
   	// private_data* pd = (private_data*)p;
   	return o << X.dGetCoef(1) << " " << X.dGetCoef(2)
     		<< " " << XP.dGetCoef(1) << " " << XP.dGetCoef(2);
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
	init,
	size,
	grad,
	func,
	out,
	destroy
};

void* ff = (void *)&_ff;

