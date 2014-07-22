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

#include "myassert.h"
#include "solman.h"
#include "harwrap.h"

struct private_data {
   int i;
};

int read(void** pp, const char* user_defined)
{
   *pp = (void*)new private_data;
   return 0;
}

int init(void* p, VectorHandler& X, VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   X.Reset(0.);
   XP.Reset(0.);
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 1;
}

int jac(void* p, MatrixHandler& JP, MatrixHandler& J,
	const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   J.PutCoef(1, 1, 1.);
   return 0;
}

int res(void* p, VectorHandler& R, const doublereal&, 
	const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   R.PutCoef(1, -1.*X(1)-1.*XP(1));
   return 0;
}

ostream& out(void* p, ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   // private_data* pd = (private_data*)p;
   return o;
}

int destroy(void** p)
{
   // private_data* pd = (private_data*)p;
   delete *p;
   *p = NULL;
   return 0;
}
