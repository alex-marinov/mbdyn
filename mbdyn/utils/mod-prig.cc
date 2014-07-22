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

#include <iostream>
#include <fstream>
#include <cmath>

#include <myassert.h>
#include <solman.h>

struct private_data {
   doublereal m;
   doublereal l;
   doublereal g;
   doublereal x[2];
};

int read(void** pp, const char* user_defined)
{
   *pp = (void*)new private_data;
   private_data* pd = (private_data*)*pp;
   
   if (user_defined != NULL) {
      // cerr << "opening file \"" << user_defined << "\"" << endl;
      std::ifstream in(user_defined);
      if (!in) {
	 std::cerr << "unable to open file \"" << user_defined << "\"" << std::endl;
	 exit(EXIT_FAILURE);
      }
      in >> pd->m >> pd->l >> pd->g
	>> pd->x[0] >> pd->x[1];
   } else {
      pd->m = 1.;     
      pd->l = 1.;
      pd->g = 9.81;     
      pd->x[0] = 0.;
      pd->x[1] = 0.;
   }
   
   std::cerr << "m=" << pd->m << ", l=" << pd->l << ", g=" << pd->g << std::endl
     << "x={" << pd->x[0] << "," << pd->x[1] << "}" << std::endl;
   
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 2;
}

int init(void* p, VectorHandler& X)
{
   private_data* pd = (private_data*)p;
   X.Reset();
   for (int i = 1; i <= size(p); i++) {      
      X.PutCoef(i, pd->x[i-1]); /* posiz. iniziale */
   }
   return 0;
}

int grad(void* p, MatrixHandler& J, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal theta = X(1); 
//    doublereal phi = X(2);
   
   doublereal ctheta = cos(theta);
//    doublereal stheta = sin(theta);
//    doublereal m = pd->m;
   doublereal l = pd->l;
   doublereal g = pd->g;
   
   J.PutCoef(1, 2, 1.);
   J.PutCoef(2, 1, -g*ctheta/l);

   return 0;
}

int func(void* p, VectorHandler& R, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal theta = X(1);
   doublereal phi = X(2);
   
//    doublereal ctheta = cos(theta);
   doublereal stheta = sin(theta);
//    doublereal m = pd->m;
   doublereal l = pd->l;
   doublereal g = pd->g;

   R.PutCoef(1, phi);
   R.PutCoef(2, -g*stheta/l);

   return 0;
}

std::ostream& out(void* p, std::ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   private_data* pd = (private_data*)p;
  
   doublereal theta = X(1);
   doublereal ctheta = cos(theta);
   doublereal stheta = sin(theta);
   doublereal phi = X(2);
   doublereal m = pd->m;   
   doublereal l = pd->l;
   doublereal g = pd->g;
   doublereal x = l*stheta;
   doublereal y = -l*ctheta;
   doublereal xp = l*ctheta*phi;
   doublereal yp = l*stheta*phi;
   
   doublereal E = .5*m*(xp*xp+yp*yp)+m*g*y;
  
   
   return o << theta
     << " " << X(2)
     << " " << XP(1) << " " << XP(2)
     << " " << x << " " << y << " " << E;
}

int destroy(void** p)
{
   private_data* pd = (private_data*)(*p);
   delete pd;
   *p = NULL;
   return 0;
}
