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

struct private_data {
   doublereal m;
   doublereal c;
   doublereal k;
   doublereal l;
   doublereal g;
   doublereal x[4];
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
   
   doublereal theta = pd->x[0];
   doublereal uu = pd->x[1];
   doublereal thetap = pd->x[2];
   doublereal uup = pd->x[3];
   
   pd->x[0] = (pd->l+uu)*sin(theta);
   pd->x[1] = -(pd->l+uu)*cos(theta);
   pd->x[2] = thetap*(pd->l+uu)*cos(theta)+uup*sin(theta);
   pd->x[3] = thetap*(pd->l+uu)*sin(theta)-uup*cos(theta);   
   
   std::cerr << "m=" << pd->m << ", c=" << pd->c << ", k=" << pd->k << std::endl
     << "l=" << pd->l << ", g=" << pd->g << std::endl
     << "x={" << pd->x[0] << "," << pd->x[1] << "," 
     << pd->x[2] << "," << pd->x[3] << "}" << std::endl;
   
   return 0;
}

int size(void* p)
{
   // private_data* pd = (private_data*)p;
   return 4;
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
   
   doublereal x = X(1);
   doublereal y = X(2);
//    doublereal u = X(3);
//    doublereal v = X(4);
   
   doublereal l = sqrt(x*x+y*y);
   doublereal m = pd->m;
   // doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l; 
//    doublereal g = pd->g;
   
   J.PutCoef(1, 3, 1.);
   J.PutCoef(2, 4, 1.);
   J.PutCoef(3, 1, -k/m*(1.-l0/l*(1.-(x*x)/(l*l))));
   J.PutCoef(3, 2, -k/m*x*y*l0/(l*l*l));   
   J.PutCoef(4, 1, -k/m*x*y*l0/(l*l*l));
   J.PutCoef(4, 2, -k/m*(1.-l0/l*(1.-(y*y)/(l*l))));

   return 0;
}

int func(void* p, VectorHandler& R, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal x = X(1);
   doublereal y = X(2);
   doublereal u = X(3);
   doublereal v = X(4);
   
   doublereal m = pd->m;
   // doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = sqrt(x*x+y*y);
   doublereal g = pd->g;

   R.PutCoef(1, u);
   R.PutCoef(2, v);
   R.PutCoef(3, -k/m*x*(1.-l0/l));
   R.PutCoef(4, -k/m*y*(1.-l0/l)-g);

   return 0;
}

std::ostream& out(void* p, std::ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   private_data* pd = (private_data*)p;
  
   doublereal x = X(1);
   doublereal y = X(2);
   doublereal u = X(3);
   doublereal v = X(4);   
   doublereal m = pd->m;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = sqrt(x*x+y*y);
   doublereal theta = atan2(x, -y);
   doublereal g = pd->g;
   
   doublereal E = .5*m*(u*u+v*v)+m*g*y+.5*k*(l-l0)*(l-l0);
  
   
   return o << theta << " " << (l-l0)
     << " " << X(3) << " " << X(4)
     << " " << XP(1) << " " << XP(2)
     << " " << XP(3) << " " << XP(4)
     << " " << x << " " << y << " " << E;
}

int destroy(void** p)
{
   private_data* pd = (private_data*)(*p);
   delete pd;
   *p = NULL;
   return 0;
}
