/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#include <myassert.h>
#include <solman.h>
#include <fstream.h>
#include <math.h>

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
      ifstream in(user_defined);
      if (!in) {
	 cerr << "unable to open file \"" << user_defined << "\"" << endl;
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
   
   cerr << "m=" << pd->m << ", c=" << pd->c << ", k=" << pd->k << endl
     << "l=" << pd->l << ", g=" << pd->g << endl
     << "x={" << pd->x[0] << "," << pd->x[1] << "," 
     << pd->x[2] << "," << pd->x[3] << "}" << endl;
   
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
   X.Reset(0.);
   for (int i = 1; i <= size(p); i++) {      
      X.fPutCoef(i, pd->x[i-1]); /* posiz. iniziale */
   }
   return 0;
}

int grad(void* p, MatrixHandler& J, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal theta = X.dGetCoef(1);
   doublereal u = X.dGetCoef(2);
   doublereal phi = X.dGetCoef(3);
   doublereal w = X.dGetCoef(4);
   
   doublereal ctheta = cos(theta);
   doublereal stheta = sin(theta);
   doublereal m = pd->m;
   doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = l0+u;
   doublereal g = pd->g;
   
   J.fPutCoef(1, 3, 1.);
   J.fPutCoef(2, 4, 1.);
   J.fPutCoef(3, 1, -g*ctheta/l);
   J.fPutCoef(3, 2, (2.*phi*w+g*stheta)/(l*l));
   J.fPutCoef(3, 3, -2.*w/l);
   J.fPutCoef(3, 4, -2.*phi/l);
   J.fPutCoef(4, 1, -g*stheta);
   J.fPutCoef(4, 2, phi*phi-k/m);
   J.fPutCoef(4, 3, 2*l*phi);
   J.fPutCoef(4, 4, -c/m);

   return 0;
}

int func(void* p, VectorHandler& R, const VectorHandler& X, const doublereal& t)
{
   private_data* pd = (private_data*)p;
   
   doublereal theta = X.dGetCoef(1);
   doublereal u = X.dGetCoef(2);
   doublereal phi = X.dGetCoef(3);
   doublereal w = X.dGetCoef(4);
   
   doublereal ctheta = cos(theta);
   doublereal stheta = sin(theta);
   doublereal m = pd->m;
   doublereal c = pd->c;
   doublereal k = pd->k;
   doublereal l0 = pd->l;
   doublereal l = l0+u;
   doublereal g = pd->g;

   R.fPutCoef(1, phi);
   R.fPutCoef(2, w);
   R.fPutCoef(3, -(2.*phi*w+g*stheta)/l);
   R.fPutCoef(4, phi*phi*l-k/m*u-c/m*w+g*ctheta);

   return 0;
}

ostream& out(void* p, ostream& o, 
	     const VectorHandler& X, const VectorHandler& XP)
{
   private_data* pd = (private_data*)p;
  
   doublereal theta = X.dGetCoef(1);
   doublereal ctheta = cos(theta);
   doublereal stheta = sin(theta);
   doublereal u = X.dGetCoef(2);
   doublereal phi = X.dGetCoef(3);
   doublereal w = X.dGetCoef(4);
   doublereal m = pd->m;
   doublereal k = pd->k;
   doublereal l = pd->l+u;
   doublereal g = pd->g;
   doublereal x = l*stheta;
   doublereal y = -l*ctheta;
   doublereal xp = w*stheta+l*ctheta*phi;
   doublereal yp = -w*ctheta+l*stheta*phi;
   
   doublereal E = .5*m*(xp*xp+yp*yp)+m*g*y+.5*k*u*u;
  
   
   return o << theta << " " << u
     << " " << X.dGetCoef(3) << " " << X.dGetCoef(4)
     << " " << XP.dGetCoef(1) << " " << XP.dGetCoef(2)
     << " " << XP.dGetCoef(3) << " " << XP.dGetCoef(4)
     << " " << x << " " << y << " " << E;
}

int destroy(void** p)
{
   // private_data* pd = (private_data*)p;
   delete *p;
   *p = NULL;
   return 0;
}
