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

#include <matvec3.h>
#include <string.h>


int main(int argn, const char* const argv[])
{
   if (argn > 1) {
      if (!strcasecmp(argv[1], "-?")
	  || !strcasecmp(argv[1], "-h") 
	  || !strcasecmp(argv[1], "--help")) {
	 cerr << endl << "usage: " << argv[0] << endl << endl
	   << "    reads a rotation matrix (row-oriented) from stdin;" << endl
	   << "    writes the rotation vector {magniture (in degs), direction} on standard output" << endl << endl
	   << "part of MBDyn package (Copyright (C) Pierangelo Masarati, 1996)" << endl << endl;
	 exit(EXIT_SUCCESS);
      }
   }   

   static doublereal d[9];
   while (1) {
      cin >> d[0];
      if (cin) {
	 cin >> d[3] >> d[6] >> d[1] >> d[4] >> d[7] >> d[2] >> d[5] >> d[8];
	 Vec3 g(gparam(Mat3x3(d, 3)));
	 doublereal d(g.Norm());
	 if (d > 0.) {
	    Vec3 phi(g*(2./d*atan(d/2.)));
	    d = phi.Norm();
	    if (d != 0.) {
	       cout << d*180./M_PI << " " << phi/d << endl;
	    } else {
	       cout << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
	    }
	 } else {
	    cout << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
	 }
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}
