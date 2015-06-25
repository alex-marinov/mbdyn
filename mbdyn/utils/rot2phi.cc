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

#include <cstring>
#include <iostream>

#include "matvec3.h"
#include "matvecexp.h"
#include "Rot.hh"

int main(int argn, const char* const argv[])
{
   if (argn > 1) {
      if (!strcasecmp(argv[1], "-?")
	  || !strcasecmp(argv[1], "-h") 
	  || !strcasecmp(argv[1], "--help")) {
	 std::cerr << std::endl << "usage: " << argv[0] << std::endl << std::endl
	   << "    reads a rotation matrix (row-oriented) from stdin;" << std::endl
	   << "    writes the rotation vector {magniture (in degs), direction} on standard output" << std::endl << std::endl
	   << "part of MBDyn package (Copyright (C) Pierangelo Masarati, 1996)" << std::endl << std::endl;
	 exit(EXIT_SUCCESS);
      }
   }   

	std::cout.precision(16);

   static doublereal d[9];
   while (true) {
      std::cin >> d[M11];
      if (std::cin) {
	 std::cin >> d[M12] >> d[M13] 
	 	>> d[M21] >> d[M22] >> d[M23] >> d[M31] >> d[M32] >> d[M33];
	 Vec3 phi(RotManip::VecRot(Mat3x3(d, 3)));
	 doublereal d = phi.Norm();
	 if (d != 0.) {
	    std::cout << d*180./M_PI << " " << phi/d << std::endl;
	 } else {
	    std::cout << 0. << " " << 0. << " " << 0. << " " << 0. << std::endl;
	 }
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}

