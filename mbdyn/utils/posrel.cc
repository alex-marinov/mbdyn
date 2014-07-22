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

#include <cstring>
#include <iostream>

#include "matvec3.h"


int
main(int argn, const char* const argv[])
{
   static doublereal d[3];

   while (true) {
      std::cin >> d[0];
      if (std::cin) {
	 std::cin >> d[1] >> d[2];
	 Vec3 x1(d);
	 std::cin >> d[0] >> d[1] >> d[2];
	 Mat3x3 R1(EulerAngles2MatR(Vec3(d)/180.*M_PI));
	 std::cin >> d[0] >> d[1] >> d[2];
	 Vec3 x2(d);

	 std::cout << R1.Transpose()*(x2-x1) << std::endl;
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}

