/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "matvec3.h"

int 
main(int argn, const char* const argv[])
{ 
   	if (argn > 1) {
      		if (!strcasecmp(argv[1], "-?")
	  	    || !strcasecmp(argv[1], "-h") 
		    || !strcasecmp(argv[1], "--help")) {
	 		std::cerr << std::endl 
				<< "usage: " << argv[0] << std::endl 
				<< std::endl
	   			<< "    reads the Euler angles (in degs)"
				" from stdin;" << std::endl
	   			<< "    writes the rotation matrix"
				" on standard output" << std::endl
				<< "    (m11, m12, m13,"
				" m21, m22, m23,"
				" m31, m32, m33)" << std::endl
				<< std::endl
	   			<< "part of MBDyn package (Copyright (C)"
				" Pierangelo Masarati, 1996-2006)" << std::endl
				<< std::endl;
	 		exit(EXIT_SUCCESS);
      		}
   	}

	std::cout.precision(16);

   	static doublereal d[3];
   	while (true) {
      		std::cin >> d[0];
      		if (std::cin) {
	 		std::cin >> d[1] >> d[2];
	 		std::cout << EulerAngles2MatR(Vec3(d)/(180./M_PI)) << std::endl;
      		} else {
	 		break;
      		}
   	}
   
   	return EXIT_SUCCESS;
}

