/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
 *
 * This file is copyrighted by:
 * 	Alessandro Fumagalli	< alessandro.fumagalli@polimi.it >
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <matvec3.h>
#include <matvecexp.h>
#include <Rot.hh>


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
				" from command line;" << std::endl
	   			<< " writes the rotation vector {magniture (in degs), direction} on standard output" 
				<< std::endl
				<< std::endl
	   			<< "part of MBDyn package (Copyright (C)"
				" Pierangelo Masarati, 1996-2006)" << std::endl
				<< " Author: Alessandro Fumagalli" << std::endl
				<< std::endl;
	 		exit(EXIT_SUCCESS);
      		}
   	}

	std::cout.precision(16);

   	static doublereal d[3];
   	if ( argn > 2)	{
		d[0] = atof(argv[1]);
		d[1] = atof(argv[2]);
		d[2] = atof(argv[3]);


		Mat3x3 R(EulerAngles2MatR(Vec3(d)/(180./M_PI)));

		Vec3 phi(RotManip::VecRot(R));
		 doublereal D = phi.Norm();
		 if (D != 0.) {
		 	std::cout << D*180./M_PI << " " << phi/D << std::endl;
	 	} else {
	    		std::cout << 0. << " " << 0. << " " << 0. << " " << 0. << std::endl;
	 	}
	}
   
   	return EXIT_SUCCESS;
}

