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

#include <matvec3.h>
#include <string.h>


int 
main(int argn, const char* const argv[])
{
   	flag f(0);
	
   	if (argn > 1) {
      		if (!strcasecmp(argv[1], "-?")
	  	    || !strcasecmp(argv[1], "-h") 
	  	    || !strcasecmp(argv[1], "--help")) {
	 		cerr << endl 
				<< "usage: " << argv[0] 
				<< " [mat|euler]" << endl 
				<< endl
	   			<< "    reads the Euler angles (in degs)"
				" of bodies 1 and 2 from stdin;" << endl
	   			<< "    writes, on standard output:" << endl
	   			<< "        default|\"euler\", the relative"
				" Euler angles," << endl
	   			<< "        \"mat\", the relative rotation"
				" matrix (column-oriented)" << endl 
				<< endl
	   			<< "part of MBDyn package (Copyright (C)"
				" Pierangelo Masarati, 1996-2000)" << endl 
				<< endl;
	 		exit(EXIT_SUCCESS);
      		} else if (!strcasecmp(argv[1], "mat")) {
	 		f = flag(1);
      		} else if (!strcasecmp(argv[1], "euler")) {
	 		f = flag(0);
      		} else {
	 		f = flag(0);
      		}
   	}
   
   	static doublereal d[3];
   	while (1) {
      		cin >> d[0];
      		if (cin) {
	 		cin >> d[1] >> d[2];
	 		Mat3x3 R1(RFromEulerAngles(Vec3(d)/180.*M_PI));
	 		cin >> d[0] >> d[1] >> d[2];
	 		Mat3x3 R2(RFromEulerAngles(Vec3(d)/180.*M_PI));
			
	 		if (f) {
	    			cout << R1.Transpose()*R2 << endl;
	 		} else {
	    			cout << EulerAngles(R1.Transpose()*R2) << endl;
	 		}
      		} else {
	 		break;
      		}
   	}
   
   	return (EXIT_SUCCESS);
}

