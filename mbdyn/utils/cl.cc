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
#include <sstream>

#include "input.h"
#include "mathp.h"

/* plugins per math parser */
#include "dummypgin.h"

int 
main(int argc, const char* const argv[])
{
	int verbose = 1;

	if (argc > 1) {
		if (strcmp(argv[1], "-?") == 0 
		    || strcmp(argv[1], "-h") == 0
		    || strcmp(argv[1], "--help") == 0) {
		    	std::cerr 
				<< "usage: " << argv[0] 
				<< "     reads from stdin" << std::endl
		   		<< "       " << argv[0] << " {-?|-h|--help}"
		   		" prints this message" << std::endl
	   			<< "       " << argv[0] 
				<< " <arg1> [<arg2> ...]"
				" evaluates the expressions" << std::endl;
	 		exit(EXIT_SUCCESS);
      		} else if (strcmp(argv[1], "-s") == 0) {
			verbose = 0;
			argv++;
			argc--;
		}
	}
	
#ifdef USE_TABLE
	//Table t(1);
	Table t(true);
#endif /* USE_TABLE */

	if (argc > 1) {
	 	for (int i = 1; i < argc; i++) {
			std::istringstream in(argv[i]);
	    		InputStream In(in);
#ifdef USE_TABLE
	    		MathParser mp(In, t);
#else /* !USE_TABLE */
	    		MathParser mp(In);
#endif /* !USE_TABLE */

			if (verbose) {
	    			std::cout << "argv[" << i << "] = ";
			}

			try {
	    			mp.GetForever(std::cout, "; ");
	    			std::cout << std::endl;
			} catch (...) {
      				std::cout << std::endl;
				exit(EXIT_FAILURE);
			}
	 	}
	 	exit(EXIT_SUCCESS);
      	}

	InputStream In(std::cin);
#ifdef USE_TABLE
      	MathParser mp(In, t);
#else /* !USE_TABLE */
      	MathParser mp(In);
#endif /* !USE_TABLE */

	try {	
		mp.GetForever(std::cout, "\n");
      		std::cout << std::endl;
      		exit(EXIT_SUCCESS);
	} catch (...) {
      		std::cout << std::endl;
		exit(EXIT_FAILURE);
	}
}

