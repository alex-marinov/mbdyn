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
#include <mbconfig.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>

#include <ac/iostream>
#include <ac/fstream>

#include <myassert.h>
#include <mynewmem.h>

#if defined(USE_AERODYNAMIC_ELEMS)
#include <aerodata.h>
#include <aerodc81.h>
#include <c81data.h>
#endif /* USE_AERODYNAMIC_ELEMS */

/*
 * uso: testc81 file alpha mach
 */
int
main(int argc, char *argv[])
{
#if defined(USE_AERODYNAMIC_ELEMS)
	if (argc != 4) {
		std::cerr << "usage: testc81 <file> <alpha> <mach>" << std::endl;
		exit(EXIT_SUCCESS);
	}

	std::ifstream in(argv[1]);
	if (!in) {
		std::cerr << "unable to open file '" << argv[1] 
			<< "'" << std::endl;
		exit(EXIT_FAILURE);
	}

	c81_data* data = new C81Data(1);

	if (read_c81_data(in, data)) {
		std::cerr << "unable to read c81 data from file '" 
			<< argv[1] << "'" << std::endl;
		exit(EXIT_FAILURE);
	}

	double alpha(atof(argv[2])), mach(atof(argv[3]));

	std::cout 
		<< "alpha: " << alpha << std::endl
		<< "mach: " << mach << std::endl
		<< "cl: " << get_c81_coef(data->NML, data->ml, data->NAL, 
				data->al, alpha, mach) << std::endl
		<< "cd: " << get_c81_coef(data->NMD, data->md, data->NAD, 
				data->ad, alpha, mach) << std::endl
		<< "cm: " << get_c81_coef(data->NMM, data->mm, data->NAM, 
				data->am, alpha, mach) << std::endl;

#else /* !USE_AERODYNAMIC_ELEMS */
	std::cerr << "compile with --with-aero to enable aerodynamic stuff" 
		<< std::endl;
#endif /* !USE_AERODYNAMIC_ELEMS */

	return 0;
}
