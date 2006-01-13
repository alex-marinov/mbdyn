/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#include <string.h>
#include <ac/iostream>

#include <matvec3.h>
#include <matvecexp.h>
#include <Rot.hh>

int main(int argn, const char* const argv[])
{
	bool normalize(false);

	if (argn > 1) {
		if (strcasecmp(argv[1], "--normalize") == 0) {
			normalize = true;
		}
	}

	static doublereal d1[3], d2[3];
	while (true) {
		std::cin >> d1[V1];
		if (!std::cin) {
			break;
		}
		std::cin >> d1[V2] >> d1[V3] 
			>> d2[V1] >> d2[V2] >> d2[V3];
		Vec3 v1(d1), v2(d2);
		if (normalize) {
			v1 /= v1.Norm();
			v2 /= v2.Norm();
		}

		Vec3 phi(v1.Cross(v2));
		doublereal d = phi.Norm();
		if (d != 0.) {
			if (normalize) {
				std::cout << std::asin(d)*180./M_PI;

			} else {
				std::cout << d;
			}
			std::cout << " " << phi/d << std::endl;
		} else {
			std::cout << 0. << " "
				<< 0. << " " << 0. << " " << 0.  << std::endl;
		}
	}      

	return (EXIT_SUCCESS);
}

