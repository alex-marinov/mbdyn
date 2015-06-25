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

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "matvec3.h"
#include "matvec6.h"
#include "matvecexp.h"
#include "Rot.hh"

int 
main(int argn, const char* const argv[])
{
   	static doublereal d[9];
   	while (true) {
      		std::cin >> d[0];
      		if (std::cin) {
	 		std::cin >> d[1] >> d[2];

#if 0
			Mat3x3 Gamma_m1 = RotManip::DRot_I(Vec3(d));
			Mat3x3 Gamma = Gamma_m1.Inv();
			doublereal dDet = Gamma.dDet();
	 		std::cout << "Gamma=" << Gamma << " det(Gamma)=" << dDet << std::endl;
#endif

#if 0
			Mat3x3 R(RotManip::Rot(d));
			Vec3 Theta(RotManip::VecRot(R));
			Mat3x3 Rcheck(RotManip::Rot(Theta));
			std::cout << "R={" << R << "}" << std::endl
				<< "  V={" << Theta << "}" << std::endl
				<< "  R={" << Rcheck << "}" << std::endl
				<< "  I={" << R.MulMT(Rcheck) << "}"
				<< std::endl;
#endif

	 		std::cin >> d[3] >> d[4] >> d[5] >> d[6] >> d[7] >> d[8];
			Mat3x3 M(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);
			doublereal b[3];
			std::cin >> b[0] >> b[1] >> b[2];
			Vec3 x = M.LDLSolve(Vec3(b[0], b[1], b[2]));
			std::cout << x << std::endl;
      		} else {
	 		break;
      		}
   	}
   
   	return EXIT_SUCCESS;
}

