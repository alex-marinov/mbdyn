/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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
#include <iomanip>

#include "dgeequ.h"

#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "naivemh.h"

static doublereal mat[5][5] = {
	{ 11.,  0., 13.,  0., 15. },
	{  0., 22.,  0., 24.,  0. },
	{ 31.,  0., 33.,  0., 35. },
	{  0., 42.,  0., 44.,  0. },
	{ 51.,  0., 53.,  0., 55. }
};

int
main(void)
{
	std::vector<integer> perm(5), invperm(5);
	perm[0] = 4;
	perm[1] = 3;
	perm[2] = 2;
	perm[3] = 1;
	perm[4] = 0;
	for (int i = 0; i < 5; i++) {
		invperm[perm[i]] = i;
	}

	NaiveMatrixHandler nm(5);
	NaivePermMatrixHandler npm(5, perm, invperm);
	FullMatrixHandler fm(5);
	SpMapMatrixHandler spm(5, 5);

	nm.Reset();
	npm.Reset();
	fm.Reset();
	spm.Reset();
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			if (mat[ir][ic] != 0.) {
				nm(ir + 1, ic + 1) = mat[ir][ic];
				npm(ir + 1, ic + 1) = mat[ir][ic];
				fm(ir + 1, ic + 1) = mat[ir][ic];
				spm(ir + 1, ic + 1) = mat[ir][ic];
			}
		}
	}

	std::vector<doublereal> Ax0, Ax1;
	std::vector<integer> Ai0, Ai1, Ap0, Ap1;

	spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);
	spm.MakeCompressedColumnForm(Ax1, Ai1, Ap1, 1);

	CColMatrixHandler<0> ccm0(Ax0, Ai0, Ap0);
	CColMatrixHandler<1> ccm1(Ax1, Ai1, Ap1);
	DirCColMatrixHandler<0> dirccm0(Ax0, Ai0, Ap0);
	DirCColMatrixHandler<1> dirccm1(Ax1, Ai1, Ap1);

	std::vector<doublereal> r, c;
	doublereal amax, rowcnd, colcnd;

	dgeequ<NaiveMatrixHandler, NaiveMatrixHandler::const_iterator>(nm, r, c, rowcnd, colcnd, amax);
	std::cout << "naive: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif
	dgeequ_scale<NaiveMatrixHandler, NaiveMatrixHandler::const_iterator>(nm, r, c);
	std::cout << nm << std::endl;

	dgeequ<NaivePermMatrixHandler, NaivePermMatrixHandler::const_iterator>(npm, r, c, rowcnd, colcnd, amax);
	std::cout << "naive permuted: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif
	dgeequ_scale<NaivePermMatrixHandler, NaivePermMatrixHandler::const_iterator>(npm, r, c);
	std::cout << npm << std::endl;

	dgeequ(fm, r, c, rowcnd, colcnd, amax);
	std::cout << "full: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif

	dgeequ_scale(fm, r, c);
	std::cout << fm << std::endl;

	dgeequ<CColMatrixHandler<0>, CColMatrixHandler<0>::const_iterator>(ccm0, r, c, rowcnd, colcnd, amax);
	std::cout << "ccol0: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif

	dgeequ_scale<CColMatrixHandler<0>, CColMatrixHandler<0>::const_iterator>(ccm0, r, c);
	std::cout << ccm0 << std::endl;

	dgeequ<CColMatrixHandler<1>, CColMatrixHandler<1>::const_iterator>(ccm1, r, c, rowcnd, colcnd, amax);
	std::cout << "ccol1: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif

	dgeequ_scale<CColMatrixHandler<1>, CColMatrixHandler<1>::const_iterator>(ccm1, r, c);
	std::cout << ccm1 << std::endl;

	dgeequ<DirCColMatrixHandler<0>, DirCColMatrixHandler<0>::const_iterator>(dirccm0, r, c, rowcnd, colcnd, amax);
	std::cout << "dirccol0: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif

	dgeequ_scale<DirCColMatrixHandler<0>, DirCColMatrixHandler<0>::const_iterator>(dirccm0, r, c);
	std::cout << dirccm0 << std::endl;

	dgeequ<DirCColMatrixHandler<1>, DirCColMatrixHandler<1>::const_iterator>(dirccm1, r, c, rowcnd, colcnd, amax);
	std::cout << "dirccol1: amax=" << amax << ", rowcnd=" << rowcnd << ", colcnd=" << colcnd << std::endl;
	for (unsigned ir = 0; ir < 5; ir++) {
		std::cout
			<< "   r[" << ir << "]=" << std::setw(12) << r[ir]
			<< "       c[" << ir << "]=" << std::setw(12) << c[ir]
			<< std::endl;
	}

#if 0
	for (unsigned ir = 0; ir < 5; ir++) {
		for (unsigned ic = 0; ic < 5; ic++) {
			std::cout << std::setw(12) << r[ir]*c[ic]*mat[ir][ic];
		}
		std::cout << std::endl;
	}
#endif

	dgeequ_scale<DirCColMatrixHandler<1>, DirCColMatrixHandler<1>::const_iterator>(dirccm1, r, c);
	std::cout << dirccm1 << std::endl;

	return 0;
}

