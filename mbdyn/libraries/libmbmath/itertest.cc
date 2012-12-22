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

#include "fullmh.h"
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

	FullMatrixHandler fm(5);
	SpMapMatrixHandler spm(5, 5);
	NaiveMatrixHandler nm(5);
	NaivePermMatrixHandler npm(5, perm, invperm);

	for (int r = 0; r < 5; r++) {
		for (int c = 0; c < 5; c++) {
			if (mat[r][c] != 0.) {
				fm(r + 1, c + 1) = mat[r][c];
				spm(r + 1, c + 1) = mat[r][c];
				nm(r + 1, c + 1) = mat[r][c];
				npm(r + 1, c + 1) = mat[r][c];
			}
		}
	}

	std::cout << "matrix in full form: " << std::endl
		<< fm << std::endl;

	std::cout << "matrix in sparse form: " << std::endl
		<< spm << std::endl;

	std::cout << "matrix in naive form: " << std::endl
		<< nm << std::endl;

	std::cout << "matrix in naive permuted form: " << std::endl
		<< npm << std::endl;

	std::vector<doublereal> Ax0;
	std::vector<integer> Ai0, Ap0;
	spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);

	CColMatrixHandler<0> ccm0(Ax0, Ai0, Ap0);

	std::cout << "matrix in cc<0> form: " << std::endl
		<< ccm0 << std::endl;
	std::cout << "matrix in cc<0> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			// NOTE: need to const_cast because ostream::operator << (const double) does not exist,
			// thus non-const double& CColMatrixHandler<0>::operator()(int, int) would be used
			std::cout << std::setw(16) << const_cast<const CColMatrixHandler<0>& >(ccm0)(ir, ic);
		}
		std::cout << std::endl;
	}

	std::vector<doublereal> Ax1;
	std::vector<integer> Ai1, Ap1;
	spm.MakeCompressedColumnForm(Ax1, Ai1, Ap1, 1);

	CColMatrixHandler<1> ccm1(Ax1, Ai1, Ap1);

	std::cout << "matrix in cc<1> form: " << std::endl
		<< ccm1 << std::endl;
	std::cout << "matrix in cc<1> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			// NOTE: need to const_cast because ostream::operator << (const double) does not exist,
			// thus non-const double& CColMatrixHandler<1>::operator()(int, int) would be used
			std::cout << std::setw(16) << const_cast<const CColMatrixHandler<1>& >(ccm1)(ir, ic);
		}
		std::cout << std::endl;
	}

	DirCColMatrixHandler<0> dirm0(Ax0, Ai0, Ap0);
	std::cout << "matrix in dir<0> form: " << std::endl
		<< dirm0 << std::endl;
	std::cout << "matrix in dir<0> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			// NOTE: need to const_cast because ostream::operator << (const double) does not exist,
			// thus non-const double& DirCColMatrixHandler<0>::operator()(int, int) would be used
			std::cout << std::setw(16) << const_cast<const DirCColMatrixHandler<0>& >(dirm0)(ir, ic);
		}
		std::cout << std::endl;
	}

	DirCColMatrixHandler<1> dirm1(Ax1, Ai1, Ap1);
	std::cout << "matrix in dir<1> form: " << std::endl
		<< dirm1 << std::endl;
	std::cout << "matrix in dir<1> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			// NOTE: need to const_cast because ostream::operator << (const double) does not exist,
			// thus non-const double& DirCColMatrixHandler<1>::operator()(int, int) would be used
			std::cout << std::setw(16) << const_cast<const DirCColMatrixHandler<1>& >(dirm1)(ir, ic);
		}
		std::cout << std::endl;
	}

	std::cout << "***************************" << std::endl
		<< "full matrix handler:" << std::endl;

	for (FullMatrixHandler::const_iterator i = fm.begin();
		i != fm.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "sparse map matrix handler:" << std::endl;

	for (SpMapMatrixHandler::const_iterator i = spm.begin();
		i != spm.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "naive sparse matrix handler:" << std::endl;
	for (NaiveMatrixHandler::const_iterator i = nm.begin();
		i != nm.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "naive permuted sparse matrix handler:" << std::endl;
	for (NaivePermMatrixHandler::const_iterator i = npm.begin();
		i != npm.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "column compressed <0> sparse matrix handler:" << std::endl;
	for (CColMatrixHandler<0>::const_iterator i = ccm0.begin();
		i != ccm0.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "column compressed <1> sparse matrix handler:" << std::endl;
	for (CColMatrixHandler<1>::const_iterator i = ccm1.begin();
		i != ccm1.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "dir column compressed <0> sparse matrix handler:" << std::endl;
	for (DirCColMatrixHandler<0>::const_iterator i = dirm0.begin();
		i != dirm0.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	std::cout << "***************************" << std::endl
		<< "dir column compressed <1> sparse matrix handler:" << std::endl;
	for (DirCColMatrixHandler<1>::const_iterator i = dirm1.begin();
		i != dirm1.end(); ++i)
	{
		std::cout << "(" << i->iRow << ", " << i->iCol << ", " << i->dCoef << ")" << std::endl;
		if (mat[i->iRow][i->iCol] != i->dCoef) {
			std::cout << "==> failed!" << std::endl;
		}
	}

	return 0;
}

