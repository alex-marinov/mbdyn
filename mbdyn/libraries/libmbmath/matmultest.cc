/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/*
 * Umfpack is used by permission; please read its Copyright,
 * License and Availability note.
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <iomanip>

#include "fullmh.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"

static doublereal mat[5][5] = {
	{ 11.,  0., 13.,  0., 15. },
	{  0., 22.,  0., 24.,  0. },
	{ 31.,  0., 33.,  0., 35. },
	{  0., 42.,  0., 44.,  0. },
	{ 51.,  0., 53.,  0., 55. }
};

static int
check_mat(MatrixHandler& mh)
{
	for (int r = 0; r < 5; r++) {
		for (int c = 0; c < 5; c++) {
			if (mat[r][c] != mh(r + 1, c + 1)) {
				return -1;
			}
		}
	}

	return 0;
}

static int
check_mat_transpose(MatrixHandler& mh)
{
	for (int r = 0; r < 5; r++) {
		for (int c = 0; c < 5; c++) {
			if (mat[c][r] != mh(r + 1, c + 1)) {
				return -1;
			}
		}
	}

	return 0;
}

static int
check_vec(VectorHandler& vh, unsigned c)
{
	c--;
	for (int r = 0; r < 5; r++) {
		if (mat[r][c] != vh(r + 1)) {
			return -1;
		}
	}

	return 0;
}

static int
check_vec_transpose(VectorHandler& vh, unsigned r)
{
	r--;
	for (int c = 0; c < 5; c++) {
		if (mat[r][c] != vh(c + 1)) {
			return -1;
		}
	}

	return 0;
}

int
main(void)
{
	SpMapMatrixHandler spm(5, 5);

	for (int r = 0; r < 5; r++) {
		for (int c = 0; c < 5; c++) {
			if (mat[r][c] != 0.) {
				spm(r + 1, c + 1) = mat[r][c];
			}
		}
	}

	std::cout << "matrix in sparse form: " << std::endl
		<< spm << std::endl;

	std::vector<doublereal> Ax0, Ax1;
	std::vector<integer> Ai0, Ai1, Ap0, Ap1;

	spm.MakeCompressedColumnForm(Ax0, Ai0, Ap0, 0);
	spm.MakeCompressedColumnForm(Ax1, Ai1, Ap1, 1);

	CColMatrixHandler<0> ccm0(Ax0, Ai0, Ap0);

	std::cout << "matrix in cc<0> form: " << std::endl
		<< ccm0 << std::endl;
	std::cout << "matrix in cc<0> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << ccm0.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	CColMatrixHandler<1> ccm1(Ax1, Ai1, Ap1);

	std::cout << "matrix in cc<1> form: " << std::endl
		<< ccm1 << std::endl;
	std::cout << "matrix in cc<1> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << ccm1.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	DirCColMatrixHandler<0> dirm0(Ax0, Ai0, Ap0);
	std::cout << "matrix in dir<0> form: " << std::endl
		<< dirm0 << std::endl;
	std::cout << "matrix in dir<0> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << dirm0.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	DirCColMatrixHandler<1> dirm1(Ax1, Ai1, Ap1);
	std::cout << "matrix in dir<1> form: " << std::endl
		<< dirm1 << std::endl;
	std::cout << "matrix in dir<1> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << dirm1.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	MyVectorHandler v(5), out(5);
	v.Reset();
	for (int i = 1; i <= 5; i++) {
		v(i) = 1.;

		out.Reset();
		spm.MatVecMul(out, v);
		std::cout << "sp*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		spm.MatTVecMul(out, v);
		std::cout << "sp^T*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec_transpose(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		ccm0.MatVecMul(out, v);
		std::cout << "cc<0>*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		ccm0.MatTVecMul(out, v);
		std::cout << "cc<0>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec_transpose(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		ccm1.MatVecMul(out, v);
		std::cout << "cc<1>*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		ccm1.MatTVecMul(out, v);
		std::cout << "cc<1>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec_transpose(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		dirm0.MatVecMul(out, v);
		std::cout << "dir<0>*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		dirm0.MatTVecMul(out, v);
		std::cout << "dir<0>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec_transpose(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		dirm1.MatVecMul(out, v);
		std::cout << "dir<1>*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		out.Reset();
		dirm1.MatTVecMul(out, v);
		std::cout << "dir<1>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;
		if (check_vec_transpose(out, i)) {
			std::cerr << "*** failed!" << std::endl;
		}

		v(i) = 0.;
	}

	FullMatrixHandler fmin(5, 5), fmout(5, 5);
	fmin.Reset();

	fmin(1, 1) = 1.;
	fmin(2, 2) = 1.;
	fmin(3, 3) = 1.;
	fmin(4, 4) = 1.;
	fmin(5, 5) = 1.;

	fmout.Reset();
	spm.MatMatMul(&fmout, fmin);
	std::cout << "sp*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	spm.MatTMatMul(&fmout, fmin);
	std::cout << "sp^T*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat_transpose(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	ccm0.MatMatMul(&fmout, fmin);
	std::cout << "cc<0>*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	ccm0.MatTMatMul(&fmout, fmin);
	std::cout << "cc<0>^T*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat_transpose(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	ccm1.MatMatMul(&fmout, fmin);
	std::cout << "cc<1>*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	ccm1.MatTMatMul(&fmout, fmin);
	std::cout << "cc<1>^T*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat_transpose(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}

	fmout.Reset();
	dirm0.MatMatMul(&fmout, fmin);
	std::cout << "dir<0>*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	dirm0.MatTMatMul(&fmout, fmin);
	std::cout << "dir<0>^T*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat_transpose(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	dirm1.MatMatMul(&fmout, fmin);
	std::cout << "dir<1>*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}
	
	fmout.Reset();
	dirm1.MatTMatMul(&fmout, fmin);
	std::cout << "dir<1>^T*eye=" << std::endl
		<< fmout << std::endl;
	if (check_mat_transpose(fmout)) {
		std::cerr << "*** failed!" << std::endl;
	}

	return 0;
}

