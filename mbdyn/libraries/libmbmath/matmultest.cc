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

#include "spmapmh.h"
#include "ccmh.h"
#include "fullmh.h"

int
main(void)
{
	SpMapMatrixHandler spm(5, 5);

	spm(1, 1) = 11.;
	spm(1, 3) = 13.;
	spm(1, 5) = 15.;
	spm(2, 2) = 22.;
	spm(2, 4) = 24.;
	spm(3, 1) = 31.;
	spm(3, 3) = 33.;
	spm(3, 5) = 35.;
	spm(4, 2) = 42.;
	spm(4, 4) = 44.;
	spm(5, 1) = 51.;
	spm(5, 3) = 53.;
	spm(5, 5) = 55.;

	std::cout << "matrix in sparse form: " << std::endl
		<< spm << std::endl;

	std::vector<doublereal> Ax;
	std::vector<integer> Ai;
	std::vector<integer> Ap;

	spm.MakeCompressedColumnForm(Ax, Ai, Ap, 0);

	CColMatrixHandler<0> ccm0(Ax, Ai, Ap);

	std::cout << "matrix in cc<0> form: " << std::endl
		<< ccm0 << std::endl;
	std::cout << "matrix in cc<0> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << ccm0.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	CColMatrixHandler<1> ccm1(Ax, Ai, Ap);

	std::cout << "matrix in cc<1> form: " << std::endl
		<< ccm1 << std::endl;
	std::cout << "matrix in cc<1> form again: " << std::endl;
	for (int ir = 1; ir <= 5; ir++) {
		for (int ic = 1; ic <= 5; ic++) {
			std::cout << std::setw(16) << ccm1.dGetCoef(ir, ic);
		}
		std::cout << std::endl;
	}

	MyVectorHandler v(5), out(5);
	v.Reset();
	out.Reset();
	for (int i = 1; i <= 5; i++) {
		v(i) = 1.;

		spm.MatVecMul(out, v);
		std::cout << "sp*v(" << i << ")=" << std::endl
			<< out << std::endl;

		spm.MatTVecMul(out, v);
		std::cout << "sp^T*v(" << i << ")=" << std::endl
			<< out << std::endl;

		ccm0.MatVecMul(out, v);
		std::cout << "cc<0>*v(" << i << ")=" << std::endl
			<< out << std::endl;

		ccm0.MatTVecMul(out, v);
		std::cout << "cc<0>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;

		ccm1.MatVecMul(out, v);
		std::cout << "cc<1>*v(" << i << ")=" << std::endl
			<< out << std::endl;

		ccm1.MatTVecMul(out, v);
		std::cout << "cc<1>^T*v(" << i << ")=" << std::endl
			<< out << std::endl;

		v(i) = 0.;
	}

	FullMatrixHandler fmin(5, 5), fmout(5, 5);
	fmin.Reset();
	fmout.Reset();

	fmin(1, 1) = 1.;
	fmin(2, 2) = 1.;
	fmin(3, 3) = 1.;
	fmin(4, 4) = 1.;
	fmin(5, 5) = 1.;

	spm.MatMatMul(&fmout, fmin);
	std::cout << "sp*eye=" << std::endl
		<< fmout << std::endl;
	
	spm.MatTMatMul(&fmout, fmin);
	std::cout << "sp^T*eye=" << std::endl
		<< fmout << std::endl;
	
	ccm0.MatMatMul(&fmout, fmin);
	std::cout << "cc<0>*eye=" << std::endl
		<< fmout << std::endl;
	
	ccm0.MatTMatMul(&fmout, fmin);
	std::cout << "cc<0>^T*eye=" << std::endl
		<< fmout << std::endl;
	
	ccm1.MatMatMul(&fmout, fmin);
	std::cout << "cc<1>*eye=" << std::endl
		<< fmout << std::endl;
	
	ccm1.MatTMatMul(&fmout, fmin);
	std::cout << "cc<1>^T*eye=" << std::endl
		<< fmout << std::endl;

	return 0;
}

