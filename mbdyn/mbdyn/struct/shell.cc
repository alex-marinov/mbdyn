/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2010
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
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
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "constltp_impl.h"
#include "tpldrive_impl.h"
#include "shell.h"
#include "mynewmem.h"


static inline void
InsertMatrix(FullMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	start_row--;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= nr+start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer i = 1; i <= nr; i++) {
		for (integer j = 1; j <= nc; j++) {
			dest(i + start_row, j + start_col) = source(i, j);
		}
	}
	return;
}

static inline void
InsertMatrix(FullMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const Mat3x3 & source)
{
	integer nr = 3;
	integer nc = 3;
	start_row--;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= nr+start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer i = 1; i <= nr; i++) {
		for (integer j = 1; j <= nc; j++) {
			dest(i + start_row, j + start_col) = source(i, j);
		}
	}
	return;
}

static inline void
InsertMatrixT(FullMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const Mat3x3 & source)
{
	integer nr = 3;
	integer nc = 3;
	start_row--;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= nr+start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer i = 1; i <= nr; i++) {
		for (integer j = 1; j <= nc; j++) {
			dest(i + start_row, j + start_col) = source(j, i);
		}
	}
	return;
}

static inline void
InsertRowVector(FullMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const Vec3 & source)
{
	integer nc = 3;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer j = 1; j <= nc; j++) {
		dest(start_row, j + start_col) = source(j);
	}
	return;
}

static inline void
CopyMatrixRow(FullMatrixHandler & dest, integer dest_row, 
	const FullMatrixHandler & source, integer source_row)
{
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumCols() == source.iGetNumCols()
	//check da fare: 1 <= dest_row <= dest.iGetNumRows();
	//check da fare: 1 <= source_row <= source.iGetNumRows();
	integer nc = dest.iGetNumCols();
	for (integer i = 1; i <= nc; i++) {
		dest(dest_row, i) = source(source_row, i);
	}
}

static inline void
CopyMatrixBlock(FullMatrixHandler & dest, integer dest_row, integer dest_col,
	const FullMatrixHandler & source, 
	integer source_start_row, integer source_end_row,
	integer source_start_col, integer source_end_col)
{
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() >= dest_row + (source_end_row - source_start_row)
	//check da fare: dest.iGetNumCols() >= dest_col + (source_end_col - source_start_col)
	//check da fare: 1 <= source_start_row <= source.iGetNumRows();
	//check da fare: 1 <= source_start_col <= source.iGetNumRows();
	//check da fare: 1 <= source_end_row <= source.iGetNumRows();
	//check da fare: 1 <= source_end_col <= source.iGetNumRows();
	//check da fare: source_start_row <= source_end_row:
	//check da fare: source_start_col <= source_end_col;
	for (integer i = source_start_row; i <= source_end_row; i++) {
		integer row = dest_row + (i - source_start_row);
		for (integer ii = source_start_col; ii <= source_end_col; ii++) {
			integer col = dest_col + (ii - source_start_col);
			dest(row, col) = source(i, ii);
		}
	}
}

static inline void
InsertVector(MyVectorHandler & dest, 
	integer start_row, 
	const Vec3 & source)
{
	integer nr = 3;
	start_row--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetSize() <= nr+start_row
	for (integer i = 1; i <= nr; i++) {
		dest(i + start_row) = source(i);
	}
	return;
}

static inline void
AssembleVector(SubVectorHandler & dest, 
	integer start_row, 
	const MyVectorHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetSize();
	start_row--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetSize() <= nr+start_row
	for (integer i = 1; i <= nr; i++) {
		dest(i + start_row) += source(i);
	}
	return;
}

static inline void
AssembleMatrix(FullSubMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	start_row--;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= nr+start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer i = 1; i <= nr; i++) {
		for (integer ii = 1; ii <= nc; ii++) {
			dest.IncCoef(i + start_row, ii + start_col, source(i, ii));
		}
	}
	return;
}

static inline void
AssembleTransposeMatrix(FullSubMatrixHandler & dest, 
	integer start_row, 
	integer start_col, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	start_row--;
	start_col--;
	//TODO. error checking: check dimensions!
	//check da fare: dest.iGetNumRows() <= nr+start_row
	//check da fare: dest.iGetNumCols() <= nc+start_col
	for (integer i = 1; i <= nr; i++) {
		for (integer ii = 1; ii <= nc; ii++) {
			dest.IncCoef(i + start_row, ii + start_col, source(ii, i));
		}
	}
	return;
}

static inline void
ExtractVec3(Vec3& dest, const MyVectorHandler & source, integer start_row)
{
	//TODO. error checking: check dimensions!
	//check da fare: source.iGetNumRows() >= start_row + 2
	dest = Vec3(source(start_row), source(start_row + 1), source(start_row + 2));
}

static inline void
ExtractMat3x3(Mat3x3& dest, const FullMatrixHandler & source, 
	integer start_row, integer start_col)
{
	//TODO. error checking: check dimensions!
	//check da fare: source.iGetNumRows() >= start_row + 2
	//check da fare: source.iGetNumCols() >= start_col + 2
	dest = Mat3x3(source(start_row, start_col), source(start_row + 1, start_col), source(start_row + 2, start_col),
		source(start_row, start_col + 1), source(start_row + 1, start_col + 1), source(start_row + 2, start_col + 1),
		source(start_row, start_col + 2), source(start_row + 1, start_col + 2), source(start_row + 2, start_col + 2)
	);
}

static inline void
RotateForward(MyVectorHandler & e, const Mat3x3 R)
{
	//TODO. error checking: check dimensions!
	//check da fare: e.iGetNumRows() == 12
	Vec3 t, t1;	
	for (integer b = 0; b < 4; b++) {
		ExtractVec3(t, e, 1 + b * 3);
		t1 = R * t;
		InsertVector(e, 1 + b * 3, t1);
	}
}

static inline void
RotateBackward(MyVectorHandler & e, const Mat3x3 R)
{
	//TODO. error checking: check dimensions!
	//check da fare: e.iGetNumRows() == 12
	Vec3 t, t1;	
	for (integer b = 0; b < 4; b++) {
		ExtractVec3(t, e, 1 + b * 3);
		t1 = R.MulTV(t);
		InsertVector(e, 1 + b * 3, t1);
	}
}

static inline void
RotateForward(FullMatrixHandler & C, const Mat3x3 R)
{
	//TODO. error checking: check dimensions!
	//check da fare: C.iGetNumRows() == 12
	//check da fare: C.iGetNumCols() == 12
	Mat3x3 m, m1;	
	for (integer rb = 0; rb < 4; rb++) {
		for (integer rc = 0; rc < 4; rc++) {
			ExtractMat3x3(m, C, 1 + rb * 3, 1 + rc * 3);
			m1 = R * m.MulMT(R);
			InsertMatrix(C, 1 + rb * 3, 1 + rc * 3, m1);
		}
	}
}

static inline doublereal L1(const doublereal xi[2]) {
	return 0.25 * (1. + xi[0]) * (1 + xi[1]);
};

static inline doublereal L2(const doublereal xi[2]) {
	return 0.25 * (1. - xi[0]) * (1 + xi[1]);
};

static inline doublereal L3(const doublereal xi[2]) {
	return 0.25 * (1. - xi[0]) * (1 - xi[1]);
};

static inline doublereal L4(const doublereal xi[2]) {
	return 0.25 * (1. + xi[0]) * (1 - xi[1]);
};

typedef doublereal (*LI_Type)(const doublereal xi[2]);
LI_Type LI[4] = {&L1, &L2, &L3, &L4};



static inline doublereal
L1_1(const doublereal xi[2])
{
	return 0.25 * (1. + xi[1]);
}

static inline doublereal
L1_2(const doublereal xi[2])
{
	return 0.25 * (1. + xi[0]);
}

static inline doublereal
L2_1(const doublereal xi[2])
{
	return -0.25 * (1. + xi[1]);
}

static inline doublereal
L2_2(const doublereal xi[2])
{
	return 0.25 * (1. - xi[0]);
}

static inline doublereal
L3_1(const doublereal xi[2])
{
	return -0.25 * (1. - xi[1]);
}

static inline doublereal
L3_2(const doublereal xi[2])
{
	return -0.25 * (1. - xi[0]);
}

static inline doublereal
L4_1(const doublereal xi[2])
{
	return 0.25 * (1. - xi[1]);
}

static inline doublereal
L4_2(const doublereal xi[2])
{
	return -0.25 * (1. + xi[0]);
}

typedef doublereal (*LI_J_Type)(const doublereal xi[2]);
LI_J_Type LI_J[4][2] = {
	{&L1_1, &L1_2},
	{&L2_1, &L2_2},
	{&L3_1, &L3_2},
	{&L4_1, &L4_2},	
};

static Vec3
Interp(const Vec3*const v, const doublereal xi[2])
{
	Vec3 r = v[0] * L1(xi) + 
		v[1] * L2(xi) + 
		v[2] * L3(xi) + 
		v[3] * L4(xi); 
	return r;
}

static Vec3
InterpDeriv1(const Vec3*const v, const FullMatrixHandler & der_mat)
{
	Vec3 r = v[0] * der_mat(1, 1) + 
		v[1] * der_mat(2, 1) + 
		v[2] * der_mat(3, 1) + 
		v[3] * der_mat(4, 1); 
	return r;
}

static Vec3
InterpDeriv2(const Vec3*const v, const FullMatrixHandler & der_mat)
{
	Vec3 r = v[0] * der_mat(1, 2) + 
		v[1] * der_mat(2, 2) + 
		v[2] * der_mat(3, 2) + 
		v[3] * der_mat(4, 2); 
	return r;
}

static void
InterpDeriv(const Vec3*const v, 
	const FullMatrixHandler & der_mat,
	Vec3 & der1, 
	Vec3 & der2)
{
	der1 = InterpDeriv1(v, der_mat);
	der2 = InterpDeriv2(v, der_mat);
	return;
}

static Vec3
InterpDeriv_xi1(const Vec3*const v, const doublereal xi[2])
{
	Vec3 r = v[0] * L1_1(xi) + 
		v[1] * L2_1(xi) + 
		v[2] * L3_1(xi) + 
		v[3] * L4_1(xi); 
	return r;
}

static Vec3
InterpDeriv_xi2(const Vec3*const v, const doublereal xi[2])
{
	Vec3 r = v[0] * L1_2(xi) + 
		v[1] * L2_2(xi) + 
		v[2] * L3_2(xi) + 
		v[3] * L4_2(xi); 
	return r;
}

static void
Inv2x2(const FullMatrixHandler& a, FullMatrixHandler & out)
{
	//FIXME: Mettere controlli dimensioni matrici In/Out
	//FIXME: mettere controllo det \neq 0
	doublereal det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1);

	out(1, 1) =  a(2, 2) / det;
	out(1, 2) = -a(2, 1) / det;
	out(2, 1) = -a(1, 2) / det;
	out(2, 2) =  a(1, 1) / det;

	return;
}

static void
Inv4x4(const FullMatrixHandler& a, FullMatrixHandler & out)
{
	//FIXME: Mettere controlli dimensioni matrici In/Out
	//FIXME: mettere controllo det \neq 0
	doublereal det =
		a(1, 1)*a(2, 2)*a(3, 3)*a(4, 4)-a(1, 1)*a(2, 2)*a(3, 4)*a(4, 3)-a(1, 1)*a(3, 2)*a(2, 3)*a(4, 4)+a(1, 1)*a(3, 2)*a(2, 4)*a(4, 3)+a(1, 1)*a(4, 2)*a(2, 3)*a(3, 4)-a(1, 1)*a(4, 2)*a(2, 4)*a(3, 3)-a(2, 1)*a(1, 2)*a(3, 3)*a(4, 4)+a(2, 1)*a(1, 2)*a(3, 4)*a(4, 3)+a(2, 1)*a(3, 2)*a(1, 3)*a(4, 4)-a(2, 1)*a(3, 2)*a(1, 4)*a(4, 3)-a(2, 1)*a(4, 2)*a(1, 3)*a(3, 4)+a(2, 1)*a(4, 2)*a(1, 4)*a(3, 3)+a(3, 1)*a(1, 2)*a(2, 3)*a(4, 4)-a(3, 1)*a(1, 2)*a(2, 4)*a(4, 3)-a(3, 1)*a(2, 2)*a(1, 3)*a(4, 4)+a(3, 1)*a(2, 2)*a(1, 4)*a(4, 3)+a(3, 1)*a(4, 2)*a(1, 3)*a(2, 4)-a(3, 1)*a(4, 2)*a(1, 4)*a(2, 3)-a(4, 1)*a(1, 2)*a(2, 3)*a(3, 4)+a(4, 1)*a(1, 2)*a(2, 4)*a(3, 3)+a(4, 1)*a(2, 2)*a(1, 3)*a(3, 4)-a(4, 1)*a(2, 2)*a(1, 4)*a(3, 3)-a(4, 1)*a(3, 2)*a(1, 3)*a(2, 4)+a(4, 1)*a(3, 2)*a(1, 4)*a(2, 3);

	out(1, 1) =  (a(2, 2)*a(3, 3)*a(4, 4)-a(2, 2)*a(3, 4)*a(4, 3)-a(3, 2)*a(2, 3)*a(4, 4)+a(3, 2)*a(2, 4)*a(4, 3)+a(4, 2)*a(2, 3)*a(3, 4)-a(4, 2)*a(2, 4)*a(3, 3))/det;
	out(1, 2) = -(a(1, 2)*a(3, 3)*a(4, 4)-a(1, 2)*a(3, 4)*a(4, 3)-a(3, 2)*a(1, 3)*a(4, 4)+a(3, 2)*a(1, 4)*a(4, 3)+a(4, 2)*a(1, 3)*a(3, 4)-a(4, 2)*a(1, 4)*a(3, 3))/det;
	out(1, 3) =  (a(1, 2)*a(2, 3)*a(4, 4)-a(1, 2)*a(2, 4)*a(4, 3)-a(2, 2)*a(1, 3)*a(4, 4)+a(2, 2)*a(1, 4)*a(4, 3)+a(4, 2)*a(1, 3)*a(2, 4)-a(4, 2)*a(1, 4)*a(2, 3))/det;
	out(1, 4) = -(a(1, 2)*a(2, 3)*a(3, 4)-a(1, 2)*a(2, 4)*a(3, 3)-a(2, 2)*a(1, 3)*a(3, 4)+a(2, 2)*a(1, 4)*a(3, 3)+a(3, 2)*a(1, 3)*a(2, 4)-a(3, 2)*a(1, 4)*a(2, 3))/det;

	out(2, 1) = -(a(2, 1)*a(3, 3)*a(4, 4)-a(2, 1)*a(3, 4)*a(4, 3)-a(3, 1)*a(2, 3)*a(4, 4)+a(3, 1)*a(2, 4)*a(4, 3)+a(4, 1)*a(2, 3)*a(3, 4)-a(4, 1)*a(2, 4)*a(3, 3))/det;
	out(2, 2) =  (a(1, 1)*a(3, 3)*a(4, 4)-a(1, 1)*a(3, 4)*a(4, 3)-a(3, 1)*a(1, 3)*a(4, 4)+a(3, 1)*a(1, 4)*a(4, 3)+a(4, 1)*a(1, 3)*a(3, 4)-a(4, 1)*a(1, 4)*a(3, 3))/det;
	out(2, 3) = -(a(1, 1)*a(2, 3)*a(4, 4)-a(1, 1)*a(2, 4)*a(4, 3)-a(2, 1)*a(1, 3)*a(4, 4)+a(2, 1)*a(1, 4)*a(4, 3)+a(4, 1)*a(1, 3)*a(2, 4)-a(4, 1)*a(1, 4)*a(2, 3))/det;
	out(2, 4) =  (a(1, 1)*a(2, 3)*a(3, 4)-a(1, 1)*a(2, 4)*a(3, 3)-a(2, 1)*a(1, 3)*a(3, 4)+a(2, 1)*a(1, 4)*a(3, 3)+a(3, 1)*a(1, 3)*a(2, 4)-a(3, 1)*a(1, 4)*a(2, 3))/det;

	out(3, 1) =  (a(2, 1)*a(3, 2)*a(4, 4)-a(2, 1)*a(3, 4)*a(4, 2)-a(3, 1)*a(2, 2)*a(4, 4)+a(3, 1)*a(2, 4)*a(4, 2)+a(4, 1)*a(2, 2)*a(3, 4)-a(4, 1)*a(2, 4)*a(3, 2))/det;
	out(3, 2) = -(a(1, 1)*a(3, 2)*a(4, 4)-a(1, 1)*a(3, 4)*a(4, 2)-a(3, 1)*a(1, 2)*a(4, 4)+a(3, 1)*a(1, 4)*a(4, 2)+a(4, 1)*a(1, 2)*a(3, 4)-a(4, 1)*a(1, 4)*a(3, 2))/det;
	out(3, 3) =  (a(1, 1)*a(2, 2)*a(4, 4)-a(1, 1)*a(2, 4)*a(4, 2)-a(2, 1)*a(1, 2)*a(4, 4)+a(2, 1)*a(1, 4)*a(4, 2)+a(4, 1)*a(1, 2)*a(2, 4)-a(4, 1)*a(1, 4)*a(2, 2))/det;
	out(3, 4) = -(a(1, 1)*a(2, 2)*a(3, 4)-a(1, 1)*a(2, 4)*a(3, 2)-a(2, 1)*a(1, 2)*a(3, 4)+a(2, 1)*a(1, 4)*a(3, 2)+a(3, 1)*a(1, 2)*a(2, 4)-a(3, 1)*a(1, 4)*a(2, 2))/det;

	out(4, 1) = -(a(2, 1)*a(3, 2)*a(4, 3)-a(2, 1)*a(3, 3)*a(4, 2)-a(3, 1)*a(2, 2)*a(4, 3)+a(3, 1)*a(2, 3)*a(4, 2)+a(4, 1)*a(2, 2)*a(3, 3)-a(4, 1)*a(2, 3)*a(3, 2))/det;
	out(4, 2) =  (a(1, 1)*a(3, 2)*a(4, 3)-a(1, 1)*a(3, 3)*a(4, 2)-a(3, 1)*a(1, 2)*a(4, 3)+a(3, 1)*a(1, 3)*a(4, 2)+a(4, 1)*a(1, 2)*a(3, 3)-a(4, 1)*a(1, 3)*a(3, 2))/det;
	out(4, 3) = -(a(1, 1)*a(2, 2)*a(4, 3)-a(1, 1)*a(2, 3)*a(4, 2)-a(2, 1)*a(1, 2)*a(4, 3)+a(2, 1)*a(1, 3)*a(4, 2)+a(4, 1)*a(1, 2)*a(2, 3)-a(4, 1)*a(1, 3)*a(2, 2))/det;
	out(4, 4) =  (a(1, 1)*a(2, 2)*a(3, 3)-a(1, 1)*a(2, 3)*a(3, 2)-a(2, 1)*a(1, 2)*a(3, 3)+a(2, 1)*a(1, 3)*a(3, 2)+a(3, 1)*a(1, 2)*a(2, 3)-a(3, 1)*a(1, 3)*a(2, 2))/det;

	return;
}

Shell::Shell(unsigned uLabel, const DofOwner* pDO, flag fOut)
: Elem(uLabel, fOut),
ElemGravityOwner(uLabel, fOut),
ElemWithDofs(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut)
{
	NO_OP;
}

Shell::~Shell(void)
{
	NO_OP;
}


/***************************************
* 
*
****************************************/
 
doublereal Shell4::xi_i[Shell4::NUMIP][2] = {
	{-1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.), -1. / std::sqrt(3.)},
	{ 1. / std::sqrt(3.),  1. / std::sqrt(3.)},
	{-1. / std::sqrt(3.),  1. / std::sqrt(3.)}
};

doublereal Shell4::xi_A[Shell4::NUMSSEP][2] =  {
	{ 0.,  1.},
	{-1.,  0.},
	{ 0., -1.},
	{ 1.,  0.}
};

doublereal Shell4::xi_n[Shell4::NUMNODES][2] = {
	{-1., -1.},
	{ 1., -1.},
	{ 1.,  1.},
	{-1.,  1.}
};

doublereal Shell4::xi_0[2] = {0., 0.};

void
Shell4::UpdateNodalAndAveragePosAndOrientation(void)
{
	Mat3x3 T_a_avg(0.);
	for (integer i = 0; i < NUMNODES; i++) {
// 		xa[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
		xa[i] = pNode[i]->GetXCurr();
		Ta[i] = pNode[i]->GetRCurr();
		T_a_avg += Ta[i];
	}
	T_a_avg /= 4.;
	//FIXME; alternative solution: polar decomposition of T_a_avg = T0*U
	T_overline = RotManip::Rot(RotManip::VecRot(T_a_avg));
	for (integer i = 0; i < NUMNODES; i++) {
		Ra[i] = T_overline.MulTM(Ta[i]);
		phia[i] = RotManip::VecRot(Ra[i]);
	}
}
	
void
Shell4::InterpolateOrientation(void)
{
	for (integer i = 0; i < NUMIP; i++) {
		phi_tilde_i[i] = Interp(phia, xi_i[i]);
		T_i[i] = T_overline * RotManip::Rot(phi_tilde_i[i]);
	}
	for (integer i = 0; i < NUMSSEP; i++) {
		phi_tilde_A[i] = Interp(phia, xi_A[i]);
		T_A[i] = T_overline * RotManip::Rot(phi_tilde_A[i]);
	}
	phi_tilde_0 = Interp(phia, xi_0);
	T_0 = T_overline * RotManip::Rot(phi_tilde_0);
}

void
Shell4::ComputeIPSEPRotations(void)
{
	for (integer i = 0; i < NUMIP; i++) {
		Q_i[i] = T_i[i].MulMT(T_0_i[i]);
	}
	for (integer i = 0; i < NUMSSEP; i++) {
		Q_A[i] = T_A[i].MulMT(T_0_A[i]);
	}
}

void
Shell4::ComputeIPCurvature(void)
{
	for (integer i = 0; i < NUMIP; i++) {
		InterpDeriv(phia, L_alpha_beta_i[i], phi_tilde_1_i[i], phi_tilde_2_i[i]);
		AxT_1_i[i] = T_overline * RotManip::DRot(phi_tilde_i[i]) * phi_tilde_1_i[i];
		AxT_2_i[i] = T_overline * RotManip::DRot(phi_tilde_i[i]) * phi_tilde_2_i[i];
		T_1_i[i] = AxT_1_i[i] *  T_overline.MulMT(T_i[i]);
		T_2_i[i] = AxT_2_i[i] *  T_overline.MulMT(T_i[i]);
	}
}

Shell4::Shell4(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pN1, const StructNode* pN2,
	const StructNode* pN3, const StructNode* pN4,
#if 0 // TODO: offset
 	const Vec3& f1, const Vec3& f2,
 	const Vec3& f3, const Vec3& f4,
#endif
	const Mat3x3& R1, const Mat3x3& R2,
	const Mat3x3& R3, const Mat3x3& R4,
#ifdef USE_CL_IN_SHELL
	const ConstitutiveLaw<vh, fmh>** pDTmp, 
#else // ! USE_CL_IN_SHELL
	const fmh& pDTmp,
#endif // ! USE_CL_IN_SHELL
	flag fOut)
: 
Elem(uL, fOut), 
Shell(uL, pDO, fOut),
S_alpha_beta_0(2, 2),
S_alpha_beta_i(NUMIP, fmh(2, 2) ),
S_alpha_beta_A(NUMSSEP, fmh(2, 2) ),
// doublereal alpha_0,
// doublereal alpha_i[NUMIP],
// doublereal alpha_A[NUMSSEP],
L_alpha_B_i(NUMIP, fmh(4, 2) ),
L_alpha_B_A(NUMSSEP, fmh(4, 2) ),
L_alpha_beta_i(NUMIP, fmh(4, 2) ),
L_alpha_beta_A(NUMSSEP, fmh(4, 2) ),

B_overline_i(NUMIP, fmh(12, 24) ),
B_overline_m_i(NUMIP, fmh(4, 24) ),
B_overline_S_i(NUMIP, fmh(2, 24) ),
D_overline_i(NUMIP, fmh(15, 24) ),
G_i(NUMIP, fmh(15, 15) ),

K_M_i(NUMIP, fmh(24, 24) ),
K_G_i(NUMIP, fmh(24, 24) ),
K_alpha_beta_i(NUMIP, fmh(8, 24) ),

p_i(NUMIP, vh(24) ),
r_d(NUMIP, vh(24) ),
r_beta(NUMIP, vh(8) ),

K_T_i(NUMIP, fmh(24, 24) ),
r_i(NUMIP, vh(24) ),

// Vec3 phi_tilde_1_i[NUMIP],
// Vec3 phi_tilde_2_i[NUMIP],
// Mat3x3 T0_1_i[NUMIP],
// Mat3x3 T0_2_i[NUMIP],
M_0(4, 4),
M_0_Inv(4, 4),
M_2_i(NUMIP, fmh(2, 2) ),
M_2_A(NUMSSEP, fmh(2, 2) ),
H_i(NUMIP, fmh(4, 8) ),
H_A(NUMSSEP, fmh(4, 8) ),

P_i(NUMIP, fmh(4, 8) ),
P_A(NUMSSEP, fmh(4, 8) ),

K_beta_beta_i(NUMIP, fmh(8, 8) ),

S_i(NUMIP, vh(12) ),

L_i(NUMIP, fmh(6, 24)),

beta(8),
epsilon_hat(4),
epsilon(12),

DRef(NUMIP, fmh(12, 12)),
stress_i(NUMIP, vh(12))
{
#ifdef USE_CL_IN_SHELL
	for (integer i = 0; i < NUMIP; i++) {
		pD[i] = 0;
		SAFENEWWITHCONSTRUCTOR(pD[i],
			ConstitutiveLawOwnerType,
			ConstitutiveLawOwnerType(pDTmp[i]));	
	}
#else // ! USE_CL_IN_SHELL
	for (unsigned i = 0; i < NUMIP; i++) {
		DRef[i] = pDTmp;
	}
#endif // ! USE_CL_IN_SHELL

	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	pNode[NODE3] = pN3;
	pNode[NODE4] = pN4;
#warning "TODO: orientazione iniziale con normale messa a posto"
// 	f[NODE1] = f1;
// 	f[NODE2] = f2;
// 	f[NODE3] = f3;
// 	f[NODE4] = f4;

	UpdateNodalAndAveragePosAndOrientation();
	InterpolateOrientation();
	// copy ref values
	T0_overline = T_overline;
	T_0_0 = T_0;
	for (integer i = 0; i < NUMIP; i++) {
		xa_0[i] = xa[i];
		T_0_i[i] = T_i[i];
	}
	for (integer i = 0; i < NUMSSEP; i++) {
		T_0_A[i] = T_A[i];
	}

	for (integer i = 0; i < NUMIP; i++) {
		phi_tilde_1_i[i] = Vec3(0.);
		phi_tilde_2_i[i] = Vec3(0.);
	}


	{
		Vec3 x_1 = InterpDeriv_xi1(xa, xi_0);
		Vec3 x_2 = InterpDeriv_xi2(xa, xi_0);
		S_alpha_beta_0(1, 1) = T_0_0.GetCol(1) * x_1;
		S_alpha_beta_0(2, 1) = T_0_0.GetCol(2) * x_1;
		S_alpha_beta_0(1, 2) = T_0_0.GetCol(1) * x_2;
		S_alpha_beta_0(2, 2) = T_0_0.GetCol(2) * x_2;
		alpha_0 = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2) -
			S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);

		M_0(1, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(1, 1);
		M_0(1, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(1, 2);
		M_0(1, 3) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(1, 1);
		M_0(1, 4) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(1, 2);	

		M_0(2, 1) = S_alpha_beta_0(2, 1) * S_alpha_beta_0(2, 1);
		M_0(2, 2) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 2);
		M_0(2, 3) = S_alpha_beta_0(2, 2) * S_alpha_beta_0(2, 1);
		M_0(2, 4) = S_alpha_beta_0(2, 1) * S_alpha_beta_0(2, 2);
	
		M_0(3, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		M_0(3, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		M_0(3, 3) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2);
		M_0(3, 4) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);

		M_0(4, 1) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 1);
		M_0(4, 2) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 2);
		M_0(4, 3) = S_alpha_beta_0(1, 2) * S_alpha_beta_0(2, 1);
		M_0(4, 4) = S_alpha_beta_0(1, 1) * S_alpha_beta_0(2, 2);
		
		Inv4x4(M_0, M_0_Inv);
	}

	for (integer i = 0; i < NUMIP; i++) {
		Vec3 x_1 = InterpDeriv_xi1(xa, xi_i[i]);
		Vec3 x_2 = InterpDeriv_xi2(xa, xi_i[i]);
		S_alpha_beta_i[i](1, 1) = T_0_i[i].GetCol(1) * x_1;
		S_alpha_beta_i[i](2, 1) = T_0_i[i].GetCol(2) * x_1;
		S_alpha_beta_i[i](1, 2) = T_0_i[i].GetCol(1) * x_2;
		S_alpha_beta_i[i](2, 2) = T_0_i[i].GetCol(2) * x_2;
		//alpha_i = det(S_alpha_beta_i)
		alpha_i[i] = S_alpha_beta_i[i](1, 1) * S_alpha_beta_i[i](2, 2) -
			S_alpha_beta_i[i](1, 2) * S_alpha_beta_i[i](2, 1);

		//xi_i_i = S_alpha_beta_i^{-1}
		FullMatrixHandler xi_i_i(2, 2);
		Inv2x2(S_alpha_beta_i[i], xi_i_i);
		
		for (integer n = 0; n < NUMNODES; n++) {
			for (integer ii = 0; ii < 2; ii++) {
				L_alpha_B_i[i](n + 1, ii + 1) = LI_J[n][ii](xi_i[i]);
			}
		}
		
		L_alpha_B_i[i].MatMatMul(L_alpha_beta_i[i], xi_i_i); 

		//calcola derivate di phia
		InterpDeriv(phia, L_alpha_beta_i[i], phi_tilde_1_i[i], phi_tilde_2_i[i]);

		Mat3x3 R_1_i = 
			Mat3x3(RotManip::DRot(phi_tilde_i[i]) * phi_tilde_1_i[i]) * 
			RotManip::Rot(phi_tilde_i[i]);
		Mat3x3 R_2_i = 
			Mat3x3(RotManip::DRot(phi_tilde_i[i]) * phi_tilde_2_i[i]) * 
			RotManip::Rot(phi_tilde_i[i]);
		T0_1_i[i] = T0_overline * R_1_i;
		T0_2_i[i] = T0_overline * R_2_i;
		
		doublereal t = xi_i[i][0] * xi_i[i][1];
		H_i[i](1, 1) = xi_i[i][0];
		H_i[i](1, 2) = t;

		H_i[i](2, 3) = xi_i[i][1];
		H_i[i](2, 4) = t;
		
		H_i[i](3, 5) = xi_i[i][0];
		H_i[i](3, 6) = t;
		
		H_i[i](4, 7) = xi_i[i][1];
		H_i[i](4, 8) = t;
		
		M_0_Inv.MatTMatMul(P_i[i], H_i[i]);
		P_i[i].ScalarMul(alpha_0 / alpha_i[i]);
	}
	//save initial axial values
	ComputeIPCurvature();
	for (integer i = 0; i < NUMIP; i++) {
		AxT_1_0_i[i] = AxT_1_i[i];
		AxT_2_0_i[i] = AxT_2_i[i];
	}
	for (integer i = 0; i < NUMSSEP; i++) {
		Vec3 x_1 = InterpDeriv_xi1(xa, xi_A[i]);
		Vec3 x_2 = InterpDeriv_xi2(xa, xi_A[i]);
		S_alpha_beta_A[i](1, 1) = T_0_A[i].GetCol(1) * x_1;
		S_alpha_beta_A[i](2, 1) = T_0_A[i].GetCol(2) * x_1;
		S_alpha_beta_A[i](1, 2) = T_0_A[i].GetCol(1) * x_2;
		S_alpha_beta_A[i](2, 2) = T_0_A[i].GetCol(2) * x_2;
		alpha_A[i] = S_alpha_beta_A[i](1, 1) * S_alpha_beta_A[i](2, 2) -
			S_alpha_beta_A[i](1, 2) * S_alpha_beta_A[i](2, 1);

		//xi_A_i = S_alpha_beta_A^{-1}
		FullMatrixHandler xi_A_i(2, 2);
		Inv2x2(S_alpha_beta_A[i], xi_A_i);

		for (integer n = 0; n < NUMNODES; n++) {
			for (integer ii = 0; ii < 2; ii++) {
				L_alpha_B_A[i](n + 1, ii + 1) = LI_J[n][ii](xi_i[i]);
			}
		}
		
		L_alpha_B_A[i].MatMatMul(L_alpha_beta_A[i], xi_A_i); 

		doublereal t = xi_A[i][0] * xi_A[i][1];
		H_A[i](1, 1) = xi_A[i][0];
		H_A[i](1, 2) = t;

		H_A[i](2, 3) = xi_A[i][1];
		H_A[i](2, 4) = t;
		
		H_A[i](3, 5) = xi_A[i][0];
		H_A[i](3, 6) = t;
		
		H_A[i](4, 7) = xi_A[i][1];
		H_A[i](4, 8) = t;	

		M_0_Inv.MatTMatMul(P_A[i], H_A[i]);
		P_A[i].ScalarMul(alpha_0 / alpha_A[i]);
	}
	//parte costante di L_i
	for (integer i = 0; i < NUMIP; i++) {
		for (integer n = 0; n < NUMNODES; n++) {
			InsertMatrix(L_i[i], 1,  1 + 6 * n, Mat3x3(LI[n](xi_i[i])));
		}
	}
	//parte costante di B_overline_i e D_overline_i
	for (integer i = 0; i < NUMIP; i++) {
		for (integer n = 0; n < NUMNODES; n++) {
			InsertMatrix(B_overline_i[i], 1, 1 + n * 6, Mat3x3(L_alpha_beta_i[i](n + 1, 1)));
			InsertMatrix(B_overline_i[i], 4, 1 + n * 6, Mat3x3(L_alpha_beta_i[i](n + 1, 2)));

			InsertMatrix(D_overline_i[i], 1, 1 + n * 6, Mat3x3(L_alpha_beta_i[i](n + 1, 1)));
			InsertMatrix(D_overline_i[i], 4, 1 + n * 6, Mat3x3(L_alpha_beta_i[i](n + 1, 2)));
		}
	}
}

Shell4::~Shell4(void)
{
#ifdef USE_CL_IN_SHELL
	for (integer i = 0; i < NUMIP; i++) {
		ASSERT(pD[i] != 0);
		if (pD[i] != 0) {
			SAFEDELETE(pD[i]);
		}
	}
#endif // USE_CL_IN_SHELL
}

SubVectorHandler& Shell4::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera e setta gli indici */
	for (integer i = 0; i < 4; i++) {
		integer iNodeFirstMomIndex = pNode[i]->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WorkVec.PutRowIndex(iCnt + 6 * i, iNodeFirstMomIndex + iCnt);
		}
	}

	integer iFirstReactionIndex = iGetFirstIndex();
	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
	}

	UpdateNodalAndAveragePosAndOrientation();
	Vec3 u[NUMNODES];
	for (integer i = 0; i < NUMNODES; i++) {
		u[i] = xa[i] - xa_0[i];
	}
	InterpolateOrientation();
	ComputeIPSEPRotations();
	for (integer i = 1; i <= 8; i++) {
		beta(i) = XCurr(iFirstReactionIndex + i);
	}

	// parte variabile di L_i
	for (integer i = 0; i < NUMIP; i++) {
		for (integer n = 0; n < NUMNODES; n++) {
			InsertMatrixT(L_i[i], 4,  4 + 6 * n, T_i[i]*LI[n](xi_i[i]));
		}
	}
	ComputeIPCurvature();
	for (integer i = 0; i < NUMIP; i++) {
		Vec3 u_der[2];
		InterpDeriv(u, L_alpha_beta_i[i], u_der[0], u_der[1]);
		eps_tilde_1_i[i] = u_der[1] + (Mat3x3(1.) - Q_i[i]) * T_0_i[i].GetVec(1);
		eps_tilde_2_i[i] = u_der[2] + (Mat3x3(1.) - Q_i[i]) * T_0_i[i].GetVec(2);
		k_tilde_1_i[i] = AxT_1_i[i] - Q_i[i] * AxT_1_0_i[i];
		k_tilde_2_i[i] = AxT_2_i[i] - Q_i[i] * AxT_2_0_i[i];
	}
	for (integer i = 0; i < NUMSSEP; i++) {
		Vec3 u_der[2];
		InterpDeriv(u, L_alpha_beta_A[i], u_der[0], u_der[1]);
		eps_tilde_1_A[i] = u_der[1] + (Mat3x3(1.) - Q_A[i]) * T_0_A[i].GetVec(1);
		eps_tilde_2_A[i] = u_der[2] + (Mat3x3(1.) - Q_A[i]) * T_0_A[i].GetVec(2);
	}
	// parte variabile di B_overline_i
	for (integer i = 0; i < NUMIP; i++) {
		for (integer n = 0; n < NUMNODES; n++) {
			InsertMatrix(B_overline_i[i],  1, 4 + 6 * n, 
				Mat3x3(T_i[i].GetVec(1) + eps_tilde_1_i[i]).MulMT(T_i[i]) * LI[n](xi_i[i]));
			InsertMatrix(B_overline_i[i],  4, 4 + 6 * n, 
				Mat3x3(T_i[i].GetVec(2) + eps_tilde_2_i[i]).MulMT(T_i[i]) * LI[n](xi_i[i]));
			InsertMatrix(B_overline_i[i],  7, 4 + 6 * n,
				T_1_i[i].Transpose() * LI[n](xi_i[i]) + 
					T_i[i].Transpose() * L_alpha_beta_i[i](n + 1, 1));
			InsertMatrix(B_overline_i[i], 10, 4 + 6 * n,
				T_2_i[i].Transpose() * LI[n](xi_i[i]) + 
					T_i[i].Transpose() * L_alpha_beta_i[i](n + 1, 2));
		}
	}
	{
		integer tmpidx1[5] = {0, 1, 4, 2, 5};
		for (integer i = 0; i < NUMIP; i++) {
			for (integer n = 1; n <=4; n++) {
				CopyMatrixRow(B_overline_m_i[i], n, B_overline_i[i], tmpidx1[n]);
			}
		}
	}
	vfmh B_tilde_S_A(NUMSSEP, fmh(2, 24));
	for (integer i = 0; i < NUMSSEP; i++) {
		Vec3 b1 = (T_A[i].MulMT(Mat3x3(T_A[i].GetVec(1) + eps_tilde_1_A[i]))).GetVec(3);
		Vec3 b2 = (T_A[i].MulMT(Mat3x3(T_A[i].GetVec(2) + eps_tilde_2_A[i]))).GetVec(3);
		
		for (integer n = 0; n < NUMNODES; n++) {
			B_tilde_S_A[i](1, n * 6 + 3) = L_alpha_beta_A[i](n + 1, 1);
			B_tilde_S_A[i](2, n * 6 + 3) = L_alpha_beta_A[i](n + 1, 2);
		}
		for (integer n = 0; n < NUMNODES; n++) {
			InsertRowVector(B_tilde_S_A[i], 1,  4 + 6 * n, b1 * LI[n](xi_A[i]));
			InsertRowVector(B_tilde_S_A[i], 2,  4 + 6 * n, b2 * LI[n](xi_A[i]));
		}
	}
	{
		FullMatrixHandler B_arrow_ABCD(4, 24);
		
		FullMatrixHandler tempmat(2, 24);
		
		integer tmpidx1[4] = {0, 2, 3, 1};
		integer tmpidx2[4] = {1, 1, 2, 2};
		for (integer i = 0; i < NUMSSEP; i++) {
				S_alpha_beta_A[i].MatTMatMul(tempmat, B_tilde_S_A[tmpidx1[i]]);
			CopyMatrixRow(B_arrow_ABCD, i + 1, tempmat, tmpidx2[i]);
		}

		B_arrow_ABCD.ScalarMul(0.5);
		
		
		for (integer i = 0; i < NUMIP; i++) {
			FullMatrixHandler sh(2, 4);
				sh(1, 1) = 1. + xi_i[i][1];
				sh(1, 2) = 1. - xi_i[i][1];
				sh(2, 3) = 1. + xi_i[i][0];
				sh(2, 4) = 1. - xi_i[i][0];
			FullMatrixHandler B_cup_S_i(2, 24);
				sh.MatMatMul(B_cup_S_i, B_arrow_ABCD);
			FullMatrixHandler Sabi(2, 2); 
				Inv2x2(S_alpha_beta_i[i], Sabi);
			Sabi.MatTMatMul(B_overline_S_i[i], B_cup_S_i);

			CopyMatrixRow(B_overline_i[i], 3, B_overline_S_i[i], 1);
			CopyMatrixRow(B_overline_i[i], 6, B_overline_S_i[i], 2);
		}
	}
	for (integer i = 0; i < NUMIP; i++) {
		for (integer n = 0; n < NUMNODES; n++) {
			CopyMatrixBlock(D_overline_i[i], 7, 4 + 6 * n, 
				B_overline_i[i],
				7, 12, 4 + 6 * n, 6 + 6 * n);

			Mat3x3 t;
			t = T_i[i].Transpose()* LI[n](xi_i[i]);
			InsertMatrix(D_overline_i[i], 13, 4 + 6 * n, t);
		}	
	}
	/* Calcola le azioni interne */
	for (integer i = 0; i < NUMIP; i++) {
		InsertVector(epsilon,  1, eps_tilde_1_i[i]);
		InsertVector(epsilon,  4, eps_tilde_2_i[i]);
		InsertVector(epsilon,  7, k_tilde_1_i[i]);
		InsertVector(epsilon, 10, k_tilde_1_i[i]);
		// TODO: recupera epsilon_hat con l'ordine giusto per qua
		P_i[i].MatVecMul(epsilon_hat, beta);
		// FIXME: set indexes as appropriate
		epsilon(1) += epsilon_hat(1);
		epsilon(2) += epsilon_hat(2);
		epsilon(3) += epsilon_hat(3);
		epsilon(4) += epsilon_hat(4);
		RotateBackward(epsilon, Q_i[i]);
#ifdef USE_CL_IN_SHELL
		pD[i]->Update(epsilon);
		stress_i[i] = pD[i]->GetF();
#else // ! USE_CL_IN_SHELL
		DRef[i].MatVecMul(stress_i[i], epsilon);
#endif // ! USE_CL_IN_SHELL
		RotateForward(stress_i[i], Q_i[i]);
		
		Vec3 n1, n2, m1, m2;
		ExtractVec3(n1, stress_i[i],  1);
		ExtractVec3(n2, stress_i[i],  4);
		ExtractVec3(m1, stress_i[i],  7);
		ExtractVec3(m2, stress_i[i], 10);
		
		Mat3x3 Hh;
		Vec3 t1 = T_i[i].GetCol(1) + eps_tilde_1_i[i];
		Vec3 t2 = T_i[i].GetCol(2) + eps_tilde_2_i[i];
		Hh = n1.Tens(t1) - Mat3x3(n1.Dot(t1)) + n2.Tens(t2) - Mat3x3(n2.Dot(t2));
		
		InsertMatrix(G_i[i],  1, 13, Mat3x3(-n1));
		InsertMatrix(G_i[i],  4, 13, Mat3x3(-n2));
		InsertMatrix(G_i[i],  7, 13, Mat3x3(-m1));
		InsertMatrix(G_i[i], 10, 13, Mat3x3(-m2));
		InsertMatrix(G_i[i], 13,  1, Mat3x3(n1));
		InsertMatrix(G_i[i], 13,  4, Mat3x3(n2));
		InsertMatrix(G_i[i], 13, 13, Hh);
	}
	
	//Residuo
	//forze di volume
	MyVectorHandler pe(24);
	MyVectorHandler rd(24);
	MyVectorHandler membranal(4);
	MyVectorHandler rbeta(8);
	for (integer i = 0; i < NUMIP; i++) {
		//L_i[i].MatTMul(pe, p); //p: forze di volume
		B_overline_i[i].MatTVecMul(rd, stress_i[i]);
		membranal(1) = stress_i[i](1);
		membranal(2) = stress_i[i](4);
		membranal(3) = stress_i[i](2);
		membranal(4) = stress_i[i](5);
		P_i[i].MatTVecMul(rbeta, membranal);
		
		AssembleVector(WorkVec,  1, rd, -alpha_i[i]);
		AssembleVector(WorkVec, 25, rbeta, -alpha_i[i]/dCoef);
	}

	return WorkVec;
}

// Jacobian matrix assembly
VariableSubMatrixHandler&
Shell4::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera e setta gli indici */
	integer iFirstReactionIndex = iGetFirstIndex();
	for (integer i = 0; i < 4; i++) {
		integer iNodeFirstMomIndex = pNode[i]->iGetFirstMomentumIndex();
		integer iNodeFirstPosIndex = pNode[i]->iGetFirstPositionIndex();
		for (int iCnt = 1; iCnt <= 6; iCnt++) {
			WM.PutRowIndex(iCnt + 6 * i, iNodeFirstMomIndex + iCnt);
			WM.PutColIndex(iCnt + 6 * i, iNodeFirstPosIndex + iCnt);
		}
	}
	for (int iCnt = 1; iCnt <= 8; iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
	}

	// tangente
#if 0
	doublereal Kg_pd[24*24], *Kg_ppd[24];
	FullMatrixHandler Kg(Kg_pd, Kg_ppd, 24*24, 24, 24);
#endif

	FullMatrixHandler Kg(24, 24);
	FullMatrixHandler Km(24, 24);
	FullMatrixHandler K_beta_q(8, 24);
	FullMatrixHandler K_beta_beta(8, 8);

	FullMatrixHandler Ktg(24, 15);
	FullMatrixHandler Ktm(24, 12);
	
	FullMatrixHandler Ktbetaq(8, 4);
	FullMatrixHandler Cm(4, 4);
	
	FullMatrixHandler C(12, 12);
	for (integer i = 0; i < NUMIP; i++) {
		D_overline_i[i].MatTMatMul(Ktg, G_i[i]);
		Ktg.MatMatMul(Kg, D_overline_i[i]);

#ifdef USE_CL_IN_SHELL
		C = pD[i]->GetFDE();
#else // ! USE_CL_IN_SHELL
		C = DRef[i];
#endif // ! USE_CL_IN_SHELL

		RotateForward(C, Q_i[i]);
		B_overline_i[i].MatTMatMul(Ktm, C);
		Ktm.MatMatMul(Km, B_overline_i[i]);
		
		// extract Cm matrix from C
		{
			static const integer perm[5] = {0, 1, 4, 2, 5};
			for (integer ir = 1; ir <= 4; ir++) {
				for (integer ic = 1; ic <= 4; ic++) {
					Cm(ir, ic) = C(perm[ir], perm[ic]);
				}
			}
		}
		
		P_i[i].MatTMatMul(Ktbetaq, Cm);
		Ktbetaq.MatMatMul(K_beta_q, B_overline_m_i[i]);
		
		Ktbetaq.MatMatMul(K_beta_beta, P_i[i]);
		
#warning	"FIXME: check dCoef!!!"
		AssembleMatrix(WM, 1, 1, Kg, alpha_i[i]*dCoef);
		AssembleMatrix(WM, 1, 1, Km, alpha_i[i]*dCoef);
		AssembleTransposeMatrix(WM, 1, 25, K_beta_q, alpha_i[i]);
		AssembleMatrix(WM, 25, 1, K_beta_q, alpha_i[i]);
		AssembleMatrix(WM, 25, 25, K_beta_beta, alpha_i[i]/dCoef);
	}
	
	return WorkMat;

}

// Contribution to restart file
std::ostream&
Shell4::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

// Initial settings
void
Shell4::SetValue(DataManager *pDM,
	VectorHandler& /* X */ , VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

VariableSubMatrixHandler&
Shell4::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
Shell4::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
	WorkVec.Resize(0);
	return WorkVec;
}

// Access to nodes
const StructNode*
Shell4::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 4);
	switch (i) {
	case 1:
	case 2:
	case 3:
	case 4:
		return pNode[i - 1];
	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
Shell4::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		if (OH.UseText(OutputHandler::BEAMS)) {
			OH.Beams() << std::setw(8) << GetLabel()
				// TODO: complete
				<< std::endl;
		}
	}
}


typedef LinearElasticGenericConstitutiveLaw<Shell4::vh, Shell4::fmh> LEGCLShell;

Elem*
ReadShell4(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
	StructNode* pN[4];
	Mat3x3 R[4];
	for (unsigned i = 0; i < 4; i++) {
		pN[i] = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		ReferenceFrame RF(pN[i]);
		if (HP.IsKeyWord("orientation")) {
			R[i] = HP.GetRotRel(RF);
		} else {
			R[i] = Eye3;
		}
	}

#ifdef USE_CL_IN_SHELL
	const ConstitutiveLaw<Shell4::vh, Shell4::fmh> *pD[4];

	/* Prestress and prestrain */
	Shell4::vh PreStress(12);
	PreStress.Reset();

	TplDriveCaller<Shell4::vh>* pTplDC = new ZeroTplDriveCaller<Shell4::vh>;

	Shell4::fmh S(12, 12);
	for (unsigned ir = 1; ir <= 12; ir++) {
		for (unsigned ic = 1; ic <= 12; ic++) {
			S(ir, ic) = HP.GetReal();
		}
	}

	pD[0] = 0;
	SAFENEWWITHCONSTRUCTOR(pD[0], LEGCLShell, LEGCLShell(pTplDC, PreStress, S));

	for (unsigned i = 1; i < 4; i++) {
		pD[i] = pD[0]->Copy();
	}
#else // ! USE_CL_IN_SHELL
	Shell4::fmh pD(12, 12);
	if (HP.IsKeyWord("diag")) {
		pD.Reset();
		for (unsigned ir = 1; ir <= 12; ir++) {
			pD(ir, ir) = HP.GetReal();
		}

	} else if (HP.IsKeyWord("sym")) {
		for (unsigned ir = 1; ir <= 12; ir++) {
			pD(ir, ir) = HP.GetReal();
			for (unsigned ic = ir + 1; ic <= 12; ic++) {
				doublereal d = HP.GetReal();
				pD(ir, ic) = d;
				pD(ic, ir) = d;
			}
		}

	} else {
		for (unsigned ir = 1; ir <= 12; ir++) {
			for (unsigned ic = 1; ic <= 12; ic++) {
				pD(ir, ic) = HP.GetReal();
			}
		}
	}
#endif // ! USE_CL_IN_SHELL

	flag fOut = pDM->fReadOutput(HP, Elem::PLATE);

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, Shell4,
		Shell4(uLabel, pDO,
			pN[0], pN[1], pN[2], pN[3],
			R[0], R[1], R[2], R[3],
			pD, fOut));

	return pEl;
}
