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

#ifndef SHELL_HC
#define SHELL_HC

#if 0
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
#endif

#if 0
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
#endif

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
		dest(i + start_row) += source(i) * dCoef;
	}
	return;
}

#if 0
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
	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			dest.IncCoef(ir + start_row, ic + start_col, source(ir, ic) * dCoef);
			// dest.ppdColsm1[ic + start_col][ir + start_row] += dCoef*source.ppdColsm1[ic][ir];
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
	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			dest.IncCoef(ir + start_row, ic + start_col, source(ic, ir) * dCoef);
			// dest.ppdColsm1[ic + start_col][ir + start_row] += dCoef*source.ppdColsm1[ir][ic];
		}
	}

	return;
}
#endif

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
RotateForward(MyVectorHandler & e, const Mat3x3& R)
{
	//TODO. error checking: check dimensions!
	//check da fare: e.iGetNumRows() == 12
	Vec3 t, t1;	
	for (integer b = 0; b < 4; b++) {
		ExtractVec3(t, e, 1 + b * 3);
		t1 = R * t;
		// InsertVector(e, 1 + b * 3, t1);
		e.Put(1 + b * 3, t1);
	}
}

static inline void
RotateBackward(MyVectorHandler & e, const Mat3x3& R)
{
	//TODO. error checking: check dimensions!
	//check da fare: e.iGetNumRows() == 12
	Vec3 t, t1;	
	for (integer b = 0; b < 4; b++) {
		ExtractVec3(t, e, 1 + b * 3);
		t1 = R.MulTV(t);
		// InsertVector(e, 1 + b * 3, t1);
		e.Put(1 + b * 3, t1);
	}
}

static inline void
RotateForward(FullMatrixHandler & C, const Mat3x3& R)
{
	//TODO. error checking: check dimensions!
	//check da fare: C.iGetNumRows() == 12
	//check da fare: C.iGetNumCols() == 12
	Mat3x3 m, m1;	
	for (integer rb = 0; rb < 4; rb++) {
		for (integer rc = 0; rc < 4; rc++) {
			ExtractMat3x3(m, C, 1 + rb * 3, 1 + rc * 3);
			m1 = R * m.MulMT(R);
			// InsertMatrix(C, 1 + rb * 3, 1 + rc * 3, m1);
			C.Put(1 + rb * 3, 1 + rc * 3, m1);
		}
	}
}

static inline doublereal L1(const doublereal xi[2]) {
	return 0.25 * (1. + xi[0]) * (1. + xi[1]);
};

static inline doublereal L2(const doublereal xi[2]) {
	return 0.25 * (1. - xi[0]) * (1. + xi[1]);
};

static inline doublereal L3(const doublereal xi[2]) {
	return 0.25 * (1. - xi[0]) * (1. - xi[1]);
};

static inline doublereal L4(const doublereal xi[2]) {
	return 0.25 * (1. + xi[0]) * (1. - xi[1]);
};

typedef doublereal (*LI_Type)(const doublereal xi[2]);
static LI_Type LI[4] = {&L1, &L2, &L3, &L4};



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
static LI_J_Type LI_J[4][2] = {
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
	out(1, 2) = -a(1, 2) / det;
	out(2, 1) = -a(2, 1) / det;
	out(2, 2) =  a(1, 1) / det;

	return;
}

static void
Inv3x3(const FullMatrixHandler& a, FullMatrixHandler & out)
{
	//FIXME: Mettere controlli dimensioni matrici In/Out
	//FIXME: mettere controllo det \neq 0
	doublereal det =
		a(1, 1)*a(2, 2)*a(3, 3)-a(1, 1)*a(2, 3)*a(3, 2)-a(2, 1)*a(1, 2)*a(3, 3)+a(2, 1)*a(1, 3)*a(3, 2)+a(3, 1)*a(1, 2)*a(2, 3)-a(3, 1)*a(1, 3)*a(2, 2);

	out(1, 1) =  (a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2))/det;
	out(1, 2) = -(a(1, 2)*a(3, 3)+a(1, 3)*a(3, 2))/det;
	out(1, 3) =  (a(1, 2)*a(2, 3)-a(1, 3)*a(2, 2))/det;

	out(2, 1) = -(a(2, 1)*a(3, 3)+a(2, 3)*a(3, 1))/det;
	out(2, 2) =  (a(1, 1)*a(3, 3)-a(1, 3)*a(3, 1))/det;
	out(2, 3) = -(a(1, 1)*a(2, 3)+a(1, 3)*a(2, 1))/det;

	out(3, 1) =  (a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1))/det;
	out(3, 2) = -(a(1, 1)*a(3, 2)+a(1, 2)*a(3, 1))/det;
	out(3, 3) =  (a(1, 1)*a(2, 2)-a(1, 2)*a(2, 1))/det;

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

static void
InvBlockDiagonal3_2x3_2(const FullMatrixHandler& a, FullMatrixHandler & out)
{
	//FIXME: Mettere controlli dimensioni matrici In/Out
	//FIXME: mettere controllo det \neq 0
	
	// 3x3
	doublereal det =
		a(1, 1)*a(2, 2)*a(3, 3)-a(1, 1)*a(2, 3)*a(3, 2)-a(2, 1)*a(1, 2)*a(3, 3)+a(2, 1)*a(1, 3)*a(3, 2)+a(3, 1)*a(1, 2)*a(2, 3)-a(3, 1)*a(1, 3)*a(2, 2);

	out(1, 1) =  (a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2))/det;
	out(1, 2) = -(a(1, 2)*a(3, 3)+a(1, 3)*a(3, 2))/det;
	out(1, 3) =  (a(1, 2)*a(2, 3)-a(1, 3)*a(2, 2))/det;

	out(2, 1) = -(a(2, 1)*a(3, 3)+a(2, 3)*a(3, 1))/det;
	out(2, 2) =  (a(1, 1)*a(3, 3)-a(1, 3)*a(3, 1))/det;
	out(2, 3) = -(a(1, 1)*a(2, 3)+a(1, 3)*a(2, 1))/det;

	out(3, 1) =  (a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1))/det;
	out(3, 2) = -(a(1, 1)*a(3, 2)+a(1, 2)*a(3, 1))/det;
	out(3, 3) =  (a(1, 1)*a(2, 2)-a(1, 2)*a(2, 1))/det;


	// 2x2
	det = a(4, 4)*a(5, 5) - a(4, 5)*a(5, 4);

	out(4, 4) =  a(5, 5) / det;
	out(4, 5) = -a(4, 5) / det;
	out(5, 4) = -a(5, 4) / det;
	out(5, 5) =  a(4, 4) / det;

	return;
}

static void
InvBlockDiagonal4_2x4_2(const FullMatrixHandler& a, FullMatrixHandler & out)
{
	//FIXME: Mettere controlli dimensioni matrici In/Out
	//FIXME: mettere controllo det \neq 0
	
	// 4x4
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


	// 2x2
	det = a(5, 5)*a(6, 6) - a(5, 6)*a(6, 5);

	out(5, 5) =  a(6, 6) / det;
	out(5, 6) = -a(5, 6) / det;
	out(6, 5) = -a(6, 5) / det;
	out(6, 6) =  a(5, 5) / det;

	return;
}

#endif // SHELL_HC

// vim:ft=c
