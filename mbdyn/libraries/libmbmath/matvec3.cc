/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* vettori 3 e matrici 3x3 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>
#include <cfloat>
#include <limits>

#include "matvec3.h"

/* noteworthy constant */
const doublereal Zero1(0.);
const Vec3 Zero3(0., 0., 0.);
const Mat3x3 Eye3(1., 0., 0., 0., 1., 0., 0., 0., 1.);
const Mat3x3 Zero3x3(0., 0., 0., 0., 0., 0., 0., 0., 0.);


/* Vec3 - begin */

Mat3x3
Vec3::Tens(const Vec3& v) const
{
   return Mat3x3(pdVec[V1]*v.pdVec[V1],
		 pdVec[V2]*v.pdVec[V1],
		 pdVec[V3]*v.pdVec[V1],
		 pdVec[V1]*v.pdVec[V2],
		 pdVec[V2]*v.pdVec[V2],
		 pdVec[V3]*v.pdVec[V2],
		 pdVec[V1]*v.pdVec[V3],
		 pdVec[V2]*v.pdVec[V3],
		 pdVec[V3]*v.pdVec[V3]);
}

/* Prodotto "tensore".  Restituisce se stesso per se stesso */
Mat3x3
Vec3::Tens(void) const
{
	return Tens(*this);
}     

/* Prodotto vettore per matrice */
Mat3x3
Vec3::Cross(const Mat3x3& m) const {
   return Mat3x3(pdVec[V2]*m.pdMat[M31]-pdVec[V3]*m.pdMat[M21],
		 pdVec[V3]*m.pdMat[M11]-pdVec[V1]*m.pdMat[M31],
		 pdVec[V1]*m.pdMat[M21]-pdVec[V2]*m.pdMat[M11],
		 pdVec[V2]*m.pdMat[M32]-pdVec[V3]*m.pdMat[M22],
		 pdVec[V3]*m.pdMat[M12]-pdVec[V1]*m.pdMat[M32],
		 pdVec[V1]*m.pdMat[M22]-pdVec[V2]*m.pdMat[M12],
		 pdVec[V2]*m.pdMat[M33]-pdVec[V3]*m.pdMat[M23],
		 pdVec[V3]*m.pdMat[M13]-pdVec[V1]*m.pdMat[M33],
		 pdVec[V1]*m.pdMat[M23]-pdVec[V2]*m.pdMat[M13]);
}


/**
 * Scalar product
 * multiplies self by matrix m; equivalent to m.Transpose() * this.
 */
Vec3
Vec3::operator * (const Mat3x3& m) const
{
	return
		Vec3(
			m.pdMat[M11]*pdVec[V1]
				+ m.pdMat[M21]*pdVec[V2]
				+ m.pdMat[M31]*pdVec[V3],
			m.pdMat[M12]*pdVec[V1]
				+ m.pdMat[M22]*pdVec[V2]
				+ m.pdMat[M32]*pdVec[V3],
			m.pdMat[M13]*pdVec[V1]
				+ m.pdMat[M23]*pdVec[V2]
				+ m.pdMat[M33]*pdVec[V3]);
}

void
Vec3::Reset(void)
{
	pdVec[V1] = 0.;
	pdVec[V2] = 0.;
	pdVec[V3] = 0.;
}

/* Vec3 - end */

/* Mat3x3 - begin */

/* inversione */
doublereal
Mat3x3::dDet(void) const
{
   doublereal* p = (doublereal*)pdMat;

   return p[M11]*(p[M22]*p[M33]-p[M23]*p[M32])
     +p[M12]*(p[M23]*p[M31]-p[M21]*p[M33])
     +p[M13]*(p[M21]*p[M32]-p[M22]*p[M31]);
}
   
/* inversione */
Mat3x3
Mat3x3::Inv(const doublereal &d) const
{
   ASSERT(fabs(d) > std::numeric_limits<doublereal>::epsilon());

   doublereal* p = (doublereal*)pdMat;

   return Mat3x3((p[M22]*p[M33]-p[M23]*p[M32])/d,
		 (p[M23]*p[M31]-p[M21]*p[M33])/d,
		 (p[M21]*p[M32]-p[M22]*p[M31])/d,
		 (p[M13]*p[M32]-p[M12]*p[M33])/d,
		 (p[M11]*p[M33]-p[M13]*p[M31])/d,
		 (p[M12]*p[M31]-p[M11]*p[M32])/d,
		 (p[M12]*p[M23]-p[M13]*p[M22])/d,
		 (p[M13]*p[M21]-p[M11]*p[M23])/d,
		 (p[M11]*p[M22]-p[M12]*p[M21])/d);
}

/* inversione */
Mat3x3
Mat3x3::Inv(void) const
{
   doublereal d = dDet();
   if (fabs(d) < std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("matrix is singular" << std::endl);
      throw MatrixHandler::ErrMatrixIsSingular(MBDYN_EXCEPT_ARGS);
   }
   
   return Inv(d);
}


/* soluzione */
Vec3
Mat3x3::Solve(const doublereal& d, const Vec3& v) const
{
   const doublereal* p = pdMat;
   const doublereal* pv = v.pGetVec();

   ASSERT(fabs(d) > std::numeric_limits<doublereal>::epsilon());
   
   return Vec3((pv[V1]*(p[M22]*p[M33] - p[M23]*p[M32])
		+pv[V2]*(p[M13]*p[M32] - p[M12]*p[M33])
		+pv[V3]*(p[M12]*p[M23] - p[M13]*p[M22]))/d,
	       (pv[V1]*(p[M23]*p[M31] - p[M21]*p[M33])
		+pv[V2]*(p[M11]*p[M33] - p[M13]*p[M31])
		+pv[V3]*(p[M13]*p[M21] - p[M11]*p[M23]))/d,
	       (pv[V1]*(p[M21]*p[M32] - p[M22]*p[M31])
		+pv[V2]*(p[M12]*p[M31] - p[M11]*p[M32])
		+pv[V3]*(p[M11]*p[M22] - p[M12]*p[M21]))/d);
}

/* soluzione */
Vec3
Mat3x3::Solve(const Vec3& v) const
{
   doublereal d = dDet();
   
   if (fabs(d) < std::numeric_limits<doublereal>::epsilon()) {
      silent_cerr("matrix is singular" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }

   return Solve(d, v);
}

Vec3
Mat3x3::LDLSolve(const Vec3& v) const
{
	doublereal d1 = 0., d2 = 0., d3 = 0., l21 = 0., l31 = 0., l32 = 0.;

	d1 = pdMat[M11];
	ASSERT(d1 >= 0.);
	if (d1 > std::numeric_limits<doublereal>::epsilon()) {
		l21 = (pdMat[M21] + pdMat[M12])/(2.*d1);
		l31 = (pdMat[M31] + pdMat[M13])/(2.*d1);
	}

	d2 = pdMat[M22] - l21*l21*d1;
	ASSERT(d2 >= 0.);
	if (d2 > std::numeric_limits<doublereal>::epsilon()) {
		l32 = ((pdMat[M32] + pdMat[M23])/2. - l21*l31*d1)/d2;
	}

	d3 = pdMat[M33] - l31*l31*d1 - l32*l32*d2;

	// L * D * L^T * x = v
	// L^T * x = y
	// D * y = z
	// L * z = v

	// z = L^-1 * v
	doublereal z1 = v(1);
	doublereal z2 = v(2) - l21*z1;
	doublereal z3 = v(3) - l31*z1 - l32*z2;

	// y = D^-1 * z
	if (d1 > std::numeric_limits<doublereal>::epsilon()) {
		z1 /= d1;
	} else {
		z1 = 0.;
	}

	if (d2 > std::numeric_limits<doublereal>::epsilon()) {
		z2 /= d2;
	} else {
		z2 = 0.;
	}

	if (d3 > std::numeric_limits<doublereal>::epsilon()) {
		z3 /= d3;
	} else {
		z3 = 0.;
	}

	// x = L^-T * y
	z2 -= l32*z3;
	z1 -= l21*z2 + l31*z3;

	return Vec3(z1, z2, z3);
}

bool
Mat3x3::EigSym(Vec3& EigenValues) const
{
	Mat3x3 EigenVectors;

	return EigSym(EigenValues, EigenVectors);
}

#if 0
bool
Mat3x3::EigSym(Vec3& EigenValues,
	Mat3x3& EigenVectors,
	const Mat3x3& Hints) const
{
	Mat3x3 Tmp;

	if (!EigSym(EigenValues, Tmp)) {
		return false;
	}

	/* TBC */
}
#endif

bool
Mat3x3::EigSym(Vec3& EigenValues, Mat3x3& EigenVectors) const
{
	// From:
	// W.M. Scherzinger, C.R. Dohrmann,
	// `A robust algorithm for finding the eigenvalues and eigenvectors
	// of 3x3 symmetric matrices'
	// Comput. Methods Appl. Mech. Engrg. 2008
	// doi:10.1016/j.cma.2008.03.031

	if (!IsSymmetric(std::numeric_limits<doublereal>::epsilon())) {
		return false;
	}

	if (IsDiag(std::numeric_limits<doublereal>::epsilon())) {
		EigenVectors = Eye3;
		EigenValues = Vec3(pdMat[M11], pdMat[M22], pdMat[M33]);

		return true;
	}

	Mat3x3 AA = *this;

	doublereal trA_3 = Trace()/3;

	AA(1, 1) -= trA_3;
	AA(2, 2) -= trA_3;
	AA(3, 3) -= trA_3;

	doublereal J2 = (AA*AA).Trace()/2;

	if (fabs(J2) < std::numeric_limits<doublereal>::epsilon()) {
		EigenVectors = Eye3;
		EigenValues = Vec3(pdMat[M11], pdMat[M22], pdMat[M33]);

		return true;
	}

	doublereal J2_dmy = sqrt(J2/3.);
	doublereal J3 = AA.dDet();
	doublereal dmy = J3/2/(J2_dmy*J2_dmy*J2_dmy);
	doublereal alpha;

	// NOTE: we want real eigenvalues; this requires the matrix to be
	// positive definite or semi-definite
	if (dmy < -1.) {
		dmy = -1.;
		alpha = M_PI;

	} else if (dmy > 1.) {
		dmy = 1.;
		alpha = 0.;

	} else {
		alpha = acos(dmy)/3;
	}

	int idx1;
	if (alpha < M_PI/6.) {
		idx1 = 1;

	} else {
		idx1 = 3;
	}

	doublereal eta1 = 2*J2_dmy*cos(alpha + 2./3.*M_PI*(idx1 - 1));
	EigenValues(idx1) = eta1 + trA_3;

	// NOTE: there's a typo in the original paper;
	// AA must be used instead of A
	Vec3 r1 = AA.GetVec(1);
	r1(1) -= eta1;
	Vec3 r2 = AA.GetVec(2);
	r2(2) -= eta1;
	Vec3 r3 = AA.GetVec(3);
	r3(3) -= eta1;

	doublereal
		nr1 = r1.Norm(),
		nr2 = r2.Norm(),
		nr3 = r3.Norm();

	int irmax = 1;
	doublereal nrmax = nr1;

	if (nr2 > nrmax) {
		irmax = 2;
		nrmax = nr2;
	}

	if (nr3 > nrmax) {
		irmax = 3;
		nrmax = nr3;
	}

	if (irmax == 2) {
		Vec3 rtmp = r2;
		r2 = r1;
		nr2 = nr1;
		r1 = rtmp;
		nr1 = nrmax;

	} else if (irmax == 3) {
		Vec3 rtmp = r3;
		r3 = r1;
		nr3 = nr1;
		r1 = rtmp;
		nr1 = nrmax;
	}

	Vec3 s1, s2;
	s1 = r1/nr1;
	Vec3 t2 = r2 - s1*(s1*r2);
	doublereal nt2 = t2.Norm();
	Vec3 t3 = r3 - s1*(s1*r3);
	doublereal nt3 = t3.Norm();

	if (nt2 > nt3) {
		s2 = t2/nt2;

	} else {
		s2 = t3/nt3;
	}

	Vec3 v1(s1.Cross(s2));
	EigenVectors.PutVec(idx1, v1);

	Mat3x3 AAA(eta1, 0., 0.,
		0., s1*(AA*s1), s2*(AA*s1),
		0., s1*(AA*s2), s2*(AA*s2));

	int idx2 = idx1%3 + 1;
	int idx3 = (idx1 + 1)%3 + 1;

	doublereal AAA22p33 = AAA(2, 2) + AAA(3, 3);
	doublereal AAA22m33 = AAA(2, 2) - AAA(3, 3);
	doublereal eta2 = (AAA22p33 - copysign(1., AAA22m33)*sqrt(AAA22m33*AAA22m33 + 4*AAA(2, 3)*AAA(3, 2)))/2;
	doublereal eta3 = AAA22p33 - eta2;

	EigenValues(idx2) = eta2 + trA_3;
	EigenValues(idx3) = eta3 + trA_3;

	Vec3 u1 = AA*s1 - s1*eta2;
	doublereal nu1 = u1.Norm();
	Vec3 u2 = AA*s2 - s2*eta2;
	doublereal nu2 = u2.Norm();

	Vec3 w1;
	if (nu1 > nu2) {
		w1 = u1/nu1;

	} else {
		w1 = u2/nu2;
	}

	Vec3 v2(w1.Cross(v1));
	EigenVectors.PutVec(idx2, v2);
	EigenVectors.PutVec(idx3, v1.Cross(v2));

	return true;
}

/**
 * multiply by another matrix, transposed: this * m^T
 */
Mat3x3
Mat3x3::MulMT(const Mat3x3& m) const
{
	return Mat3x3(
		pdMat[M11]*m.pdMat[M11]
			+ pdMat[M12]*m.pdMat[M12]
			+ pdMat[M13]*m.pdMat[M13],
		pdMat[M21]*m.pdMat[M11]
			+ pdMat[M22]*m.pdMat[M12]
			+ pdMat[M23]*m.pdMat[M13],
		pdMat[M31]*m.pdMat[M11]
			+ pdMat[M32]*m.pdMat[M12]
			+ pdMat[M33]*m.pdMat[M13],

		pdMat[M11]*m.pdMat[M21]
			+ pdMat[M12]*m.pdMat[M22]
			+ pdMat[M13]*m.pdMat[M23],
		pdMat[M21]*m.pdMat[M21]
			+ pdMat[M22]*m.pdMat[M22]
			+ pdMat[M23]*m.pdMat[M23],
		pdMat[M31]*m.pdMat[M21]
			+ pdMat[M32]*m.pdMat[M22]
			+ pdMat[M33]*m.pdMat[M23],

		pdMat[M11]*m.pdMat[M31]
			+ pdMat[M12]*m.pdMat[M32]
			+ pdMat[M13]*m.pdMat[M33],
		pdMat[M21]*m.pdMat[M31]
			+ pdMat[M22]*m.pdMat[M32]
			+ pdMat[M23]*m.pdMat[M33],
		pdMat[M31]*m.pdMat[M31]
			+ pdMat[M32]*m.pdMat[M32]
			+ pdMat[M33]*m.pdMat[M33]);
}

/**
 * multiply self transposed by a vector: this^T * v
 */
Vec3
Mat3x3::MulTV(const Vec3& v) const
{
	return Vec3(
		pdMat[M11]*v.pdVec[V1]
			+ pdMat[M21]*v.pdVec[V2]
			+ pdMat[M31]*v.pdVec[V3],
		pdMat[M12]*v.pdVec[V1]
			+ pdMat[M22]*v.pdVec[V2]
			+ pdMat[M32]*v.pdVec[V3],
		pdMat[M13]*v.pdVec[V1]
			+ pdMat[M23]*v.pdVec[V2]
			+ pdMat[M33]*v.pdVec[V3]);
}

/**
 * multiply self transposed by another matrix: this^T * m
 */
Mat3x3
Mat3x3::MulTM(const Mat3x3& m) const
{
	return Mat3x3(
		pdMat[M11]*m.pdMat[M11]
			+ pdMat[M21]*m.pdMat[M21]
			+ pdMat[M31]*m.pdMat[M31],
		pdMat[M12]*m.pdMat[M11]
			+ pdMat[M22]*m.pdMat[M21]
			+ pdMat[M32]*m.pdMat[M31],
		pdMat[M13]*m.pdMat[M11]
			+ pdMat[M23]*m.pdMat[M21]
			+ pdMat[M33]*m.pdMat[M31],

		pdMat[M11]*m.pdMat[M12]
			+ pdMat[M21]*m.pdMat[M22]
			+ pdMat[M31]*m.pdMat[M32],
		pdMat[M12]*m.pdMat[M12]
			+ pdMat[M22]*m.pdMat[M22]
			+ pdMat[M32]*m.pdMat[M32],
		pdMat[M13]*m.pdMat[M12]
			+ pdMat[M23]*m.pdMat[M22]
			+ pdMat[M33]*m.pdMat[M32],

		pdMat[M11]*m.pdMat[M13]
			+ pdMat[M21]*m.pdMat[M23]
			+ pdMat[M31]*m.pdMat[M33],
		pdMat[M12]*m.pdMat[M13]
			+ pdMat[M22]*m.pdMat[M23]
			+ pdMat[M32]*m.pdMat[M33],
		pdMat[M13]*m.pdMat[M13]
			+ pdMat[M23]*m.pdMat[M23]
			+ pdMat[M33]*m.pdMat[M33]);
}

/**
 * multiply self transposed by another matrix, transposed: this^T * m^T
 */
Mat3x3
Mat3x3::MulTMT(const Mat3x3& m) const
{
	return Mat3x3(
		pdMat[M11]*m.pdMat[M11]
			+ pdMat[M21]*m.pdMat[M12]
			+ pdMat[M31]*m.pdMat[M13],
		pdMat[M12]*m.pdMat[M11]
			+ pdMat[M22]*m.pdMat[M12]
			+ pdMat[M32]*m.pdMat[M13],
		pdMat[M13]*m.pdMat[M11]
			+ pdMat[M23]*m.pdMat[M12]
			+ pdMat[M33]*m.pdMat[M13],

		pdMat[M11]*m.pdMat[M21]
			+ pdMat[M21]*m.pdMat[M22]
			+ pdMat[M31]*m.pdMat[M23],
		pdMat[M12]*m.pdMat[M21]
			+ pdMat[M22]*m.pdMat[M22]
			+ pdMat[M32]*m.pdMat[M23],
		pdMat[M13]*m.pdMat[M21]
			+ pdMat[M23]*m.pdMat[M22]
			+ pdMat[M33]*m.pdMat[M23],

		pdMat[M11]*m.pdMat[M31]
			+ pdMat[M21]*m.pdMat[M32]
			+ pdMat[M31]*m.pdMat[M33],
		pdMat[M12]*m.pdMat[M31]
			+ pdMat[M22]*m.pdMat[M32]
			+ pdMat[M32]*m.pdMat[M33],
		pdMat[M13]*m.pdMat[M31]
			+ pdMat[M23]*m.pdMat[M32]
			+ pdMat[M33]*m.pdMat[M33]);
}

/**
 * multiply self times vCross
 */
Mat3x3
Mat3x3::MulVCross(const Vec3& v) const
{
	return Mat3x3(
		pdMat[M12]*v.pdVec[V3] - pdMat[M13]*v.pdVec[V2],
		pdMat[M22]*v.pdVec[V3] - pdMat[M23]*v.pdVec[V2],
		pdMat[M32]*v.pdVec[V3] - pdMat[M33]*v.pdVec[V2],

		pdMat[M13]*v.pdVec[V1] - pdMat[M11]*v.pdVec[V3],
		pdMat[M23]*v.pdVec[V1] - pdMat[M21]*v.pdVec[V3],
		pdMat[M33]*v.pdVec[V1] - pdMat[M31]*v.pdVec[V3],

		pdMat[M11]*v.pdVec[V2] - pdMat[M12]*v.pdVec[V1],
		pdMat[M21]*v.pdVec[V2] - pdMat[M22]*v.pdVec[V1],
		pdMat[M31]*v.pdVec[V2] - pdMat[M32]*v.pdVec[V1]);
}

/**
 * multiply self transposed times vCross
 */
Mat3x3
Mat3x3::MulTVCross(const Vec3& v) const
{
	return Mat3x3(
		pdMat[M21]*v.pdVec[V3] - pdMat[M31]*v.pdVec[V2],
		pdMat[M22]*v.pdVec[V3] - pdMat[M32]*v.pdVec[V2],
		pdMat[M23]*v.pdVec[V3] - pdMat[M33]*v.pdVec[V2],

		pdMat[M31]*v.pdVec[V1] - pdMat[M11]*v.pdVec[V3],
		pdMat[M32]*v.pdVec[V1] - pdMat[M12]*v.pdVec[V3],
		pdMat[M33]*v.pdVec[V1] - pdMat[M13]*v.pdVec[V3],

		pdMat[M11]*v.pdVec[V2] - pdMat[M21]*v.pdVec[V1],
		pdMat[M12]*v.pdVec[V2] - pdMat[M22]*v.pdVec[V1],
		pdMat[M13]*v.pdVec[V2] - pdMat[M23]*v.pdVec[V1]);
}

void
Mat3x3::Reset(void)
{
	pdMat[M11] = 0.;
	pdMat[M12] = 0.;
	pdMat[M13] = 0.;

	pdMat[M21] = 0.;
	pdMat[M22] = 0.;
	pdMat[M23] = 0.;

	pdMat[M31] = 0.;
	pdMat[M32] = 0.;
	pdMat[M33] = 0.;
}

Mat3x3
Mat3x3::Skew(void) const
{
	return Mat3x3(MatCross, this->Ax());
}

/* Mat3x3 - end */

const Mat3x3Zero_Manip Mat3x3Zero;
const Mat3x3DEye_Manip Mat3x3DEye;
const Mat3x3Diag_Manip Mat3x3Diag;
const MatCross_Manip MatCross;
const MatCrossCross_Manip MatCrossCross;

/* Manipolatori */
namespace CGR_Rot {
const Param_Manip Param;
const MatR_Manip MatR;
const MatG_Manip MatG;
const MatGm1_Manip MatGm1;
};


Vec3 operator - (const Vec3& v)
{
   return Vec3(-v.pdVec[V1],
	       -v.pdVec[V2],
	       -v.pdVec[V3]);
}


Mat3x3 operator - (const Mat3x3& m)
{
   const doublereal* pdMat = m.pGetMat();
   return Mat3x3(-pdMat[M11],
		 -pdMat[M21],
		 -pdMat[M31],
		 -pdMat[M12],
		 -pdMat[M22],
		 -pdMat[M32],
		 -pdMat[M13],
		 -pdMat[M23],
		 -pdMat[M33]);
}


const char sForm[] = "%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e";
const char sDefFill[] = " ";

/* output di matrici */

std::ostream& 
operator << (std::ostream& out, const Mat3x3& m)
{
   const doublereal* pd = m.pGetMat();

   out
     << pd[M11] << sDefFill << pd[M12] << sDefFill << pd[M13] << sDefFill
     << pd[M21] << sDefFill << pd[M22] << sDefFill << pd[M23] << sDefFill
     << pd[M31] << sDefFill << pd[M32] << sDefFill << pd[M33];
   
   return out;
}


std::ostream& 
Write(std::ostream& out, const Mat3x3& m, const char* s, const char* s2)
{
   return m.Write(out, s, s2);
}


/* output di vettori */

std::ostream&
operator << (std::ostream& out, const Vec3& v)
{
   const doublereal* pd = v.pGetVec();
   
   out << pd[0] << sDefFill << pd[1] << sDefFill << pd[2];
   
   return out;
}


std::ostream& 
Write(std::ostream& out, const Vec3& v, const char* s)
{
   return v.Write(out, s);
}


/* Output di matrici */
std::ostream&
Mat3x3::Write(std::ostream& out, const char* sFill, const char* sFill2) const
{
   char* sF2 = (char*)sFill2; 
   if (sFill2 == NULL) {
      sF2 = (char*)sFill;
   }
   out 
     << pdMat[M11] << sFill << pdMat[M12] << sFill << pdMat[M13] << sF2
     << pdMat[M21] << sFill << pdMat[M22] << sFill << pdMat[M23] << sF2
     << pdMat[M31] << sFill << pdMat[M32] << sFill << pdMat[M33];
   
   return out;
}


std::ostream&
Vec3::Write(std::ostream& out, const char* sFill) const
{
   out << pdVec[V1] << sFill << pdVec[V2] << sFill << pdVec[V3];
   
   return out;
}


std::ostream&
Write(std::ostream& out, const doublereal& d, const char*)
{
   return out << d;
}


/* calcolo dei parametri di rotazione a partire dalla matrice R */

Vec3
MatR2gparam(const Mat3x3& m)
{
	/* test di singolarita' */
	doublereal d = 1. + m.Trace();
   
	if (fabs(d) < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("MatR2gparam(): divide by zero, "
		"probably due to singularity in rotation parameters" << std::endl);
		throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
	}
   
	return m.Ax()*(4./d);
}

/* So-called linear parametrization, consisting in ax(R - R') */
Vec3
MatR2LinParam(const Mat3x3& m)
{
	return m.Ax();
}

/* Calcolo della matrice R a partire da due vettori sghembi */

Mat3x3 MatR2vec(unsigned short int ia, const Vec3& va, 
		unsigned short int ib, const Vec3& vb)
{
   ASSERT(ia >= 1 && ia <= 3); 
   ASSERT(ib >= 1 && ib <= 3); 
   
   Vec3 r[3];

   DEBUGCOUT("MatR2vec: ia = " << ia << " (" << va << "),"
	     << " ib = " << ib << " (" << vb << ")" << std::endl);
   
   if (ia < 1 || ia > 3) {
      silent_cerr("MatR2vec: first index is illegal" 
	      << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   int i1 = ia-1;
   int i2 = ia%3;
   int i3 = (ia+1)%3;
   
   if (ib == (ia%3)+1) {
      doublereal d = va.Norm();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: first vector must be non-null" << std::endl );
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      r[i1] = va/d;
      d = vb.Norm();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: second vector must be non-null" << std::endl );
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      r[i3] = r[i1].Cross(vb);
      d = r[i3].Dot();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: vectors must be distinct" 
		 << std::endl);
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }	
      d = sqrt(d);
      r[i3] /= d;
      r[i2] = r[i3].Cross(r[i1]);
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << std::endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else if (ib == ((ia+1)%3+1)) {
      doublereal d = va.Norm();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: first vector must be non-null" << std::endl );
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      r[i1] = va/d;
      d = vb.Norm();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: second vector must be non-null" << std::endl );
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      r[i2] = vb.Cross(r[i1]);
      d = r[i2].Dot();
      if (d <= std::numeric_limits<doublereal>::epsilon()) {
	 silent_cerr("MatR2vec: vectors must be distinct" 
		 << std::endl);
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }	
      d = sqrt(d);
      r[i2] /= d;
      r[i3] = r[i1].Cross(r[i2]);  
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << std::endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else {
      silent_cerr("MatR2vec: second index is illegal" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   return ::Zero3x3; // phony call, not reachable
}


/* Subroutine che calcola gli angoli di Eulero a partire dalla matrice R
      SUBROUTINE RPYPAR(ALFA, PHI)
      
      IMPLICIT NONE
      REAL*8 PHI(3), ALFA(3,3)

      REAL*8 RADEGR
      PARAMETER (RADEGR = 180.D0/3.141592653589793D0)

      REAL*8 CA, SA

      PHI(1) = ATAN2(- ALFA(2,3), ALFA(3,3))
      CA = COS(PHI(1))
      SA = SIN(PHI(1))
      PHI(1) = PHI(1)*RADEGR
      PHI(2) = ATAN2(ALFA(1,3), - SA*ALFA(2,3) + CA*ALFA(3,3))*RADEGR
      PHI(3) = ATAN2(CA*ALFA(2,1) + SA*ALFA(3,1),
     &               CA*ALFA(2,2) + SA*ALFA(3,2))*RADEGR
      RETURN
      END

 */

const doublereal dRaDegr = 180./M_PI;

Vec3
MatR2EulerAngles(const Mat3x3& R)
{
	return MatR2EulerAngles123(R);
}

Vec3
MatR2EulerAngles123(const Mat3x3& R)
{  
	doublereal dAlpha = atan2(-R(2, 3), R(3, 3));
	doublereal dCosAlpha = cos(dAlpha);
	doublereal dSinAlpha = sin(dAlpha);
      
	return Vec3(dAlpha,
	       atan2(R(1, 3), 
			dCosAlpha*R(3, 3) - dSinAlpha*R(2, 3)),
	       atan2(dCosAlpha*R(2, 1) + dSinAlpha*R(3, 1),
			dCosAlpha*R(2, 2) + dSinAlpha*R(3, 2)));   
}

/*
 R_313 = R_a * R_b * R_c

       [  ca -sa   0 ]
 R_a = [  sa  ca   0 ]
       [   0   0   1 ]

       [   1   0   0 ]
 R_a = [   0  cb -sb ]
       [   0  sb  cb ]

       [  cc -sc   0 ]
 R_c = [  sc  cc   0 ]
       [   0   0   1 ]

         [  ca*cc-cb*sa*sc -ca*sc-cb*cc*sa  sa*sb ]
 R_313 = [  cb*sa+ca*cb*sc  ca*cb*cc-sa*sc -ca*sb ]
         [           sb*sc           cc*sb     cb ]
 */

Vec3
MatR2EulerAngles313(const Mat3x3& R)
{
	doublereal dAlpha = atan2(R(1, 3), -R(2, 3));
	doublereal dCosAlpha = cos(dAlpha);
	doublereal dSinAlpha = sin(dAlpha);
      
	return Vec3(dAlpha,
	       atan2(dSinAlpha*R(1, 3) - dCosAlpha*R(2, 3),
			R(3, 3)),
	       atan2(-dCosAlpha*R(1, 2) - dSinAlpha*R(2, 2),
			dCosAlpha*R(1, 1) + dSinAlpha*R(2, 1)));
}

Vec3
MatR2EulerAngles321(const Mat3x3& R)
{  
	doublereal dAlpha = atan2(R(2, 1), R(1, 1));
	doublereal dCosAlpha = cos(dAlpha);
	doublereal dSinAlpha = sin(dAlpha);
      
	return Vec3(dAlpha,
	       atan2(-R(3, 1),
			dCosAlpha*R(1, 1) + dSinAlpha*R(2, 1)),
	       atan2(dSinAlpha*R(1, 3) - dCosAlpha*R(2, 3),
			-dSinAlpha*R(1, 2) + dCosAlpha*R(2, 2)));
}

void MatR2EulerParams(const Mat3x3& R, doublereal& e0, Vec3& e)
{
   doublereal t = R.Tr();
   doublereal T[4];
   T[0] = 1. + t;
   T[1] = 1. - t + 2.*R.dGet(1, 1);
   T[2] = 1. - t + 2.*R.dGet(2, 2);
   T[3] = 1. - t + 2.*R.dGet(3, 3);
   
   int k = 0;
   for (int i = 1; i <= 3; i++) {
      if (fabs(T[i]) >= fabs(T[k])) {
	 k = i;       
      }
   }
   
   switch (k) {
    case 0:   
      e0 = .5*sqrt(T[0]);
      e = Vec3(R.dGet(3, 2)-R.dGet(2, 3),
	       R.dGet(1, 3)-R.dGet(3, 1),
	       R.dGet(2, 1)-R.dGet(1, 2))/(4.*e0);
      break;
      
    case 1: {
       doublereal e1 = .5*sqrt(T[1]);       
       e0 = (R.dGet(3, 2) - R.dGet(2, 3))/(4.*e1);
       e = Vec3(T[1],
		R.dGet(2, 1)+R.dGet(1, 2),
		R.dGet(3, 1)+R.dGet(1, 3))/(4.*e1);
       break;
    }
      
    case 2: {
       doublereal e2 = .5*sqrt(T[2]);
       e0 = (R.dGet(1, 3) - R.dGet(3, 1))/(4.*e2);
       e = Vec3(R.dGet(1, 2)+R.dGet(2, 1),
		T[2],
		R.dGet(3, 2)+R.dGet(2, 3))/(4.*e2);
       break;
    }
      
    case 3: {
       doublereal e3 = .5*sqrt(T[3]);
       e0 = (R.dGet(2, 1) - R.dGet(1, 2))/(4.*e3);
       e = Vec3(R.dGet(1, 3)+R.dGet(3, 1),
		R.dGet(2, 3)+R.dGet(3, 2),
		T[3])/(4.*e3);
       break;
    }
   }
}

Mat3x3
EulerAngles2MatR(const Vec3& v)
{
	return EulerAngles123_2MatR(v);
}

Mat3x3
EulerAngles123_2MatR(const Vec3& v)
{
	doublereal d = v(1);
	doublereal dCosAlpha(cos(d));
	doublereal dSinAlpha(sin(d));
	d = v(2);
	doublereal dCosBeta(cos(d));
	doublereal dSinBeta(sin(d));
	d = v(3);
	doublereal dCosGamma(cos(d));
	doublereal dSinGamma(sin(d));
   
	return Mat3x3(
		dCosBeta*dCosGamma,
		dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma,
		dSinAlpha*dSinGamma - dCosAlpha*dSinBeta*dCosGamma,
		-dCosBeta*dSinGamma,
		dCosAlpha*dCosGamma - dSinAlpha*dSinBeta*dSinGamma,
		dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma,
		dSinBeta,
		-dSinAlpha*dCosBeta,
		dCosAlpha*dCosBeta);
};

Mat3x3
EulerAngles313_2MatR(const Vec3& v)
{
	doublereal d = v(1);
	doublereal dCosAlpha(cos(d));
	doublereal dSinAlpha(sin(d));
	d = v(2);
	doublereal dCosBeta(cos(d));
	doublereal dSinBeta(sin(d));
	d = v(3);
	doublereal dCosGamma(cos(d));
	doublereal dSinGamma(sin(d));
   
	return Mat3x3(
		dCosAlpha*dCosGamma - dSinAlpha*dCosBeta*dSinGamma,
		dSinAlpha*dCosGamma + dCosAlpha*dCosBeta*dSinGamma,
		dSinBeta*dSinGamma,
		-dCosAlpha*dSinGamma - dSinAlpha*dCosBeta*dCosGamma,
		-dSinAlpha*dSinGamma + dCosAlpha*dCosBeta*dCosGamma,
		dSinBeta*dCosGamma,
		dSinAlpha*dSinBeta,
		-dCosAlpha*dSinBeta,
		dCosBeta);
}

Mat3x3
EulerAngles321_2MatR(const Vec3& v)
{
	doublereal d = v(1);
	doublereal dCosAlpha(cos(d));
	doublereal dSinAlpha(sin(d));
	d = v(2);
	doublereal dCosBeta(cos(d));
	doublereal dSinBeta(sin(d));
	d = v(3);
	doublereal dCosGamma(cos(d));
	doublereal dSinGamma(sin(d));
   
	return Mat3x3(
		dCosAlpha*dCosBeta,
		dSinAlpha*dCosBeta,
		-dSinBeta,
		-dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma,
		dCosAlpha*dCosGamma + dSinAlpha*dSinBeta*dSinGamma,
		dCosBeta*dSinGamma,
		dSinAlpha*dSinGamma + dCosAlpha*dSinBeta*dCosGamma,
		-dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma,
		dCosBeta*dCosGamma);
}

Vec3
Unwrap(const Vec3& vPrev, const Vec3& v)
{
	doublereal dTheta = v.Norm();
	if (dTheta > 0.) {
		doublereal dThetaPrev = vPrev.Norm();
		if (dThetaPrev > std::numeric_limits<doublereal>::epsilon()) {
			bool b(false);

			if (vPrev*v < 0) {
				dTheta = -dTheta;
			}

			doublereal dThetaOld = dTheta;

			while (dTheta - dThetaPrev > M_PI) {
				dTheta -= 2.*M_PI;
				b = true;
			}

			while (dThetaPrev - dTheta > M_PI) {
				dTheta += 2.*M_PI;
				b = true;
			}

			if (b) {
				return v*(dTheta/dThetaOld);
			}
		}
	}

	return v;
}

template <>
bool
IsNull(const doublereal& d)
{
	return d == 0.;
}
 
template <>
bool
IsExactlySame(const doublereal& d1, const doublereal& d2)
{
	return d1 == d2;
}

template <>
bool
IsSame(const doublereal& d1, const doublereal& d2, const doublereal& dTol)
{
	return fabs(d1 - d2) <= dTol;
}

Vec3
MultRV(const Vec3& v, const Mat3x3& R)
{
	return R*v;
}

Mat3x3
MultRM(const Mat3x3& m, const Mat3x3& R)
{
	return R*m;
}

Mat3x3
MultMRt(const Mat3x3& m, const Mat3x3& R)
{
	return m.MulMT(R);
}

Mat3x3
MultRMRt(const Mat3x3& m, const Mat3x3& R)
{
	return R*m.MulMT(R);
}

