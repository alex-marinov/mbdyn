/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include <matvec3n.h>

/* VecN - begin */

#ifdef DEBUG
void VecN::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdVec != NULL);
}
#endif /* DEBUG */

void VecN::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdVec != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdVec, doublereal, ns);
}
 

void VecN::Destroy_(void)
{
   if (pdVec != NULL) {
      SAFEDELETEARR(pdVec);
   }
}


VecN::VecN(void)
: iMaxRows(0), iNumRows(0), pdVec(NULL)
{
   NO_OP;
}


VecN::VecN(integer ns)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   Create_(ns);
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


VecN::VecN(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   Create_(ns);
   Reset(d);
}

VecN::VecN(const VecN& v)
: iMaxRows(v.iNumRows), iNumRows(v.iNumRows), pdVec(NULL) 
{
   Create_(v.iNumRows);
   for (integer iCnt = 0; iCnt < iNumRows; iCnt++) {
      pdVec[iCnt] = v.pdVec[iCnt];
   }
}

/*
  Costruttore da VectorHandler. Prende i valori da iFirstIndex 
  a iFirstIndex+ns. Nota: gli indici del VectorHandler partono da 1, 
  in stile FORTRAN.
*/

VecN::VecN(const VectorHandler& vh, integer ns, integer iFirstIndex)
:iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   ASSERT(iFirstIndex > 0 && iFirstIndex + ns - 1 <= vh.iGetSize());
   Create_(ns);
   for(integer iCnt = 0; iCnt < ns; iCnt++) {
      pdVec[iCnt] = vh(iFirstIndex+iCnt);
   }
}


const VecN& 
VecN::Copy(const VectorHandler& vh, integer iFirstIndex)
{
   ASSERT(iFirstIndex > 0 && iFirstIndex + iNumRows - 1 <= vh.iGetSize());
   for (integer iCnt = 0; iCnt < iNumRows; iCnt++) {
      pdVec[iCnt] = vh(iFirstIndex+iCnt);
   }

   return *this;
}


VecN::~VecN(void)
{
   Destroy_();
}


void VecN::Resize(integer ns)
{
   ASSERT(ns > 0);
   if (pdVec != NULL) {
      Destroy_();
   }
   Create_(ns);
   iMaxRows = iNumRows = ns;
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


void VecN::Reset(const doublereal d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer i = iNumRows; i-- > 0; ) {
      pdVec[i] = d;
   }
}

void VecN::RightMult(const MatNx3& n, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer i = iNumRows; i-- > 0; ) {
      pdVec[i] = n.pdCols[0][i]*v.pdVec[V1]
                 + n.pdCols[1][i]*v.pdVec[V2]
                 + n.pdCols[2][i]*v.pdVec[V3];
   }
}

const VecN& VecN::operator += (const VecN& m)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
    for (integer i = m.iNumRows; i-- > 0; ) {
	 pdVec[i] += m.pdVec[i];
      }
   return *this;
}

const VecN& 
VecN::operator *= (const doublereal& d)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	for (integer i = iNumRows; i-- > 0; ) {
		pdVec[i] *= d;
	}
	return *this;
}

const VecN& VecN::Mult(const MatNxN& m, const VecN& n)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(iNumRows == m.iNumRows);
   ASSERT(m.iNumRows == n.iNumRows);

   for (integer i = 0; i < iNumRows; i++ ) {
      pdVec[i] = 0.;
      for (integer j = 0; j < iNumRows; j++ ) { 
         pdVec[i] += m.pdMat[j][i]*n.pdVec[j];
      }
   }

   return *this;
}

const VecN& 
VecN::Mult(const MatNxN& m, const ArrayView& vm, 
		const VecN& n, const ArrayView& vn)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(iNumRows == m.iNumRows);
   ASSERT(m.iNumRows >= vm.Last());
   ASSERT(n.iNumRows >= vn.Last());
   ASSERT(vm.Number() == vm.Number());

   for(integer i = 0; i < iNumRows; i++) {
      integer jm = vm.Start() - 1;
      integer jn = vn.Start() - 1;

      pdVec[i] = 0.;
      for(integer j = 0; j < vm.Number(); j++) {
         pdVec[i] += m.pdMat[jm][i]*n.pdVec[jn];

	 jm += vm.Offset();
	 jn += vn.Offset();
      }
   }

   return *this;
}

/* VecN - end */


/* Mat3xN - begin */
#ifdef DEBUG
void Mat3xN::IsValid(void) const
{
   ASSERT(iNumCols > 0);
   ASSERT(iMaxCols > 0);
   ASSERT(pdRows[0] != NULL);
}
#endif /* DEBUG */


void Mat3xN::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdRows[0] != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdRows[0], doublereal, 3*ns);
   pdRows[1] = pdRows[0]+ns;
   pdRows[2] = pdRows[1]+ns;
}
  

void Mat3xN::Destroy_(void)
{
   if (pdRows[0] != NULL) {
      SAFEDELETEARR(pdRows[0]);
   }
}


Mat3xN::Mat3xN(void)
: iMaxCols(0), iNumCols(0)
{
   pdRows[0] = NULL;
   pdRows[1] = NULL;
   pdRows[2] = NULL;
}


Mat3xN::Mat3xN(integer nc)
: iMaxCols(nc), iNumCols(nc)
{
   ASSERT(iNumCols > 0);
   pdRows[0] = NULL;
   pdRows[1] = NULL;
   pdRows[2] = NULL;
   
   Create_(iNumCols);
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


Mat3xN::Mat3xN(integer nc, const doublereal& d)
: iMaxCols(nc), iNumCols(nc) 
{
   ASSERT(iNumCols > 0);
   pdRows[0] = NULL;
   pdRows[1] = NULL;
   pdRows[2] = NULL;

   Create_(iNumCols);
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   Reset(d);
}


Mat3xN::~Mat3xN(void)
{
   Destroy_();
}


void Mat3xN::Resize(integer ns)
{
   ASSERT(ns > 0);
   if (ns <= iMaxCols) {
      iNumCols = ns;
   } else {
      Destroy_();
      Create_(ns);
      iMaxCols = iNumCols = ns;
   }
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


void Mat3xN::Reset(const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   for (integer j = iNumCols; j-- > 0; ) {
      pdRows[0][j] = d;
      pdRows[1][j] = d;
      pdRows[2][j] = d;
   }
}


const Mat3xN& Mat3xN::LeftMult(const Mat3x3& m) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer i = iNumCols; i-- > 0; ) {
      doublereal d[3] = { 0., 0., 0. };

      d[0] =
	m.pdMat[M11] * pdRows[0][i]
	+ m.pdMat[M12] * pdRows[1][i]
	+ m.pdMat[M13] * pdRows[2][i];
      d[1] =
	m.pdMat[M21] * pdRows[0][i]
	+ m.pdMat[M22] * pdRows[1][i]
	+ m.pdMat[M23] * pdRows[2][i];
      d[2] =
	m.pdMat[M31] * pdRows[0][i]
	+ m.pdMat[M32] * pdRows[1][i]
	+ m.pdMat[M33] * pdRows[2][i];

      pdRows[0][i] = d[0];
      pdRows[1][i] = d[1];
      pdRows[2][i] = d[2];
   }
   return *this;
}


const Mat3xN& Mat3xN::LeftMult(const Mat3x3& m, const Mat3xN& n)
{
#ifdef DEBUG
   IsValid();
   n.IsValid();
#endif /* DEBUG */

   if (iNumCols != n.iNumCols) {
      Resize(n.iNumCols);
      /* FIXME: sicuri di voler fare resize? non si azzera? */
   }   
   
   for (integer i = iNumCols; i-- > 0; ) {
      pdRows[0][i] =
	m.pdMat[M11] * n.pdRows[0][i]
	+ m.pdMat[M12] * n.pdRows[1][i]
	+ m.pdMat[M13] * n.pdRows[2][i];
      pdRows[1][i] =
	m.pdMat[M21] * n.pdRows[0][i]
	+ m.pdMat[M22] * n.pdRows[1][i]
	+ m.pdMat[M23] * n.pdRows[2][i];
      pdRows[2][i] =
	m.pdMat[M31] * n.pdRows[0][i]
	+ m.pdMat[M32] * n.pdRows[1][i]
	+ m.pdMat[M33] * n.pdRows[2][i];
   }
   return *this;
}


const Mat3xN& Mat3xN::Copy(const Mat3xN& m) 
{
#ifdef DEBUG
   m.IsValid();
#endif /* DEBUG */
   Resize(m.iNumCols);
   for (integer j = m.iNumCols; j-- > 0; ) {
      pdRows[0][j] = m.pdRows[0][j];
      pdRows[1][j] = m.pdRows[1][j];
      pdRows[2][j] = m.pdRows[2][j];
   }
   return *this;
}


const Mat3xN& Mat3xN::operator *= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   if (d != 1.) {
      for (integer j = iNumCols; j-- > 0; ) {
         pdRows[0][j] *= d;
         pdRows[1][j] *= d;
         pdRows[2][j] *= d;
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator /= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   if (d == 0.) {
      silent_cerr("division by zero" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   if (d != 1.) {
      for (integer j = iNumCols; j-- > 0; ) {
         pdRows[0][j] /= d;
         pdRows[1][j] /= d;
         pdRows[2][j] /= d;
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator += (const Mat3xN& m)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer j = m.iNumCols; j-- > 0; ) {
      pdRows[0][j] += m.pdRows[0][j];
      pdRows[1][j] += m.pdRows[1][j];
      pdRows[2][j] += m.pdRows[2][j];
   }
   return *this;
}


const Mat3xN& Mat3xN::operator -= (const Mat3xN& m)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer j = m.iNumCols; j-- > 0; ) {
      pdRows[0][j] -= m.pdRows[0][j];
      pdRows[1][j] -= m.pdRows[1][j];
      pdRows[2][j] -= m.pdRows[2][j];
   }
   return *this;
}

  
Vec3 
Mat3xN::operator * (const VecN& v) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iNumCols == v.iNumRows);

	doublereal d[3] = { 0., 0., 0. };
	for (integer j = iNumCols; j-- > 0; ) {
		d[0] += pdRows[0][j]*v.pdVec[j];
		d[1] += pdRows[1][j]*v.pdVec[j];
		d[2] += pdRows[2][j]*v.pdVec[j];
	}

	return Vec3(d);
}

Vec3 
Mat3xN::Mult(const ArrayView& vm, const VecN& v) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iNumCols >= vm.Last());
	ASSERT(vm.Number() == v.iNumRows);

	doublereal d[3] = { 0., 0., 0. };
	integer jm = vm.Start() - 1;
	for (integer j = 0; j < vm.Number(); j++) {
		d[0] += pdRows[0][jm]*v.pdVec[j];
		d[1] += pdRows[1][jm]*v.pdVec[j];
		d[2] += pdRows[2][jm]*v.pdVec[j];

		jm += vm.Offset();
	}

	return Vec3(d);
}

Vec3 
Mat3xN::Mult(const ArrayView& vm, const VecN& v, const ArrayView& vv) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iNumCols >= vm.Last());
	ASSERT(v.iNumRows >= vv.Last());
	ASSERT(vm.Number() == vv.Number());

	doublereal d[3] = { 0., 0., 0. };
	integer jm = vm.Start() - 1;
	integer jv = vv.Start() - 1;
	for (integer j = 0; j < vm.Number(); j++) {
		d[0] += pdRows[0][jm]*v.pdVec[jv];
		d[1] += pdRows[1][jm]*v.pdVec[jv];
		d[2] += pdRows[2][jm]*v.pdVec[jv];

		jm += vm.Offset();
		jv += vv.Offset();
	}

	return Vec3(d);
}

Vec3 
Mat3xN::GetVec(integer iCol) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iCol > 0 && iCol <= iNumCols);

	--iCol;
	return Vec3(pdRows[0][iCol], pdRows[1][iCol], pdRows[2][iCol]);
}

void 
Mat3xN::PutVec(integer iCol, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iCol > 0 && iCol <= iNumCols);

	--iCol;
	pdRows[0][iCol] = v.pdVec[0];
	pdRows[1][iCol] = v.pdVec[1];
	pdRows[2][iCol] = v.pdVec[2];
}

void 
Mat3xN::AddVec(integer iCol, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iCol > 0 && iCol <= iNumCols);

	--iCol;
	pdRows[0][iCol] += v.pdVec[0];
	pdRows[1][iCol] += v.pdVec[1];
	pdRows[2][iCol] += v.pdVec[2];
}

void 
Mat3xN::SubVec(integer iCol, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iCol > 0 && iCol <= iNumCols);

	--iCol;
	pdRows[0][iCol] -= v.pdVec[0];
	pdRows[1][iCol] -= v.pdVec[1];
	pdRows[2][iCol] -= v.pdVec[2];
}

Mat3x3 
Mat3xN::GetMat3x3(integer iFirstCol) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iFirstCol >= 1 && iFirstCol <= iNumCols-2);
	
	--iFirstCol;
	return Mat3x3(
			pdRows[0][iFirstCol],
			pdRows[1][iFirstCol],
			pdRows[2][iFirstCol],
			pdRows[0][iFirstCol+1],
			pdRows[1][iFirstCol+1],
			pdRows[2][iFirstCol+1],
			pdRows[0][iFirstCol+2],
			pdRows[1][iFirstCol+2],
			pdRows[2][iFirstCol+2]);
}

void
Mat3xN::PutMat3x3(integer iFirstCol, const Mat3x3& m)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iFirstCol >= 1 && iFirstCol <= iNumCols-2);
	
	--iFirstCol;
	pdRows[V1][iFirstCol] = m.pdMat[M11];
	pdRows[V2][iFirstCol] = m.pdMat[M21];
	pdRows[V3][iFirstCol] = m.pdMat[M31];

	iFirstCol++;
	pdRows[V1][iFirstCol] = m.pdMat[M12];
	pdRows[V2][iFirstCol] = m.pdMat[M22];
	pdRows[V3][iFirstCol] = m.pdMat[M32];

	iFirstCol++;
	pdRows[V1][iFirstCol] = m.pdMat[M13];
	pdRows[V2][iFirstCol] = m.pdMat[M23];
	pdRows[V3][iFirstCol] = m.pdMat[M33];
}

void
Mat3xN::AddMat3x3(integer iFirstCol, const Mat3x3& m)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iFirstCol >= 1 && iFirstCol <= iNumCols-2);
	
	--iFirstCol;
	pdRows[V1][iFirstCol] += m.pdMat[M11];
	pdRows[V2][iFirstCol] += m.pdMat[M21];
	pdRows[V3][iFirstCol] += m.pdMat[M31];

	iFirstCol++;
	pdRows[V1][iFirstCol] += m.pdMat[M12];
	pdRows[V2][iFirstCol] += m.pdMat[M22];
	pdRows[V3][iFirstCol] += m.pdMat[M32];

	iFirstCol++;
	pdRows[V1][iFirstCol] += m.pdMat[M13];
	pdRows[V2][iFirstCol] += m.pdMat[M23];
	pdRows[V3][iFirstCol] += m.pdMat[M33];
}

void
Mat3xN::SubMat3x3(integer iFirstCol, const Mat3x3& m)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iFirstCol >= 1 && iFirstCol <= iNumCols-2);
	
	--iFirstCol;
	pdRows[V1][iFirstCol] -= m.pdMat[M11];
	pdRows[V2][iFirstCol] -= m.pdMat[M21];
	pdRows[V3][iFirstCol] -= m.pdMat[M31];

	iFirstCol++;
	pdRows[V1][iFirstCol] -= m.pdMat[M12];
	pdRows[V2][iFirstCol] -= m.pdMat[M22];
	pdRows[V3][iFirstCol] -= m.pdMat[M32];

	iFirstCol++;
	pdRows[V1][iFirstCol] -= m.pdMat[M13];
	pdRows[V2][iFirstCol] -= m.pdMat[M23];
	pdRows[V3][iFirstCol] -= m.pdMat[M33];
}

Mat3x3
Mat3xN::GetMat3x3ScalarMult(integer iFirstCol, const doublereal& d) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iFirstCol >= 1 && iFirstCol <= iNumCols-2);
	
	--iFirstCol;
	return Mat3x3(
			pdRows[0][iFirstCol]*d,
			pdRows[1][iFirstCol]*d,
			pdRows[2][iFirstCol]*d,
			pdRows[0][iFirstCol+1]*d,
			pdRows[1][iFirstCol+1]*d,
			pdRows[2][iFirstCol+1]*d,
			pdRows[0][iFirstCol+2]*d,
			pdRows[1][iFirstCol+2]*d,
			pdRows[2][iFirstCol+2]*d);

}

/* Mat3xN - end */


/* MatNx3 - begin */
#ifdef DEBUG
void MatNx3::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdCols[0] != NULL);
   ASSERT(pdCols[1] != NULL);
   ASSERT(pdCols[2] != NULL);
}
#endif /* DEBUG */


void MatNx3::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdCols[0] != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdCols[0], doublereal, 3*ns);
   pdCols[1] = pdCols[0] + ns;
   pdCols[2] = pdCols[1] + ns;
}
 

void MatNx3::Destroy_(void)
{
   if (pdCols[0] != NULL) {
      SAFEDELETEARR(pdCols[0]);
   }
}


MatNx3::MatNx3(void)
: iMaxRows(0), iNumRows(0)
{
   pdCols[0] = NULL;
}


MatNx3::MatNx3(integer ns)
: iMaxRows(ns), iNumRows(ns)
{
   pdCols[0] = NULL;
   Create_(ns);
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


MatNx3::MatNx3(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns)
{
   pdCols[0] = NULL;
   Create_(ns);
   Reset(d);
}


MatNx3::~MatNx3(void)
{
   Destroy_();
}

void MatNx3::Resize(integer ns)
{
   ASSERT(ns > 0);
   if (pdCols[0] != NULL) {
      Destroy_();
   }
   Create_(ns);
   iMaxRows = iNumRows = ns;
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}

void MatNx3::Reset(const doublereal d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer i = iNumRows; i-- > 0; ) {
      pdCols[0][i] = d;
      pdCols[1][i] = d;
      pdCols[2][i] = d;
   }
}

const MatNx3& MatNx3::Copy(const MatNx3& m) 
{
#ifdef DEBUG
   m.IsValid();
#endif /* DEBUG */
   Resize(m.iNumRows);
   for (integer j = m.iNumRows; j-- > 0; ) {
      pdCols[0][j] = m.pdCols[0][j];
      pdCols[1][j] = m.pdCols[1][j];
      pdCols[2][j] = m.pdCols[2][j];
   }
   return *this;
}

const MatNx3& MatNx3::operator *= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   if (d != 1.) {
      for (integer j = iNumRows; j-- > 0; ) {
         pdCols[0][j] *= d;
         pdCols[1][j] *= d;
         pdCols[2][j] *= d;
      }
   }
   return *this;
}


const MatNx3& MatNx3::operator /= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   if (d == 0.) {
      silent_cerr("division by zero" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   if (d != 1.) {
      for (integer j = iNumRows; j-- > 0; ) {
         pdCols[0][j] /= d;
         pdCols[1][j] /= d;
         pdCols[2][j] /= d;
      }
   }
   return *this;
}

const MatNx3& MatNx3::RightMult(const MatNx3& n, const Mat3x3& m)
{
#ifdef DEBUG
   IsValid();
   n.IsValid();
#endif /* DEBUG */

   if (iNumRows != n.iNumRows) {
      Resize(n.iNumRows);
   }   
   
   for (integer i = iNumRows; i-- > 0; ) {
      pdCols[0][i] = 
	n.pdCols[0][i]*m.pdMat[M11]
	+n.pdCols[1][i]*m.pdMat[M21]
	+n.pdCols[2][i]*m.pdMat[M31];
      pdCols[1][i] = 
	n.pdCols[0][i]*m.pdMat[M12]
	+n.pdCols[1][i]*m.pdMat[M22]
	+n.pdCols[2][i]*m.pdMat[M32];
      pdCols[2][i] =
	n.pdCols[0][i]*m.pdMat[M13]
	+n.pdCols[1][i]*m.pdMat[M23]
	+n.pdCols[2][i]*m.pdMat[M33];
   }

   return *this;
}

Vec3 
MatNx3::GetVec(integer iRow) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iRow > 0 && iRow <= iNumRows);

	--iRow;
	return Vec3(pdCols[0][iRow], pdCols[1][iRow], pdCols[2][iRow]);
}

void 
MatNx3::PutVec(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iRow > 0 && iRow <= iNumRows);

	--iRow;
	pdCols[0][iRow] = v.pdVec[0];
	pdCols[1][iRow] = v.pdVec[1];
	pdCols[2][iRow] = v.pdVec[2];
}

void 
MatNx3::AddVec(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iRow > 0 && iRow <= iNumRows);

	--iRow;
	pdCols[0][iRow] += v.pdVec[0];
	pdCols[1][iRow] += v.pdVec[1];
	pdCols[2][iRow] += v.pdVec[2];
}

void 
MatNx3::SubVec(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	ASSERT(iRow > 0 && iRow <= iNumRows);

	--iRow;
	pdCols[0][iRow] -= v.pdVec[0];
	pdCols[1][iRow] -= v.pdVec[1];
	pdCols[2][iRow] -= v.pdVec[2];
}


const MatNx3& MatNx3::Transpose(const Mat3xN& n)
{
#ifdef DEBUG
  IsValid();
  n.IsValid();
#endif /* DEBUG */

  for (integer i = iNumRows; i-- > 0; ) {
     pdCols[0][i] = n.pdRows[0][i];
     pdCols[1][i] = n.pdRows[1][i];
     pdCols[2][i] = n.pdRows[2][i];
  }
  return *this;
}
/* MatNx3 - end */


/* MatNxN - begin */

#ifdef DEBUG
void MatNxN::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdVec != NULL);
   ASSERT(pdMat != NULL);
}
#endif /* DEBUG */


void MatNxN::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdVec != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdVec, doublereal, ns*ns);
   SAFENEWARR(pdMat, doublereal*, ns);

   pdMat[0] = pdVec;
   for (integer i = 1; i < ns; i++) {
      pdMat[i] = pdMat[i-1]+ns;
   }
   iNumRows = iMaxRows = ns;
}
 

void MatNxN::Destroy_(void)
{
   if (pdVec != NULL) {
      SAFEDELETEARR(pdVec);
   }

   if (pdMat != NULL) {
      SAFEDELETEARR(pdMat);
   }
   iNumRows = iMaxRows = 0;
}


MatNxN::MatNxN(void)
: iMaxRows(0), iNumRows(0), pdVec(NULL), pdMat(NULL)
{
   NO_OP;
}


MatNxN::MatNxN(integer ns)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL), pdMat(NULL)
{
   Create_(ns);
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}


MatNxN::MatNxN(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL), pdMat(NULL)
{
   Create_(ns);
   Reset(d);
}


MatNxN::~MatNxN(void)
{
   Destroy_();
}

void MatNxN::Reset(const doublereal d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   for (integer i = iNumRows*iNumRows; i-- > 0; ) {
      pdVec[i] = d;
   }
}

const MatNxN&
MatNxN::Copy(const MatNxN& m) 
{
#ifdef DEBUG
   m.IsValid();
#endif /* DEBUG */
   if (iMaxRows < m.iNumRows) Create_(m.iNumRows);
   else iNumRows = m.iNumRows;
	   
   for (integer c = 0; c < iNumRows; c++) {
#ifdef HAVE_MEMCPY
      memcpy(pdMat[c], m.pdMat[c], sizeof(doublereal)*iNumRows);
#else // ! HAVE_MEMCPY
      for (integer r = 0; r < iNumRows; r++) {
         pdMat[c][r] = m.pdMat[c][r];
      }
#endif // ! HAVE_MEMCPY
   }
   return *this;
}

const MatNxN& MatNxN::operator *= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   if (d != 1.) {
      for (integer j = iNumRows * iNumRows; j-- > 0; ) {
         pdVec[j] *= d;
      }
   }
   return *this;
}


const MatNxN& MatNxN::operator /= (const doublereal& d) 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   if (d == 0.) {
      silent_cerr("division by zero" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   if (d != 1.) {
      for (integer j = iNumRows * iNumRows; j-- > 0; ) {
         pdVec[j] /= d;
      }
   }
   return *this;
}

const MatNxN& MatNxN::Mult(const MatNx3& m, const Mat3xN& n)
{
#ifdef DEBUG
   IsValid();
   n.IsValid();
   m.IsValid();
#endif /* DEBUG */

   ASSERT(m.iNumRows == iNumRows);  
   ASSERT(m.iNumRows == n.iNumCols); 

   for (integer i = 0; i < iNumRows; i++) {
      for (integer j = 0; j < iNumRows; j++) {
         pdMat[j][i] = 
		 m.pdCols[0][i]*n.pdRows[0][j]
		 + m.pdCols[1][i]*n.pdRows[1][j]
		 + m.pdCols[2][i]*n.pdRows[2][j];
      }
   }

   return *this;
}

std::ostream&
operator << (std::ostream& out, const MatNxN& m)
{
	for (integer r = 0; r < m.iNumRows; r++) {
		for (integer c = 0; c < m.iNumRows; c++) {
			out << " " << m.pdMat[c][r];
		}
		out << std::endl;
	}

	return out;
}

/* MatNxN - end */

