/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <matvec3n.h>

void VecN::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdVec != NULL);
}


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
   IsValid();
}


VecN::VecN(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   Create_(ns);
   Reset(d);
}

/*
  Costruttore da VectorHandler. Prende i valori da iFirstIndex 
  a iFirstIndex+ns. Nota: gli indici del VectorHandler partono da 1, 
  in stile FORTRAN.
*/

 VecN::VecN(const VectorHandler& vh, integer ns, integer iFirstIndex)
   :iMaxRows(ns), iNumRows(ns), pdVec(NULL)
 {
      ASSERT(iFirstIndex > 0 && iFirstIndex <= vh.iGetSize()-(ns-1));
      Create_(ns);
      for(integer iCnt = 0; iCnt < ns; iCnt++) {
        pdVec[iCnt] = vh.dGetCoef(iFirstIndex+iCnt);
      }
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
   IsValid();
}


void VecN::Reset(const doublereal& d)
{
   IsValid();
   for (integer i = iNumRows; i-- > 0; ) {
      pdVec[i] = d;
   }
}

void VecN::RightMult(const MatNx3&n, const Vec3& v)
{
   IsValid();
   for (integer i = iNumRows; i-- > 0; ) {
      pdVec[i] = n.pdVec1[i]*v.dGet(1)+
                 n.pdVec2[i]*v.dGet(2)+
                 n.pdVec3[i]*v.dGet(3);
   }
}

const VecN& VecN::operator += (const VecN& m)
{
   IsValid();
    for (integer i = m.iNumRows; i-- > 0; ) {
	 pdVec[i] += m.pdVec[i];
      }
   return *this;
}

const VecN& VecN::operator * (const doublereal& d)
{
   IsValid();
    for (integer i = iNumRows; i-- > 0; ) {
	 pdVec[i] = pdVec[i]*d;
      }
   return *this;
}

const VecN& VecN::Mult(const MatNxN& m, const VecN& n)
{
  IsValid();
  ASSERT(iNumRows=n.iNumRows);
  for(integer i = iNumRows; i-- > 0; ) {
  double temp=0;
   for(integer j = iNumRows; j-- > 0; ) { 
      temp+=m.pdVec[i*iNumRows+j]*n.pdVec[j];
     }
  pdVec[i] = temp; 
  }
  return *this;
}


/* VecN - end */


/* Mat3xN - begin */

void Mat3xN::IsValid(void) const
{
   ASSERT(iNumCols > 0);
   ASSERT(iMaxCols > 0);
   ASSERT(pvdMat[0] != NULL);
}


void Mat3xN::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pvdMat[0] != NULL) {
      Destroy_();
   }
   SAFENEWARR(pvdMat[0], doublereal, 3*ns);
   pvdMat[1] = pvdMat[0]+ns;
   pvdMat[2] = pvdMat[1]+ns;
}
  

void Mat3xN::Destroy_(void)
{
   if (pvdMat[0] != NULL) {
      SAFEDELETEARR(pvdMat[0]);
   }
}


Mat3xN::Mat3xN(void)
: iMaxCols(0), iNumCols(0)
{
   pvdMat[0] = NULL;
   pvdMat[1] = NULL;
   pvdMat[2] = NULL;
}


Mat3xN::Mat3xN(integer nc)
: iMaxCols(nc), iNumCols(nc)
{
   ASSERT(iNumCols > 0);
   pvdMat[0] = NULL;
   pvdMat[1] = NULL;
   pvdMat[2] = NULL;
   
   Create_(iNumCols);
   IsValid();
}


Mat3xN::Mat3xN(integer nc, const doublereal& d)
: iMaxCols(nc), iNumCols(nc) 
{
   ASSERT(iNumCols > 0);
   pvdMat[0] = NULL;
   pvdMat[1] = NULL;
   pvdMat[2] = NULL;

   Create_(iNumCols);
   IsValid();
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
   IsValid();
}


void Mat3xN::Reset(const doublereal& d)
{
   IsValid();
   for (int i = 3; i-- > 0; ) {
      for (integer j = iNumCols; j-- > 0; ) {
	 pvdMat[i][j] = d;
      }
   }
}


const Mat3xN& Mat3xN::LeftMult(const Mat3x3& m) 
{
   IsValid();
   for (integer i = iNumCols; i-- > 0; ) {
      doublereal d1 =
	pvdMat[0][i]*m.dGet(1, 1)
	+pvdMat[1][i]*m.dGet(1, 2)
	+pvdMat[2][i]*m.dGet(1, 3);
      doublereal d2 =
	pvdMat[0][i]*m.dGet(2, 1)
	+pvdMat[1][i]*m.dGet(2, 2)
	+pvdMat[2][i]*m.dGet(2, 3);
      doublereal d3 =
	pvdMat[0][i]*m.dGet(3, 1)
	+pvdMat[1][i]*m.dGet(3, 2)
	+pvdMat[2][i]*m.dGet(3, 3);
      pvdMat[0][i] = d1;
      pvdMat[1][i] = d2;
      pvdMat[2][i] = d3;
   }
   return *this;
}


const Mat3xN& Mat3xN::LeftMult(const Mat3x3& m, const Mat3xN& n)
{
   IsValid();
   n.IsValid();
   if (iNumCols != n.iNumCols) {
      Resize(n.iNumCols);
   }   
   
   for (integer i = iNumCols; i-- > 0; ) {
      pvdMat[0][i] =
	n.pvdMat[0][i]*m.dGet(1, 1)
	+n.pvdMat[1][i]*m.dGet(1, 2)
	+n.pvdMat[2][i]*m.dGet(1, 3);
      pvdMat[1][i] =
	n.pvdMat[0][i]*m.dGet(2, 1)
	+n.pvdMat[1][i]*m.dGet(2, 2)
	+n.pvdMat[2][i]*m.dGet(2, 3);
      pvdMat[2][i] =
	n.pvdMat[0][i]*m.dGet(3, 1)
	+n.pvdMat[1][i]*m.dGet(3, 2)
	+n.pvdMat[2][i]*m.dGet(3, 3);
   }
   return *this;
}


const Mat3xN& Mat3xN::Copy(const Mat3xN& m) 
{
   m.IsValid();
   Resize(m.iNumCols);
   for (integer i = 3; i-- > 0; ) {
      for (integer j = m.iNumCols; j-- > 0; ) {
	 pvdMat[i][j] = m.pvdMat[i][j];
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator *= (const doublereal& d) 
{
   IsValid();
   if (d != 1.) {
      for (int i = 3; i-- > 0; ) {
	 for (integer j = iNumCols; j-- > 0; ) {
	    pvdMat[i][j] *= d;
	 }
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator /= (const doublereal& d) 
{
   IsValid();
   if (d == 0.) {
      std::cerr << "division by zero" << std::endl;
      THROW(ErrGeneric());
   }
   if (d != 1.) {
      for (int i = 3; i-- > 0; ) {
	 for (integer j = iNumCols; j-- > 0; ) {
	    pvdMat[i][j] /= d;
	 }
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator += (const Mat3xN& m)
{
   IsValid();
   for (int i = 3; i-- > 0; ) {
      for (integer j = m.iNumCols; j-- > 0; ) {
	 pvdMat[i][j] += m.pvdMat[i][j];
      }
   }
   return *this;
}


const Mat3xN& Mat3xN::operator -= (const Mat3xN& m)
{
   IsValid();
   for (int i = 3; i-- > 0; ) {
      for (integer j = m.iNumCols; j-- > 0; ) {
	 pvdMat[i][j] -= m.pvdMat[i][j];
      }
   }
   return *this;
}

  
Vec3 Mat3xN::operator * (const VecN& v)
{
   IsValid();
   ASSERT(iNumCols == v.iNumRows);
   doublereal d[3];
   for (int i = 3; i-- > 0; ) {
      d[i] = 0.;
      for (integer j = iNumCols; j-- > 0; ) {
	 d[i] += pvdMat[i][j]*v.pdVec[j];
      }
   }
   return Vec3(d);
}

/* Mat3xN - end */

/* MatNx3 - begin */

void MatNx3::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdVec1 != NULL);
   ASSERT(pdVec2 != NULL);
   ASSERT(pdVec3 != NULL);
}


void MatNx3::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdVec1 != NULL) {
      Destroy_();
   }
   if (pdVec2 != NULL) {
      Destroy_();
   }
   if (pdVec3 != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdVec1, doublereal, ns);
   SAFENEWARR(pdVec2, doublereal, ns);
   SAFENEWARR(pdVec3, doublereal, ns);

}
 

void MatNx3::Destroy_(void)
{
   if (pdVec1 != NULL) {
      SAFEDELETEARR(pdVec1);
   }
    if (pdVec2 != NULL) {
      SAFEDELETEARR(pdVec2);
   }
    if (pdVec3 != NULL) {
      SAFEDELETEARR(pdVec3);
  }
}


MatNx3::MatNx3(void)
: iMaxRows(0), iNumRows(0), pdVec1(NULL), pdVec2(NULL), pdVec3(NULL)
{
   NO_OP;
}


MatNx3::MatNx3(integer ns)
: iMaxRows(ns), iNumRows(ns), pdVec1(NULL), pdVec2(NULL), pdVec3(NULL)
{
   Create_(ns);
   IsValid();
}


MatNx3::MatNx3(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns), pdVec1(NULL), pdVec2(NULL), pdVec3(NULL)
{
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
   if (pdVec1 != NULL) {
      Destroy_();
   }
   if (pdVec2 != NULL) {
      Destroy_();
   }
   if (pdVec3 != NULL) {
      Destroy_();
   }
   Create_(ns);
   iMaxRows = iNumRows = ns;
   IsValid();
}

void MatNx3::Reset(const doublereal& d)
{
   IsValid();
   for (integer i = iNumRows; i-- > 0; ) {
      pdVec1[i] = d;
      pdVec2[i] = d;
      pdVec3[i] = d;
   }
}

const MatNx3& MatNx3::RightMult(const MatNx3& n, const Mat3x3& m)
{  
   doublereal temp1=0, temp2=0, temp3=0;
   IsValid();
   n.IsValid();
   if (iNumRows != n.iNumRows) {
      Resize(n.iNumRows);
   }   
   
   for (integer i = iNumRows; i-- > 0; ) {
      temp1 =
	n.pdVec1[i]*m.dGet(1, 1)
	+n.pdVec2[i]*m.dGet(2, 1)
	+n.pdVec3[i]*m.dGet(3, 1);
      temp2 =
	n.pdVec1[i]*m.dGet(1, 2)
	+n.pdVec2[i]*m.dGet(2, 2)
	+n.pdVec3[i]*m.dGet(3, 2);
      temp3 =
	n.pdVec1[i]*m.dGet(1, 3)
	+n.pdVec2[i]*m.dGet(2, 3)
	+n.pdVec3[i]*m.dGet(3, 3);
      pdVec1[i]=temp1;
      pdVec2[i]=temp2;
      pdVec3[i]=temp3;
   }
   return *this;
}

const MatNx3& MatNx3::Transpose(const Mat3xN& n)
{
  IsValid();
  n.IsValid();
  for (integer i = iNumRows; i-- > 0; ) {
   pdVec1[i]=n.pvdMat[0][i];
   pdVec2[i]=n.pvdMat[1][i];
   pdVec3[i]=n.pvdMat[2][i];
  }
  return *this;
}
/* MatNx3 - end */

/* MatNxN - begin */

void MatNxN::IsValid(void) const
{      
   ASSERT(iMaxRows > 0);
   ASSERT(iNumRows > 0);
   ASSERT(pdVec != NULL);
}


void MatNxN::Create_(integer ns)
{
   ASSERT(ns > 0);
   if (pdVec != NULL) {
      Destroy_();
   }
   SAFENEWARR(pdVec, doublereal, ns*ns);
}
 

void MatNxN::Destroy_(void)
{
   if (pdVec != NULL) {
      SAFEDELETEARR(pdVec);
   }
}


MatNxN::MatNxN(void)
: iMaxRows(0), iNumRows(0), pdVec(NULL)
{
   NO_OP;
}


MatNxN::MatNxN(integer ns)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   Create_(ns*ns);
   IsValid();
}


MatNxN::MatNxN(integer ns, const doublereal& d)
: iMaxRows(ns), iNumRows(ns), pdVec(NULL)
{
   Create_(ns*ns);
   Reset(d);
}


MatNxN::~MatNxN(void)
{
   Destroy_();
}

void MatNxN::Reset(const doublereal& d)
{
   IsValid();
   for (integer i = iNumRows*iNumRows; i-- > 0; ) {
      pdVec[i] = d;
   }
}

const MatNxN& MatNxN::Mult(const MatNx3& m, const Mat3xN& n)
{
   IsValid();
   n.IsValid();
   m.IsValid();
   ASSERT(m.iNumRows==iNumRows);  
   ASSERT(m.iNumRows==n.iNumCols); 
   for (integer i = 0; i < iNumRows; i++) {
     for (integer j = 0; j < iNumRows; j++) {
       double temp = 0;
       for(integer iCnt = 0; iCnt < 3; iCnt++) {
         temp += m.dGet(i+1, iCnt+1)*n.dGet(iCnt+1, j+1);
       } 
      pdVec[i*iNumRows+j] = temp;
     }
   }      
   return *this;
}



/* MatNxN - end */
