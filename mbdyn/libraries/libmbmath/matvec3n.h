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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifndef MATVEC3N_H
#define MATVEC3N_H

#include "matvec3.h"

/* VecN - begin */

class VecN {
   friend class Mat3xN;
   friend class MatNx3;
   friend class MatNxN;
   
 private:
   
   /* not defined */
   VecN(const VecN&);
   const VecN& operator = (const VecN&);
   
 protected:
   integer iMaxRows;
   integer iNumRows;
   doublereal* pdVec;
   
   void IsValid(void) const;   
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   VecN(void);
   VecN(integer ns);
   VecN(integer ns, const doublereal& d);
 
  /* costruttore da VectorHandler (Aggiunta) */
   VecN(const VectorHandler& vh, integer ns, integer iFirstIndex);

   virtual ~VecN(void);
   
   inline integer iGetNumRows(void) const;
   
   void Resize(integer ns);   
   void Reset(const doublereal& d = 0.);
   
   inline void Put(integer i, const doublereal& d);
   inline void Add(integer i, const doublereal& d);
   inline void Sub(integer i, const doublereal& d);   
   inline const doublereal& dGet(integer i) const;
  
   void RightMult(const MatNx3& n, const Vec3& v);
   const VecN& Mult(const MatNxN& m, const VecN& n); 
   const VecN& operator += (const VecN& m);
 

   /* prodotto per scalare */
   const VecN& operator * (const doublereal & d);
 
};


inline integer VecN::iGetNumRows(void) const 
{
   IsValid();
   return iNumRows;
}


inline void VecN::Put(integer i, const doublereal& d)
{   
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] = d;
}


inline void VecN::Add(integer i, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] += d;
}

inline void VecN::Sub(integer i, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] -= d;
}


inline const doublereal& VecN::dGet(integer i) const
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   return pdVec[--i];
}

/* VecN - end */


/* Mat3xN - begin */

class Mat3xN {

 friend class MatNx3;
 friend class MatNxN;

 private:
   
   /* not defined */
   Mat3xN(const Mat3xN&);
   const Mat3xN& operator = (const Mat3xN&);
   
 protected:
   integer iMaxCols;
   integer iNumCols;
   doublereal* pvdMat[3];

   void IsValid(void) const;
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   Mat3xN(void); /* to allow arrays of Mat3xN */
   Mat3xN(integer nc);
   Mat3xN(integer nc, const doublereal& d);
   virtual ~Mat3xN(void);
   
   void Resize(integer ns);
   void Reset(const doublereal& d = 0.);
   
   inline integer iGetNumCols(void) const;
   
   inline void Put(int i, integer j, const doublereal& d);
   inline void Add(int i, integer j, const doublereal& d);
   inline void Sub(int i, integer j, const doublereal& d);
   inline const doublereal& dGet(int i, integer j) const;

   /* *this = m x *this */
   const Mat3xN& LeftMult(const Mat3x3& m);

   /* *this = m x *this */
   const Mat3xN& LeftMult(const Mat3x3& m, const Mat3xN& n);

   const Mat3xN& Mult(const Mat3xN& m, const MatNxN& n);
   const Mat3xN& Copy(const Mat3xN& m);
   
   const Mat3xN& operator *= (const doublereal& d);
   const Mat3xN& operator /= (const doublereal& d);

   const Mat3xN& operator += (const Mat3xN& m);
   const Mat3xN& operator -= (const Mat3xN& m);
   
   Vec3 operator * (const VecN& v);
};


inline integer Mat3xN::iGetNumCols(void) const
{
   IsValid();
   return iNumCols;
}
   
inline void Mat3xN::Put(int i, integer j, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pvdMat[--i][--j] = d;
}

inline void Mat3xN::Add(int i, integer j, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pvdMat[--i][--j] += d;
}

inline void Mat3xN::Sub(int i, integer j, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pvdMat[--i][--j] -= d;
};   

inline const doublereal& Mat3xN::dGet(int i, integer j) const
{
   IsValid();
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   return pvdMat[--i][--j];
}
      
/* Mat3xN - end */

/* MatNx3 - begin */

 /* classe aggiunta per gestire operazioni con matrici Nx3. Per adesso
    sono memorizzate come tre array Nx1  */

class MatNx3 {

 friend class VecN;
 friend class Mat3xN;
 friend class MatNxN;
 private:
   
   /* not defined */
   MatNx3(const MatNx3&);
   const MatNx3& operator = (const MatNx3&);
   
 protected:
   integer iMaxRows;
   integer iNumRows;
   doublereal* pdVec1;
   doublereal* pdVec2;
   doublereal* pdVec3;
   
   void IsValid(void) const;   
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   MatNx3(void);
   MatNx3(integer ns);
   MatNx3(integer ns, const doublereal& d);   
   virtual ~MatNx3(void);
   
   inline integer iGetNumRows(void) const;
   void Resize(integer ns); 
   void Reset(const doublereal& d = 0.);
   
   inline void Put(integer i, integer j, const doublereal& d);
   inline void Add(integer i, integer j, const doublereal& d);
   inline void Sub(integer i, integer j, const doublereal& d);   
   inline const doublereal& dGet(integer i, integer j) const;
   /* *this = *this x m */
   const MatNx3& RightMult(const MatNx3& n, const Mat3x3& m);
   /* *this = [3xN]T    */
   const MatNx3& Transpose(const Mat3xN& n);
   const MatNx3& Mult(const MatNxN& m, const MatNx3& n);
};


inline integer MatNx3::iGetNumRows(void) const 
{
   IsValid();
   return iNumRows;
}


inline void MatNx3::Put(integer i, integer j, const doublereal& d)
{   
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   if(j==1)
   pdVec1[--i] = d;
   else if(j==2)
   pdVec2[--i] = d;
   else
   pdVec3[--i] = d;
}


inline void MatNx3::Add(integer i, integer j, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   if(j==1)
   pdVec1[--i] += d;
   else if(j==2)
   pdVec2[--i] += d;
   else
   pdVec3[--i] += d;
}

inline const doublereal& MatNx3::dGet(integer i, integer j) const
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   if(j==1)
   return pdVec1[--i];
   else if(j==2)
   return pdVec2[--i];
   else
   return pdVec3[--i];
}

/* Mat Nx3 - end */

/* MatNxN - begin */
/* Questa classe memorizza matrici quadrate di dimensioni n x n in array (n^2) x 1.
   Nota: in origine l'ho messa per memorizzare le matrici di massa e rigidezza modale.
   In verita' queste matrici devono comunque essere diagonali (o non si riesce a escludere i
   modi rigidi quando si esportano i modi nel multi-corpo) quindi non dovrebbe servire piu' */

class MatNxN {
   friend class Mat3xN;
   friend class MatNx3;
   friend class VecN;
   
 private:
   
   /* not defined */
   MatNxN(const MatNxN&);
   const MatNxN& operator = (const MatNxN&);
   
 protected:
   integer iMaxRows;
   integer iNumRows;
   doublereal* pdVec;
   
   void IsValid(void) const;   
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   MatNxN(void);
   MatNxN(integer ns);
   MatNxN(integer ns, const doublereal& d);   
   virtual ~MatNxN(void);
   
   inline integer iGetNumRows(void) const;
   void Reset(const doublereal& d = 0.);
   inline void Put(integer i, integer j, const doublereal& d);
   inline void Add(integer i, integer j, const doublereal& d);
   inline const doublereal& dGet(integer i, integer j) const;

   const MatNxN& Mult(const MatNx3&m, const Mat3xN& n);
};


inline integer MatNxN::iGetNumRows(void) const 
{
   IsValid();
   return iNumRows;
}


inline void MatNxN::Put(integer i, integer j, const doublereal& d)
{   
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   pdVec[--i*iNumRows+(--j)] = d;
}


inline void MatNxN::Add(integer i, integer j, const doublereal& d)
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   pdVec[--i*iNumRows+(--j)] += d;
}


inline const doublereal& MatNxN::dGet(integer i, integer j) const
{
   IsValid();
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   return pdVec[--i*iNumRows+(--j)];
}

/* MatNxN - end */


#endif // MATVEC3N_H
