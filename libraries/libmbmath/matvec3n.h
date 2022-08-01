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

#ifndef MATVEC3N_H
#define MATVEC3N_H

#include "matvec3.h"
class MatNx3;
class MatNxN;

/* VecN - begin */

/*
 * FIXME: why not derived from VectorHandler, or viceversa? 
 */

class ArrayView {
private:
	integer start, offset, number;

public:
	ArrayView(integer s, integer o, integer n)
	: start(s), offset(o), number(n) { 
		NO_OP; 
	};

	~ArrayView(void) {
		NO_OP;
	};

	inline integer Start(void) const {
		return start;
	};

	inline integer Offset(void) const {
		return offset;
	};

	inline integer Number(void) const {
		return number;
	};

	inline integer Last(void) const {
		return start + (number-1)*offset;
	};
};

class VecN {
   friend class Mat3xN;
   friend class MatNx3;
   friend class MatNxN;
   
 private:
   
   /* not defined */
   const VecN& operator = (const VecN&);
   
 protected:
   integer iMaxRows;
   integer iNumRows;
   doublereal* pdVec;

#ifdef DEBUG 
   void IsValid(void) const;   
#endif /* DEBUG */
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   VecN(void);
   VecN(integer ns);
   VecN(integer ns, const doublereal& d);

   /* costruttore copia */
   VecN(const VecN&);
 
   /* costruttore da VectorHandler (Aggiunta) */
   VecN(const VectorHandler& vh, integer ns, integer iFirstIndex);

   ~VecN(void);
   
   inline integer iGetNumRows(void) const;
   
   void Resize(integer ns);   
   void Reset(const doublereal d = 0.);
   
   inline void Put(integer i, const doublereal& d);
   inline void Add(integer i, const doublereal& d);
   inline void Sub(integer i, const doublereal& d);   
   inline const doublereal& dGet(integer i) const;

   const VecN& Copy(const VectorHandler& vh, integer iFirstIndex = 1);
 
   /*
    Dirty job: restituisce il puntatore al vettore (deprecato).
    */
   const doublereal* pGetVec(void) const { 
      return pdVec;
   };
      
   doublereal* pGetVec(void) { 
      return pdVec;
   };
      
   /* *this = n * v */ 
   void RightMult(const MatNx3& n, const Vec3& v);

   /* *this = m * n */
   const VecN& Mult(const MatNxN& m, const VecN& n); 
   const VecN& Mult(const MatNxN& m, const ArrayView& vm, 
		   const VecN& n, const ArrayView& vn);
   const VecN& operator += (const VecN& m);

   /* prodotto per scalare */
   const VecN& operator *= (const doublereal & d);

   inline doublereal& operator () (integer i);
   inline const doublereal& operator () (integer i) const;
};


inline integer VecN::iGetNumRows(void) const 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   return iNumRows;
}


inline void VecN::Put(integer i, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] = d;
}


inline void VecN::Add(integer i, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] += d;
}

inline void VecN::Sub(integer i, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   pdVec[--i] -= d;
}


inline const doublereal& VecN::dGet(integer i) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   return pdVec[--i];
}

inline doublereal &
VecN::operator () (integer i)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   return pdVec[--i];
}

inline const doublereal &
VecN::operator () (integer i) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   return pdVec[--i];
}

/* VecN - end */


/* Mat3xN - begin */

class Mat3xN: public sp_grad::SpConstMatElemAdapter<Mat3xN>
{

 friend class MatNx3;
 friend class MatNxN;

 private:
   
   /* not defined */
   Mat3xN(const Mat3xN&);
   const Mat3xN& operator = (const Mat3xN&);
   
 protected:
   integer iMaxCols;
   integer iNumCols;
   doublereal* pdRows[3];

#ifdef DEBUG
   void IsValid(void) const;
#endif /* DEBUG */
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   Mat3xN(void); /* to allow arrays of Mat3xN */
   explicit Mat3xN(integer nc); /* Attention: Mat3xN could be constructed from a doublereal if this is not explicit! */
   Mat3xN(integer nc, const doublereal& d);
   ~Mat3xN(void);
   
   void Resize(integer ns);
   void Reset(const doublereal& d = 0.);
   
   inline integer iGetNumCols(void) const;
   inline integer iGetNumRows(void) const { return 3; }
   
   inline void Put(int i, integer j, const doublereal& d);
   inline void Add(int i, integer j, const doublereal& d);
   inline void Sub(int i, integer j, const doublereal& d);
   inline const doublereal& dGet(int i, integer j) const;

   /* *this = m x *this */
   const Mat3xN& LeftMult(const Mat3x3& m);

   /* *this = m x n */
   const Mat3xN& LeftMult(const Mat3x3& m, const Mat3xN& n);

   /* *this = m * n */
   const Mat3xN& Mult(const Mat3xN& m, const MatNxN& n);
   const Mat3xN& Copy(const Mat3xN& m);
   
   const Mat3xN& operator *= (const doublereal& d);
   const Mat3xN& operator /= (const doublereal& d);

   const Mat3xN& operator += (const Mat3xN& m);
   const Mat3xN& operator -= (const Mat3xN& m);
   
   Vec3 operator * (const VecN& v) const;

   Vec3 Mult(const ArrayView& vm, const VecN& v) const;
   Vec3 Mult(const ArrayView& vm, const VecN& v, const ArrayView& vv) const;

   Vec3 GetVec(integer iCol) const;
   void PutVec(integer iCol, const Vec3& v);
   void AddVec(integer iCol, const Vec3& v);
   void SubVec(integer iCol, const Vec3& v);

   Mat3x3 GetMat3x3(integer iFirstCol) const;
   void PutMat3x3(integer iCol, const Mat3x3& m);
   void AddMat3x3(integer iCol, const Mat3x3& m);
   void SubMat3x3(integer iCol, const Mat3x3& m);
   Mat3x3 GetMat3x3ScalarMult(integer iFirstCol, const doublereal& d) const;

   inline doublereal & operator () (integer i, integer j);
   inline const doublereal & operator () (integer i, integer j) const;

     static constexpr sp_grad::index_type iNumRowsStatic = 3;
     static constexpr sp_grad::index_type iNumColsStatic = sp_grad::SpMatrixSize::DYNAMIC;
     inline sp_grad::index_type iGetRowOffset() const noexcept { return iGetNumCols(); }
     inline constexpr sp_grad::index_type iGetColOffset() const noexcept { return 1; }
     inline const doublereal* begin() const noexcept { return &pdRows[0][0]; }
     inline const doublereal* end() const noexcept { return &pdRows[0][iGetNumRows() * iGetNumCols()]; }
     doublereal inline dGetValue(sp_grad::index_type i, sp_grad::index_type j) const noexcept { return (*this)(i, j); }
};


inline integer Mat3xN::iGetNumCols(void) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   return iNumCols;
}
   
inline void Mat3xN::Put(int i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pdRows[--i][--j] = d;
}

inline void Mat3xN::Add(int i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pdRows[--i][--j] += d;
}

inline void Mat3xN::Sub(int i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   pdRows[--i][--j] -= d;
}

inline const doublereal& Mat3xN::dGet(int i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   return pdRows[--i][--j];
}

inline doublereal &
Mat3xN::operator () (integer i, integer j)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   return pdRows[--i][--j];
}
      
inline const doublereal &
Mat3xN::operator () (integer i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= 3);
   ASSERT(j > 0 && j <= iNumCols);
   return pdRows[--i][--j];
}
      
/* Mat3xN - end */

/* MatNx3 - begin */

 /* classe aggiunta per gestire operazioni con matrici Nx3. Per adesso
    sono memorizzate come tre array Nx1  */

class MatNx3: public sp_grad::SpConstMatElemAdapter<MatNx3>
{

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
   doublereal* pdCols[3];

#ifdef DEBUG 
   void IsValid(void) const;   
#endif /* DEBUG */
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   MatNx3(void);
   MatNx3(integer ns);
   MatNx3(integer ns, const doublereal& d);   
   ~MatNx3(void);
   
   const MatNx3& Copy(const MatNx3& m);

   inline integer iGetNumRows(void) const;
   void Resize(integer ns); 
   void Reset(const doublereal d = 0.);
   
   inline void Put(integer i, integer j, const doublereal& d);
   inline void Add(integer i, integer j, const doublereal& d);
   inline void Sub(integer i, integer j, const doublereal& d);   
   inline const doublereal& dGet(integer i, integer j) const;
   /* *this = n x m */
   const MatNx3& RightMult(const MatNx3& n, const Mat3x3& m);

   /* *this = [3xN]T */
   const MatNx3& Transpose(const Mat3xN& n);

   const MatNx3& operator *= (const doublereal& d);
   const MatNx3& operator /= (const doublereal& d);

   /* *this = m * n */
   const MatNx3& Mult(const MatNxN& m, const MatNx3& n);

   Vec3 GetVec(integer iRow) const;
   void PutVec(integer iRow, const Vec3& v);
   void AddVec(integer iRow, const Vec3& v);
   void SubVec(integer iRow, const Vec3& v);

   inline doublereal & operator () (integer i, integer j);
   inline const doublereal & operator () (integer i, integer j) const;

     static constexpr sp_grad::index_type iNumRowsStatic = sp_grad::SpMatrixSize::DYNAMIC;
     static constexpr sp_grad::index_type iNumColsStatic = 3;
     inline constexpr sp_grad::index_type iGetNumCols() const noexcept { return 3; }
     inline constexpr sp_grad::index_type iGetRowOffset() const noexcept { return 1; }
     inline sp_grad::index_type iGetColOffset() const noexcept { return iGetNumRows(); }
     inline const doublereal* begin() const noexcept { return &pdCols[0][0]; }
     inline const doublereal* end() const noexcept { return &pdCols[0][iGetNumRows() * iGetNumCols()]; }
     doublereal inline dGetValue(sp_grad::index_type i, sp_grad::index_type j) const noexcept { return (*this)(i, j); }
};


inline integer 
MatNx3::iGetNumRows(void) const 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   return iNumRows;
}


inline void 
MatNx3::Put(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   pdCols[--j][--i] = d;
}

inline void 
MatNx3::Add(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   pdCols[--j][--i] += d;
}

inline void 
MatNx3::Sub(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   pdCols[--j][--i] -= d;
}

inline const doublereal& 
MatNx3::dGet(integer i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   return pdCols[--j][--i];
}

inline doublereal &
MatNx3::operator () (integer i, integer j)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   return pdCols[--j][--i];
}

inline const doublereal &
MatNx3::operator () (integer i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= 3);
   return pdCols[--j][--i];
}

/* Mat Nx3 - end */

/* MatNxN - begin */
/* Questa classe memorizza matrici quadrate di dimensioni n x n in array (n^2) x 1.
   Nota: in origine l'ho messa per memorizzare le matrici di massa e rigidezza modale.
   In verita' queste matrici devono comunque essere diagonali (o non si riesce a escludere i
   modi rigidi quando si esportano i modi nel multi-corpo) quindi non dovrebbe servire piu' */

class MatNxN: public sp_grad::SpConstMatElemAdapter<MatNxN>
{
   friend class Mat3xN;
   friend class MatNx3;
   friend class VecN;
   friend std::ostream& operator << (std::ostream&, const MatNxN&);
   
 private:
   
   /* not defined */
   MatNxN(const MatNxN&);
   const MatNxN& operator = (const MatNxN&);
   
 protected:
   integer iMaxRows;
   integer iNumRows;
   doublereal* pdVec;
   doublereal** pdMat;

#ifdef DEBUG 
   void IsValid(void) const;   
#endif /* DEBUG */
   void Create_(integer ns);
   void Destroy_(void);
   
 public:
   MatNxN(void);
   MatNxN(integer ns);
   MatNxN(integer ns, const doublereal& d);   
   ~MatNxN(void);
   
   void Resize(integer ns) {Create_(ns);};
   inline integer iGetNumRows(void) const;
   inline integer iGetNumCols(void) const;
   void Reset(const doublereal d = 0.);
   inline void Put(integer i, integer j, const doublereal& d);
   inline void Add(integer i, integer j, const doublereal& d);
   inline void Sub(integer i, integer j, const doublereal& d);
   inline const doublereal& dGet(integer i, integer j) const;

   const MatNxN& operator *= (const doublereal& d);
   const MatNxN& operator /= (const doublereal& d);

   const MatNxN& Copy(const MatNxN& m);

   /* *this = m * n */
   const MatNxN& Mult(const MatNx3& m, const Mat3xN& n);

   inline doublereal& operator () (integer i, integer j);
   inline const doublereal& operator () (integer i, integer j) const;

     static constexpr sp_grad::index_type iNumRowsStatic = sp_grad::SpMatrixSize::DYNAMIC;
     static constexpr sp_grad::index_type iNumColsStatic = sp_grad::SpMatrixSize::DYNAMIC;
     inline constexpr sp_grad::index_type iGetRowOffset() const noexcept { return 1; }
     inline sp_grad::index_type iGetColOffset() const noexcept { return iGetNumRows(); }
     inline const doublereal* begin() const noexcept { return &pdMat[0][0]; }
     inline const doublereal* end() const noexcept { return &pdMat[0][iGetNumRows() * iGetNumCols()]; }
     doublereal inline dGetValue(sp_grad::index_type i, sp_grad::index_type j) const noexcept { return (*this)(i, j); }
};


inline integer MatNxN::iGetNumRows(void) const 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   return iNumRows;
}

inline integer MatNxN::iGetNumCols(void) const 
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   return iNumRows;
}

inline void MatNxN::Put(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   pdMat[--j][--i] = d;
}


inline void MatNxN::Add(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   pdMat[--j][--i] += d;
}


inline void MatNxN::Sub(integer i, integer j, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   pdMat[--j][--i] -= d;
}


inline const doublereal& 
MatNxN::dGet(integer i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   return pdMat[--j][--i];
}

inline doublereal&
MatNxN::operator () (integer i, integer j)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   return pdMat[--j][--i];
}

inline const doublereal&
MatNxN::operator () (integer i, integer j) const
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
   ASSERT(i > 0 && i <= iNumRows);
   ASSERT(j > 0 && j <= iNumRows);
   return pdMat[--j][--i];
}

std::ostream& operator << (std::ostream&, const MatNxN&);

/* MatNxN - end */


#endif /* MATVEC3N_H */

