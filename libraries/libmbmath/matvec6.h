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

#ifndef MATVEC6_H
#define MATVEC6_H

#include "matvec3.h"

class Vec6
#ifdef USE_SPARSE_AUTODIFF
     :public sp_grad::SpConstMatElemAdapter<Vec6>
#endif
{
 protected:
   Vec3 v[2];
   
 public:
   Vec6(void) { 
      NO_OP; 
   };
   
   ~Vec6(void) { 
      NO_OP;
   };
 
   Vec6(const Vec6& vin) {
      v[0] = vin.v[0];
      v[1] = vin.v[1];
   };

   Vec6(const doublereal& d1, const doublereal& d2, const doublereal& d3,
	const doublereal& d4, const doublereal& d5, const doublereal& d6)
   {
      v[0] = Vec3(d1, d2, d3);
      v[1] = Vec3(d4, d5, d6);
   };
   
   Vec6(const Vec3& v1, const Vec3& v2) {
      v[0] = v1;
      v[1] = v2;
   };

   Vec6(const doublereal *pd) {
      v[0] = Vec3(pd);
      v[1] = Vec3(&pd[3]);
   };
   
   inline const Vec3& GetVec1(void) const {
      return v[0];
   };
   
   inline const Vec3& GetVec2(void) const {
      return v[1];
   };
   
   inline const Vec3& GetVec(integer i) const {
      ASSERT(i == 0 || i == 1);
      return v[i];
   };
   
   Vec3 GetVec1(void) {
      return v[0];
   };
   
   Vec3 GetVec2(void) {
      return v[1];
   };

   Vec3 GetVec(integer i) {
      ASSERT(i == 0 || i == 1);
      return v[i];
   };
   
   inline const doublereal* pGetVec(integer i) const {
      ASSERT(i == 0 || i == 1);
      return v[i].pGetVec();
   };
   
   inline Vec6& operator = (const Vec6& x) {
      v[0] = x.GetVec1();
      v[1] = x.GetVec2();
      return *this;
   };

#ifdef USE_SPARSE_AUTODIFF
     template <typename DERIVED>
     Vec6& operator = (const sp_grad::SpMatElemExprBase<doublereal, DERIVED>& v) {
          using namespace sp_grad;

          static_assert(v.iNumRowsStatic == iNumRowsStatic);
          static_assert(v.iNumColsStatic == iNumColsStatic);
          
          for (index_type i = 1; i <= iNumRowsStatic; ++i) {
               (*this)(i) = v.dGetValue(i, 1);
          }
          
          return *this;
     }
#endif
     
   inline Vec6& operator += (const Vec6& x) {
      v[0] += x.GetVec1();
      v[1] += x.GetVec2();
      return *this;
   };
   
   inline Vec6& operator -= (const Vec6& x) {
      v[0] -= x.GetVec1();
      v[1] -= x.GetVec2();
      return *this;
   };
   
   Vec6& operator *= (const doublereal& d) {
      if (d == 1.) {
	 return *this; /* No operations */
      }
      if (d == 0.) {
	 v[0].Reset(); /* Reset vector */
	 v[1].Reset();
	 return *this;
      }
      /* else */
      v[0] *= d; /* Multiply */
      v[1] *= d;
      return *this;
   };   
   
   Vec6& operator /= (const doublereal& d) {
      if (d == 1.) {
	 return *this; /* No operations */
      }
      if (d == 0.) {
	 throw ErrDivideByZero(MBDYN_EXCEPT_ARGS); /* error */
      }
      /* else */
      v[0] /= d; /* divide */
      v[1] /= d;
      return *this;
   };
   
   inline Vec6 operator + (const Vec6& x) const {
      return Vec6(v[0] + x.GetVec1(), v[1] + x.GetVec2());
   };
   
   inline Vec6 operator - (const Vec6& x) const {
      return Vec6(v[0] - x.GetVec1(), v[1] - x.GetVec2());
   };
   
   inline Vec6 operator * (const doublereal& d) const {
      return Vec6(v[0]*d, v[1]*d);
   };   
   
   inline Vec6 operator / (const doublereal& d) const {
      ASSERT(d != 0.);
      return Vec6(v[0]/d, v[1]/d);
   };   
  
   inline doublereal operator * (const Vec6& x) const {
      return v[0]*(x.v[0]) + v[1]*(x.v[1]);
   };  

   inline doublereal Dot(const Vec6& x) const {
      return v[0].Dot(x.GetVec1()) + v[1].Dot(x.GetVec2());
   };
   
   inline doublereal Dot(void) const {
      return v[0].Dot() + v[1].Dot();
   };
   
   inline doublereal Norm(void) const {
      return sqrt(v[0].Dot() + v[1].Dot());
   };

   inline const doublereal& dGet(integer i) const {
      ASSERT(i > 0 && i < 7);
      if (i < 1 || i > 6) {
	 throw ErrRowIndexOutOfRange(i, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      integer j = (i - 1)/3;
      return v[j].dGet(i - 3*j);
   };

   inline const doublereal& operator ()(integer i) const {
      ASSERT(i > 0 && i < 7);
      integer j = (i - 1)/3;
      return v[j](i - 3*j);
   };
   
   inline doublereal& operator ()(integer i) {
      ASSERT(i > 0 && i < 7);
      integer j = (i - 1)/3;
      return v[j](i - 3*j);
   };
   
   inline void Put(integer i, const doublereal& d) {
      ASSERT(i > 0 && i < 7);
      if (i < 1 || i > 6) {
	 throw ErrRowIndexOutOfRange(i, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      integer j = (i-1)/3;
      v[j].Put(i-3*j, d);
   };
   
   /*
    Scrive se stesso sull'array pd.
    Si assume che l'array pd sia lungo almeno 6 
    */
   void PutTo(doublereal* pd) const {
      ASSERT(pd != NULL);
      v[0].PutTo(pd);
      v[1].PutTo(&pd[3]);
   };   
   
   void Reset(void);

   std::ostream& Write(std::ostream& out, const char* sFill = " ") const;
     
#ifdef USE_SPARSE_AUTODIFF
     static constexpr sp_grad::index_type iNumRowsStatic = 6;
     static constexpr sp_grad::index_type iNumColsStatic = 1;
     inline constexpr sp_grad::index_type iGetNumRows() const noexcept { return iNumRowsStatic; }
     inline constexpr sp_grad::index_type iGetNumCols() const noexcept { return iNumColsStatic; }
     doublereal inline dGetValue(sp_grad::index_type i, sp_grad::index_type j) const noexcept { return (*this)(i); }     
#endif     
};


extern Vec6 operator + (const Vec6& v);
extern Vec6 operator - (const Vec6& v);
extern std::ostream& operator << (std::ostream& out, const Vec6& m);
extern std::ostream& Write(std::ostream& out, const Vec6& v, const char* sFill = " ");


class Mat6x6
#ifdef USE_SPARSE_AUTODIFF
     :public sp_grad::SpConstMatElemAdapter<Mat6x6>
#endif
{
 protected:
   Mat3x3 m[2][2];
   
 public: 
   Mat6x6(void) {
      NO_OP;
   };
   
   ~Mat6x6(void) {
      NO_OP;
   };

   Mat6x6(const Mat6x6& min) {
      m[0][0] = min.m[0][0];
      m[1][0] = min.m[1][0];
      m[0][1] = min.m[0][1];
      m[1][1] = min.m[1][1];
   };

   Mat6x6(const doublereal& d11,
	  const doublereal& d21,
	  const doublereal& d31,
	  const doublereal& d41,
	  const doublereal& d51,
	  const doublereal& d61,
	  const doublereal& d12,
	  const doublereal& d22,
	  const doublereal& d32,
	  const doublereal& d42,
	  const doublereal& d52,
	  const doublereal& d62,
	  const doublereal& d13,
	  const doublereal& d23,
	  const doublereal& d33,
	  const doublereal& d43,
	  const doublereal& d53,
	  const doublereal& d63,
	  const doublereal& d14,
	  const doublereal& d24,
	  const doublereal& d34,
	  const doublereal& d44,
	  const doublereal& d54,
	  const doublereal& d64,
	  const doublereal& d15,
	  const doublereal& d25,
	  const doublereal& d35,
	  const doublereal& d45,
	  const doublereal& d55,
	  const doublereal& d65,
	  const doublereal& d16,
	  const doublereal& d26,
	  const doublereal& d36,
	  const doublereal& d46,
	  const doublereal& d56,
	  const doublereal& d66) {
      m[0][0] = Mat3x3(d11, d21, d31, d12, d22, d32, d13, d23, d33);
      m[0][1] = Mat3x3(d14, d24, d34, d15, d25, d35, d16, d26, d36);
      m[1][0] = Mat3x3(d41, d51, d61, d42, d52, d62, d43, d53, d63);
      m[1][1] = Mat3x3(d44, d54, d64, d45, d55, d65, d46, d56, d66);
   };
   
   Mat6x6(const doublereal* pd, unsigned int i = 6) {
      ASSERT(i >= 6);
      m[0][0] = Mat3x3(*(pd+0),*(pd+1),*(pd+2),
		       *(pd+i+0),*(pd+i+1),*(pd+i+2),
		       *(pd+2*i+0),*(pd+2*i+1),*(pd+2*i+2));
      m[0][1] = Mat3x3(*(pd+3*i+0),*(pd+3*i+1),*(pd+3*i+2),
		       *(pd+4*i+0),*(pd+4*i+1),*(pd+4*i+2),
		       *(pd+5*i+0),*(pd+5*i+1),*(pd+5*i+2));
      m[1][0] = Mat3x3(*(pd+3),*(pd+4),*(pd+5),
		       *(pd+i+3),*(pd+i+4),*(pd+i+5),
		       *(pd+2*i+3),*(pd+2*i+4),*(pd+2*i+5));
      m[1][1] = Mat3x3(*(pd+3*i+3),*(pd+3*i+4),*(pd+3*i+5),
		       *(pd+4*i+3),*(pd+4*i+4),*(pd+4*i+5),
		       *(pd+5*i+3),*(pd+5*i+4),*(pd+5*i+5));
   };   
   
   Mat6x6(const Mat3x3& m11, const Mat3x3& m21, 
	  const Mat3x3& m12, const Mat3x3& m22) {
      m[0][0] = m11;
      m[1][0] = m21;
      m[0][1] = m12;
      m[1][1] = m22;
   };
   
   Mat3x3 GetMat11(void) {
      return m[0][0];
   };
   
   Mat3x3 GetMat21(void) {
      return m[1][0];
   };
   
   Mat3x3 GetMat12(void) {
      return m[0][1];
   };
   
   Mat3x3 GetMat22(void) {
      return m[1][1];
   };
   
   Mat3x3 GetMat(integer i, integer j) {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      return m[i][j];
   };
   
   inline const Mat3x3& GetMat11(void) const {
      return m[0][0];
   };
   
   inline const Mat3x3& GetMat21(void) const {
      return m[1][0];
   };
   
   inline const Mat3x3& GetMat12(void) const {
      return m[0][1];
   };
   
   inline const Mat3x3& GetMat22(void) const {
      return m[1][1];
   };
   
   inline const Mat3x3& GetMat(integer i, 
			       integer j) const {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      return m[i][j];
   };

   
   void PutMat11(const Mat3x3& x) {
      m[0][0] = x;
   };
   
   void PutMat21(const Mat3x3& x) {
      m[1][0] = x;
   };
   
   void PutMat12(const Mat3x3& x) {
      m[0][1] = x;
   };
   
   void PutMat22(const Mat3x3& x) {
      m[1][1] = x;
   };
   
   void PutMat(integer i, integer j, const Mat3x3& x) {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      m[i][j] = x;
   };
   
   
   void AddMat11(const Mat3x3& x) {
      m[0][0] += x;
   };
   
   void AddMat21(const Mat3x3& x) {
      m[1][0] += x;
   };
   
   void AddMat12(const Mat3x3& x) {
      m[0][1] += x;
   };
   
   void AddMat22(const Mat3x3& x) {
      m[1][1] += x;
   };
   
   void AddMat(integer i, integer j, const Mat3x3& x) {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      m[i][j] += x;
   };
   
   
   void SubMat11(const Mat3x3& x) {
      m[0][0] -= x;
   };
   
   void SubMat21(const Mat3x3& x) {
      m[1][0] -= x;
   };
   
   void SubMat12(const Mat3x3& x) {
      m[0][1] -= x;
   };
   
   void SubMat22(const Mat3x3& x) {
      m[1][1] -= x;
   };
   
   void SubMat(integer i, integer j, const Mat3x3& x) {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      m[i][j] -= x;
   };
   
   

   
   inline const doublereal* pGetMat(integer i, 
				    integer j) const {
      ASSERT((i == 0 || i == 1) && (j == 0 || j == 1));
      return m[i][j].pGetMat();
   };

   inline Mat6x6& operator = (const Mat6x6& x) {
      m[0][0] = x.GetMat11();
      m[1][0] = x.GetMat21();
      m[0][1] = x.GetMat12();
      m[1][1] = x.GetMat22();
      return *this;
   };

#ifdef USE_SPARSE_AUTODIFF
     template <typename DERIVED>
     Mat6x6& operator = (const sp_grad::SpMatElemExprBase<doublereal, DERIVED>& m) {
          using namespace sp_grad;

          static_assert(m.iNumRowsStatic == iNumRowsStatic);
          static_assert(m.iNumColsStatic == iNumColsStatic);

          for (index_type j = 1; j <= iNumColsStatic; ++j) {
               for (index_type i = 1; i <= iNumRowsStatic; ++i) {
                    (*this)(i, j) = m.dGetValue(i, j);
               }
          }
          
          return *this;
     }
#endif
     
   inline Mat6x6& operator += (const Mat6x6& x) {
      m[0][0] += x.GetMat11();
      m[1][0] += x.GetMat21();
      m[0][1] += x.GetMat12();
      m[1][1] += x.GetMat22();
      return *this;
   };
   
   inline Mat6x6& operator -= (const Mat6x6& x) {
      m[0][0] -= x.GetMat11();
      m[1][0] -= x.GetMat21();
      m[0][1] -= x.GetMat12();
      m[1][1] -= x.GetMat22();
      return *this;
   };
   
   Mat6x6 operator + (const Mat6x6& x) const {
      return Mat6x6(m[0][0]+x.GetMat11(), m[1][0]+x.GetMat21(),
		    m[0][1]+x.GetMat12(), m[1][1]+x.GetMat22());
   };
   
   Mat6x6 operator - (const Mat6x6& x) const {
      return Mat6x6(m[0][0]-x.GetMat11(), m[1][0]-x.GetMat21(),
		    m[0][1]-x.GetMat12(), m[1][1]-x.GetMat22());
   };

   Mat6x6 operator * (const doublereal& d) const {
      return Mat6x6(m[0][0]*d, m[1][0]*d, m[0][1]*d, m[1][1]*d);
   };   
   
   Mat6x6 operator / (const doublereal& d) const {
      ASSERT(d != 0.);
      return Mat6x6(m[0][0]/d, m[1][0]/d, m[0][1]/d, m[1][1]/d);
   };   
              
   Mat6x6& operator *= (const doublereal& d) {
      m[0][0] *= d;
      m[1][0] *= d;
      m[0][1] *= d;
      m[1][1] *= d;
      return *this;
   };   
   
   Mat6x6& operator /= (const doublereal& d) {
      ASSERT(d != 0.);
      m[0][0] /= d;
      m[1][0] /= d;
      m[0][1] /= d;
      m[1][1] /= d;
      return *this;
   };   
              
   Vec6 operator * (const Vec6& v) const {
      return Vec6(m[0][0]*v.GetVec1()+m[0][1]*v.GetVec2(),
		  m[1][0]*v.GetVec1()+m[1][1]*v.GetVec2());
   };
   
   Mat6x6 operator * (const Mat6x6& x) const {
      return Mat6x6(m[0][0]*x.GetMat11()+m[0][1]*x.GetMat21(),
		    m[1][0]*x.GetMat11()+m[1][1]*x.GetMat21(),
		    m[0][0]*x.GetMat12()+m[0][1]*x.GetMat22(),
		    m[1][0]*x.GetMat12()+m[1][1]*x.GetMat22());
   };  

   bool IsNull(void) const {
      return (m[0][0].IsNull()
           && m[0][1].IsNull()
           && m[1][0].IsNull()
           && m[1][1].IsNull());
   };
     
   bool IsExactlySame(const Mat6x6& x) const {
      return (m[0][0].IsExactlySame(x.GetMat11())
           && m[0][1].IsExactlySame(x.GetMat12())
           && m[1][0].IsExactlySame(x.GetMat21())
           && m[1][1].IsExactlySame(x.GetMat22()));
   };
 
   bool IsSame(const Mat6x6& x, const doublereal& dTol) const {
      return (m[0][0].IsSame(x.GetMat11(), dTol)
           && m[0][1].IsSame(x.GetMat12(), dTol)
           && m[1][0].IsSame(x.GetMat21(), dTol)
           && m[1][1].IsSame(x.GetMat22(), dTol));
   };

   Mat6x6 Transpose(void) {
      return Mat6x6(m[0][0].Transpose(),
		    m[0][1].Transpose(),
		    m[1][0].Transpose(),
		    m[1][1].Transpose());
   };   

   const doublereal& dGet(integer ir, integer ic) const {
      ASSERT((ir > 0 && ir < 7) && (ic > 0 && ic < 7));
      if (ir < 1 || ir > 6) {
	 throw ErrRowIndexOutOfRange(ir, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      if (ic < 1 || ic > 6) {
	 throw ErrColIndexOutOfRange(ic, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      integer jr = (ir-1)/3;
      integer jc = (ic-1)/3;
      return m[jr][jc].dGet(ir-3*jr, ic-3*jc);
   };      
   
   const doublereal& operator ()(integer ir, integer ic) const {
      ASSERT((ir > 0 && ir < 7) && (ic > 0 && ic < 7));
      integer jr = (ir - 1)/3;
      integer jc = (ic - 1)/3;
      return m[jr][jc](ir - 3*jr, ic - 3*jc);
   };
   
   doublereal& operator ()(integer ir, integer ic) {
      ASSERT((ir > 0 && ir < 7) && (ic > 0 && ic < 7));
      integer jr = (ir - 1)/3;
      integer jc = (ic - 1)/3;
      return m[jr][jc](ir - 3*jr, ic - 3*jc);
   };
   
   void Put(integer ir, integer ic, const doublereal& d) {
      ASSERT((ir > 0 && ir < 7) && (ic > 0 && ic < 7));
      if (ir < 1 || ir > 6) {
	 throw ErrRowIndexOutOfRange(ir, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      if (ic < 1 || ic > 6) {
	 throw ErrColIndexOutOfRange(ic, 1, 6, MBDYN_EXCEPT_ARGS);
      }
      integer jr = (ir-1)/3;
      integer jc = (ic-1)/3;
      m[jr][jc].Put(ir-3*jr, ic-3*jc, d);
   };      
   
   void Reset(void);

   /* Scrittura su ostream della matrice */
   std::ostream& Write(std::ostream& out, 
		  const char* sFill = " ", 
		  const char* sFill2 = NULL) const;

#ifdef USE_SPARSE_AUTODIFF
     static constexpr sp_grad::index_type iNumRowsStatic = 6;
     static constexpr sp_grad::index_type iNumColsStatic = 6;
     inline constexpr sp_grad::index_type iGetNumRows() const noexcept { return iNumRowsStatic; }
     inline constexpr sp_grad::index_type iGetNumCols() const noexcept { return iNumColsStatic; }
     doublereal inline dGetValue(sp_grad::index_type i, sp_grad::index_type j) const noexcept { return (*this)(i, j); }
#endif     
};

extern std::ostream& operator << (std::ostream& out, const Mat6x6& m);
extern std::ostream& Write(std::ostream& out,
		      const Mat6x6& m,
		      const char* sFill = " ", 
		      const char* sFill2 = NULL);


extern Vec6 MultRV(const Vec6& v, const Mat3x3& R);

extern Mat6x6 MultRM(const Mat6x6& m, const Mat3x3& R);
extern Mat6x6 MultMRt(const Mat6x6& m, const Mat3x3& R);
extern Mat6x6 MultRMRt(const Mat6x6& m, const Mat3x3& R);
extern Mat6x6 MultRMRt(const Mat6x6& m, const Mat3x3& R, const doublereal& c);


/* esegue l'operazione:
 * [I   0] [     ]
 * [     ] [  m  ]
 * [vx  I] [     ] */
extern Mat6x6 MultVCrossMat(const Mat6x6& m, const Vec3& v); 

/* esegue l'operazione:
 * [I vxT] [     ]
 * [     ] [  m  ]
 * [0   I] [     ] */
extern Mat6x6 MultVCrossTMat(const Mat6x6& m, const Vec3& v); 

/* esegue l'operazione:
 * [     ] [I  vx] 
 * [  m  ] [     ] 
 * [     ] [0   I] */
extern Mat6x6 MultMatVCross(const Mat6x6& m, const Vec3& v); 

/* esegue l'operazione:
 * [     ] [I   0] 
 * [  m  ] [     ] 
 * [     ] [vx  I] */
extern Mat6x6 MultMatVCrossT(const Mat6x6& m, const Vec3& v); 


extern const Vec6 Zero6;
extern const Mat6x6 Zero6x6;
extern const Mat6x6 Eye6;

template <>
inline const Vec6& mb_zero<Vec6>(void)
{
	return Zero6;
}

template <>
inline const Mat6x6& mb_zero<Mat6x6>(void)
{
	return Zero6x6;
}

template <>
inline Mat6x6 mb_deye<Mat6x6>(const doublereal d)
{
	// TODO: optimize
	return Mat6x6(mb_deye<Mat3x3>(d), Zero3x3, Zero3x3, mb_deye<Mat3x3>(d));
}

template <>
inline Mat6x6& mb_deye<Mat6x6>(Mat6x6& out, const doublereal d)
{
	out.PutMat11(mb_deye<Mat3x3>(d));
	out.PutMat12(Zero3x3);
	out.PutMat21(Zero3x3);
	out.PutMat22(mb_deye<Mat3x3>(d));

	return out;
}

#endif // MATVEC6_H
