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

#ifndef MATVECEXP_H
#define MATVECEXP_H

#include "matvec3.h"
#include "matvec6.h"

class VecExp {
 protected:
   Vec3 x;
   Vec3 g;
   
 public:
   VecExp(void) { 
      NO_OP; 
   };
   
   ~VecExp(void) { 
      NO_OP;
   };
 
   VecExp(const VecExp& vin) {
      x = vin.GetX();
      g = vin.GetG();
   };

   VecExp(const Vec6& vin) {
      x = vin.GetVec1();
      g = vin.GetVec2();
   };

   VecExp(const doublereal& d1, 
	const doublereal& d2 = 0.,
	const doublereal& d3 = 0.,
	const doublereal& d4 = 0., 
	const doublereal& d5 = 0., 
	const doublereal& d6 = 0.) {
      x = Vec3(d1, d2, d3);
      g = Vec3(d4, d5, d6);
   };
   
   VecExp(const Vec3& v1, const Vec3& v2) {
      x = v1;
      g = v2;
   };
   
   inline const Vec3& GetX(void) const {
      return x;
   };
   
   inline const Vec3& GetG(void) const {
      return g;
   };
   
   inline const Vec3& GetVec(unsigned short int i) const {
      ASSERT(i == 0 || i == 1);
      if (i == 0) {
	 return x;
      } else if (i == 1) {
	 return g;
      }
      THROW(ErrGeneric());
   };
   
   Vec3 GetX(void) {
      return x;
   };
   
   Vec3 GetG(void) {
      return g;
   };

   Vec3 GetVec(unsigned short int i) {
      ASSERT(i == 0 || i == 1);
      if (i == 0) {
	 return x;
      } else if (i == 1) {
	 return g;
      }
      THROW(ErrGeneric());
   };
   
   inline const doublereal* pGetVec(unsigned short int i) const {
      ASSERT(i == 0 || i == 1);
      if (i == 0) {
	 return x.pGetVec();
      } else if (i == 1) {
	 return g.pGetVec();
      }
      THROW(ErrGeneric());
   };
   
   inline const VecExp& operator = (const VecExp& v) {
      x = v.GetX();
      g = v.GetG();
      return *this;
   };
   
   inline const VecExp& operator += (const VecExp& v) {
      x += v.GetX();
      g += v.GetG();
      return *this;
   };
   
   inline const VecExp& operator -= (const VecExp& v) {
      x -= v.GetX();
      g -= v.GetG();
      return *this;
   };
   
   const VecExp& operator *= (const doublereal& d) {
      if (d == 1.) {
	 return *this; /* No operations */
      }
      if (d == 0.) {
	 x = Vec3(0.); /* Reset vector */
	 g = Vec3(0.);
	 return *this;
      }
      /* else */
      x *= d; /* Multiply */
      g *= d;
      return *this;
   };   
   
   const VecExp& operator /= (const doublereal& d) {
      if (d == 1.) {
	 return *this; /* No operations */
      }
      if (d == 0.) {
	 THROW(ErrDivideByZero()); /* error */
#if 0	 
	 exit(1);
#endif /* 0 */
      }
      /* else */
      x /= d; /* divide */
      g /= d;
      return *this;
   };
   
   inline VecExp operator + (const VecExp& v) const {
      return VecExp(x+v.GetX(), g+v.GetG());
   };
   
   inline VecExp operator - (const VecExp& v) const {
      return VecExp(x-v.GetX(), g-v.GetG());
   };
   
   inline VecExp operator * (const doublereal& d) const {
      return VecExp(x*d, g*d);
   };   
   
   inline VecExp operator / (const doublereal& d) const {
      ASSERT(d != 0.);
      return VecExp(x/d, g/d);
   };   
  
   inline const doublereal& dGet(unsigned short int i) const {
      ASSERT(i > 0 && i < 7);
      if (i < 1 || i > 6) {
	 THROW(ErrOutOfRange());
      }
      unsigned short int j = (i-1)/3;
      if (j == 0) {
	 return x.dGet(i-3*j);
      } else if (j == 1) {
	 return g.dGet(i-3*j);
      }
      THROW(ErrGeneric());
   };
   
   inline void Put(unsigned short int i, const doublereal& d) {
      ASSERT(i > 0 && i < 7);
      if (i < 1 || i > 6) {
	 THROW(ErrOutOfRange());
      }
      unsigned short int j = (i-1)/3;
      if (j == 0) {
	 x.Put(i-3*j, d);
      } else if (j == 1) {
	 g.Put(i-3*j, d);
      }
      THROW(ErrGeneric());
   };
   
   ostream& Write(ostream& out, const char* sFill = " ") const;   
};


extern VecExp operator + (const VecExp& v);
extern VecExp operator - (const VecExp& v);
extern ostream& operator << (ostream& out, const VecExp& v);
extern ostream& Write(ostream& out, const VecExp& v, const char* sFill = " ");


class MatExp {
 protected:
   Mat3x3 a;
   Mat3x3 xa;

 public: 
   MatExp(void) {
      NO_OP;
   };
   
   ~MatExp(void) {
      NO_OP;
   };

   MatExp(const MatExp& min) {
      a = min.a;
      xa = min.xa;
   };

   MatExp(const doublereal& d1, const doublereal& d2 = 0.) {
      a = Mat3x3(d1);
      xa = Mat3x3(d2);
   };

   MatExp(const Vec3& vx) {
      a = Eye3;
      xa = Mat3x3(vx);
   };
   
   MatExp(const Mat3x3& ma, const Mat3x3& mxa) {
      a = ma;
      xa = mxa*a;
   };
   
   MatExp(const Vec3& vx, const Mat3x3& ma) {
      a = ma;
      xa = Mat3x3(vx)*a;
   };
   
   MatExp(const Vec3& vx, const Vec3& vg) {
      a = Mat3x3(MatR, vg);
      xa = Mat3x3(vx)*a;
   };
   
   Mat3x3 GetA(void) {
      return a;
   };
   
   Mat3x3 GetXA(void) {
      return xa;
   };
   
   inline const Mat3x3& GetA(void) const {
      return a;
   };
   
   inline const Mat3x3& GetXA(void) const {
      return xa;
   };
   
   void PutA(const Mat3x3& ma) {
      a = ma;
   };
   
   void PutXA(const Vec3& vx) {
      xa = Mat3x3(vx)*a;
   };
   
   void PutXA(const Mat3x3& mxa) {
      xa = mxa;
   };
   
   void PutXA(const Vec3& vx, const Mat3x3& ma) {
      xa = Mat3x3(vx)*ma;
   };
   
   void PutXA(const Vec3& vx, const Vec3& vg) {
      a = Mat3x3(MatR, vg);
      xa = Mat3x3(vx)*a;
   };
   
   inline const MatExp& operator = (const MatExp& m) {
      a = m.GetA();
      xa = m.GetXA();
      return *this;
   };
   
   MatExp operator * (const doublereal& d) const {
      return MatExp(a*d, xa*d);
   };
   
   MatExp operator / (const doublereal& d) const {
      ASSERT(d != 0.);
      return MatExp(a/d, xa/d);
   };   
              
   VecExp operator * (const VecExp& v) const {
      return VecExp(a*v.GetX()+xa*v.GetG(), a*v.GetG());
   };
   
   MatExp operator * (const MatExp& m) const {
      return MatExp(a*m.GetA(), a*m.GetXA()+xa*m.GetA());
   };  
   
   MatExp Transpose(void) {
      return MatExp(a.Transpose(), xa.Transpose());
   };   

   /* Scrittura su ostream della matrice */
   ostream& Write(ostream& out, 
		  const char* sFill = " ", 
		  const char* sFill2 = NULL) const;
};

extern ostream& operator << (ostream& out, const MatExp& m);
extern ostream& Write(ostream& out,
		      const MatExp& m,
		      const char* sFill = " ", 
		      const char* sFill2 = NULL);


extern VecExp MultRV(const VecExp& v, const Mat3x3& R);

extern MatExp MultRM(const MatExp& m, const Mat3x3& R);
extern MatExp MultMRt(const MatExp& m, const Mat3x3& R);
extern MatExp MultRMRt(const MatExp& m, const Mat3x3& R);

extern const VecExp ZeroExp;
extern const MatExp EyeExp;

#endif /* MATVECEXP_H */

