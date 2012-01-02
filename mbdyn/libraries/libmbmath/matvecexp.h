/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * This code is a partial merge of HmFe and MBDyn.
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
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2012
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 * Teodoro Merlini  <merlini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef MATVECEXP_H
#define MATVECEXP_H

#include "matvec3.h"
#include "matvec6.h"

/* forward declaration */
class ScalExp;
class VecExp;
class MatExp;

class ScalExp {
protected:
	doublereal vec;
	doublereal mom;

public:
	ScalExp(void) { 
		NO_OP; 
	};
   
	~ScalExp(void) { 
		NO_OP;
	};

	ScalExp(const ScalExp& vin) {
		vec = vin.vec;
		mom = vin.mom;
	};

	ScalExp(
			const doublereal& d1, 
			const doublereal& d2 = 0.
	) {
		vec = d1;
		mom = d2;
	};
   
	inline const doublereal& GetVec(void) const {
		return vec;
	};

	inline const doublereal& GetMom(void) const {
		return mom;
	};
   
	inline void PutVec(const doublereal& x) {
		vec = x;
	};

	inline void PutMom(const doublereal& x) {
		mom = x;
	};
	
	inline const ScalExp& operator = (const ScalExp& v) {
		vec = v.vec;
		mom = v.mom;
		return *this;
	};

	inline const ScalExp& operator += (const ScalExp& v) {
		vec += v.vec;
		mom += v.mom;
		return *this;
	};

	inline const ScalExp& operator -= (const ScalExp& v) {
		vec -= v.vec;
		mom -= v.mom;
		return *this;
	};

	inline ScalExp operator + (void) const {
		return *this;
	};

	inline ScalExp operator - (void) const {
		return ScalExp(-vec, -mom);
	};
	
	inline ScalExp operator + (const ScalExp& v) const {
		return ScalExp(vec+v.vec, mom+v.mom);
	};

	inline ScalExp operator - (const ScalExp& v) const {
		return ScalExp(vec-v.vec, mom-v.mom);
	};
	
	inline ScalExp operator * (const ScalExp& v) const {
		return ScalExp(vec*v.vec,mom*v.vec+vec*v.mom);;
	};

	inline ScalExp operator / (const ScalExp& v) const {
		return ScalExp(vec/v.vec,(mom*v.vec-vec*v.mom)/(v.vec*v.vec));
	};
	

	std::ostream& 
		Write(std::ostream& out, const char* sFill = " ") const;   
};

extern ScalExp pow(const ScalExp &d, const doublereal &e);
extern ScalExp sqrt(const ScalExp &d);
extern ScalExp sin(const ScalExp &d);
extern ScalExp cos(const ScalExp &d);
extern ScalExp exp(const ScalExp &d);
//extern ScalExp operator + (const ScalExp& v);
//extern ScalExp operator - (const ScalExp& v);
extern std::ostream& operator << (std::ostream& out, const ScalExp& v);
extern std::ostream& Write(std::ostream& out, const ScalExp& v, const char* sFill = " ");



class VecExp {
protected:
	Vec3 vec;
	Vec3 mom;

public:
	VecExp(void) { 
		NO_OP; 
	};
   
	~VecExp(void) { 
		NO_OP;
	};

	VecExp(const VecExp& vin) {
		vec = vin.vec;
		mom = vin.mom;
	};

	VecExp(const Vec6& vin) {
		vec = vin.GetVec1();
		mom = vin.GetVec2();
	};

	VecExp(const doublereal& d1, const doublereal& d2,
		const doublereal& d3, const doublereal& d4, 
		const doublereal& d5, const doublereal& d6)
	{
		vec = Vec3(d1, d2, d3);
		mom = Vec3(d4, d5, d6);
	};
   
	VecExp(const Vec3& v1, const Vec3& v2) {
		vec = v1;
		mom = v2;
	};
   
	inline const Vec3& GetVec(void) const {
		return vec;
	};

	inline const Vec3& GetMom(void) const {
		return mom;
	};
   
	inline void PutVec(const Vec3& x) {
		vec = x;
	};

	inline void PutMom(const Vec3& x) {
		mom = x;
	};
	
	inline const VecExp& operator = (const VecExp& v) {
		vec = v.vec;
		mom = v.mom;
		return *this;
	};

	inline const VecExp& operator += (const VecExp& v) {
		vec += v.vec;
		mom += v.mom;
		return *this;
	};

	inline const VecExp& operator -= (const VecExp& v) {
		vec -= v.vec;
		mom -= v.mom;
		return *this;
	};

	VecExp& operator *= (const doublereal& d) {
		//M> ma si guadagna? spero di no, chi e' cosi' fessacchiotto
		//M> da moltiplicare spesso
		//M> per 1. o per zero di questi tempi? 
		//P> tutte le volte che moltiplichi per una variabile
		//P> di cui non conosci il valore (forse sono paranoico ...)
#ifdef __MBDYN_PARANOID__
		if (d == 1.) {
			return *this; /* No operations */
		}
		if (d == 0.) {
			vec = Vec3(0.); /* Reset vector */
			mom = Vec3(0.);
			return *this;
		}
		/* else */
#endif /* __MBDYN_PARANOID__ */
		vec *= d; /* Multiply */
		mom *= d;
		return *this;
	};   

	VecExp& operator /= (const doublereal& d) {
#ifdef __MBDYN_PARANOID__
		if (d == 1.) {
			return *this; /* No operations */
		}
		if (d == 0.) {
			throw ErrDivideByZero(MBDYN_EXCEPT_ARGS); /* error */
		}
		/* else */
#endif /* __MBDYN_PARANOID__ */
		vec /= d; /* divide */
		mom /= d;
		return *this;
	};

	inline VecExp operator + (const VecExp& v) const {
		return VecExp(vec+v.vec, mom+v.mom);
	};

	inline VecExp operator - (const VecExp& v) const {
		return VecExp(vec-v.vec, mom-v.mom);
	};
	
	inline VecExp operator * (const doublereal& d) const {
		return VecExp(vec*d, mom*d);
	};   
	
	inline VecExp operator / (const doublereal& d) const {
		ASSERT(d != 0.);
		return VecExp(vec/d, mom/d);
	};

	inline VecExp operator * (const ScalExp& d) const {
		return VecExp(vec*d.GetVec(), mom*d.GetVec()+vec*d.GetMom());
	};   
	
	inline VecExp operator / (const ScalExp& d) const {
		ASSERT(d.GetVec() != 0.);
		return VecExp(vec/d.GetVec(),
			(mom*d.GetVec()-vec*d.GetMom())/(d.GetVec()*d.GetVec()));
	};

	inline ScalExp operator * (const VecExp& v) const {
		return ScalExp(vec*v.vec, mom*v.vec+vec*v.mom);
	};   

	inline VecExp Cross(const VecExp &v) const {
		return VecExp(vec.Cross(v.vec),
				vec.Cross(v.mom)+mom.Cross(v.vec));
	};

	inline MatExp Cross(void) const;
	
	inline MatExp Tens(const VecExp &v) const;

	std::ostream& Write(std::ostream& out, const char* sFill = " ") const;
};

extern VecExp operator + (const VecExp& v);
extern VecExp operator - (const VecExp& v);
extern std::ostream& operator << (std::ostream& out, const VecExp& v);
extern std::ostream& Write(std::ostream& out, const VecExp& v, const char* sFill = " ");


class MatExp {
protected:
	Mat3x3 vec;
	Mat3x3 mom;

public: 
	MatExp(void) {
		NO_OP;
	};

	~MatExp(void) {
		NO_OP;
	};

	MatExp(const MatExp& min) {
		vec = min.vec;
		mom = min.mom;
	};

	MatExp(const VecExp& vin) {
		vec = Mat3x3(MatCross, vin.GetVec());
		mom = Mat3x3(MatCross, vin.GetMom());
	};

	MatExp(const VecExp& v1, const VecExp& v2) {
		vec = Mat3x3(MatCrossCross, v1.GetVec(), v2.GetVec());
		mom = Mat3x3(MatCrossCross, v1.GetMom(), v2.GetVec());
		mom += Mat3x3(MatCrossCross, v1.GetVec(), v2.GetMom());
	};

	MatExp(const doublereal& d, const VecExp& v2) {
		vec = Mat3x3(d,v2.GetVec());
		mom = Mat3x3(MatCross, v2.GetMom());
	};

	MatExp(const Mat3x3& ma, const Mat3x3& mxa) {
		vec = ma;
		mom = mxa;
	};
   
	inline const Mat3x3& GetVec(void) const {
		return vec;
	};

	inline const Mat3x3& GetMom(void) const {
		return mom;
	};

	inline void PutVec(const Mat3x3& x) {
		vec = x;
	};

	inline void PutMom(const Mat3x3& x) {
		mom = x;
	};

	inline const MatExp& operator = (const MatExp& m) {
		vec = m.vec;
		mom = m.mom;
		return *this;
	};

	inline const MatExp& operator += (const MatExp& v) {
		vec += v.vec;
		mom += v.mom;
		return *this;
	};

	inline const MatExp& operator -= (const MatExp& v) {
		vec -= v.vec;
		mom -= v.mom;
		return *this;
	};

	MatExp operator * (const doublereal& d) const {
		return MatExp(vec*d, mom*d);
	};

	MatExp operator / (const doublereal& d) const {
		ASSERT(d != 0.);
		return MatExp(vec/d, mom/d);
	};

	VecExp operator * (const VecExp& v) const {
		return VecExp(vec*v.GetVec(), vec*v.GetMom()+mom*v.GetVec());
	};
	
	MatExp operator * (const MatExp& m) const {
		return MatExp(vec*m.vec, vec*m.mom+mom*m.vec);
	};

	MatExp Transpose(void) const {
		return MatExp(vec.Transpose(), mom.Transpose());
	};

	VecExp Ax(void) const {
		return VecExp(vec.Ax(),mom.Ax());
	};

	/* Scrittura su ostream della matrice */
	std::ostream& Write(std::ostream& out, const char* sFill = " ", 
			const char* sFill2 = NULL) const;
};

inline MatExp
VecExp::Cross(void) const {
	return MatExp(Mat3x3(MatCross, vec), Mat3x3(MatCross, mom));
}


extern std::ostream& operator << (std::ostream& out, const MatExp& m);
extern std::ostream& Write(std::ostream& out, const MatExp& m, 
		const char* sFill = " ", const char* sFill2 = NULL);

inline MatExp
VecExp::Tens(const VecExp &v) const {
	return MatExp(vec.Tens(v.GetVec()), mom.Tens(v.GetMom()));
};


//M> questi no ho capito a cosa servono
//P> moltiplicano destra/sinistra/entrambi per R e R^t
#if 0
extern VecExp MultRV(const VecExp& v, const Mat3x3& R);
extern MatExp MultRM(const MatExp& m, const Mat3x3& R);
extern MatExp MultMRt(const MatExp& m, const Mat3x3& R);
extern MatExp MultRMRt(const MatExp& m, const Mat3x3& R);
#endif /* 0 */

extern const VecExp ZeroExp;
extern const MatExp EyeExp;

#endif /* MATVECEXP_H */

