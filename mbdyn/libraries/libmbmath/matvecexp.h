/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
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
class MatExp;

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

	VecExp(
			const doublereal& d1, 
			const doublereal& d2 = 0.,
			const doublereal& d3 = 0.,
			const doublereal& d4 = 0., 
			const doublereal& d5 = 0., 
			const doublereal& d6 = 0.
	) {
		vec = Vec3(d1, d2, d3);
		mom = Vec3(d4, d5, d6);
	};
   
	VecExp(const Vec3& v1, const Vec3& v2) {
		vec = v1;
		mom = v2;
	};
   
	//M> accesso in solo lettura 
	inline const Vec3& GetVec(void) const {
		return vec;
	};

	inline const Vec3& GetMom(void) const {
		return mom;
	};
   
	//M> accesso in scrittura
	//P> e' pericoloso: non sai mai quale viene usato. La cosa migliore
	//P> e' fare const Vec3& get() const per sola lettura,
	//P> e put(const Vec3&) per scrittura
	//M> come vuoi
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
			THROW(ErrDivideByZero()); /* error */
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

	//M> questo mi serve.
	inline MatExp Cross(void) const;

	ostream& Write(ostream& out, const char* sFill = " ") const;   
};

extern VecExp operator + (const VecExp& v);
extern VecExp operator - (const VecExp& v);
extern ostream& operator << (ostream& out, const VecExp& v);
extern ostream& Write(ostream& out, const VecExp& v, const char* sFill = " ");


class MatExp {
protected:
	//P> Che senso hanno questi nomi? perche' chiami "vec" una matrice
	//P> e "mom" l'altra?
	//M> residui storici a cui sono affezionato. 
	//M> quando con ste robe a 6 (VecExp) indichi delle forze
	//M> (anche se non e' questo il caso), vec e' la forza,
	//M> mom il momento della forza messa in x piu' la coppia,
	//M> mom=x.Cross(forza)+coppia
	//M> sostanzialmente ne avevo bisogno per orientarmi, adesso
	//M> se preferisci posso tranquillamente tornare a a e xa.
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

	MatExp(const doublereal& d1, const doublereal& d2 = 0.) {
		vec = Mat3x3(d1);
		mom = Mat3x3(d2);
	};

	//P> tutti questi costruttori potrebbero venire comodi
	//P> in certe circostanze
	//M> dopo aver visto come hai taroccato Rot.C penso di
	//M> capire cosa intendi
	//M> aspetta un attimo di vedere cosa serve e poi aggiungi
	//M> quello che ti sembra meglio
	//M> forse il piu' carino da fare e' MatExp(VecExp) che ti caccia
	//M> fuori quello che adesso e' VecExp.Cross()

#if 0
	MatExp(const Vec3& vx) {
		vec = Eye3;
		xa = Mat3x3(vx);
	};
#endif /* 0 */

	MatExp(const Mat3x3& ma, const Mat3x3& mxa) {
		vec = ma;
		mom = mxa;
	};
   
#if 0
	MatExp(const Vec3& vx, const Mat3x3& ma) {
		vec = ma;
		mom = Mat3x3(vx)*vec;
	};

	MatExp(const Vec3& vx, const Vec3& vg) {
		vec = Mat3x3(MatR, vg);
		xa = Mat3x3(vx)*vec;
	};


#endif /* 0 */

	inline const Mat3x3& GetVec(void) const {
		return vec;
	};

	inline const Mat3x3& GetMom(void) const {
		return mom;
	};

	//P> Come sai, io non sono d'accordo con questo metodo ...
	//M> come vuoi
	inline void PutVec(const Mat3x3& x) {
		vec = x;
	};

	inline void PutMom(const Mat3x3& x) {
		mom = x;
	};

//M> metodi assurdi, dato MatExp x; Mat3x3 y; si fa x.Vec() = y;
//M> secondo me
//M> non bisogna lavorare sulla struttura di un MatExp dove mom=x.Cross(a),
//M> questo deve essere fatto da fuori
//M> anche perche' un MatExp non e' sempre con questa struttura.
#if 0
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
#endif /* 0 */

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

	//M> questo mi serve
	VecExp Ax(void) const {
		return VecExp(vec.Ax(),mom.Ax());
	};

	/* Scrittura su ostream della matrice */
	ostream& Write(
			ostream& out,
			const char* sFill = " ",
			const char* sFill2 = NULL
	) const;
};

inline MatExp
VecExp::Cross(void) const {
	return MatExp(Mat3x3(vec), Mat3x3(mom));
}


extern ostream& operator << (ostream& out, const MatExp& m);
extern ostream& Write(
		ostream& out,
		const MatExp& m,
		const char* sFill = " ", 
		const char* sFill2 = NULL
);


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

