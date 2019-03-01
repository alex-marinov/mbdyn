/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2017
 *
 * Marco Morandini  <morandini@aero.polimi.it>
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
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include <cmath>
#include <typeinfo>

#include "myassert.h"

#include "ScalarFunctionsImpl.h"
#include "interp.h"

#include "mbpar.h"
#include "dataman.h"

//ptr_cast: I use an helper stuct to avoid troubles
//to be re-written when compilers will
//support partial specialization of functions template

//notare che non ho StlOverload.C dove metto il codice dei metodi statici,
//in teoria cosi' non e' corretto
template<class T1, class T2> struct ptr_cast_helper {
        static T1 cast(T2&);
};

template<class T1, class T2> struct ptr_cast_helper<T1*, T2> {
        static T1*  cast(T2&arg) {
                T1* out = dynamic_cast<T1*>(arg);
                if (out == 0) {
#ifdef MYDEBUG
                        std::cerr << "ptr_cast error" << std::endl;
#endif
                        throw std::bad_cast();
                }
                return out;
        };
};

template<class T1, class T2> struct ptr_cast_helper<T1*const, T2> {
        static T1*const  cast(T2&arg) {
                T1*const out = dynamic_cast<T1*const>(arg);
                if (out == 0) {
#ifdef MYDEBUG
                        std::cerr << "ptr_cast error" << std::endl;
#endif
                        throw std::bad_cast();
                }
                return out;
        };
};

template<class T1, class T2> struct ptr_cast_helper<const T1*, T2> {
        static const T1* cast(T2&arg) {
                const T1* out = dynamic_cast<const T1*>(arg);
                if (out == 0) {
#ifdef MYDEBUG
                        std::cerr << "ptr_cast error" << std::endl;
#endif
                        throw std::bad_cast();
                }
                return out;
        };
};

template<class T1, class T2> struct ptr_cast_helper<const T1*const, T2> {
        static const T1*const  cast(T2&arg) {
                const T1*const out = dynamic_cast<const T1*const>(arg);
                if (out == 0) {
#ifdef MYDEBUG
                        std::cerr << "ptr_cast error" << std::endl;
#endif
                        throw std::bad_cast();
                }
                return out;
        };
};

template<class T1,class T2>
T1 ptr_cast(T2& arg) {
        return ptr_cast_helper<T1,T2>::cast(arg);
}

// BasicScalarFunction
BasicScalarFunction::~BasicScalarFunction(void)
{
	NO_OP;
}

// DifferentiableScalarFunction
DifferentiableScalarFunction::~DifferentiableScalarFunction(void)
{
	NO_OP;
}

// ConstScalarFunction
ConstScalarFunction::ConstScalarFunction(const doublereal v)
: y(v)
{
	NO_OP;
}

ConstScalarFunction::~ConstScalarFunction(void)
{
	NO_OP;
}

doublereal
ConstScalarFunction::operator()(const doublereal x) const
{
	return y;
}

doublereal
ConstScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in ConstScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	default:
		return 0;
		break;
	}
}

// ScalarFunction parsing functional object
struct ConstSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal c = HP.GetReal();
		return new ConstScalarFunction(c);
	};
};

// LinearScalarFunction
LinearScalarFunction::LinearScalarFunction(
	const doublereal t_i,
	const doublereal y_i,
	const doublereal t_f,
	const doublereal y_f)
{
	ASSERTMSGBREAK(t_i != t_f, "LinearScalarFunction error, t_i == t_f");
	m = (y_f - y_i)/(t_f - t_i);
	y0 = y_i - m*t_i;
}

LinearScalarFunction::~LinearScalarFunction(void)
{
	NO_OP;
}

doublereal
LinearScalarFunction::operator()(const doublereal x) const
{
	return y0 + m*x;
}

doublereal
LinearScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in LinearScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);

	case 1: 
		return m;

	default:
		return 0.;

	}
}

// ScalarFunction parsing functional object
struct LinearSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal t_i = HP.GetReal();
		doublereal y_i = HP.GetReal();
		doublereal t_f = HP.GetReal();
		doublereal y_f = HP.GetReal();
		return new LinearScalarFunction(t_i,y_i,t_f,y_f);
	};
};

// PowScalarFunction
PowScalarFunction::PowScalarFunction(const doublereal p)
: pw(p)
{
	NO_OP;
}

PowScalarFunction::~PowScalarFunction(void)
{
	NO_OP;
}

doublereal
PowScalarFunction::operator()(const doublereal x) const
{
	return pow(x, pw);
}

doublereal
PowScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in PowScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);

	default:
		doublereal mul = 1.;
		for (integer i = 0; i < order; i++) {
			mul *= pw - i;
		}
		return mul*pow(x, pw - order);

	}
}

// ScalarFunction parsing functional object
struct PowSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal p = HP.GetReal();
		return new PowScalarFunction(p);
	};
};

// LogScalarFunction
LogScalarFunction::LogScalarFunction(
	const doublereal& ml,
	const doublereal& b,
	const doublereal& c)
: mul_input(ml), mul_const(ml), base(b), coef(c)
{
	if (b != 1.) {
		mul_const /= log(b);
	}
}

LogScalarFunction::~LogScalarFunction()
{
	NO_OP;
}

doublereal
LogScalarFunction::operator()(const doublereal x) const
{
	doublereal xx = coef*x;
	if (xx <= 0.) {
		silent_cerr("LogScalarFunction: argument must be positive" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return log(xx)*mul_const;
}

doublereal
LogScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >= 0, "Error in LogScalarFunction::ComputeDiff, order<0");
	if (order == 0) {
		return this->operator()(x);
	}

	if (x <= 0.) {
		silent_cerr("LogScalarFunction: argument must be positive" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	 *

 d^i                                                                mul_const
 ---  ( mul_const * log( coef * x ) ) = - ( -1 ) ^ i * ( i - 1 )! * ---------
 dx^i                                                                  x^i

	 *
	 */

	doublereal d = mul_const/x;
	for (int i = 1; i < order; i++) {
		d *= -i/x;
	}

	return d;
}

// ScalarFunction parsing functional object
struct LogSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal b = 1.;
		if (HP.IsKeyWord("base")) {
			b = HP.GetReal();
			if (b <= 0.) {
				silent_cerr("LogSFR: "
					"invalid base " << b
					<< " at line "
					<< HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal c = 1.;
		if (HP.IsKeyWord("coefficient")) {
			c = HP.GetReal();
			// note: c*x must be > 0, but x could be negative
			if (c == 0.) {
				silent_cerr("LogSFR: "
					"invalid coefficient " << c
					<< " at line "
					<< HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal m = HP.GetReal();
		return new LogScalarFunction(m, b, c);
	};
};

// ExpScalarFunction
ExpScalarFunction::ExpScalarFunction(
	const doublereal& ml,
	const doublereal& b,
	const doublereal& c)
: mul(ml), base(b), coef_input(c), coef_const(c)
{
	if (base != 1.) {
		coef_const *= log(b);
	}
}

ExpScalarFunction::~ExpScalarFunction()
{
	NO_OP;
}

doublereal
ExpScalarFunction::operator()(const doublereal x) const
{
	return exp(coef_const*x)*mul;
}

doublereal
ExpScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >= 0, "Error in ExpScalarFunction::ComputeDiff, order<0");

	/*
	 *

 d^i                                                                mul_const
 ---  ( mul * log( coef * x ) ) = - ( -1 ) ^ i * ( i - 1 )! * ---------
 dx^i                                                                  x^i

	 *
	 */

	doublereal d = 1.;
	for (int i = 0; i < order; i++) {
		d *= coef_input;
	}

	return d * this->operator()(x);
}

// ScalarFunction parsing functional object
struct ExpSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal b = 1.;
		if (HP.IsKeyWord("base")) {
			b = HP.GetReal();
			if (b <= 0.) {
				silent_cerr("ExpSFR: "
					"invalid base " << b
					<< " at line "
					<< HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal c = 1.;
		if (HP.IsKeyWord("coefficient")) {
			c = HP.GetReal();
			if (c == 0.) {
				silent_cerr("ExpSFR: "
					"invalid coefficient " << c
					<< " at line "
					<< HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal m = HP.GetReal();
		return new ExpScalarFunction(m, b, c);
	};
};

// CubicSplineScalarFunction
CubicSplineScalarFunction::CubicSplineScalarFunction(
	const std::vector<doublereal>& y_i,
	const std::vector<doublereal>& x_i,
	bool doNotExtrapolate)
: Y_i(y_i), X_i(x_i),
doNotExtrapolate(doNotExtrapolate)
{
	ASSERTMSGBREAK(Y_i.size() == X_i.size(),
		"CubicSplineScalarFunction error, Y_i.size() != X_i.size()");
	std::vector<doublereal>::iterator xi, xe;
	xi = X_i.begin();
	xe = X_i.end() - 1;
	for (unsigned i = 0; xi != xe; ++xi, ++i) {
		if (*xi >= *(xi + 1)) {
			silent_cerr("CubicSplineScalarFunction error, "
				"X is not ordered: "
				"X[" << i << "]=" << *xi
				<< " is not less than "
				"X[" << i + 1 << "]=" << *(xi + 1)
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	spline(X_i, Y_i, b, c, d);
}

CubicSplineScalarFunction::~CubicSplineScalarFunction(void)
{
	NO_OP;
}

doublereal
CubicSplineScalarFunction::operator()(const doublereal x) const
{
	if (doNotExtrapolate) {
		if (x <= X_i[0]) {
			return Y_i[0];
		}

		int s = X_i.size() - 1;
		if (x >= X_i[s]) {
			return Y_i[s];
		}
	}

	return seval(x, X_i, Y_i, b, c, d);
}

doublereal
CubicSplineScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in CubicSplineScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);

	case 1: 
		return seval(x, X_i, Y_i, b, c, d, 1);

	case 2: 
		return seval(x, X_i, Y_i, b, c, d, 2);

	case 3: 
		return seval(x, X_i, Y_i, b, c, d, 3);

	default:
		return 0.;
	}
}

// ScalarFunction parsing functional object
struct CubicSplineSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		bool doNotExtrapolate(false);
		if (HP.IsKeyWord("do" "not" "extrapolate")) {
			doNotExtrapolate = true;
		}
		std::vector<doublereal> y_i;
		std::vector<doublereal> x_i;
		y_i.resize(3);
		x_i.resize(3);
		for (int i=0; i<3; i++) {
			x_i[i] = HP.GetReal();
			y_i[i] = HP.GetReal();
		}
		while (HP.IsArg() && !HP.IsKeyWord("end")) {
			int size = x_i.size();
			x_i.resize(size+1);
			y_i.resize(size+1);
			x_i[size] = HP.GetReal();
			y_i[size] = HP.GetReal();
		}
		return new CubicSplineScalarFunction(y_i, x_i,
			doNotExtrapolate);
	};
};

// MultiLinearScalarFunction
MultiLinearScalarFunction::MultiLinearScalarFunction(
	const std::vector<doublereal>& y_i,
	const std::vector<doublereal>& x_i,
	bool doNotExtrapolate)
: Y_i(y_i), X_i(x_i),
doNotExtrapolate(doNotExtrapolate)
{
	ASSERTMSGBREAK(X_i.size() == Y_i.size(),
		"MultiLinearScalarFunction error, Y_i.size() != X_i.size()");
	std::vector<doublereal>::iterator xi, xe;
	xi = X_i.begin();
	xe = X_i.end()-1;
	for (unsigned i = 0; xi != xe; ++xi, ++i) {
		if (*xi >= *(xi + 1)) {
			silent_cerr("MultiLinearScalarFunction error, "
				"X is not ordered: "
				"X[" << i << "]=" << *xi
				<< " is not less than "
				"X[" << i + 1 << "]=" << *(xi + 1)
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

MultiLinearScalarFunction::~MultiLinearScalarFunction(void)
{
	NO_OP;
}

doublereal
MultiLinearScalarFunction::operator()(const doublereal x) const
{
	if (doNotExtrapolate) {
		if (x <= X_i[0]) {
			return Y_i[0];
		}
		
		int s = X_i.size() - 1;
		if (x >= X_i[s]) {
			return Y_i[s];
		}
	}

	return leval(x, X_i, Y_i);
}

doublereal
MultiLinearScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in MultiLinearScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return operator()(x);

	case 1: 
		return leval(x, X_i, Y_i, order);

	default:
		return 0.;
	}
}

// ScalarFunction parsing functional object
struct MultiLinearSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		bool doNotExtrapolate(false);
		if (HP.IsKeyWord("do" "not" "extrapolate")) {
			doNotExtrapolate = true;
		}
		std::vector<doublereal> y_i;
		std::vector<doublereal> x_i;
		y_i.resize(2);
		x_i.resize(2);
		for (int i=0; i<2; i++) {
			x_i[i] = HP.GetReal();
			y_i[i] = HP.GetReal();
		}
		while (HP.IsArg() && !HP.IsKeyWord("end")) {
			int size = x_i.size();
			x_i.resize(size+1);
			y_i.resize(size+1);
			x_i[size] = HP.GetReal();
			y_i[size] = HP.GetReal();
		}
		return new MultiLinearScalarFunction(y_i, x_i,
			doNotExtrapolate);
	};
};

// ChebychevScalarFunction
ChebychevScalarFunction::ChebychevScalarFunction(
	const std::vector<doublereal>& v,
	const doublereal& a, const doublereal& b, bool dne)
: vCoef(v), da(a), dfa(0.), dfap(0.), db(b), dfb(0.), dfbp(0.),
doNotExtrapolate(dne)
{
	if (!doNotExtrapolate) {
		const_cast<doublereal&>(dfa) = this->operator()(a);
		const_cast<doublereal&>(dfap) = this->ComputeDiff(a);
		const_cast<doublereal&>(dfb) = this->operator()(b);
		const_cast<doublereal&>(dfbp) = this->ComputeDiff(b);
	}
}

ChebychevScalarFunction::~ChebychevScalarFunction(void)
{
	NO_OP;
}

doublereal
ChebychevScalarFunction::operator()(const doublereal x) const
{
	if (x < da) {
		if (doNotExtrapolate) {
			silent_cerr("Chebychev interpolation: "
				"x=" << x << " is out of range "
				"[" << da << "," << db << "]" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return dfa + dfap*(x - da);

	} else if (x > db) {
		if (doNotExtrapolate) {
			silent_cerr("Chebychev interpolation: "
				"x=" << x << " is out of range "
				"[" << da << "," << db << "]" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return dfb + dfbp*(x - db);
	}

	doublereal xi = (2.*x - (da + db))/(db - da);
	doublereal d[2] = { 1., xi };
	doublereal val = vCoef[0] + vCoef[1]*xi;

	for (unsigned i = 2; i < vCoef.size(); i++) {
		doublereal Tx = 2.*xi*d[1 - i%2] - d[i%2];
		val += vCoef[i]*Tx;
		d[i%2] = Tx;
	}

	return val;
}

doublereal
ChebychevScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in ChebychevScalarFunction::ComputeDiff, order<0");

	switch (order) {
	case 0: 
		return operator()(x);

	case 1:
		break;

	default:
		silent_cerr("differentiation of order " << order << " not supported yet" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (x < da) {
		if (doNotExtrapolate) {
			silent_cerr("Chebychev interpolation: "
				"x=" << x << " is out of range "
				"[" << da << "," << db << "]" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return dfap;

	} else if (x > db) {
		if (doNotExtrapolate) {
			silent_cerr("Chebychev interpolation: "
				"x=" << x << " is out of range "
				"[" << da << "," << db << "]" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return dfbp;
	}

	doublereal xi = (2.*x - (da + db))/(db - da);
	doublereal xip = 2./(db - da);
	doublereal d[2] = { 1., xi };
	doublereal dp[2] = { 0., 1. };
	doublereal val = vCoef[1];

	for (unsigned i = 2; i < vCoef.size(); i++) {
		doublereal Tx = 2.*xi*d[1 - i%2] - d[i%2];
		doublereal Txp = 2.*d[1 - i%2] + 2.*xi*dp[1 - i%2] - dp[i%2];
		val += vCoef[i]*Txp;
		d[i%2] = Tx;
		dp[i%2] = Txp;
	}

	return xip*val;
}

// ScalarFunction parsing functional object
struct ChebychevSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal a = HP.GetReal();
		doublereal b = HP.GetReal();
		bool doNotExtrapolate(false);
		if (HP.IsKeyWord("do" "not" "extrapolate")) {
			doNotExtrapolate = true;
		}
		if (b <= a) {
			silent_cerr("Upper interval bound "
				"of Chebychev series b=" << b
				<< " must be larger than lower bound a=" << a
				<< " at line" << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			
		}
		std::vector<doublereal> v;
		while (HP.IsArg() && !HP.IsKeyWord("end")) {
			int size = v.size();
			v.resize(size+1);
			v[size] = HP.GetReal();
		}
		unsigned order = v.size();
		if (order == 0) {
			silent_cerr("Need at least one Chebychev series coefficient "
				<< "at line" << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		return new ChebychevScalarFunction(v, a, b, doNotExtrapolate);
	};
};

// SumScalarFunction
SumScalarFunction::SumScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2)
: a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)), 
a2(ptr_cast<const DifferentiableScalarFunction*const>(b2))
{
	NO_OP;
}

SumScalarFunction::~SumScalarFunction(void)
{
	NO_OP;
}

doublereal
SumScalarFunction::operator()(const doublereal x) const
{
	return a1->operator()(x) + a2->operator()(x);
}

doublereal
SumScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >=0, "Error in SumScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return operator()(x);

	default:
		return a1->ComputeDiff(x, order)
			+ a2->ComputeDiff(x, order);
	}
}

// ScalarFunction parsing functional object
struct SumSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		const BasicScalarFunction *const
			f1(ParseScalarFunction(HP, pDM));
		const BasicScalarFunction *const 
			f2(ParseScalarFunction(HP, pDM));
		return new SumScalarFunction(f1,f2);
	};
};

// SubScalarFunction
SubScalarFunction::SubScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2)
: a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)), 
a2(ptr_cast<const DifferentiableScalarFunction*const>(b2))
{
	NO_OP;
}

SubScalarFunction::~SubScalarFunction(void)
{
	NO_OP;
}

doublereal
SubScalarFunction::operator()(const doublereal x) const
{
	return a1->operator()(x) - a2->operator()(x);
}

doublereal
SubScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >= 0, "Error in SubScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return operator()(x);

	default:
		return a1->ComputeDiff(x, order)
			- a2->ComputeDiff(x, order);
	}
}

// ScalarFunction parsing functional object
struct SubSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		const BasicScalarFunction *const
			f1(ParseScalarFunction(HP, pDM));
		const BasicScalarFunction *const 
			f2(ParseScalarFunction(HP, pDM));
		return new SubScalarFunction(f1,f2);
	};
};

// MulScalarFunction
MulScalarFunction::MulScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2)
: a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)),
a2(ptr_cast<const DifferentiableScalarFunction*const>(b2))
{
	NO_OP;
}

MulScalarFunction::~MulScalarFunction()
{
	NO_OP;
}

doublereal
MulScalarFunction::operator()(const doublereal x) const
{
	return a1->operator()(x)*a2->operator()(x);
}

doublereal
MulScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	ASSERTMSGBREAK(order >= 0, "Error in MulScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);

	default:
		return a1->ComputeDiff(x, order)*a2->operator()(x)
			+ a1->operator()(x)*a2->ComputeDiff(x, order);
	}
}

// ScalarFunction parsing functional object
struct MulSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		const BasicScalarFunction *const
			f1(ParseScalarFunction(HP, pDM));
		const BasicScalarFunction *const 
			f2(ParseScalarFunction(HP, pDM));
		return new MulScalarFunction(f1, f2);
	};
};

// DivScalarFunction
DivScalarFunction::DivScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2)
: a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)),
a2(ptr_cast<const DifferentiableScalarFunction*const>(b2))
{
	NO_OP;
}

DivScalarFunction::~DivScalarFunction()
{
	NO_OP;
}

doublereal
DivScalarFunction::operator()(const doublereal x) const
{
	doublereal n, d;
	d = a2->operator()(x);
	if (d == 0) {
		/* TODO: cleanup exception handling */
		silent_cerr("DivScalarFunction: division by zero" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	n = a1->operator()(x);
	return n/d;
}

doublereal
DivScalarFunction::ComputeDiff(const doublereal x, const integer order) const
{
	doublereal d;
	ASSERTMSGBREAK(order >= 0, "Error in DivScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);

	default:
		d = a2->operator()(x);
		if (d == 0.) {
			/* TODO: cleanup exception handling */
			silent_cerr("DivScalarFunction: division by zero" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		return a1->ComputeDiff(x, order)/d
			- a1->operator()(x)/(d*d)*a2->ComputeDiff(x, order);
	}
}

// ScalarFunction parsing functional object
struct DivSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		const BasicScalarFunction *const
			f1(ParseScalarFunction(HP, pDM));
		const BasicScalarFunction *const 
			f2(ParseScalarFunction(HP, pDM));
		return new DivScalarFunction(f1, f2);
	};
};

//---------------------------------------

typedef std::map<std::string, const ScalarFunctionRead *, ltstrcase> SFReadType;
static SFReadType SFRead;

struct SFWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return SFRead.find(std::string(s)) != SFRead.end();
	};
};
static SFWordSetType SFWordSet;

const BasicScalarFunction *const
ParseScalarFunction(MBDynParser& HP, DataManager* const pDM)
{
	std::string func_name(HP.GetStringWithDelims());
	
	const BasicScalarFunction *sf = HP.GetScalarFunction(func_name);
	if (sf == 0) {
		const char *s = HP.IsWord(SFWordSet);
		if (s == 0) {
			s = "const";
		}

		SFReadType::iterator func = SFRead.find(std::string(s));
		if (func == SFRead.end()) {
			silent_cerr("unknown scalar function type \"" << s << "\" "
				"for function \"" << func_name << "\" "
				"at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		try {
			sf = func->second->Read(pDM, HP);
		} catch (...) {
			silent_cerr("Unable to parse "
				"ScalarFunction(\"" << func_name << "\") "
				"at line " << HP.GetLineData() << std::endl);
			throw;
		}
		if (!HP.SetScalarFunction(func_name, sf)) {
			silent_cerr("scalar function \"" << func_name << "\" "
				"already defined at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (HP.IsWord(SFWordSet)) {
		silent_cerr("Error: redefinition of "
			"\"" << func_name << "\" scalar function "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	return sf;
}

//---------------------------------------

/* Scalar FunctionDrive - begin */

class ScalarFunctionDriveCaller : public DriveCaller {
private:
	const BasicScalarFunction *const sc;

public:
	ScalarFunctionDriveCaller(const DriveHandler* pDH,
    		const BasicScalarFunction *const f)
	: DriveCaller(pDH), sc(f)
	{ NO_OP; };

	virtual ~ScalarFunctionDriveCaller(void) { NO_OP; };

	/* Copia */
	virtual DriveCaller*
	pCopy(void) const {
		DriveCaller* pDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pDC, ScalarFunctionDriveCaller,
			ScalarFunctionDriveCaller(pDrvHdl, sc));
		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */	
	virtual std::ostream& Restart(std::ostream& out) const {
		silent_cerr("ScalarFunctionDriveCaller: Restart not implemented"
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

	virtual doublereal dGet(const doublereal& dVar) const {
		return (*sc)(dVar);
	};

	virtual bool bIsDifferentiable(void) const {
		return (ptr_cast<const DifferentiableScalarFunction*const>(sc) != 0);
	};

	virtual doublereal dGetP(const doublereal& dVar) const {
		return ptr_cast<const DifferentiableScalarFunction *const>(sc)->ComputeDiff(dVar);
	};
};

struct ScalarFunctionDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ScalarFunctionDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "scalar function");

	/* driver legato alle scalar function */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no driver dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const DriveHandler* pDrvHdl = pDM->pGetDrvHdl();

	const BasicScalarFunction *const sc = ParseScalarFunction(HP, (DataManager *const)pDM);

	DriveCaller *pDC = 0;

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		ScalarFunctionDriveCaller,
		ScalarFunctionDriveCaller(pDrvHdl, sc));

	return pDC;
}

// ScalarFunction constitutive laws

// ScalarFunction elastic isotropic constitutive law

template <class T, class Tder>
class ScalarFunctionIsotropicCL
: public ConstitutiveLaw<T, Tder> {
private:
	const DifferentiableScalarFunction * pSF;
	int n;

public:
	ScalarFunctionIsotropicCL(const DifferentiableScalarFunction * psf)
	: pSF(psf) {
		if (typeid(T) == typeid(Vec3)) {
			n = 3;

		} else if (typeid(T) == typeid(Vec6)) {
			n = 6;

		} else {
			silent_cerr("ScalarFunctionIsotropicCL<" << typeid(T).name() << ", " << typeid(Tder).name() << "> not implemented" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	virtual ~ScalarFunctionIsotropicCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef ScalarFunctionIsotropicCL<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(pSF));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "# not implemented!";
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */  = 0.) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		for (int i = 1; i <= n; i++) {
#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
			ConstitutiveLaw<T, Tder>::F.Put(i, (*pSF)(Eps(i)));
			ConstitutiveLaw<T, Tder>::FDE.Put(i, i, pSF->ComputeDiff(Eps(i)));
#else // !MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
			ConstitutiveLaw<T, Tder>::F(i) = (*pSF)(Eps(i));
			ConstitutiveLaw<T, Tder>::FDE(i, i) = pSF->ComputeDiff(Eps(i));
#endif // !MBDYN_X_WORKAROUND_GCC_3_3 && ! MBDYN_X_WORKAROUND_GCC_3_3
		}
	};
};

template <>
class ScalarFunctionIsotropicCL<doublereal, doublereal>
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	const DifferentiableScalarFunction *pSF;

public:
	ScalarFunctionIsotropicCL(const DifferentiableScalarFunction * psf)
	: pSF(psf) {
		NO_OP;
	};

	virtual ~ScalarFunctionIsotropicCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef ScalarFunctionIsotropicCL<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(pSF));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "# not implemented!";
	};

	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */  = 0.) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;
		ConstitutiveLaw<doublereal, doublereal>::F = (*pSF)(Eps);
		ConstitutiveLaw<doublereal, doublereal>::FDE = pSF->ComputeDiff(Eps);
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct ScalarFunctionIsotropicCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
};

template <class T, class Tder>
ConstitutiveLaw<T, Tder> *
ScalarFunctionIsotropicCLR<T, Tder>::Read(const DataManager* pDM,
	MBDynParser& HP,
	ConstLawType::Type& CLType)
{
	ConstitutiveLaw<T, Tder>* pCL = 0;

	CLType = ConstLawType::ELASTIC;

	int n = 0;
	if (typeid(T) == typeid(doublereal)) {
		n = 1;

	} else if (typeid(T) == typeid(Vec3)) {
		n = 3;

	} else if (typeid(T) == typeid(Vec6)) {
		n = 6;

	} else {
		silent_cerr("ScalarFunctionIsotropicCL"
			"<" << typeid(T).name() << ", " << typeid(Tder).name() << "> "
			"not implemented" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	const BasicScalarFunction *pSF = ParseScalarFunction(HP, (DataManager *const)pDM);
	const DifferentiableScalarFunction *psf = dynamic_cast<const DifferentiableScalarFunction *>(pSF);
	if (psf == 0) {
		silent_cerr("ScalarFunction must be differentiable "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	typedef ScalarFunctionIsotropicCL<T, Tder> L;
	SAFENEWWITHCONSTRUCTOR(pCL, L, L(psf));

	return pCL;
}

// ScalarFunction elastic orthotropic constitutive law

template <class T, class Tder>
class ScalarFunctionOrthotropicCL
: public ConstitutiveLaw<T, Tder> {
private:
	std::vector<const DifferentiableScalarFunction *> SF;
	unsigned n;

public:
	ScalarFunctionOrthotropicCL(const std::vector<const DifferentiableScalarFunction *>& sf)
	{
		if (typeid(T) == typeid(Vec3)) {
			n = 3;

		} else if (typeid(T) == typeid(Vec6)) {
			n = 6;

		} else {
			silent_cerr("ScalarFunctionOrthotropicCL<" << typeid(T).name() << ", " << typeid(Tder).name() << "> not implemented" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(sf.size() == n);
		SF.resize(n);
		for (unsigned i = 0; i < n; i++) {
			SF[i] = sf[i];
		}
	};

	virtual ~ScalarFunctionOrthotropicCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef ScalarFunctionOrthotropicCL<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(SF));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "# not implemented!";
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */  = 0.) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		for (unsigned i = 1; i <= n; i++) {
			/* skip null scalar functions */
			if (SF[i - 1] == 0) {
				continue;
			}

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
			ConstitutiveLaw<T, Tder>::F.Put(i, (*SF[i - 1])(Eps(i)));
			ConstitutiveLaw<T, Tder>::FDE.Put(i, i, SF[i - 1]->ComputeDiff(Eps(i)));
#else // !MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
			ConstitutiveLaw<T, Tder>::F(i) = (*SF[i - 1])(Eps(i));
			ConstitutiveLaw<T, Tder>::FDE(i, i) = SF[i - 1]->ComputeDiff(Eps(i));
#endif // !MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
		}
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct ScalarFunctionOrthotropicCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
};

template <class T, class Tder>
ConstitutiveLaw<T, Tder> *
ScalarFunctionOrthotropicCLR<T, Tder>::Read(const DataManager* pDM,
	MBDynParser& HP,
	ConstLawType::Type& CLType)
{
	ConstitutiveLaw<T, Tder>* pCL = 0;

	CLType = ConstLawType::ELASTIC;

	int n = 1;
	if (typeid(T) == typeid(Vec3)) {
		n = 3;

	} else if (typeid(T) == typeid(Vec6)) {
		n = 6;

	} else if (typeid(T) != typeid(doublereal)) {
		silent_cerr("ScalarFunctionOrthotropicCL<" << typeid(T).name() << ", " << typeid(Tder).name() << "> not implemented" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	std::vector<const DifferentiableScalarFunction *> SF(n);
	for (int i = 0; i < n; i++) {
		if (HP.IsKeyWord("null")) {
			SF[i] = 0;

		} else {
			const BasicScalarFunction *pSF = ParseScalarFunction(HP, (DataManager *const)pDM);
			SF[i] = dynamic_cast<const DifferentiableScalarFunction *>(pSF);
			if (SF[i] == 0) {
				silent_cerr("ScalarFunction #" << i + 1 << " must be differentiable "
					"at line " << HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	typedef ScalarFunctionOrthotropicCL<T, Tder> L;
	SAFENEWWITHCONSTRUCTOR(pCL, L, L(SF));

	return pCL;
}

bool
SetSF(const std::string &s, const ScalarFunctionRead *rf)
{
	pedantic_cout("registering scalar function \"" << s << "\""
		<< std::endl );
	return SFRead.insert(SFReadType::value_type(s, rf)).second;
}

static unsigned done = 0;

void
InitSF(void)
{
	if (::done++ > 0) {
		return;
	}

	bool b;

	b = SetSF("const", new ConstSFR);
	ASSERT(b);
	b = SetSF("linear", new LinearSFR);
	ASSERT(b);
	b = SetSF("pow", new PowSFR);
	ASSERT(b);
	b = SetSF("log", new LogSFR);
	ASSERT(b);
	b = SetSF("exp", new ExpSFR);
	ASSERT(b);
	b = SetSF("sum", new SumSFR);
	ASSERT(b);
	b = SetSF("sub", new SubSFR);
	ASSERT(b);
	b = SetSF("mul", new MulSFR);
	ASSERT(b);
	b = SetSF("div", new DivSFR);
	ASSERT(b);
	b = SetSF("cubic" "spline", new CubicSplineSFR);
	ASSERT(b);
	b = SetSF("multilinear", new MultiLinearSFR);
	ASSERT(b);
	b = SetSF("chebychev", new ChebychevSFR);
	ASSERT(b);

	/* this is about initializing the scalar function drive */
	ScalarFunctionDCR *rf = new ScalarFunctionDCR;
	if (!SetDriveCallerData("scalar" "function", rf)) {
		delete rf;

		silent_cerr("unable to register scalar function drive caller"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* this is about initializing the scalar function constitutive law(s) */
	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new ScalarFunctionIsotropicCLR<doublereal, doublereal>;
	if (!SetCL1D("scalar" "function" "elastic" "isotropic", rf1D)) {
		delete rf1D;

		silent_cerr("unable to register scalar function isotropic 1D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rf1D = new ScalarFunctionIsotropicCLR<doublereal, doublereal>;
	if (!SetCL1D("scalar" "function" "elastic" "orthotropic", rf1D)) {
		delete rf1D;

		silent_cerr("unable to register scalar function orthotropic 1D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rf1D = new ScalarFunctionIsotropicCLR<doublereal, doublereal>;
	if (!SetCL1D("scalar" "function" "elastic", rf1D)) {
		delete rf1D;

		silent_cerr("unable to register scalar function 1D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new ScalarFunctionIsotropicCLR<Vec3, Mat3x3>;
	if (!SetCL3D("scalar" "function" "elastic" "isotropic", rf3D)) {
		delete rf3D;

		silent_cerr("unable to register scalar function isotropic 3D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rf3D = new ScalarFunctionOrthotropicCLR<Vec3, Mat3x3>;
	if (!SetCL3D("scalar" "function" "elastic" "orthotropic", rf3D)) {
		delete rf3D;

		silent_cerr("unable to register scalar function orthotropic 3D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ConstitutiveLawRead<Vec6, Mat6x6> *rf6D = new ScalarFunctionIsotropicCLR<Vec6, Mat6x6>;
	if (!SetCL6D("scalar" "function" "elastic" "isotropic", rf6D)) {
		delete rf6D;

		silent_cerr("unable to register scalar function isotropic 6D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	rf6D = new ScalarFunctionOrthotropicCLR<Vec6, Mat6x6>;
	if (!SetCL6D("scalar" "function" "elastic" "orthotropic", rf6D)) {
		delete rf6D;

		silent_cerr("unable to register scalar function orthotropic 6D constitutive law"
			<< std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
DestroySF(void)
{
	if (::done == 0) {
		silent_cerr("DestroySF() called once too many" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--done > 0) {
		return;
	}

	for (SFReadType::iterator i = SFRead.begin(); i != SFRead.end(); ++i) {
		delete i->second;
	}
	SFRead.clear();
}

