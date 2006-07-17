/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2006
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
 * Copyright (C) 1996-2006
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
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include "myassert.h"

#include "ScalarFunctionsImpl.h"
#include "interp.h"

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
};


BasicScalarFunction::~BasicScalarFunction() {};

DifferentiableScalarFunction::~DifferentiableScalarFunction() {};

ConstScalarFunction::ConstScalarFunction(const doublereal v) : y(v) {
};

ConstScalarFunction::~ConstScalarFunction() {};
doublereal ConstScalarFunction::operator()(const doublereal x) const {return y;};
doublereal ConstScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in ConstScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	default:
		return 0;
		break;
	}
};

LinearScalarFunction::LinearScalarFunction(
	const doublereal y_i,
	const doublereal y_f,
	const doublereal t_i,
	const doublereal t_f) {
	ASSERTMSGBREAK(t_i != t_f, "LinearScalarFunction error, t_i == t_f");
	m=(y_f-y_i)/(t_f-t_i);
	y0=y_i-m*t_i;
};
LinearScalarFunction::~LinearScalarFunction() {};
doublereal LinearScalarFunction::operator()(const doublereal x) const {
	return (y0+m*x);
};
doublereal LinearScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in LinearScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	case 1: 
		return m;
		break;
	default:
		return 0;
		break;
	}
};

PowScalarFunction::PowScalarFunction(const doublereal p) : pw(p) {
};
PowScalarFunction::~PowScalarFunction() {};
doublereal PowScalarFunction::operator()(const doublereal x) const {
	return (pow(x,pw));
};
doublereal PowScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in PowScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	default:
		doublereal mul = 1.;
		for (integer i=0; i<order; i++) {
			mul *= (pw-i);
		}
		return mul*pow(x,pw-order);
		break;
	}
};

LogScalarFunction::LogScalarFunction(const doublereal ml) : mul_const(ml){
};
LogScalarFunction::~LogScalarFunction() {};
doublereal LogScalarFunction::operator()(const doublereal x) const {
	if (x>0.) {
		return log(x)*mul_const;
	} else {
#warning FIXME: bogus solution to the problem (avoid fp error)
		return 0.;
	}
	//will never reach this point :(
	return log(x)*mul_const;
};
doublereal LogScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in LogScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	case 1: 
		return mul_const/x;
		break;
	case 2: 
		return -mul_const/(x*x);
		break;
	case 3: 
		return 2.*mul_const/pow(x,3.);
		break;
	case 4: 
		return -6.*mul_const/pow(x,4.);
		break;
	default:
		silent_cerr("Error, LogScalarFunction::ComputeDiff "
			<< "with diff order " << order
			<< " while the maximum implemented is 4" << std::endl);
		throw ErrGeneric();
		break;
	}
};

SumScalarFunction::SumScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2
) : 
	a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)), 
	a2(ptr_cast<const DifferentiableScalarFunction*const>(b2)) {
};
SumScalarFunction::~SumScalarFunction() {};
doublereal SumScalarFunction::operator()(const doublereal x) const {
	doublereal t;
	t = a1->operator()(x);
	t += a2->operator()(x);
	return t;
};
doublereal SumScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in SumScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return operator()(x);
		break;
	default:
		return a1->ComputeDiff(x,order)+
			a2->ComputeDiff(x,order);
		break;
	}
};

MulScalarFunction::MulScalarFunction(
	const BasicScalarFunction*const b1,
	const BasicScalarFunction*const b2
) :
	a1(ptr_cast<const DifferentiableScalarFunction*const>(b1)),
	a2(ptr_cast<const DifferentiableScalarFunction*const>(b2)) {
};
MulScalarFunction::~MulScalarFunction() {};
doublereal MulScalarFunction::operator()(const doublereal x) const {
	doublereal t;
	t = a1->operator()(x);
	t *= a2->operator()(x);
	return t;
};
doublereal MulScalarFunction::ComputeDiff(const doublereal x, const integer order) const {
	ASSERTMSGBREAK(order >=0, "Error in MulScalarFunction::ComputeDiff, order<0");
	switch (order) {
	case 0: 
		return this->operator()(x);
		break;
	default:
		return a1->ComputeDiff(x,order)*a2->operator()(x)+
			a1->operator()(x)*a2->ComputeDiff(x,order);
		break;
	}
};

CubicSplineScalarFunction::CubicSplineScalarFunction(
	const std::vector<doublereal> y_i,
	const std::vector<doublereal> x_i,
	bool doNotExtrapolate)
:doNotExtrapolate(doNotExtrapolate)
{
	Y_i = y_i;
	X_i = x_i;
	ASSERTMSGBREAK(Y_i.size() == X_i.size(), "CubicSplineScalarFunction error, Y_i.size() != X_i.size()");
	std::vector<doublereal>::iterator xi,xe;
	xi = X_i.begin();
	xe = X_i.end()-1;
	for (; xi != xe; xi++) {
		if (*xi >= *(xi+1)) {
			silent_cerr("CubicSplineScalarFunction error, X_i is not ordered" << std::endl);
			throw ErrGeneric();
		}
	}
	spline(X_i,Y_i,b,c,d);
}

CubicSplineScalarFunction::~CubicSplineScalarFunction()
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
		break;
	case 1: 
		return seval(x, X_i, Y_i, b, c, d, 1);
		break;
	case 2: 
		return seval(x, X_i, Y_i, b, c, d, 2);
		break;
	case 3: 
		return seval(x, X_i, Y_i, b, c, d, 3);
		break;
	default:
		return 0.;
		break;
	}
}

MultiLinearScalarFunction::MultiLinearScalarFunction(
	const std::vector<doublereal> y_i,
	const std::vector<doublereal> x_i,
	bool doNotExtrapolate)
: doNotExtrapolate(doNotExtrapolate)
{
	Y_i = y_i;
	X_i = x_i;
	ASSERTMSGBREAK(X_i.size() == Y_i.size(), "MultiLinearScalarFunction error, Y_i.size() != X_i.size()");
	std::vector<doublereal>::iterator xi,xe;
	xi = X_i.begin();
	xe = X_i.end()-1;
	for (; xi != xe; xi++) {
		if (*xi >= *(xi+1)) {
			silent_cerr("MultiLinearScalarFunction error, X_i is not ordered" << std::endl);
			throw ErrGeneric();
		}
	}
}

MultiLinearScalarFunction::~MultiLinearScalarFunction()
{
	NO_OP;
}

doublereal MultiLinearScalarFunction::operator()(const doublereal x) const
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
		break;
	case 1: 
		return leval(x,X_i,Y_i,order);
		break;
	default:
		return 0.;
		break;
	}
}

//---------------------------------------

#include "mbpar.h"
#include "dataman.h"

typedef std::map<std::string, const ScalarFunctionRead *, ltstrcase> SFReadType;
static SFReadType SFRead;

struct SFWordSetType : public HighParser::WordSet {
	bool IsWord(const char *s) const {
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
			throw DataManager::ErrGeneric();
		}

		sf = func->second->Read(pDM, HP);
		if (!HP.SetScalarFunction(func_name, sf)) {
			silent_cerr("scalar function \"" << func_name << "\" "
				"already defined at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}

	} else if (HP.IsWord(SFWordSet)) {
		silent_cerr("Error: redefinition of "
			"\"" << func_name << "\" scalar function "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric();
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
		throw ErrGeneric();
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
			<< "aborting ..." << std::endl);
		throw DataManager::ErrGeneric();
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

struct ConstSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal c = HP.GetReal();
		return new ConstScalarFunction(c);
	};
};

struct PowSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal p = HP.GetReal();
		return new PowScalarFunction(p);
	};
};

struct LogSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal m = HP.GetReal();
		return new LogScalarFunction(m);
	};
};

struct LinearSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		doublereal y_i = HP.GetReal();
		doublereal y_f = HP.GetReal();
		doublereal t_i = HP.GetReal();
		doublereal t_f = HP.GetReal();
		return new LinearScalarFunction(y_i,y_f,t_i,t_f);
	};
};

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
		while (HP.IsArg()) {
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
		while (HP.IsArg()) {
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

struct MulSFR: public ScalarFunctionRead {
	virtual const BasicScalarFunction *
	Read(DataManager* const pDM, MBDynParser& HP) const {
		const BasicScalarFunction *const
			f1(ParseScalarFunction(HP, pDM));
		const BasicScalarFunction *const 
			f2(ParseScalarFunction(HP, pDM));
		return new MulScalarFunction(f1,f2);
	};
};

bool
SetSF(const std::string &s, const ScalarFunctionRead *rf)
{
	return SFRead.insert(SFReadType::value_type(s, rf)).second;
}

static bool done = false;

void
InitSF(void)
{
	if (::done) {
		return;
	}

	::done = true;

	SetSF("const", new ConstSFR);
	SetSF("linear", new LinearSFR);
	SetSF("pow", new PowSFR);
	SetSF("log", new LogSFR);
	SetSF("sum", new SumSFR);
	SetSF("mul", new MulSFR);
	SetSF("cubic" "spline", new CubicSplineSFR);
	SetSF("multilinear", new MultiLinearSFR);

	/* this is about initializing the scalar function drive */
	ScalarFunctionDCR *rf = new ScalarFunctionDCR;
	if (!SetDriveData("scalar" "function", rf)) {
		delete rf;

		silent_cerr("unable to register scalar function drive caller"
			<< std::endl);

		throw ErrGeneric();
	}
}

void
DestroySF(void)
{
	if (!::done) {
		return;
	}

	done = false;

	for (SFReadType::iterator i = SFRead.begin(); i != SFRead.end(); i++) {
		delete i->second;
	}
}

