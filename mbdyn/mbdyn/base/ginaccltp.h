/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* GiNaC constitutive law */

#ifndef GINACCLTP_H
#define GINACCLTP_H

#include <typeinfo>

#include <ginac/ginac.h>

#include "symcltp.h"

/* GiNaCElasticConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCElasticConstitutiveLaw 
: public SymbolicElasticConstitutiveLaw<T, Tder> {
private:
	unsigned dim;

	std::vector<GiNaC::symbol *> gEps;		/* parameter symbols */

	std::vector<GiNaC::ex> gExpr;			/* expressions */
	std::vector<std::vector<GiNaC::ex> > gExprDEps;	/* derivatives */

public:
	GiNaCElasticConstitutiveLaw(
		const TplDriveCaller<T>* pDC,
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& expression);
     	virtual ~GiNaCElasticConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
};

typedef GiNaCElasticConstitutiveLaw<doublereal, doublereal> 
	GiNaCElasticConstitutiveLaw1D;
typedef GiNaCElasticConstitutiveLaw<Vec3, Mat3x3>
	GiNaCElasticConstitutiveLaw3D;
typedef GiNaCElasticConstitutiveLaw<Vec6, Mat6x6> 
	GiNaCElasticConstitutiveLaw6D;

template <class T, class Tder>
GiNaCElasticConstitutiveLaw<T, Tder>::GiNaCElasticConstitutiveLaw(
	const TplDriveCaller<T>* pDC,
	const T& PStress, 
	std::vector<std::string>& epsilon,
	std::vector<std::string>& expression)
: SymbolicElasticConstitutiveLaw<T, Tder>(pDC, PStress, epsilon, expression)
{
	if (typeid(T) == typeid(Vec3)) {
		dim = 3;

	} else if (typeid(T) == typeid(Vec6)) {
		dim = 6;

	} else {
		throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(MBDYN_EXCEPT_ARGS)); 
	}

	gEps.resize(dim);
	gExpr.resize(dim);

	gExprDEps.resize(dim);
	for (unsigned row = 0; row < dim; row++) {
		gExprDEps[row].resize(dim);
	}

	ConstitutiveLaw<T, Tder>::FDE = mb_zero<Tder>();

	GiNaC::lst l; 

	for (unsigned row = 0; row < dim; row++) {
		gEps[row] = new GiNaC::symbol(epsilon[row]);
		l.append(*gEps[row]);
	}

	for (unsigned row = 0; row < dim; row++) {
		try {
			gExpr[row] = GiNaC::ex(expression[row], l);

		} catch (std::exception e) {
			silent_cerr("GiNaCElasticConstitutiveLaw<T, Tder>: expression #" << row << " parsing "
				"failed: " << e.what() << std::endl);
			throw e;
		}

		for (unsigned col = 0; col < dim; col++) {
			try {
				gExprDEps[row][col] = gExpr[row].diff(*gEps[col]);

			} catch (std::exception e) {
				silent_cerr("GiNaCElasticConstitutiveLaw<T, Tder>: expression #" << row << " differentiation "
					"wrt/ Eps #" << col << "failed: "
					<< e.what() << std::endl);
				throw e;
			}
		}
	}

	silent_cout("\tGiNaCElasticConstitutiveLaw:" << std::endl);
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tEps[" << row << "]:              \"" << *gEps[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tConstitutive law[" << row << "]: \"" << gExpr[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		for (unsigned col = 0; col < dim; col++) {
			silent_cout("\t\tDer/Eps[" << row << "][" << col << "]:          \"" << gExprDEps[row][col] << "\"" << std::endl);
		}
	}

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<T, Tder>::Epsilon, ElasticConstitutiveLaw<T, Tder>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCElasticConstitutiveLaw<T, Tder>::GiNaCElasticConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
}	
 
template <class T, class Tder>
GiNaCElasticConstitutiveLaw<T, Tder>::~GiNaCElasticConstitutiveLaw(void)
{
	for (unsigned row = 0; row < dim; row++) {
		delete gEps[row];
	}
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCElasticConstitutiveLaw<T, Tder>::pCopy(void) const
{
	ConstitutiveLaw<T, Tder>* pCL = 0;

	std::vector<std::string> epsilon(dim);
	std::vector<std::string> expression(dim);

	for (unsigned row = 0; row < dim; row++) {
		std::ostringstream	eps;
		std::ostringstream	expr;

		eps << *gEps[row];
		expr << gExpr[row];

		epsilon[row] = eps.str();
		expression[row] = expr.str();
	}

	typedef GiNaCElasticConstitutiveLaw<T, Tder> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<T, Tder>::PreStress,
			epsilon, expression));
      
	return pCL;
}

template <class T, class Tder> std::ostream& 
GiNaCElasticConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
	out << "symbolic elastic, epsilon";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << *gEps[row] << "\"";
	}
	out << "\", expression";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << gExpr[row] << "\"";
	}

  	return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
}

template <class T, class Tder> void
GiNaCElasticConstitutiveLaw<T, Tder>::Update(const T& Eps, 
	const T& /* EpsPrime */ )
{
	GiNaC::lst l;

	ElasticConstitutiveLaw<T, Tder>::Epsilon = Eps;

	T e = Eps - ElasticConstitutiveLaw<T, Tder>::Get();

	for (unsigned row = 0; row < dim; row++) {
		l.append(*gEps[row] == e(row + 1));
	}

	ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress;

	for (unsigned row = 0; row < dim; row++) {
		GiNaC::ex f_expr = gExpr[row].subs(l);

		ConstitutiveLaw<T, Tder>::F(row + 1)
			+= GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

		for (unsigned col = 0; col < dim; col++) {
			GiNaC::ex f_derEps = gExprDEps[row][col].subs(l);

			ConstitutiveLaw<T, Tder>::FDE(row + 1, col + 1)
				= GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();
		}
	}
}

/* specialize for scalar constitutive law */

template <>
class GiNaCElasticConstitutiveLaw<doublereal, doublereal>
: public SymbolicElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gExprDEps;		/* derivative */

public:
	GiNaCElasticConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC,
		const doublereal& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& expression);
     	virtual ~GiNaCElasticConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.);
};

GiNaCElasticConstitutiveLaw<doublereal, doublereal>::GiNaCElasticConstitutiveLaw(
	const TplDriveCaller<doublereal>* pDC,
	const doublereal& PStress, 
	std::vector<std::string>& epsilon,
	std::vector<std::string>& expression)
: SymbolicElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress, epsilon, expression),
gEps(epsilon[0])
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;

	GiNaC::lst l; 

	l.append(gEps);

	try {
		gExpr = GiNaC::ex(expression[0], l);

	} catch (std::exception &e) {
		silent_cerr("GiNaCElasticConstitutiveLaw<doublereal, doublereal>: expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gExprDEps = gExpr.diff(gEps);

	} catch (std::exception &e) {
		silent_cerr("GiNaCElasticConstitutiveLaw<doublereal, doublereal>: expression differentiation wrt/ Eps failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCElasticConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gExprDEps << "\"" << std::endl);

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon, ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCElasticConstitutiveLaw<doublereal, doublereal>::GiNaCElasticConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
}
 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::~GiNaCElasticConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	std::ostringstream	eps;
	std::ostringstream	expr;

	eps << gEps;
	expr << gExpr;

	std::vector<std::string> epsilon;
	std::vector<std::string> expression;

	epsilon.push_back(eps.str());
	expression.push_back(expr.str());

	typedef GiNaCElasticConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			epsilon, expression));
      
	return pCL;
}

std::ostream& 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic elastic isotropic, epsilon, \"" << gEps
		<< "\", expression, \"" << gExpr << "\"";

  	return Restart_int(out);
}

void
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::Update(const doublereal& Eps, 
	const doublereal& /* EpsPrime */ )
{
	GiNaC::lst l;

	ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

	doublereal e = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();

	l.append(gEps == e);

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_derEps = gExprDEps.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();
}

/* GiNaCElasticConstitutiveLaw - end */

/* GiNaCViscousConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCViscousConstitutiveLaw 
: public SymbolicViscousConstitutiveLaw<T, Tder> {
private:
	unsigned dim;

	std::vector<GiNaC::symbol *> gEpsPrime;			/* parameter symbols */

	std::vector<GiNaC::ex> gExpr;				/* expressions */
	std::vector<std::vector<GiNaC::ex> > gExprDEpsPrime;	/* derivatives */

public:
	GiNaCViscousConstitutiveLaw(
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscousConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
};

typedef GiNaCViscousConstitutiveLaw<doublereal, doublereal> 
	GiNaCViscousConstitutiveLaw1D;
typedef GiNaCViscousConstitutiveLaw<Vec3, Mat3x3>
	GiNaCViscousConstitutiveLaw3D;
typedef GiNaCViscousConstitutiveLaw<Vec6, Mat6x6> 
	GiNaCViscousConstitutiveLaw6D;

template <class T, class Tder>
GiNaCViscousConstitutiveLaw<T, Tder>::GiNaCViscousConstitutiveLaw(
	const T& PStress, 
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscousConstitutiveLaw<T, Tder>(PStress, epsilonPrime, expression)
{ 
	if (typeid(T) == typeid(Vec3)) {
		dim = 3;

	} else if (typeid(T) == typeid(Vec6)) {
		dim = 6;

	} else {
		throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(MBDYN_EXCEPT_ARGS)); 
	}

	gEpsPrime.resize(dim);
	gExpr.resize(dim);

	gExprDEpsPrime.resize(dim);
	for (unsigned row = 0; row < dim; row++) {
		gExprDEpsPrime[row].resize(dim);
	}

	ConstitutiveLaw<T, Tder>::FDE = mb_zero<Tder>();
	ConstitutiveLaw<T, Tder>::FDEPrime = mb_zero<Tder>();

	GiNaC::lst l; 

	for (unsigned row = 0; row < dim; row++) {
		gEpsPrime[row] = new GiNaC::symbol(epsilonPrime[row]);
		l.append(*gEpsPrime[row]);
	}

	for (unsigned row = 0; row < dim; row++) {
		try {
			gExpr[row] = GiNaC::ex(expression[row], l);

		} catch (std::exception e) {
			silent_cerr("GiNaCViscousConstitutiveLaw<T, Tder>: expression #" << row << " parsing "
				"failed: " << e.what() << std::endl);
			throw e;
		}

		for (unsigned col = 0; col < dim; col++) {
			try {
				gExprDEpsPrime[row][col] = gExpr[row].diff(*gEpsPrime[col]);

			} catch (std::exception e) {
				silent_cerr("GiNaCViscousConstitutiveLaw<T, Tder>: expression #" << row << " differentiation "
					"wrt/ EpsPrime #" << col << " failed: " << e.what()
					<< std::endl);
				throw e;
			}
		}
	}

	silent_cout("\tGiNaCViscousConstitutiveLaw:" << std::endl);
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tEpsPrime[" << row << "]:         \"" << *gEpsPrime[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tConstitutive law[" << row << "]: \"" << gExpr[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		for (unsigned col = 0; col < dim; col++) {
			silent_cout("\t\tDer/EpsPrime[" << row << "][" << col << "]:  \"" << gExprDEpsPrime[row][col] << "\"" << std::endl);
		}
	}

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<T, Tder>::Epsilon, ElasticConstitutiveLaw<T, Tder>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCViscousConstitutiveLaw<T, Tder>::GiNaCViscousConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
};

template <class T, class Tder>
GiNaCViscousConstitutiveLaw<T, Tder>::~GiNaCViscousConstitutiveLaw(void)
{
	for (unsigned row = 0; row < dim; row++) {
		delete gEpsPrime[row];
	}
}

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCViscousConstitutiveLaw<T, Tder>::pCopy(void) const
{
	ConstitutiveLaw<T, Tder>* pCL = 0;

	std::vector<std::string> epsilonPrime(dim);
	std::vector<std::string> expression(dim);

	for (unsigned row = 0; row < dim; row++) {
		std::ostringstream	epsPrime;
		std::ostringstream	expr;

		epsPrime << *gEpsPrime[row];
		expr << gExpr[row];

		epsilonPrime[row] = epsPrime.str();
		expression[row] = expr.str();
	}

	typedef GiNaCViscousConstitutiveLaw<T, Tder> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<T, Tder>::PreStress,
			epsilonPrime, expression));
      
	return pCL;
}

template <class T, class Tder> std::ostream& 
GiNaCViscousConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
	out << "symbolic viscous, epsilonPrime";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << *gEpsPrime[row] << "\"";
	}
	out << "\", expression";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << gExpr[row] << "\"";
	}

  	return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
}

template <class T, class Tder> void
GiNaCViscousConstitutiveLaw<T, Tder>::Update(const T& Eps, 
	const T& EpsPrime)
{
	GiNaC::lst l;

	ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

	for (unsigned row = 0; row < dim; row++) {
		l.append(*gEpsPrime[row] == EpsPrime(row + 1));
	}

	for (unsigned row = 0; row < dim; row++) {
		GiNaC::ex f_expr = gExpr[row].subs(l);

		ConstitutiveLaw<T, Tder>::F(row + 1)
			= GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

		for (unsigned col = 0; col < dim; col++) {
			GiNaC::ex f_derEpsPrime = gExprDEpsPrime[row][col].subs(l);

			ConstitutiveLaw<T, Tder>::FDEPrime(row + 1, col + 1)
				= GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
		}
	}
}

/* specialize for scalar constitutive law */

template <>
class GiNaCViscousConstitutiveLaw<doublereal, doublereal>
: public SymbolicViscousConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gExprDEpsPrime;		/* derivative */

public:
	GiNaCViscousConstitutiveLaw(
		const doublereal& PStress,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscousConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.);
};

GiNaCViscousConstitutiveLaw<doublereal, doublereal>::GiNaCViscousConstitutiveLaw(
	const doublereal& PStress, 
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscousConstitutiveLaw<doublereal, doublereal>(PStress, epsilonPrime, expression),
gEpsPrime(epsilonPrime[0])
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

	GiNaC::lst l; 

	l.append(gEpsPrime);

	try {
		gExpr = GiNaC::ex(expression[0], l);
	} catch (std::exception e) {
		silent_cerr("GiNaCViscousConstitutiveLaw<doublereal, doublereal>: expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gExprDEpsPrime = gExpr.diff(gEpsPrime);
	} catch (std::exception e) {
		silent_cerr("GiNaCViscousConstitutiveLaw<doublereal, doublereal>: expression differentiation wrt/ EpsPrime failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCViscousConstitutiveLaw:" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gExprDEpsPrime << "\"" << std::endl);

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon, ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCViscousConstitutiveLaw<doublereal, doublereal>::GiNaCViscousConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
}
 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::~GiNaCViscousConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	std::ostringstream	epsPrime;
	std::ostringstream	expr;

	epsPrime << gEpsPrime;
	expr << gExpr;

	std::vector<std::string> epsilonPrime;
	std::vector<std::string> expression;

	epsilonPrime.push_back(epsPrime.str());
	expression.push_back(expr.str());

	typedef GiNaCViscousConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			epsilonPrime, expression));
      
	return pCL;
}

std::ostream& 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic viscous isotropic, "
		<< "epsilon prime, \"" << gEpsPrime
		<< "\", expression, \"" << gExpr << "\"";

  	return Restart_int(out);
}

void
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::Update(
	const doublereal& /* Eps */ , 
	const doublereal& EpsPrime)
{
	GiNaC::lst l;

	ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

	l.append(gEpsPrime == EpsPrime);

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_derEpsPrime = gExprDEpsPrime.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
}

/* GiNaCViscousConstitutiveLaw - end */

/* GiNaCViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCViscoElasticConstitutiveLaw 
: public SymbolicViscoElasticConstitutiveLaw<T, Tder> {
private:
	unsigned dim;

	std::vector<GiNaC::symbol *> gEps;			/* parameter symbols */
	std::vector<GiNaC::symbol *> gEpsPrime;			/* parameter symbols */

	std::vector<GiNaC::ex> gExpr;				/* expressions */
	std::vector<std::vector<GiNaC::ex> > gExprDEps;		/* derivatives */
	std::vector<std::vector<GiNaC::ex> > gExprDEpsPrime;	/* derivatives */

public:
	GiNaCViscoElasticConstitutiveLaw(
		const TplDriveCaller<T>* pDC,
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscoElasticConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = mb_zero<T>());
};

typedef GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal> 
	GiNaCViscoElasticConstitutiveLaw1D;
typedef GiNaCViscoElasticConstitutiveLaw<Vec3, Mat3x3>
	GiNaCViscoElasticConstitutiveLaw3D;
typedef GiNaCViscoElasticConstitutiveLaw<Vec6, Mat6x6> 
	GiNaCViscoElasticConstitutiveLaw6D;

template <class T, class Tder>
GiNaCViscoElasticConstitutiveLaw<T, Tder>::GiNaCViscoElasticConstitutiveLaw(
	const TplDriveCaller<T>* pDC,
	const T& PStress, 
	std::vector<std::string>& epsilon,
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscoElasticConstitutiveLaw<T, Tder>(pDC, PStress, epsilon, epsilonPrime, expression)
{ 
	if (typeid(T) == typeid(Vec3)) {
		dim = 3;

	} else if (typeid(T) == typeid(Vec6)) {
		dim = 6;

	} else {
		throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(MBDYN_EXCEPT_ARGS)); 
	}

	gEps.resize(dim);
	gEpsPrime.resize(dim);
	gExpr.resize(dim);

	gExprDEps.resize(dim);
	gExprDEpsPrime.resize(dim);
	for (unsigned row = 0; row < dim; row++) {
		gExprDEps[row].resize(dim);
		gExprDEpsPrime[row].resize(dim);
	}

	ConstitutiveLaw<T, Tder>::FDE = mb_zero<Tder>();
	ConstitutiveLaw<T, Tder>::FDEPrime = mb_zero<Tder>();

	GiNaC::lst l; 

	for (unsigned row = 0; row < dim; row++) {
		gEps[row] = new GiNaC::symbol(epsilon[row]);
		l.append(*gEps[row]);
		gEpsPrime[row] = new GiNaC::symbol(epsilonPrime[row]);
		l.append(*gEpsPrime[row]);
	}

	for (unsigned row = 0; row < dim; row++) {
		try {
			gExpr[row] = GiNaC::ex(expression[row], l);

		} catch (std::exception e) {
			silent_cerr("GiNaCViscoElasticConstitutiveLaw<T, Tder>: expression #" << row << " parsing "
				"failed: " << e.what() << std::endl);
			throw e;
		}

		for (unsigned col = 0; col < dim; col++) {
			try {
				gExprDEps[row][col] = gExpr[row].diff(*gEps[col]);

			} catch (std::exception e) {
				silent_cerr("GiNaCViscoElasticConstitutiveLaw<T, Tder>: expression #" << row << " differentiation "
					"wrt/ Eps #" << col << "failed: "
					<< e.what() << std::endl);
				throw e;
			}

			try {
				gExprDEpsPrime[row][col] = gExpr[row].diff(*gEpsPrime[col]);

			} catch (std::exception e) {
				silent_cerr("GiNaCViscoElasticConstitutiveLaw<T, Tder>: expression #" << row << " differentiation "
					"wrt/ EpsPrime #" << col << "failed: "
					<< e.what() << std::endl);
				throw e;
			}
		}
	}

	silent_cout("\tGiNaCViscoElasticConstitutiveLaw:" << std::endl);
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tEps[" << row << "]:              \"" << *gEps[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tEpsPrime[" << row << "]:              \"" << *gEpsPrime[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		silent_cout("\t\tConstitutive law[" << row << "]: \"" << gExpr[row] << "\"" << std::endl);
	}
	for (unsigned row = 0; row < dim; row++) {
		for (unsigned col = 0; col < dim; col++) {
			silent_cout("\t\tDer/Eps[" << row << "][" << col << "]:          \"" << gExprDEps[row][col] << "\"" << std::endl);
		}
	}
	for (unsigned row = 0; row < dim; row++) {
		for (unsigned col = 0; col < dim; col++) {
			silent_cout("\t\tDer/EpsPrime[" << row << "][" << col << "]:          \"" << gExprDEpsPrime[row][col] << "\"" << std::endl);
		}
	}

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<T, Tder>::Epsilon, ElasticConstitutiveLaw<T, Tder>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCViscoElasticConstitutiveLaw<T, Tder>::GiNaCViscoElasticConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
}
 
template <class T, class Tder>
GiNaCViscoElasticConstitutiveLaw<T, Tder>::~GiNaCViscoElasticConstitutiveLaw(void)
{
	for (unsigned row = 0; row < dim; row++) {
		delete gEps[row];
		delete gEpsPrime[row];
	}
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCViscoElasticConstitutiveLaw<T, Tder>::pCopy(void) const
{
	ConstitutiveLaw<T, Tder>* pCL = 0;

	std::vector<std::string> epsilon(dim);
	std::vector<std::string> epsilonPrime(dim);
	std::vector<std::string> expression(dim);

	for (unsigned row = 0; row < dim; row++) {
		std::ostringstream	eps;
		std::ostringstream	epsPrime;
		std::ostringstream	expr;

		eps << *gEps[row];
		epsPrime << *gEpsPrime[row];
		expr << gExpr[row];

		epsilon[row] = eps.str();
		epsilonPrime[row] = epsPrime.str();
		expression[row] = expr.str();
	}

	typedef GiNaCViscoElasticConstitutiveLaw<T, Tder> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<T, Tder>::PreStress,
			epsilon, epsilonPrime, expression));
      
	return pCL;
}

template <class T, class Tder> std::ostream& 
GiNaCViscoElasticConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
	out << "symbolic viscoelastic, epsilon";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << *gEps[row] << "\"";
	}
	out << "\", epsilon prime";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << *gEpsPrime[row] << "\"";
	}
	out << "\", expression";
	for (unsigned row = 0; row < dim; row++) {
		out << ", \"" << gExpr[row] << "\"";
	}

  	return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
}

template <class T, class Tder> void
GiNaCViscoElasticConstitutiveLaw<T, Tder>::Update(const T& Eps, 
	const T& EpsPrime)
{
	GiNaC::lst l;

	ElasticConstitutiveLaw<T, Tder>::Epsilon = Eps;
	ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

	T e = Eps - ElasticConstitutiveLaw<T, Tder>::Get();

	for (unsigned row = 0; row < dim; row++) {
		l.append(*gEps[row] == e(row + 1));
		l.append(*gEpsPrime[row] == EpsPrime(row + 1));
	}

	ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress;

	for (unsigned row = 0; row < dim; row++) {
		GiNaC::ex f_expr = gExpr[row].subs(l);

		ConstitutiveLaw<T, Tder>::F(row + 1)
			+= GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

		for (unsigned col = 0; col < dim; col++) {
			GiNaC::ex f_derEps = gExprDEps[row][col].subs(l);

			ConstitutiveLaw<T, Tder>::FDE(row + 1, col + 1)
				= GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();

			GiNaC::ex f_derEpsPrime = gExprDEpsPrime[row][col].subs(l);

			ConstitutiveLaw<T, Tder>::FDEPrime(row + 1, col + 1)
				= GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
		}
	}
}

/* specialize for scalar constitutive law */

template <>
class GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>
: public SymbolicViscoElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gExprDEps;		/* derivative */
	GiNaC::ex gExprDEpsPrime;		/* derivative */

public:
	GiNaCViscoElasticConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC,
		const doublereal& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscoElasticConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.);
};

GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::GiNaCViscoElasticConstitutiveLaw(
	const TplDriveCaller<doublereal>* pDC,
	const doublereal& PStress, 
	std::vector<std::string>& epsilon,
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscoElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress, epsilon, epsilonPrime, expression),
gEps(epsilon[0]), gEpsPrime(epsilonPrime[0])
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

	GiNaC::lst l; 

	l.append(gEps);
	l.append(gEpsPrime);

	try {
		gExpr = GiNaC::ex(expression[0], l);
	} catch (std::exception e) {
		silent_cerr("GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>: expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gExprDEps = gExpr.diff(gEps);
	} catch (std::exception e) {
		silent_cerr("GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>: expression differentiation wrt/ Eps failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gExprDEpsPrime = gExpr.diff(gEpsPrime);
	} catch (std::exception e) {
		silent_cerr("GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>: expression differentiation wrt/ EpsPrime failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCViscoElasticConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gExprDEps << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gExprDEpsPrime << "\"" << std::endl);

	// try and evaluate the constitutive law
	try {
		Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon, ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon);
	}
	catch (std::exception e) {
		silent_cerr("GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::GiNaCViscoElasticConstitutiveLaw: Update() failed (" << e.what() << ")" << std::endl);
		throw e;
	}
}
 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::~GiNaCViscoElasticConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	std::ostringstream	eps;
	std::ostringstream	epsPrime;
	std::ostringstream	expr;

	eps << gEps;
	epsPrime << gEpsPrime;
	expr << gExpr;

	std::vector<std::string> epsilon;
	std::vector<std::string> epsilonPrime;
	std::vector<std::string> expression;

	epsilon.push_back(eps.str());
	epsilonPrime.push_back(epsPrime.str());
	expression.push_back(expr.str());

	typedef GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			epsilon, epsilonPrime, expression));
      
	return pCL;
}

std::ostream& 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic elastic isotropic, "
		 << "epsilon, \"" << gEps
		<< "epsilon prime, \"" << gEpsPrime
		<< "\", expression, \"" << gExpr << "\"";

  	return Restart_int(out);
}

void
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::Update(
	const doublereal& Eps, 
	const doublereal& EpsPrime)
{
	GiNaC::lst l;

	ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

	doublereal e = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();

	l.append(gEps == e);

	ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

	l.append(gEpsPrime == EpsPrime);

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_derEps = gExprDEps.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();

	GiNaC::ex f_derEpsPrime = gExprDEpsPrime.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
}

/* GiNaCViscoElasticConstitutiveLaw - end */

#endif /* GINACCLTP_H */

