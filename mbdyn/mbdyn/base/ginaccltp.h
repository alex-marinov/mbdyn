/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifdef HAVE_GINAC

#include "symcltp.h"
#include <ginac/ginac.h>
#include <typeinfo>

/* GiNaCElasticConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCElasticConstitutiveLaw 
: public SymbolicElasticConstitutiveLaw<T, Tder> {
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
	virtual void IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ = 0.);
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
	throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(std::cerr,
		"GiNaCElasticConstitutiveLaw() is not defined "
		"for the requested dimensionality")); 
}
 
template <class T, class Tder>
GiNaCElasticConstitutiveLaw<T, Tder>::~GiNaCElasticConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCElasticConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return 0;
}

template <class T, class Tder> std::ostream& 
GiNaCElasticConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
GiNaCElasticConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
GiNaCElasticConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

/* specialize for scalar constitutive law */

template <>
class GiNaCElasticConstitutiveLaw<doublereal, doublereal>
: public SymbolicElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEps;		/* derivative */

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
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& /* EpsPrime */ = 0.);
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
	} catch (std::exception e) {
		silent_cerr("expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gDerEps = gExpr.diff(gEps);
	} catch (std::exception e) {
		silent_cerr("expression differentiation wrt/ Eps failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCElasticConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gDerEps << "\"" << std::endl);
}
 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::~GiNaCElasticConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

#if defined(HAVE_SSTREAM)
	std::ostringstream	eps;
	std::ostringstream	expr;
#else /* HAVE_STRSTREAM_H */
	ostrstream		eps;
	ostrstream		expr;
#endif /* HAVE_STRSTREAM_H */

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

	GiNaC::ex f_derEps = gDerEps.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();
}

void 
GiNaCElasticConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(
	const doublereal& DeltaEps, 
	const doublereal& /* EpsPrime */ )
{
	Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon + DeltaEps);
}

/* GiNaCElasticConstitutiveLaw - end */

/* GiNaCViscousConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCViscousConstitutiveLaw 
: public SymbolicViscousConstitutiveLaw<T, Tder> {
public:
	GiNaCViscousConstitutiveLaw(
		const TplDriveCaller<T>* pDC,
		const T& PStress,
		std::vector<std::string>& epsilon,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscousConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ = 0.);
};

typedef GiNaCViscousConstitutiveLaw<doublereal, doublereal> 
	GiNaCViscousConstitutiveLaw1D;
typedef GiNaCViscousConstitutiveLaw<Vec3, Mat3x3>
	GiNaCViscousConstitutiveLaw3D;
typedef GiNaCViscousConstitutiveLaw<Vec6, Mat6x6> 
	GiNaCViscousConstitutiveLaw6D;

template <class T, class Tder>
GiNaCViscousConstitutiveLaw<T, Tder>::GiNaCViscousConstitutiveLaw(
	const TplDriveCaller<T>* pDC,
	const T& PStress, 
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscousConstitutiveLaw<T, Tder>(pDC, PStress, epsilonPrime, expression)
{ 
	throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(std::cerr,
		"GiNaCElasticConstitutiveLaw() is not defined "
		"for the requested dimensionality")); 
}
 
template <class T, class Tder>
GiNaCViscousConstitutiveLaw<T, Tder>::~GiNaCViscousConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCViscousConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return 0;
}

template <class T, class Tder> std::ostream& 
GiNaCViscousConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
GiNaCViscousConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
GiNaCViscousConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

/* specialize for scalar constitutive law */

template <>
class GiNaCViscousConstitutiveLaw<doublereal, doublereal>
: public SymbolicViscousConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEpsPrime;		/* derivative */

public:
	GiNaCViscousConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC,
		const doublereal& PStress,
		std::vector<std::string>& epsilonPrime,
		std::vector<std::string>& expression);
     	virtual ~GiNaCViscousConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.);
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& EpsPrime = 0.);
};

GiNaCViscousConstitutiveLaw<doublereal, doublereal>::GiNaCViscousConstitutiveLaw(
	const TplDriveCaller<doublereal>* pDC,
	const doublereal& PStress, 
	std::vector<std::string>& epsilonPrime,
	std::vector<std::string>& expression)
: SymbolicViscousConstitutiveLaw<doublereal, doublereal>(pDC, PStress, epsilonPrime, expression),
gEpsPrime(epsilonPrime[0])
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

	GiNaC::lst l; 

	l.append(gEpsPrime);

	try {
		gExpr = GiNaC::ex(expression[0], l);
	} catch (std::exception e) {
		silent_cerr("expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gDerEpsPrime = gExpr.diff(gEpsPrime);
	} catch (std::exception e) {
		silent_cerr("expression differentiation wrt/ EpsPrime failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCViscousConstitutiveLaw:" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gDerEpsPrime << "\"" << std::endl);
}
 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::~GiNaCViscousConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

#if defined(HAVE_SSTREAM)
	std::ostringstream	epsPrime;
	std::ostringstream	expr;
#else /* HAVE_STRSTREAM_H */
	ostrstream		epsPrime;
	ostrstream		expr;
#endif /* HAVE_STRSTREAM_H */

	epsPrime << gEpsPrime;
	expr << gExpr;

	std::vector<std::string> epsilonPrime;
	std::vector<std::string> expression;

	epsilonPrime.push_back(epsPrime.str());
	expression.push_back(expr.str());

	typedef GiNaCViscousConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
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

	GiNaC::ex f_derEpsPrime = gDerEpsPrime.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
}

void 
GiNaCViscousConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(const doublereal& /* DeltaEps */ ,
		const doublereal& EpsPrime)
{
	Update(0, ElasticConstitutiveLaw<doublereal, doublereal>::EpsilonPrime + EpsPrime);
}

/* GiNaCViscousConstitutiveLaw - end */

/* GiNaCViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class GiNaCViscoElasticConstitutiveLaw 
: public SymbolicViscoElasticConstitutiveLaw<T, Tder> {
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
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ = 0.);
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
	throw (typename ConstitutiveLaw<T, Tder>::ErrNotAvailable(std::cerr,
		"GiNaCElasticConstitutiveLaw() is not defined "
		"for the requested dimensionality")); 
}
 
template <class T, class Tder>
GiNaCViscoElasticConstitutiveLaw<T, Tder>::~GiNaCViscoElasticConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
GiNaCViscoElasticConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return 0;
}

template <class T, class Tder> std::ostream& 
GiNaCViscoElasticConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
GiNaCViscoElasticConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
GiNaCViscoElasticConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& EpsPrime )
{
	NO_OP;
}

/* specialize for scalar constitutive law */

template <>
class GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>
: public SymbolicViscoElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEps;		/* derivative */
	GiNaC::ex gDerEpsPrime;		/* derivative */

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
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& EpsPrime = 0.);
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
		silent_cerr("expression parsing failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gDerEps = gExpr.diff(gEps);
	} catch (std::exception e) {
		silent_cerr("expression differentiation wrt/ Eps failed: " << e.what() << std::endl);
		throw e;
	}

	try {
		gDerEpsPrime = gExpr.diff(gEpsPrime);
	} catch (std::exception e) {
		silent_cerr("expression differentiation wrt/ EpsPrime failed: " << e.what() << std::endl);
		throw e;
	}

	silent_cout("\tGiNaCViscoElasticConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gDerEps << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gDerEpsPrime << "\"" << std::endl);
}
 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::~GiNaCViscoElasticConstitutiveLaw(void)
{
	NO_OP;
};

ConstitutiveLaw<doublereal, doublereal>* 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

#if defined(HAVE_SSTREAM)
	std::ostringstream	eps;
	std::ostringstream	epsPrime;
	std::ostringstream	expr;
#else /* HAVE_STRSTREAM_H */
	ostrstream		eps;
	ostrstream		epsPrime;
	ostrstream		expr;
#endif /* HAVE_STRSTREAM_H */

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

	GiNaC::ex f_derEps = gDerEps.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();

	GiNaC::ex f_derEpsPrime = gDerEpsPrime.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
}

void 
GiNaCViscoElasticConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(
	const doublereal& DeltaEps, 
	const doublereal& EpsPrime)
{
	Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon + DeltaEps,
		ElasticConstitutiveLaw<doublereal, doublereal>::EpsilonPrime + EpsPrime);
}

/* GiNaCViscoElasticConstitutiveLaw - end */

#endif /* HAVE_GINAC */

#endif /* GINACCLTP_H */

