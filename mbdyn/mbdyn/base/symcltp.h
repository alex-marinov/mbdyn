/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* Symbolic constitutive law */


#ifndef SYMCLTP_H
#define SYMCLTP_H


#include <constltp.h>
#ifdef HAVE_GINAC
#include <ginac/ginac.h>
#endif /* HAVE_GINAC */

/* SymbolicElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class SymbolicElasticIsotropicConstitutiveLaw 
: public ElasticConstitutiveLaw<T, Tder> {
public:
	SymbolicElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const char *epsilon,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicElasticIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const T& Eps, const T& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const T& DeltaEps, 
			const T& /* EpsPrime */ = 0.);
};

typedef SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal> 
	SymbolicElasticIsotropicConstitutiveLaw1D;
typedef SymbolicElasticIsotropicConstitutiveLaw<Vec3, Mat3x3>
	SymbolicElasticIsotropicConstitutiveLaw3D;
typedef SymbolicElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> 
	SymbolicElasticIsotropicConstitutiveLaw6D;

template <class T, class Tder>
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::SymbolicElasticIsotropicConstitutiveLaw(
		const TplDriveCaller<T>* pDC, const T& PStress, 
		const char *epsilon, const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<T, Tder>(pDC, PStress)
{ 

#ifndef HAVE_GINAC
	std::cerr << "SymbolicElasticIsotropicConstitutiveLaw needs "
		"GiNaC symbolic algebra package" << std::endl;
	THROW(DataManager::ErrGeneric());
#else /* HAVE_GINAC */
	THROW(Err(std::cerr, "symbolic constitutive law "
				"is allowed only for scalar")); 
#endif /* HAVE_GINAC */
}
 
template <class T, class Tder>
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::~SymbolicElasticIsotropicConstitutiveLaw(void)
{
	NO_OP;
};

template <class T, class Tder> ConstitutiveLaw<T, Tder>* 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::pCopy(void) const
{
	return NULL;
}

template <class T, class Tder> std::ostream& 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const
{
  	return out;
}

template <class T, class Tder> void
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::Update(const T& Eps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}

template <class T, class Tder> void 
SymbolicElasticIsotropicConstitutiveLaw<T, Tder>::IncrementalUpdate(const T& DeltaEps, 
		const T& /* EpsPrime */ )
{
	NO_OP;
}


#ifdef HAVE_GINAC

/* specialize for scalar constitutive law */

class SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */
	GiNaC::symbol **gSymList;	/* list of other symbols */
	double *vals;			/* values of other symbols */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDer;			/* derivative */

#if 0
	SymbolicElasticIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const GiNaC::symbol& epsilon,
			const GiNaC::ex& expression,
			const GiNaC::ex& derivative,
			const GiNaC::symbol **symbols = NULL,
			doublereal *v = NULL);
#endif

public:
	SymbolicElasticIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress, const char *epsilon,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicElasticIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.);
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
			const doublereal& /* EpsPrime */ = 0.);
};

SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::SymbolicElasticIsotropicConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC, const doublereal& PStress, 
		const char *epsilon, const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress),
gEps(epsilon), gSymList(0), vals(0)
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;

	GiNaC::lst l; 

	l.append(gEps);

	if (symbols) {
		int i;

		for (i = 0; symbols[i]; i++);

		gSymList = new GiNaC::symbol*[i+1];
		vals = new double[i];

		for (i = 0; symbols[i]; i++) {
			char *s = strdup(symbols[i]), *v = 0;

			if ((v = strchr(s, '='))) {
				*v = '\0';
			}

			gSymList[i] = new GiNaC::symbol(s);
			l.append(*gSymList[i]);

			if (v) {
				vals[i] = atof(&v[1]);
			}

			free(s);
		}

		gSymList[i] = 0;
	}

	try {
		gExpr = GiNaC::ex(expression, l);
	} catch (std::exception e) {
		std::cerr << "expression parsing failed: " << e.what() << std::endl;
		throw e;
	}

	try {
		gDer = gExpr.diff(gEps);
	} catch (std::exception e) {
		std::cerr << "expression differentiation failed: " << e.what() << std::endl;
		throw e;
	}

	silent_cout("Variable:         \"" << gEps << "\"" << std::endl
		<< "Constitutive law: \"" << gExpr << "\"" << std::endl
		<< "Derivative:       \"" << gDer << "\"" << std::endl);
}
 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::~SymbolicElasticIsotropicConstitutiveLaw(void)
{
	if (gSymList) {
		for (int i = 0; gSymList[i]; i++) {
			delete gSymList[i];
		}
		delete[] gSymList;
		delete[] vals;
	}
};

ConstitutiveLaw<doublereal, doublereal>* 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	typedef SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			/* gEps.get_name() */ 0, 0, 0)); /* FIXME! */
      
	return pCL;
}

std::ostream& 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic elastic isotropic, epsilon, \"" << gEps
		<< "\", expression, \"" << gExpr << "\"";

	if (gSymList) {
		int i;

		for (i = 0; gSymList[i]; i++);
		out << "symbols, " << i;

		for (i = 0; gSymList[i]; i++) {
			out << ", \"" << gSymList[i] << "=" << vals[i] << "\"";
		}
	}

  	return Restart_(out);
}

void
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::Update(const doublereal& Eps, 
		const doublereal& /* EpsPrime */ )
{
	GiNaC::lst l;

	ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

	doublereal e = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();

	l.append(gEps == e);

	if (gSymList) {
		for (int i = 0; gSymList[i]; i++) {
			l.append(*gSymList[i] == vals[i]);
		}
	}

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_der = gDer.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_der).to_double();
}

void 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& /* EpsPrime */ )
{
	Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon+DeltaEps);
}

/* LinearElasticIsotropicConstitutiveLaw - end */

#endif /* HAVE_GINAC */

#endif /* SYMCLTP_H */

