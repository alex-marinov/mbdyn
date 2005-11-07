/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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


#ifndef GINACCLTP_H
#define GINACCLTP_H

#ifdef HAVE_GINAC

#include <symcltp.h>
#include <ginac/ginac.h>

/* specialize for scalar constitutive law */

/* SymbolicElasticIsotropicConstitutiveLaw - begin */

template<>
class SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */
	GiNaC::symbol **gSymList;	/* list of other symbols */
	doublereal *vals;		/* values of other symbols */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEps;		/* derivative */

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
		unsigned int i;

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
				char	*next = NULL;
				vals[i] = strtod(&v[1], &next);
				if (next == NULL || next[0] != '\0') {
					silent_cerr("unable to parse value "
						<< &v[1] << std::endl);
					throw ErrGeneric();
				}
			}

			free(s);
		}

		gSymList[i] = 0;
	}

	try {
		gExpr = GiNaC::ex(expression, l);
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

	silent_cout("\tSymbolicElasticIsotropicConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gDerEps << "\"" << std::endl);
}
 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::~SymbolicElasticIsotropicConstitutiveLaw(void)
{
	if (gSymList) {
		for (unsigned int i = 0; gSymList[i]; i++) {
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
		unsigned int i;

		for (i = 0; gSymList[i]; i++);
		out << "symbols, " << i;

		for (i = 0; gSymList[i]; i++) {
			out << ", \"" << gSymList[i] << "=" << vals[i] << "\"";
		}
	}

  	return Restart_int(out);
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
		for (unsigned int i = 0; gSymList[i]; i++) {
			l.append(*gSymList[i] == vals[i]);
		}
	}

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_derEps = gDerEps.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDE = GiNaC::ex_to<GiNaC::numeric>(f_derEps).to_double();
}

void 
SymbolicElasticIsotropicConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& /* EpsPrime */ )
{
	Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon + DeltaEps);
}

/* SymbolicElasticIsotropicConstitutiveLaw - end */

/* SymbolicViscousIsotropicConstitutiveLaw - begin */

template<>
class SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */
	GiNaC::symbol **gSymList;	/* list of other symbols */
	doublereal *vals;		/* values of other symbols */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEpsPrime;		/* derivative */

#if 0
	SymbolicViscousIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const GiNaC::symbol& epsilonPrime,
			const GiNaC::ex& expression,
			const GiNaC::ex& derivative,
			const GiNaC::symbol **symbols = NULL,
			doublereal *v = NULL);
#endif

public:
	SymbolicViscousIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const char *epsilonPrime,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicViscousIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.);
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
			const doublereal& EpsPrime = 0.);
};

SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::SymbolicViscousIsotropicConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC, const doublereal& PStress, 
		const char *epsilonPrime,
		const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress),
gEpsPrime(epsilonPrime), gSymList(0), vals(0)
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

	GiNaC::lst l; 

	l.append(gEpsPrime);

	if (symbols) {
		unsigned int i;

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
				char	*next = NULL;
				vals[i] = strtod(&v[1], &next);
				if (next == NULL || next[0] != '\0') {
					silent_cerr("unable to parse value "
						<< &v[1] << std::endl);
					throw ErrGeneric();
				}
			}

			free(s);
		}

		gSymList[i] = 0;
	}

	try {
		gExpr = GiNaC::ex(expression, l);
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

	silent_cout("\tSymbolicViscousIsotropicConstitutiveLaw:" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gDerEpsPrime << "\"" << std::endl);
}
 
SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::~SymbolicViscousIsotropicConstitutiveLaw(void)
{
	if (gSymList) {
		for (unsigned int i = 0; gSymList[i]; i++) {
			delete gSymList[i];
		}
		delete[] gSymList;
		delete[] vals;
	}
};

ConstitutiveLaw<doublereal, doublereal>* 
SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	typedef SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			/* gEpsPrime.get_name() */ 0, 0, 0)); /* FIXME! */
      
	return pCL;
}

std::ostream& 
SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic viscous isotropic, "
		<< "epsilon prime, \"" << gEpsPrime
		<< "\", expression, \"" << gExpr << "\"";

	if (gSymList) {
		unsigned int i;

		for (i = 0; gSymList[i]; i++);
		out << "symbols, " << i;

		for (i = 0; gSymList[i]; i++) {
			out << ", \"" << gSymList[i] << "=" << vals[i] << "\"";
		}
	}

  	return Restart_int(out);
}

void
SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::Update(const doublereal& /* Eps */ , 
		const doublereal& EpsPrime)
{
	GiNaC::lst l;

	ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

	l.append(gEpsPrime == EpsPrime);

	if (gSymList) {
		for (unsigned int i = 0; gSymList[i]; i++) {
			l.append(*gSymList[i] == vals[i]);
		}
	}

	GiNaC::ex f_expr = gExpr.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::F = 
		ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
		+ GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();

	GiNaC::ex f_derEpsPrime = gDerEpsPrime.subs(l);
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = GiNaC::ex_to<GiNaC::numeric>(f_derEpsPrime).to_double();
}

void 
SymbolicViscousIsotropicConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(const doublereal& /* DeltaEps */ ,
		const doublereal& EpsPrime)
{
	Update(0, ElasticConstitutiveLaw<doublereal, doublereal>::EpsilonPrime + EpsPrime);
}

/* SymbolicViscousIsotropicConstitutiveLaw - end */

/* SymbolicViscoElasticIsotropicConstitutiveLaw - begin */

template<>
class SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	GiNaC::symbol gEps;		/* parameter symbol */
	GiNaC::symbol gEpsPrime;	/* parameter derivative symbol */
	GiNaC::symbol **gSymList;	/* list of other symbols */
	doublereal *vals;		/* values of other symbols */

	GiNaC::ex gExpr;		/* expression */
	GiNaC::ex gDerEps;		/* derivative */
	GiNaC::ex gDerEpsPrime;		/* derivative */

#if 0
	SymbolicViscoElasticIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const GiNaC::symbol& epsilon,
			const GiNaC::symbol& epsilonPrime,
			const GiNaC::ex& expression,
			const GiNaC::ex& derivative,
			const GiNaC::symbol **symbols = NULL,
			doublereal *v = NULL);
#endif

public:
	SymbolicViscoElasticIsotropicConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const char *epsilon, const char *epsilonPrime,
			const char *expression, const char *const *symbols);
     	virtual ~SymbolicViscoElasticIsotropicConstitutiveLaw(void);
     	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.);
	virtual void IncrementalUpdate(const doublereal& DeltaEps, 
			const doublereal& EpsPrime = 0.);
};

SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::SymbolicViscoElasticIsotropicConstitutiveLaw(
		const TplDriveCaller<doublereal>* pDC, const doublereal& PStress, 
		const char *epsilon, const char *epsilonPrime,
		const char *expression, 
		const char *const *symbols = NULL)
: ElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress),
gEps(epsilon), gEpsPrime(epsilonPrime), gSymList(0), vals(0)
{ 
	ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;
	ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

	GiNaC::lst l; 

	l.append(gEps);
	l.append(gEpsPrime);

	if (symbols) {
		unsigned int i;

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
				char	*next = NULL;
				vals[i] = strtod(&v[1], &next);
				if (next == NULL || next[0] != '\0') {
					silent_cerr("unable to parse value "
						<< &v[1] << std::endl);
					throw ErrGeneric();
				}
			}

			free(s);
		}

		gSymList[i] = 0;
	}

	try {
		gExpr = GiNaC::ex(expression, l);
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

	silent_cout("\tSymbolicViscoElasticIsotropicConstitutiveLaw:" << std::endl
		<< "\t\tEps:              \"" << gEps << "\"" << std::endl
		<< "\t\tEpsPrime:         \"" << gEpsPrime << "\"" << std::endl
		<< "\t\tConstitutive law: \"" << gExpr << "\"" << std::endl
		<< "\t\tDer/Eps:          \"" << gDerEps << "\"" << std::endl
		<< "\t\tDer/EpsPrime:     \"" << gDerEpsPrime << "\"" << std::endl);
}
 
SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::~SymbolicViscoElasticIsotropicConstitutiveLaw(void)
{
	if (gSymList) {
		for (unsigned int i = 0; gSymList[i]; i++) {
			delete gSymList[i];
		}
		delete[] gSymList;
		delete[] vals;
	}
};

ConstitutiveLaw<doublereal, doublereal>* 
SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::pCopy(void) const
{
	ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

	typedef SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal> cl;
	SAFENEWWITHCONSTRUCTOR(pCL, 
		cl, 
		cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(), 
			ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
			/* gEps.get_name() */ 0, 0, 0, 0)); /* FIXME! */
      
	return pCL;
}

std::ostream& 
SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::Restart(std::ostream& out) const
{
	out << "symbolic elastic isotropic, "
		 << "epsilon, \"" << gEps
		<< "epsilon prime, \"" << gEpsPrime
		<< "\", expression, \"" << gExpr << "\"";

	if (gSymList) {
		unsigned int i;

		for (i = 0; gSymList[i]; i++);
		out << "symbols, " << i;

		for (i = 0; gSymList[i]; i++) {
			out << ", \"" << gSymList[i] << "=" << vals[i] << "\"";
		}
	}

  	return Restart_int(out);
}

void
SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::Update(const doublereal& Eps, 
		const doublereal& EpsPrime)
{
	GiNaC::lst l;

	ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

	doublereal e = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();

	l.append(gEps == e);

	ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

	l.append(gEpsPrime == EpsPrime);

	if (gSymList) {
		for (unsigned int i = 0; gSymList[i]; i++) {
			l.append(*gSymList[i] == vals[i]);
		}
	}

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
SymbolicViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal>::IncrementalUpdate(const doublereal& DeltaEps, 
		const doublereal& EpsPrime)
{
	Update(ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon + DeltaEps,
			ElasticConstitutiveLaw<doublereal, doublereal>::EpsilonPrime + EpsPrime);
}

/* SymbolicViscoElasticIsotropicConstitutiveLaw - end */

#endif /* HAVE_GINAC */

#endif /* GINACCLTP_H */

