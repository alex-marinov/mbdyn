/* $Header$ */
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

/*
 * This family of constitutive laws was sponsored by Hutchinson CdR
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "ac/math.h"
#include "ac/float.h"

#include "dataman.h"
#include "ScalarFunctionsImpl.h"
#include "constltp_impl.h"

template <class T, class Tder>
class NLSFViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
private:
	Tder FDE0;
	Tder FDEPrime0;
	std::vector<const DifferentiableScalarFunction *> FDEsc;
	std::vector<const DifferentiableScalarFunction *> FDEPrimesc;

public:
	NLSFViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
		const T& PStress,
		const Tder& fde0,
		const std::vector<const DifferentiableScalarFunction *>& fdesc,
		const Tder& fdeprime0,
		const std::vector<const DifferentiableScalarFunction *>& fdeprimesc)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
	FDE0(fde0), FDEPrime0(fdeprime0),
	FDEsc(fdesc), FDEPrimesc(fdeprimesc)
	{
		ConstitutiveLaw<T, Tder>::FDE = FDE0;
		ConstitutiveLaw<T, Tder>::FDEPrime = FDEPrime0;
	};

	virtual ~NLSFViscoElasticConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef NLSFViscoElasticConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl,
			cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<T, Tder>::PreStress,
				FDE0, FDEsc, FDEPrime0, FDEPrimesc));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		silent_cerr("NLSFViscoElasticConstitutiveLaw: Restart not implemented"
			<< std::endl);
		throw ErrGeneric();

#if 0 /* ScalarFunction restart not implemented yet */
		out << "nlp viscoelastic, ",
			FDE0.Write(out, ", ");

		for (unsigned i = 0; i < FDEsc.size(); i++) {
			out << ", ";
			if (FDEsc[i]) {
				FDEsc[i]->Restart(out);

			} else {
				out << "null";
			}
		}

		out << ", ", FDEPrime0.Write(out, ", ");

		for (unsigned i = 0; i < FDEPrimesc.size(); i++) {
			out << ", ";
			if (FDEsc[i]) {
				FDEPrimesc[i]->Restart(out);

			} else {
				out << "null";
			}
		}

		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
#endif
	};

	virtual void Update(const T& Eps, const T& EpsPrime = 0.) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		T E = ElasticConstitutiveLaw<T, Tder>::Epsilon
			- ElasticConstitutiveLaw<T, Tder>::Get();

		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
			+ FDE0*E + FDEPrime0*EpsPrime;
		ConstitutiveLaw<T, Tder>::FDE = FDE0;
		ConstitutiveLaw<T, Tder>::FDEPrime = FDEPrime0;

		for (unsigned i = 0; i < FDEsc.size(); i++) {

			/*
			 * each diagonal coefficient is:
			 *
			 * f = f0'*eps + f'(eps) + f0''*epsPrime + f''(epsPrime)
			 *
			 * so the Jacobian matrix contributions are
			 *
			 * f/eps = f0' + f'/eps
			 * f/epsPrime = f0'' + f''/epsPrime
			 */

			doublereal dEps = E(i + 1);
			doublereal dEpsPrime = EpsPrime(i + 1);

			doublereal df1 = 0.;
			doublereal df1DE = 0.;
			doublereal df2 = 0.;
			doublereal df2DEPrime = 0.;

			if (FDEsc[i]) {
				df1 = (*FDEsc[i])(dEps);
				df1DE = FDEsc[i]->ComputeDiff(dEps);
			}

			if (FDEPrimesc[i]) {
				df2 = (*FDEPrimesc[i])(dEpsPrime);
				df2DEPrime = FDEPrimesc[i]->ComputeDiff(dEpsPrime);
			}

			ConstitutiveLaw<T, Tder>::F(i + 1) += df1 + df2;
			ConstitutiveLaw<T, Tder>::FDE(i + 1, i + 1) += df1DE;
			ConstitutiveLaw<T, Tder>::FDEPrime(i + 1, i + 1) += df2DEPrime;
		}
	};

	virtual void IncrementalUpdate(const T& DeltaEps, const T& EpsPrime = 0.) {
		Update(ConstitutiveLaw<T, Tder>::Epsilon + DeltaEps, EpsPrime);
	};
};

template <>
class NLSFViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	doublereal FDE0;
	doublereal FDEPrime0;
	std::vector<const DifferentiableScalarFunction *> FDEsc;
	std::vector<const DifferentiableScalarFunction *> FDEPrimesc;

public:
	NLSFViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
		const doublereal& PStress,
		const doublereal& fde0,
		const std::vector<const DifferentiableScalarFunction *>& fdesc,
		const doublereal& fdeprime0,
		const std::vector<const DifferentiableScalarFunction *>& fdeprimesc)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress),
	FDE0(fde0), FDEPrime0(fdeprime0),
	FDEsc(fdesc), FDEPrimesc(fdeprimesc)
	{
		ConstitutiveLaw<doublereal, doublereal>::FDE = FDE0;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = FDEPrime0;
	};

	virtual ~NLSFViscoElasticConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef NLSFViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl,
			cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
				FDE0, FDEsc, FDEPrime0, FDEPrimesc));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		silent_cerr("NLSFViscoElasticConstitutiveLaw: Restart not implemented"
			<< std::endl);
		throw ErrGeneric();

#if 0 /* ScalarFunction restart not implemented yet */
		out << "nlp viscoelastic, ",
			FDE0.Write(out, ", ");

		for (unsigned i = 0; i < FDEsc.size(); i++) {
			out << ", ";
			if (FDEsc[i]) {
				FDEsc[i]->Restart(out);

			} else {
				out << "null";
			}
		}

		out << ", ", FDEPrime0.Write(out, ", ");

		for (unsigned i = 0; i < FDEPrimesc.size(); i++) {
			out << ", ";
			if (FDEsc[i]) {
				FDEPrimesc[i]->Restart(out);

			} else {
				out << "null";
			}
		}

		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
#endif
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal dEps = ElasticConstitutiveLaw<doublereal, doublereal>::Epsilon
			- ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		doublereal dEpsPrime = EpsPrime;

		ConstitutiveLaw<doublereal, doublereal>::F = ElasticConstitutiveLaw<doublereal, doublereal>::PreStress
			+ FDE0*dEps + FDEPrime0*dEpsPrime;
		ConstitutiveLaw<doublereal, doublereal>::FDE = FDE0;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = FDEPrime0;

		/*
		 * each diagonal coefficient is:
		 *
		 * f = (f0' + f'(eps))*eps + (f0'' + f''(eps))*epsPrime
		 *
		 * so the Jacobian matrix contributions are
		 *
		 * f/eps = f0' + f' + f'/eps * eps + f''/eps * epsPrime
		 * f/epsPrime = f''
		 */

		doublereal df1 = 0.;
		doublereal df1DE = 0.;
		doublereal df2 = 0.;
		doublereal df2DEPrime = 0.;

		if (FDEsc[0]) {
			df1 = (*FDEsc[0])(dEps);
			df1DE = FDEsc[0]->ComputeDiff(dEps);
		}

		if (FDEPrimesc[0]) {
			df2 = (*FDEPrimesc[0])(dEpsPrime);
			df2DEPrime = FDEPrimesc[0]->ComputeDiff(dEpsPrime);
		}

		ConstitutiveLaw<doublereal, doublereal>::F += df1 + df2;
		ConstitutiveLaw<doublereal, doublereal>::FDE += df1DE;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime += df2DEPrime;
	};

	virtual void IncrementalUpdate(const doublereal& DeltaEps, const doublereal& EpsPrime = 0.) {
		Update(ConstitutiveLaw<doublereal, doublereal>::Epsilon + DeltaEps, EpsPrime);
	};
};

typedef NLSFViscoElasticConstitutiveLaw<doublereal, doublereal> NLSFViscoElasticConstitutiveLaw1D;
typedef NLSFViscoElasticConstitutiveLaw<Vec3, Mat3x3> NLSFViscoElasticConstitutiveLaw3D;
typedef NLSFViscoElasticConstitutiveLaw<Vec6, Mat6x6> NLSFViscoElasticConstitutiveLaw6D;

/* specific functional object(s) */
template <class T, class Tder>
struct NLSFViscoElasticCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType)
	{
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		unsigned dim;
		if (typeid(T) == typeid(doublereal)) {
			dim = 1;

		} else if (typeid(T) == typeid(Vec3)) {
			dim = 3;

		} else if (typeid(T) == typeid(Vec6)) {
			dim = 6;

		} else {
			silent_cerr("Invalid dimensionality "
				"for NLSF viscoelastic constitutive law "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric();
		}

		/* stiffness */
		Tder FDE0(0.);
		FDE0 = HP.Get(FDE0);

		std::vector<const DifferentiableScalarFunction *> FDEsc(dim);
		for (unsigned i = 0; i < dim; i++) {
			if (HP.IsKeyWord("null")) {
				FDEsc[i] = 0;

			} else {
				const BasicScalarFunction *const sc = ParseScalarFunction(HP, (DataManager *const)pDM);
				FDEsc[i] = dynamic_cast<const DifferentiableScalarFunction *>(sc);
				if (FDEsc[i] == 0) {
					silent_cerr("NLSFViscoElasticCLR: "
						"stiffness scalar function #" << i << " "
						"at line " << HP.GetLineData() << " "
						"must be differentiable" << std::endl);
					throw ErrGeneric();
				}
			}
		}

		/* damping */
		Tder FDEPrime0(0.);
		FDEPrime0 = HP.Get(FDEPrime0);

		std::vector<const BasicScalarFunction *> FDEPrimesc(dim);
		for (unsigned i = 0; i < dim; i++) {
			if (HP.IsKeyWord("null")) {
				FDEPrimesc[i] = 0;

			} else {
				const BasicScalarFunction *const sc = ParseScalarFunction(HP, (DataManager *const)pDM);
				FDEPrimesc[i] = dynamic_cast<const DifferentiableScalarFunction *>(sc);
				if (FDEPrimesc[i] == 0) {
					silent_cerr("NLSFViscoElasticCLR: "
						"damping scalar function #" << i << " "
						"at line " << HP.GetLineData() << " "
						"must be differentiable" << std::endl);
					throw ErrGeneric();
				}
			}
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef NLSFViscoElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L,
			L(pTplDC, PreStress,
				FDE0, (const std::vector<const DifferentiableScalarFunction *>&)FDEsc,
				FDEPrime0, (const std::vector<const DifferentiableScalarFunction *>&)FDEPrimesc));

		return pCL;
	};
};

int
NLSF_init(void)
{
	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new NLSFViscoElasticCLR<doublereal, doublereal>;
	if (!SetCL1D("nlsf" "viscoelastic", rf1D)) {
		delete rf1D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw1D: "
			"init failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new NLSFViscoElasticCLR<Vec3, Mat3x3>;
	if (!SetCL3D("nlsf" "viscoelastic", rf3D)) {
		delete rf3D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw3D: "
			"init failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec6, Mat6x6> *rf6D = new NLSFViscoElasticCLR<Vec6, Mat6x6>;
	if (!SetCL6D("nlsf" "viscoelastic", rf6D)) {
		delete rf6D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw6D: "
			"init failed" << std::endl);

		return -1;
	}

	return 0;
}

