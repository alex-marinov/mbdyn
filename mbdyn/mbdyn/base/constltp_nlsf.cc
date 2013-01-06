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

/*
 * This family of constitutive laws was sponsored by Hutchinson CdR
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "ScalarFunctionsImpl.h"
#include "constltp_impl.h"

template <class T, class Tder>
class NLSFViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
private:
	ConstLawType::Type CLType;
	Tder FDE0;
	Tder FDEPrime0;
	std::vector<const DifferentiableScalarFunction *> FDEsc;
	std::vector<const DifferentiableScalarFunction *> FDEPrimesc;

public:
	NLSFViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
		const T& PStress,
		const ConstLawType::Type &cltype,
		const Tder& fde0,
		const std::vector<const DifferentiableScalarFunction *>& fdesc,
		const Tder& fdeprime0,
		const std::vector<const DifferentiableScalarFunction *>& fdeprimesc)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
	CLType(cltype),
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
		return CLType;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef NLSFViscoElasticConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl,
			cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<T, Tder>::PreStress,
				CLType, FDE0, FDEsc, FDEPrime0, FDEPrimesc));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		silent_cerr("NLSFViscoElasticConstitutiveLaw: Restart not implemented"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

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

		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress;
		if (CLType & ConstLawType::ELASTIC) {
			ConstitutiveLaw<T, Tder>::F += FDE0*E;
			ConstitutiveLaw<T, Tder>::FDE = FDE0;
		}
		if (CLType & ConstLawType::VISCOUS) {
			ConstitutiveLaw<T, Tder>::F += FDEPrime0*EpsPrime;
			ConstitutiveLaw<T, Tder>::FDEPrime = FDEPrime0;
		}

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

			doublereal f = 0.;

			if ((CLType & ConstLawType::ELASTIC) && FDEsc[i]) {
				doublereal dEps = E(i + 1);
				f += (*FDEsc[i])(dEps);
				doublereal df1DE = FDEsc[i]->ComputeDiff(dEps);
				ConstitutiveLaw<T, Tder>::FDE(i + 1, i + 1) += df1DE;
			}

			if ((CLType & ConstLawType::ELASTIC) && FDEPrimesc[i]) {
				doublereal dEpsPrime = EpsPrime(i + 1);
				f += (*FDEPrimesc[i])(dEpsPrime);
				doublereal df2DEPrime = FDEPrimesc[i]->ComputeDiff(dEpsPrime);
				ConstitutiveLaw<T, Tder>::FDEPrime(i + 1, i + 1) += df2DEPrime;
			}

			ConstitutiveLaw<T, Tder>::F(i + 1) += f;
		}
	};
};

template <>
class NLSFViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	ConstLawType::Type CLType;
	doublereal FDE0;
	doublereal FDEPrime0;
	std::vector<const DifferentiableScalarFunction *> FDEsc;
	std::vector<const DifferentiableScalarFunction *> FDEPrimesc;

public:
	NLSFViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
		const doublereal& PStress,
		const ConstLawType::Type& cltype,
		const doublereal& fde0,
		const std::vector<const DifferentiableScalarFunction *>& fdesc,
		const doublereal& fdeprime0,
		const std::vector<const DifferentiableScalarFunction *>& fdeprimesc)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pDC, PStress),
	CLType(cltype),
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
		return CLType;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef NLSFViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl,
			cl(ElasticConstitutiveLaw<doublereal, doublereal>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<doublereal, doublereal>::PreStress,
				CLType, FDE0, FDEsc, FDEPrime0, FDEPrimesc));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		silent_cerr("NLSFViscoElasticConstitutiveLaw: Restart not implemented"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

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

		ConstitutiveLaw<doublereal, doublereal>::F = ElasticConstitutiveLaw<doublereal, doublereal>::PreStress;
		if (CLType & ConstLawType::ELASTIC) {
			ConstitutiveLaw<doublereal, doublereal>::F += FDE0*dEps;
			ConstitutiveLaw<doublereal, doublereal>::FDE = FDE0;
		}
		if (CLType & ConstLawType::VISCOUS) {
			ConstitutiveLaw<doublereal, doublereal>::F += FDEPrime0*dEpsPrime;
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = FDEPrime0;
		}

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

		doublereal f = 0.;

		if ((CLType & ConstLawType::ELASTIC) && FDEsc[0]) {
			f += (*FDEsc[0])(dEps);
			doublereal df1DE = FDEsc[0]->ComputeDiff(dEps);
			ConstitutiveLaw<doublereal, doublereal>::FDE += df1DE;
		}

		if ((CLType & ConstLawType::VISCOUS) && FDEPrimesc[0]) {
			f += (*FDEPrimesc[0])(dEpsPrime);
			doublereal df2DEPrime = FDEPrimesc[0]->ComputeDiff(dEpsPrime);
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime += df2DEPrime;
		}

		ConstitutiveLaw<doublereal, doublereal>::F += f;
	};
};

typedef NLSFViscoElasticConstitutiveLaw<doublereal, doublereal> NLSFViscoElasticConstitutiveLaw1D;
typedef NLSFViscoElasticConstitutiveLaw<Vec3, Mat3x3> NLSFViscoElasticConstitutiveLaw3D;
typedef NLSFViscoElasticConstitutiveLaw<Vec6, Mat6x6> NLSFViscoElasticConstitutiveLaw6D;

/* specific functional object(s) */
template <class T, class Tder, ConstLawType::Type Typ>
struct NLSFViscoElasticCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType)
	{
		ConstitutiveLaw<T, Tder>* pCL = 0;

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
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* stiffness */
		Tder FDE0(mb_zero<Tder>());
		bool bElastic(false);
		std::vector<const DifferentiableScalarFunction *> FDEsc(dim);
		for (unsigned i = 0; i < dim; i++) {
			FDEsc[i] = 0;
		}

		if (Typ & ConstLawType::ELASTIC) {
			FDE0 = HP.Get(FDE0);

			bElastic = !IsNull<Tder>(FDE0);
			for (unsigned i = 0; i < dim; i++) {
				if (!HP.IsKeyWord("null")) {
					const BasicScalarFunction *const sc
						= ParseScalarFunction(HP, (DataManager *const)pDM);
					FDEsc[i] = dynamic_cast<const DifferentiableScalarFunction *>(sc);
					if (FDEsc[i] == 0) {
						silent_cerr("NLSFViscoElasticCLR: "
							"stiffness scalar function #" << i << " "
							"at line " << HP.GetLineData() << " "
							"must be differentiable" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					bElastic = true;
				}
			}
		}

		/* damping */
		Tder FDEPrime0(mb_zero<Tder>());
		bool bViscous(false);
		std::vector<const DifferentiableScalarFunction *> FDEPrimesc(dim);
		for (unsigned i = 0; i < dim; i++) {
			FDEPrimesc[i] = 0;
		}

		if (Typ & ConstLawType::VISCOUS) {
			if ((Typ & ConstLawType::ELASTIC) && HP.IsKeyWord("proportional")) {
				FDEPrime0 = FDE0*HP.GetReal();
			} else {
				FDEPrime0 = HP.Get(FDEPrime0);
			}

			bViscous = !IsNull<Tder>(FDEPrime0);
			for (unsigned i = 0; i < dim; i++) {
				if (!HP.IsKeyWord("null")) {
					const BasicScalarFunction *const sc
						= ParseScalarFunction(HP, (DataManager *const)pDM);
					FDEPrimesc[i] = dynamic_cast<const DifferentiableScalarFunction *>(sc);
					if (FDEPrimesc[i] == 0) {
						silent_cerr("NLSFViscoElasticCLR: "
							"damping scalar function #" << i << " "
							"at line " << HP.GetLineData() << " "
							"must be differentiable" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					bViscous = true;
				}
			}
		}

		/* Prestress and prestrain */
		T PreStress(mb_zero<T>());
		GetPreStress(HP, PreStress);
		TplDriveCaller<T>* pTplDC = GetPreStrain<T>(pDM, HP);

		if (bElastic && bViscous) {
			CLType = ConstLawType::VISCOELASTIC;
		} else if (bElastic) {
			CLType = ConstLawType::ELASTIC;
		} else if (bViscous) {
			CLType = ConstLawType::VISCOUS;
		} else {
			/* needs to be at least elastic... */
			CLType = ConstLawType::ELASTIC;
		}

		typedef NLSFViscoElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L,
			L(pTplDC, PreStress,
				CLType,
				FDE0, FDEsc,
				FDEPrime0, FDEPrimesc));

		return pCL;
	};
};

int
NLSF_init(void)
{
	// 1D viscoelastic
	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new NLSFViscoElasticCLR<doublereal, doublereal, ConstLawType::VISCOELASTIC>;
	if (!SetCL1D("nlsf" "viscoelastic", rf1D)) {
		delete rf1D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw1D: "
			"init failed" << std::endl);

		return -1;
	}

	// 1D elastic
	rf1D = new NLSFViscoElasticCLR<doublereal, doublereal, ConstLawType::ELASTIC>;
	if (!SetCL1D("nlsf" "elastic", rf1D)) {
		delete rf1D;

		silent_cerr("NLSFElasticConstitutiveLaw1D: "
			"init failed" << std::endl);

		return -1;
	}

	// 1D viscous
	rf1D = new NLSFViscoElasticCLR<doublereal, doublereal, ConstLawType::VISCOUS>;
	if (!SetCL1D("nlsf" "viscous", rf1D)) {
		delete rf1D;

		silent_cerr("NLSFViscousConstitutiveLaw1D: "
			"init failed" << std::endl);

		return -1;
	}

	// 3D viscoelastic
	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D
		= new NLSFViscoElasticCLR<Vec3, Mat3x3, ConstLawType::VISCOELASTIC>;
	if (!SetCL3D("nlsf" "viscoelastic", rf3D)) {
		delete rf3D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw3D: "
			"init failed" << std::endl);

		return -1;
	}

	// 3D elastic
	rf3D = new NLSFViscoElasticCLR<Vec3, Mat3x3, ConstLawType::ELASTIC>;
	if (!SetCL3D("nlsf" "elastic", rf3D)) {
		delete rf3D;

		silent_cerr("NLSFElasticConstitutiveLaw3D: "
			"init failed" << std::endl);

		return -1;
	}

	// 3D viscous
	rf3D = new NLSFViscoElasticCLR<Vec3, Mat3x3, ConstLawType::VISCOUS>;
	if (!SetCL3D("nlsf" "viscous", rf3D)) {
		delete rf3D;

		silent_cerr("NLSFViscousConstitutiveLaw3D: "
			"init failed" << std::endl);

		return -1;
	}

	// 6D viscoelastic
	ConstitutiveLawRead<Vec6, Mat6x6> *rf6D
		= new NLSFViscoElasticCLR<Vec6, Mat6x6, ConstLawType::VISCOELASTIC>;
	if (!SetCL6D("nlsf" "viscoelastic", rf6D)) {
		delete rf6D;

		silent_cerr("NLSFViscoElasticConstitutiveLaw6D: "
			"init failed" << std::endl);

		return -1;
	}

	// 6D elastic
	rf6D = new NLSFViscoElasticCLR<Vec6, Mat6x6, ConstLawType::ELASTIC>;
	if (!SetCL6D("nlsf" "elastic", rf6D)) {
		delete rf6D;

		silent_cerr("NLSFElasticConstitutiveLaw6D: "
			"init failed" << std::endl);

		return -1;
	}

	// 6D viscous
	rf6D = new NLSFViscoElasticCLR<Vec6, Mat6x6, ConstLawType::VISCOUS>;
	if (!SetCL6D("nlsf" "viscous", rf6D)) {
		delete rf6D;

		silent_cerr("NLSFViscousConstitutiveLaw6D: "
			"init failed" << std::endl);

		return -1;
	}

	return 0;
}

