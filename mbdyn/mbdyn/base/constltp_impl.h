/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Legami costitutivi */


#ifndef CONSTLTP_IMPL_H
#define CONSTLTP_IMPL_H

#include <limits>
#include <cfloat>
#include <cmath>

#include "tpldrive_impl.h"
#include "constltp.h"
#include "hint_impl.h"
#include "elem.h"

/* ConstitutiveLawArray - begin */

template <class T, class Tder>
class ConstitutiveLawArray
: public ConstitutiveLaw<T, Tder> {
protected:
	ConstLawType::Type m_type;
	std::vector<ConstitutiveLaw<T, Tder> *> m_clv;

public:
	ConstitutiveLawArray(const std::vector<ConstitutiveLaw<T, Tder> *>& clv)
	: m_clv(clv)
	{
		unsigned type = 0;
		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::const_iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			type |= (*i)->GetConstLawType();
		}
		m_type = ConstLawType::Type(type);
	};

	virtual ~ConstitutiveLawArray(void)
	{
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return m_type;
	};

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			(*i)->SetValue(pDM, X, XP, ph);
		}
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		std::vector<ConstitutiveLaw<T, Tder> *> clv(m_clv.size());
		for (unsigned i = 0; i < m_clv.size(); i++) {
			clv[i] = m_clv[i]->pCopy();
		}

		typedef ConstitutiveLawArray<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(clv));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "array, " << m_clv.size();
		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::const_iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			out << ", ", (*i)->Restart(out);
		}
		return out;
	};

	virtual void Update(const T& Eps, const T& EpsPrime = ::mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
		ConstitutiveLaw<T, Tder>::F = ::mb_zero<T>();
		if (m_type & ConstLawType::ELASTIC) {
			ConstitutiveLaw<T, Tder>::FDE = ::mb_zero<Tder>();
		}
		if (m_type & ConstLawType::VISCOUS) {
			ConstitutiveLaw<T, Tder>::FDEPrime = ::mb_zero<Tder>();
		}

		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			(*i)->Update(Eps, EpsPrime);
			ConstitutiveLaw<T, Tder>::F += (*i)->GetF();
			if (m_type & ConstLawType::ELASTIC) {
				ConstitutiveLaw<T, Tder>::FDE += (*i)->GetFDE();
			}

			if (m_type & ConstLawType::VISCOUS) {
				ConstitutiveLaw<T, Tder>::FDEPrime += (*i)->GetFDEPrime();
			}
		}
	};

	virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			(*i)->AfterConvergence(Eps, EpsPrime);
		}
	};

	virtual std::ostream& OutputAppend(std::ostream& out) const {
		for (typename std::vector<ConstitutiveLaw<T, Tder> *>::const_iterator i = m_clv.begin(); i != m_clv.end(); i++) {
			(*i)->OutputAppend(out);
		}
		return out;
	};
};

/* ConstitutiveLawArray - end */


/* ElasticConstitutiveLaw - begin */

template <class T, class Tder>
class ElasticConstitutiveLaw
: public ConstitutiveLaw<T, Tder>, public TplDriveOwner<T> {
protected:
	T PreStress;

	virtual std::ostream&
	Restart_int(std::ostream& out) const
	{
		out << ", prestress, ",
			Write(out, PreStress /* + GetF() */ , ", ");
		if (TplDriveOwner<T>::pGetDriveCaller()) {
			out << ", prestrain, single, ",
				Write(out, -ConstitutiveLaw<T, Tder>::Epsilon, ", ") << ", one /* ",
				TplDriveOwner<T>::pGetDriveCaller()->Restart(out) << " */ ";
		}

		return out;
	};

public:
	ElasticConstitutiveLaw(const TplDriveCaller<T>* pDC, const T& PStress)
	: TplDriveOwner<T>(pDC), PreStress(PStress)
	{
		NO_OP;
	};

	virtual ~ElasticConstitutiveLaw(void)
	{
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0)
	{
		if (ph) {
			for (unsigned i = 0; i < ph->size(); i++) {
				TplVecHint<T> *pStressH = dynamic_cast<TplVecHint<T> *>((*ph)[i]);

				if (pStressH) {
					PreStress = pStressH->pCreateVec(pDM);
					continue;
				}

				TplDriveHint<T> *pStrainH = dynamic_cast<TplDriveHint<T> *>((*ph)[i]);

				if (pStrainH) {
					TplDriveCaller<T> *pDC = pStrainH->pCreateDrive(pDM);
					if (pDC == 0) {
						silent_cerr("ElasticConstitutiveLaw: "
							"unable to create prestrain drive after hint "
							"#" << i << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					TplDriveOwner<T>::Set(pDC);

					continue;
				}
			}
		}
	};

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const
	{
		if (strncasecmp(s, "prestress{" /*}*/ , STRLENOF("prestress{" /*}*/ )) == 0) {
			s += STRLENOF("prestress{" /*}*/ );
		
			size_t	len = strlen(s);

			if (s[len - 1] != '}') {
				return 0;
			}

			char *sStr = new char[len + 1];
			memcpy(sStr, s, len + 1);
			sStr[len - 1] = ';';

			return new TplVecHint<T>(sStr);

		} else if (strncasecmp(s, "prestrain{" /*}*/ , STRLENOF("prestrain{" /*}*/ )) == 0) {
			s += STRLENOF("prestrain{" /*}*/ );
		
			size_t	len = strlen(s);

			if (s[len - 1] != '}') {
				return 0;
			}

			char *sStr = new char[len + 1];
			memcpy(sStr, s, len + 1);
			sStr[len - 1] = ';';

			return new TplDriveHint<T>(sStr);
		}

		return 0;
	};
};

typedef ElasticConstitutiveLaw<doublereal, doublereal> ElasticConstitutiveLaw1D;
typedef ElasticConstitutiveLaw<Vec3, Mat3x3> ElasticConstitutiveLaw3D;
typedef ElasticConstitutiveLaw<Vec6, Mat6x6> ElasticConstitutiveLaw6D;

/* ElasticConstitutiveLaw - end */


/* LinearElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticIsotropicConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
private:
	doublereal dStiffness;

public:
	LinearElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, doublereal dStiff)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress), dStiffness(dStiff) {
		mb_deye<Tder>(ConstitutiveLaw<T, Tder>::FDE, dStiffness);
	};

	virtual ~LinearElasticIsotropicConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef LinearElasticIsotropicConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw<T, Tder>::PreStress,
					dStiffness));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "linear elastic isotropic, " << dStiffness;
		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */  = mb_zero<T>()) {
		ElasticConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
			+ (ElasticConstitutiveLaw<T, Tder>::Epsilon - ElasticConstitutiveLaw<T, Tder>::Get())*dStiffness;
	};
};

typedef LinearElasticIsotropicConstitutiveLaw<doublereal, doublereal> LinearElasticIsotropicConstitutiveLaw1D;
typedef LinearElasticIsotropicConstitutiveLaw<Vec3, Mat3x3> LinearElasticIsotropicConstitutiveLaw3D;
typedef LinearElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> LinearElasticIsotropicConstitutiveLaw6D;

/* LinearElasticIsotropicConstitutiveLaw - end */


/* LinearElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	LinearElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const Tder& Stiff)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		ConstitutiveLaw<T, Tder>::FDE = Stiff;
	};

	virtual ~LinearElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef LinearElasticGenericConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw<T, Tder>::PreStress,
					ConstitutiveLaw<T, Tder>::FDE));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "linear elastic generic, ",
			Write(out, ConstitutiveLaw<T, Tder>::FDE, ", ");
		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */ = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
			+ ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon - ElasticConstitutiveLaw<T, Tder>::Get());
	};
};

typedef LinearElasticGenericConstitutiveLaw<doublereal, doublereal> LinearElasticGenericConstitutiveLaw1D;
typedef LinearElasticGenericConstitutiveLaw<Vec3, Mat3x3> LinearElasticGenericConstitutiveLaw3D;
typedef LinearElasticGenericConstitutiveLaw<Vec6, Mat6x6> LinearElasticGenericConstitutiveLaw6D;

/* LinearElasticGenericConstitutiveLaw - end */


/* LinearElasticGenericAxialTorsionCouplingConstitutiveLaw - begin */

template <class T, class Tder>
class LinearElasticGenericAxialTorsionCouplingConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			const Tder& = mb_zero<Tder>(),
			doublereal = 0.)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ConstitutiveLaw<T, Tder>::Err(std::cerr, "axial-torsion coupling constitutive law "
						"is allowed only for beams (6x6)"));
	};

	virtual ~LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

template<>
class LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6>
: public ElasticConstitutiveLaw6D {
protected:
	Mat6x6 Stiffness;
	doublereal dRefTorsion;
	doublereal dAxialTorsionCoupling;

public:
	LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<Vec6>* pDC,
			const Vec6& PStress, const Mat6x6& Stiff,
			doublereal dAxTors)
	: ElasticConstitutiveLaw6D(pDC, PStress),
	Stiffness(Stiff),
	dRefTorsion(Stiffness.dGet(4, 4)),
	dAxialTorsionCoupling(dAxTors) {
		FDE = Stiffness;
	};

	virtual ~LinearElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy(void) const {
		ConstitutiveLaw<Vec6, Mat6x6>* pCL = 0;

		typedef LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					Stiffness,
					dAxialTorsionCoupling));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "linear elastic generic axial torsion coupling, ",
			Write(out, FDE, ", ") << ", " << dAxialTorsionCoupling;
		return Restart_int(out);
	};

	virtual void Update(const Vec6& Eps, const Vec6& /* EpsPrime */ = mb_zero<Vec6>()) {
		Epsilon = Eps;
		doublereal d = Epsilon.dGet(1);
		FDE.Put(4, 4, dRefTorsion + d*dAxialTorsionCoupling);
		F = PreStress + FDE*(Epsilon-Get());
	};
};

/* LinearElasticGenericAxialTorsionCouplingConstitutiveLaw - end */


/* CubicElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class CubicElasticGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	CubicElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const T& Stiff1, const T& Stiff2, const T& Stiff3)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ConstitutiveLaw<T, Tder>::Err(std::cerr, "cubic elastic generic constitutive law "
						"is allowed only for scalar and 3x3"));
	};

	virtual ~CubicElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

template <>
class CubicElasticGenericConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
private:
	doublereal Stiff1;
	doublereal Stiff2;
	doublereal Stiff3;

public:
	CubicElasticGenericConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress, const doublereal& Stiff1,
			const doublereal& Stiff2, const doublereal& Stiff3)
	: ElasticConstitutiveLaw1D(pDC, PStress),
		Stiff1(Stiff1), Stiff2(Stiff2), Stiff3(Stiff3)
	{
		NO_OP;
	};

	virtual ~CubicElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw1D* pCopy(void) const {
		ConstitutiveLaw1D* pCL = 0;

		typedef CubicElasticGenericConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw1D::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw1D::PreStress,
					Stiff1, Stiff2, Stiff3));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "cubic elastic generic, ",
			Write(out, Stiff1, ", ") << ", ",
			Write(out, Stiff2, ", ") << ", ",
			Write(out, Stiff3, ", ") << ", ";
		return ElasticConstitutiveLaw1D::Restart_int(out);
	};

	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.) {
		ConstitutiveLaw1D::Epsilon = Eps;
		doublereal e1 = Eps - ElasticConstitutiveLaw1D::Get();
		doublereal f1 = fabs(e1);
		doublereal e2 = e1*e1;
		doublereal f2 = f1*e1;
		doublereal e3 = e2*e1;
		ConstitutiveLaw1D::FDE = Stiff1 + 2.*Stiff2*f1 + 3.*Stiff3*e2;
		ConstitutiveLaw1D::F = ElasticConstitutiveLaw1D::PreStress
			+ Stiff1*e1 + Stiff2*f2 + Stiff3*e3;
	};
};

template <>
class CubicElasticGenericConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
private:
	Vec3 Stiff1;
	Vec3 Stiff2;
	Vec3 Stiff3;

public:
	CubicElasticGenericConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
			const Vec3& PStress, const Vec3& Stiff1,
			const Vec3& Stiff2, const Vec3& Stiff3)
	: ElasticConstitutiveLaw3D(pDC, PStress),
		Stiff1(Stiff1), Stiff2(Stiff2), Stiff3(Stiff3)
	{
		NO_OP;
	};

	virtual ~CubicElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw3D* pCopy(void) const {
		ConstitutiveLaw3D* pCL = 0;

		typedef CubicElasticGenericConstitutiveLaw<Vec3, Mat3x3> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw3D::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw3D::PreStress,
					Stiff1, Stiff2, Stiff3));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "cubic elastic generic, ",
			Write(out, Stiff1, ", ") << ", ",
			Write(out, Stiff2, ", ") << ", ",
			Write(out, Stiff3, ", ") << ", ";
		return ElasticConstitutiveLaw3D::Restart_int(out);
	};

	virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */ = Zero3) {
		ConstitutiveLaw3D::Epsilon = Eps;
		Vec3 v1 = Eps - ElasticConstitutiveLaw3D::Get();
		ConstitutiveLaw3D::F = ElasticConstitutiveLaw3D::PreStress;

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
		Vec3 FTmp;
#endif // MBDYN_X_WORKAROUND_GCC_3_2 || MBDYN_X_WORKAROUND_GCC_3_3

		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			doublereal e1 = v1(iCnt);
			doublereal f1 = fabs(e1);
			doublereal e2 = e1*e1;
			doublereal f2 = f1*e1;
			doublereal e3 = e2*e1;

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
			ConstitutiveLaw3D::FDE.Put(iCnt, iCnt,
				Stiff1(iCnt)
				+ 2.*Stiff2(iCnt)*f1
				+ 3.*Stiff3(iCnt)*e2);
			FTmp(iCnt) = Stiff1(iCnt)*e1
				+ Stiff2(iCnt)*f2
				+ Stiff3(iCnt)*e3;
#else // ! MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
			ConstitutiveLaw3D::FDE(iCnt, iCnt) = Stiff1(iCnt) + 2.*Stiff2(iCnt)*f1 + 3.*Stiff3(iCnt)*e2;
			ConstitutiveLaw3D::F(iCnt) += Stiff1(iCnt)*e1 + Stiff2(iCnt)*f2 + Stiff3(iCnt)*e3;
#endif // ! MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
		}

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
		ConstitutiveLaw3D::F += FTmp;
#endif // MBDYN_X_WORKAROUND_GCC_3_3 || MBDYN_X_WORKAROUND_GCC_3_3
	};
};

/* CubicElasticGenericConstitutiveLaw - end */


/* LogConstitutiveLaw - begin */

template <class T, class Tder>
class LogConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	LogConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			doublereal = 0.)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ElasticConstitutiveLaw<T, Tder>::Err(std::cerr, "log constitutive law is allowed only for rods"));
	};

	virtual ~LogConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

template<>
class LogConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
private:
	doublereal dStiffness;
	doublereal dCurrEps;

public:
	LogConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress, doublereal dStiff)
	: ElasticConstitutiveLaw1D(pDC, PStress),
	dStiffness(dStiff), dCurrEps(0.) {
		ASSERT(Get() < 1.);
	};

	virtual ~LogConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		typedef LogConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					dStiffness));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "log elastic, " << dStiffness;
		return Restart_int(out);
	};

	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.) {
		Epsilon = Eps;

		doublereal dPreStrain = Get();
		dCurrEps = 1. + Epsilon - dPreStrain;
		ASSERT(dCurrEps > std::numeric_limits<doublereal>::epsilon());

		if (dCurrEps < std::numeric_limits<doublereal>::epsilon()) {
			// throw ErrGeneric();
			dCurrEps = std::numeric_limits<doublereal>::epsilon();
		}

		F = PreStress + dStiffness*log(dCurrEps);
		FDE = dStiffness/dCurrEps;
	};
};

/* LogConstitutiveLaw - end */


/* DoubleLinearElasticConstitutiveLaw - begin */

template <class T, class Tder>
class DoubleLinearElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			doublereal = 0.,
			doublereal = 0.,
			doublereal = 0.,
			doublereal = 0.)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ElasticConstitutiveLaw<T, Tder>::Err(std::cerr, "double linear elastic constitutive law "
						"is allowed only for rods and 3D hinges"));
	};

	virtual ~DoubleLinearElasticConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

template<>
class DoubleLinearElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
private:
	doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
	doublereal dUpperLimitStrain;
	doublereal dLowerLimitStrain;
	doublereal dSecondStiffness;
	flag fSecondStiff;

public:
	DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			doublereal dStiff,
			doublereal dUppLimStrain,
			doublereal dLowLimStrain,
			doublereal dSecondStiff)
	: ElasticConstitutiveLaw1D(pDC, PStress),
	dStiffness(dStiff),
	dUpperLimitStrain(dUppLimStrain),
	dLowerLimitStrain(dLowLimStrain),
	dSecondStiffness(dSecondStiff) {
		FDE = dStiffness;
	};

	virtual ~DoubleLinearElasticConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		typedef DoubleLinearElasticConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					dStiffness,
					dUpperLimitStrain,
					dLowerLimitStrain,
					dSecondStiffness));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "double linear elastic, "
			<< dStiffness << ", "
			<< dUpperLimitStrain << ", "
			<< dLowerLimitStrain << ", "
			<< dSecondStiffness;
		return Restart_int(out);
	};

	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */ = 0.) {
		Epsilon = Eps;

		doublereal dPreStrain = Get();
		doublereal dCurrStrain = Epsilon-dPreStrain;
		if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
			FDE = dStiffness;
			F = PreStress + dStiffness*dCurrStrain;
		} else {
			FDE = dSecondStiffness;

			if (dCurrStrain > dUpperLimitStrain) {
				F = PreStress + dStiffness*dUpperLimitStrain
					+ dSecondStiffness*(dCurrStrain - dUpperLimitStrain);
			} else /* if (dCurrStrain < dLowerLimitStrain) */ {
				F = PreStress + dStiffness*dLowerLimitStrain
					+ dSecondStiffness*(dCurrStrain - dLowerLimitStrain);
			}
		}
	};
};


template<>
class DoubleLinearElasticConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
private:
	doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
	doublereal dUpperLimitStrain;
	doublereal dLowerLimitStrain;
	doublereal dSecondStiffness;

public:
	DoubleLinearElasticConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
			const Vec3& PStress,
			doublereal dStiff,
			doublereal dUppLimStrain,
			doublereal dLowLimStrain,
			doublereal dSecondStiff)
	: ElasticConstitutiveLaw3D(pDC, PStress),
	dStiffness(dStiff),
	dUpperLimitStrain(dUppLimStrain),
	dLowerLimitStrain(dLowLimStrain),
	dSecondStiffness(dSecondStiff) {
		Mat3x3DEye.Manipulate(FDE, dStiffness);
	};

	virtual ~DoubleLinearElasticConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		typedef DoubleLinearElasticConstitutiveLaw<Vec3, Mat3x3> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					dStiffness,
					dUpperLimitStrain,
					dLowerLimitStrain,
					dSecondStiffness));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "double linear elastic, "
			<< dStiffness << ", "
			<< dUpperLimitStrain << ", "
			<< dLowerLimitStrain << ", "
			<< dSecondStiffness;
		return Restart_int(out);
	};

	virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */ = Zero3) {
		Epsilon = Eps;

		Vec3 PreStrain = Get();
		Vec3 CurrStrain = Epsilon-PreStrain;
		doublereal dCurrStrain = CurrStrain.dGet(3);

		if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
			FDE.Put(3, 3, dStiffness);
			F = PreStress+CurrStrain*dStiffness;
		} else {
			FDE.Put(3, 3, dSecondStiffness);

			if (dCurrStrain > dUpperLimitStrain) {
				F = PreStress + Vec3(CurrStrain.dGet(1)*dStiffness,
						CurrStrain.dGet(2)*dStiffness,
						dUpperLimitStrain*dStiffness
						+ (dCurrStrain - dUpperLimitStrain)*dSecondStiffness);
			} else /* if (dCurrStrain < dLowerLimitStrain) */ {
				F = PreStress + Vec3(CurrStrain.dGet(1)*dStiffness,
						CurrStrain.dGet(2)*dStiffness,
						dLowerLimitStrain*dStiffness
						+ (dCurrStrain - dLowerLimitStrain)*dSecondStiffness);
			}
		}
	};
};

/* DoubleLinearElasticConstitutiveLaw - end */


/* IsotropicHardeningConstitutiveLaw - begin */

template <class T, class Tder>
class IsotropicHardeningConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
private:
	doublereal dStiffness;
	doublereal dAlpha;
	doublereal dBeta;

	/* Legge costitutiva:
	 *
	 *            Beta + Alpha*x^2
	 *    f = k * ---------------- * x
	 *              1 + Alpha*x^2
	 */

public:
	IsotropicHardeningConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			doublereal dStiff,
			doublereal dStiff0,
			doublereal dEpsHard)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
	dStiffness(dStiff), dAlpha(0.), dBeta(0.)
	{
		ASSERT(dEpsHard > std::numeric_limits<doublereal>::epsilon());
		ASSERT(dStiff > std::numeric_limits<doublereal>::epsilon());

		dAlpha = 3./(dEpsHard*dEpsHard);
		dBeta = dStiff0/dStiff;
	};

	virtual ~IsotropicHardeningConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef IsotropicHardeningConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw<T, Tder>::PreStress,
					dStiffness,
					dBeta*dStiffness,
					sqrt(3./dAlpha)));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "isotropic hardening elastic, " << dStiffness << ", "
			<< sqrt(3./dAlpha);
		if (dBeta != 0.) {
			out << ", linear stiffness, " << dBeta*dStiffness;
		}
		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */ = Zero3) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		T x = ConstitutiveLaw<T, Tder>::Epsilon - ElasticConstitutiveLaw<T, Tder>::Get();
		doublereal dx2 = x*x;
		doublereal dDen = 1. + dAlpha*dx2;
		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
			+ x*(dStiffness*(dBeta + dAlpha*dx2)/dDen);
		mb_deye<Tder>(ConstitutiveLaw<T, Tder>::FDE, dStiffness*(dBeta + (3. - dBeta + dAlpha*dx2)*dAlpha*dx2)/(dDen*dDen));
	};
};

typedef IsotropicHardeningConstitutiveLaw<doublereal, doublereal> IsotropicHardeningConstitutiveLaw1D;
typedef IsotropicHardeningConstitutiveLaw<Vec3, Mat3x3> IsotropicHardeningConstitutiveLaw3D;
typedef IsotropicHardeningConstitutiveLaw<Vec6, Mat6x6> IsotropicHardeningConstitutiveLaw6D;

/* IsotropicHardeningConstitutiveLaw - end */


/* ContactConstitutiveLaw - begin */

template <class T, class Tder>
class ContactConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {

public:
	ContactConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			const doublereal = 0.,
			const doublereal = 0.)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ElasticConstitutiveLaw<T, Tder>::Err(std::cerr, "contact constitutive law "
						"is allowed only for rods"));
	};

	virtual ~ContactConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

/*
 *
 * F = k * ( 1-l_0 / l ) ^ gamma
 *     ==>
 * F = k * ( eps / ( 1 + eps ) ) ^ gamma
 *
 * d F = gamma * k * ( eps / ( 1 + eps ) ) ^ ( gamma - 1 ) * 1 / ( 1 + eps ) ^ 2 d eps
 *
 */

template<>
class ContactConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
private:
	doublereal dKappa;
	doublereal dGamma;

public:
	ContactConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress,
			const doublereal& dKappa,
			const doublereal& dGamma)
	: ElasticConstitutiveLaw1D(pDC, PStress), dKappa(dKappa), dGamma(dGamma) {
		NO_OP;
	};

	virtual ~ContactConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		typedef ContactConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					dKappa,
					dGamma));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "contact elastic, "
			<< dKappa << ", "
			<< dGamma;
		return Restart_int(out);
	};

	virtual void Update(const doublereal& Eps, const doublereal& /* EpsPrime */  = 0.) {
		doublereal dE;

		Epsilon = Eps;
		dE = Epsilon-Get();
		if ( dE >= 0. ) {
			F = PreStress;
			FDE = 0.;
		} else {
			F = PreStress+dKappa*(1. - 1./pow(1. + dE, dGamma));
			FDE = dGamma*dKappa/pow(1. + dE, dGamma + 1.);
		}
	};
};


template<>
class ContactConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
private:
	doublereal dKappa;
	doublereal dGamma;

public:
	ContactConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
			const Vec3& PStress,
			const doublereal& dKappa,
			const doublereal& dGamma)
	: ElasticConstitutiveLaw3D(pDC, PStress), dKappa(dKappa), dGamma(dGamma) {
		F = Zero3;
		FDE = Zero3x3;
	};

	virtual ~ContactConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		typedef ContactConstitutiveLaw<Vec3, Mat3x3> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(pGetDriveCaller()->pCopy(),
					PreStress,
					dKappa,
					dGamma));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "contact elastic, "
			<< dKappa << ", "
			<< dGamma;
		return Restart_int(out);
	};

	virtual void Update(const Vec3& Eps, const Vec3& /* EpsPrime */  = Zero3) {
		doublereal dE;

		Epsilon = Eps;
		dE = Epsilon.dGet(3)-Get().dGet(3);
		if ( dE >= 0. ) {
			F.Put(3, 0.);
			FDE.Put(3, 3, 0.);
		} else {
			F.Put(3, dKappa*(1. - 1./pow(1. + dE, dGamma)));
			FDE.Put(3, 3, dGamma*dKappa/pow(1. + dE, dGamma + 1.));
		}
	};
};

typedef ContactConstitutiveLaw<doublereal, doublereal> ContactConstitutiveLaw1D;
typedef ContactConstitutiveLaw<Vec3, Mat3x3> ContactConstitutiveLaw3D;
typedef ContactConstitutiveLaw<Vec6, Mat6x6> ContactConstitutiveLaw6D;

/* ContactConstitutiveLaw - end */


/* LinearViscousIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscousIsotropicConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffnessPrime;  /* Isotropa: Eye*dStiffnessPrime */

 public:
   LinearViscousIsotropicConstitutiveLaw(const T& PStress,
					 doublereal dStiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(0, PStress),
     dStiffnessPrime(dStiffPrime) {
      mb_deye<Tder>(ConstitutiveLaw<T, Tder>::FDEPrime, dStiffnessPrime);
   };

   virtual ~LinearViscousIsotropicConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = 0;

      typedef LinearViscousIsotropicConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(ElasticConstitutiveLaw<T, Tder>::PreStress,
                               dStiffnessPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscous isotropic, "
        << dStiffnessPrime;
      return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
   };

   virtual void Update(const T& /* Eps */ , const T& EpsPrime = mb_zero<T>()) {
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::EpsilonPrime*dStiffnessPrime;
   };
};

/* LinearViscousIsotropicConstitutiveLaw - end */


/* LinearViscousGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscousGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   LinearViscousGenericConstitutiveLaw(const T& PStress,
				       const Tder& StiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(0, PStress) {
      ConstitutiveLaw<T, Tder>::FDEPrime = StiffPrime;
   };

   virtual ~LinearViscousGenericConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = 0;

      typedef LinearViscousGenericConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(ElasticConstitutiveLaw<T, Tder>::PreStress,
                               ConstitutiveLaw<T, Tder>::FDEPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscous generic, ",
        Write(out, ConstitutiveLaw<T, Tder>::FDEPrime, ", ");
      return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
   };

   virtual void Update(const T& /* Eps */ , const T& EpsPrime = mb_zero<T>()) {
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
   };
};

/* LinearViscousGenericConstitutiveLaw - end */


/* LinearViscoElasticIsotropicConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticIsotropicConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 private:
   doublereal dStiffness;
   doublereal dStiffnessPrime;

 public:
   LinearViscoElasticIsotropicConstitutiveLaw(const TplDriveCaller<T>* pDC,
					      const T& PStress,
					      doublereal dStiff,
					      doublereal dStiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
     dStiffness(dStiff), dStiffnessPrime(dStiffPrime) {
      mb_deye<Tder>(ConstitutiveLaw<T, Tder>::FDE, dStiffness);
      mb_deye<Tder>(ConstitutiveLaw<T, Tder>::FDEPrime, dStiffnessPrime);
   };

   virtual ~LinearViscoElasticIsotropicConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = 0;

      typedef LinearViscoElasticIsotropicConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
                               ElasticConstitutiveLaw<T, Tder>::PreStress,
                               dStiffness,
                               dStiffnessPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "linear viscoelastic isotropic, "
	<< dStiffness << ", "
	<< dStiffnessPrime;
      return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
   };

   virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get())*dStiffness+ConstitutiveLaw<T, Tder>::EpsilonPrime*dStiffnessPrime;
   };
};

typedef LinearViscoElasticIsotropicConstitutiveLaw<doublereal, doublereal> LinearViscoElasticIsotropicConstitutiveLaw1D;
typedef LinearViscoElasticIsotropicConstitutiveLaw<Vec3, Mat3x3> LinearViscoElasticIsotropicConstitutiveLaw3D;
typedef LinearViscoElasticIsotropicConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticIsotropicConstitutiveLaw6D;

/* LinearViscoElasticIsotropicConstitutiveLaw - end */


/* LinearViscoElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   LinearViscoElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
					    const T& PStress,
					    const Tder& Stiff,
					    const Tder& StiffPrime)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      ConstitutiveLaw<T, Tder>::FDE = Stiff;
      ConstitutiveLaw<T, Tder>::FDEPrime = StiffPrime;
   };

   virtual ~LinearViscoElasticGenericConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      ConstitutiveLaw<T, Tder>* pCL = 0;

      typedef LinearViscoElasticGenericConstitutiveLaw<T, Tder> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
                               ElasticConstitutiveLaw<T, Tder>::PreStress,
                               ConstitutiveLaw<T, Tder>::FDE,
                               ConstitutiveLaw<T, Tder>::FDEPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
     out << "linear viscoelastic generic, ",
       Write(out, ConstitutiveLaw<T, Tder>::FDE, ", ") << ", ",
       Write(out, ConstitutiveLaw<T, Tder>::FDEPrime, ", ");
       return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
   };

   virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
      ConstitutiveLaw<T, Tder>::Epsilon = Eps;
      ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
      ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
	+ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon-ElasticConstitutiveLaw<T, Tder>::Get())
	+ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
   };
};

typedef LinearViscoElasticGenericConstitutiveLaw<doublereal, doublereal> LinearViscoElasticGenericConstitutiveLaw1D;
typedef LinearViscoElasticGenericConstitutiveLaw<Vec3, Mat3x3> LinearViscoElasticGenericConstitutiveLaw3D;
typedef LinearViscoElasticGenericConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticGenericConstitutiveLaw6D;

/* LinearViscoElasticGenericConstitutiveLaw - end */

/* LTVViscoElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class LTVViscoElasticGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
protected:
	DriveOwner FDECoef;
	Tder FDERef;
	doublereal dPrevScaleFactor;
	DriveOwner FDEPrimeCoef;
	Tder FDEPrimeRef;
	doublereal dPrevScaleFactorPrime;

public:
	LTVViscoElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
		const T& PStress,
		const Tder& Stiff,
		const DriveCaller *pdc,
		const Tder& StiffPrime,
		const DriveCaller *pdcp)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress),
	FDECoef(pdc), FDEPrimeCoef(pdcp) {
		FDERef = Stiff;
		dPrevScaleFactor = 0.;
		FDEPrimeRef = StiffPrime;
		dPrevScaleFactorPrime = 0.;
	};

	virtual ~LTVViscoElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef LTVViscoElasticGenericConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
			cl,
			cl(ElasticConstitutiveLaw<T, Tder>::pGetDriveCaller()->pCopy(),
				ElasticConstitutiveLaw<T, Tder>::PreStress,
				FDERef,
				FDECoef.pGetDriveCaller()->pCopy(),
				FDEPrimeRef,
				FDEPrimeCoef.pGetDriveCaller()->pCopy()));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "linear time variant viscoelastic generic, ",
			Write(out, FDERef, ", ") << ", ",
			FDECoef.pGetDriveCaller()->Restart(out) << ", ",
			Write(out, FDEPrimeRef, ", ") << ", ",
			FDEPrimeCoef.pGetDriveCaller()->Restart(out);
		return ElasticConstitutiveLaw<T, Tder>::Restart_int(out);
	};

	virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
		doublereal dCurrScaleFactor = FDECoef.dGet();
		if (dCurrScaleFactor != dPrevScaleFactor) {
			dPrevScaleFactor = dCurrScaleFactor;
			ConstitutiveLaw<T, Tder>::FDE = FDERef*dCurrScaleFactor;
		}
		doublereal dCurrScaleFactorPrime = FDEPrimeCoef.dGet();
		if (dCurrScaleFactorPrime != dPrevScaleFactorPrime) {
			dPrevScaleFactorPrime = dCurrScaleFactorPrime;
			ConstitutiveLaw<T, Tder>::FDEPrime = FDEPrimeRef*dCurrScaleFactorPrime;
		}
		ConstitutiveLaw<T, Tder>::F = ElasticConstitutiveLaw<T, Tder>::PreStress
			+ ConstitutiveLaw<T, Tder>::FDE*(ConstitutiveLaw<T, Tder>::Epsilon - ElasticConstitutiveLaw<T, Tder>::Get())
			+ ConstitutiveLaw<T, Tder>::FDEPrime*ConstitutiveLaw<T, Tder>::EpsilonPrime;
	};
};

typedef LTVViscoElasticGenericConstitutiveLaw<doublereal, doublereal> LTVViscoElasticGenericConstitutiveLaw1D;
typedef LTVViscoElasticGenericConstitutiveLaw<Vec3, Mat3x3> LTVViscoElasticGenericConstitutiveLaw3D;
typedef LTVViscoElasticGenericConstitutiveLaw<Vec6, Mat6x6> LTVViscoElasticGenericConstitutiveLaw6D;

/* LTVViscoElasticGenericConstitutiveLaw - end */

/* LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw - begin */

template <class T, class Tder>
class LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress,
			const Tder& Stiff,
			const Tder& StiffPrime,
			doublereal dAxTors)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ConstitutiveLaw<T, Tder>::Err(std::cerr, "axial-torsion coupling constitutive law "
					"is allowed only for beams (6x6)"));
	};

	virtual ~LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		NO_OP;
	};
};

template <>
class LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6>
: public ElasticConstitutiveLaw<Vec6, Mat6x6> {
private:
	doublereal dRefTorsion;
	doublereal dAxialTorsionCoupling;

public:
	LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw(const TplDriveCaller<Vec6>* pDC,
			const Vec6& PStress,
			const Mat6x6& Stiff,
			const Mat6x6& StiffPrime,
			doublereal dAxTors)
	: ElasticConstitutiveLaw<Vec6, Mat6x6>(pDC, PStress),
	dRefTorsion(Stiff(4, 4)),
	dAxialTorsionCoupling(dAxTors) {
		ConstitutiveLaw<Vec6, Mat6x6>::FDE = Stiff;
		ConstitutiveLaw<Vec6, Mat6x6>::FDEPrime = StiffPrime;
	};

	virtual ~LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy(void) const {
		ConstitutiveLaw<Vec6, Mat6x6>* pCL = 0;

		typedef LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw<Vec6, Mat6x6>::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw<Vec6, Mat6x6>::PreStress,
					ConstitutiveLaw<Vec6, Mat6x6>::FDE,
					ConstitutiveLaw<Vec6, Mat6x6>::FDEPrime,
					dAxialTorsionCoupling));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		doublereal d = FDE(4, 4);
		out << "linear viscoelastic generic axial torsion coupling, ",
		Write(out, ConstitutiveLaw<Vec6, Mat6x6>::FDE, ", ") << ", ",
		Write(out, ConstitutiveLaw<Vec6, Mat6x6>::FDEPrime, ", ") << ", " << dAxialTorsionCoupling;
		ElasticConstitutiveLaw<Vec6, Mat6x6>::Restart_int(out);
		((Mat6x6&)ConstitutiveLaw<Vec6, Mat6x6>::FDE)(4, 4) = d;
		return out;
	};

	virtual void Update(const Vec6& Eps, const Vec6& EpsPrime = mb_zero<Vec6>()) {
		ConstitutiveLaw<Vec6, Mat6x6>::Epsilon = Eps;
		ConstitutiveLaw<Vec6, Mat6x6>::EpsilonPrime = EpsPrime;
		doublereal d = Epsilon.dGet(1);
		((Mat6x6&)ConstitutiveLaw<Vec6, Mat6x6>::FDE)(4, 4) = dRefTorsion + d*dAxialTorsionCoupling;
		ConstitutiveLaw<Vec6, Mat6x6>::F = ElasticConstitutiveLaw<Vec6, Mat6x6>::PreStress
			+ ConstitutiveLaw<Vec6, Mat6x6>::FDE*(ConstitutiveLaw<Vec6, Mat6x6>::Epsilon - ElasticConstitutiveLaw<Vec6, Mat6x6>::Get())
			+ ConstitutiveLaw<Vec6, Mat6x6>::FDEPrime*ConstitutiveLaw<Vec6, Mat6x6>::EpsilonPrime;
	};
};

typedef LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw<Vec6, Mat6x6> LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw6D;

/* LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw - end */


/* CubicViscoElasticGenericConstitutiveLaw - begin */

template <class T, class Tder>
class CubicViscoElasticGenericConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
public:
	CubicViscoElasticGenericConstitutiveLaw(const TplDriveCaller<T>* pDC,
			const T& PStress, const T& Stiff1, const T& Stiff2, const T& Stiff3,
			const Tder& StiffPrime)
	: ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
		throw (typename ConstitutiveLaw<T, Tder>::Err(std::cerr, "cubic viscoelastic generic constitutive law "
						"is allowed only for scalar and 3x3"));
	};

	virtual ~CubicViscoElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */ = mb_zero<T>()) {
		NO_OP;
	};
};

template <>
class CubicViscoElasticGenericConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
private:
	doublereal Stiff1;
	doublereal Stiff2;
	doublereal Stiff3;

public:
	CubicViscoElasticGenericConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
			const doublereal& PStress, const doublereal& Stiff1,
			const doublereal& Stiff2, const doublereal& Stiff3,
			const doublereal& StiffPrime)
	: ElasticConstitutiveLaw1D(pDC, PStress),
		Stiff1(Stiff1), Stiff2(Stiff2), Stiff3(Stiff3)
	{
		ConstitutiveLaw1D::FDEPrime = StiffPrime;
	};

	virtual ~CubicViscoElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw1D* pCopy(void) const {
		ConstitutiveLaw1D* pCL = 0;

		typedef CubicViscoElasticGenericConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw1D::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw1D::PreStress,
					Stiff1, Stiff2, Stiff3, ConstitutiveLaw1D::FDEPrime));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "cubic elastic generic, ",
			Write(out, Stiff1, ", ") << ", ",
			Write(out, Stiff2, ", ") << ", ",
			Write(out, Stiff3, ", ") << ", ",
			Write(out, ConstitutiveLaw1D::FDEPrime, ", ") << ", ";
		return ElasticConstitutiveLaw1D::Restart_int(out);
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		ConstitutiveLaw1D::Epsilon = Eps;
		ConstitutiveLaw1D::EpsilonPrime = EpsPrime;
		doublereal e1 = Eps - ElasticConstitutiveLaw1D::Get();
		doublereal f1 = fabs(e1);
		doublereal e2 = e1*e1;
		doublereal f2 = f1*e1;
		doublereal e3 = e2*e1;
		ConstitutiveLaw1D::FDE = Stiff1 + 2.*Stiff2*f1 + 3.*Stiff3*e2;
		ConstitutiveLaw1D::F = ElasticConstitutiveLaw1D::PreStress
			+ Stiff1*e1 + Stiff2*f2 + Stiff3*e3 + ConstitutiveLaw1D::FDEPrime*EpsPrime;
	};
};

template <>
class CubicViscoElasticGenericConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
private:
	Vec3 Stiff1;
	Vec3 Stiff2;
	Vec3 Stiff3;

public:
	CubicViscoElasticGenericConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
			const Vec3& PStress, const Vec3& Stiff1,
			const Vec3& Stiff2, const Vec3& Stiff3,
			const Mat3x3& StiffPrime)
	: ElasticConstitutiveLaw3D(pDC, PStress),
		Stiff1(Stiff1), Stiff2(Stiff2), Stiff3(Stiff3)
	{
		ConstitutiveLaw3D::FDEPrime = StiffPrime;
	};

	virtual ~CubicViscoElasticGenericConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw3D* pCopy(void) const {
		ConstitutiveLaw3D* pCL = 0;

		typedef CubicViscoElasticGenericConstitutiveLaw<Vec3, Mat3x3> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
				cl,
				cl(ElasticConstitutiveLaw3D::pGetDriveCaller()->pCopy(),
					ElasticConstitutiveLaw3D::PreStress,
					Stiff1, Stiff2, Stiff3, ConstitutiveLaw3D::FDEPrime));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "cubic elastic generic, ",
			Write(out, Stiff1, ", ") << ", ",
			Write(out, Stiff2, ", ") << ", ",
			Write(out, Stiff3, ", ") << ", ",
			Write(out, ConstitutiveLaw3D::FDEPrime, ", ") << ", ";
		return ElasticConstitutiveLaw3D::Restart_int(out);
	};

	virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = Zero3) {
		ConstitutiveLaw3D::Epsilon = Eps;
		ConstitutiveLaw3D::EpsilonPrime = EpsPrime;
		Vec3 v1 = Eps - ElasticConstitutiveLaw3D::Get();
		ConstitutiveLaw3D::F = ElasticConstitutiveLaw3D::PreStress + ConstitutiveLaw3D::FDEPrime*EpsPrime;

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
		Vec3 FTmp;
#endif // MBDYN_X_WORKAROUND_GCC_3_2 || MBDYN_X_WORKAROUND_GCC_3_3

		for (int iCnt = 1; iCnt <= 3; iCnt++) {
			doublereal e1 = v1(iCnt);
			doublereal f1 = fabs(e1);
			doublereal e2 = e1*e1;
			doublereal f2 = f1*e1;
			doublereal e3 = e2*e1;
			
#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
			ConstitutiveLaw3D::FDE.Put(iCnt, iCnt,
				Stiff1(iCnt)
				+ 2.*Stiff2(iCnt)*f1
				+ 3.*Stiff3(iCnt)*e2);
			FTmp(iCnt) = Stiff1(iCnt)*e1
				+ Stiff2(iCnt)*f2
				+ Stiff3(iCnt)*e3;
#else // ! MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
			ConstitutiveLaw3D::FDE(iCnt, iCnt) = Stiff1(iCnt)
				+ 2.*Stiff2(iCnt)*f1 + 3.*Stiff3(iCnt)*e2;
			ConstitutiveLaw3D::F(iCnt) += Stiff1(iCnt)*e1
				+ Stiff2(iCnt)*f2 + Stiff3(iCnt)*e3;
#endif // ! MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
		}

#if defined(MBDYN_X_WORKAROUND_GCC_3_2) || defined(MBDYN_X_WORKAROUND_GCC_3_3)
		ConstitutiveLaw3D::F += FTmp;
#endif // ! MBDYN_X_WORKAROUND_GCC_3_2 && ! MBDYN_X_WORKAROUND_GCC_3_3
	};
};

/* CubicViscoElasticGenericConstitutiveLaw - end */



/* DoubleLinearViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class DoubleLinearViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
					   const T& PStress,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.,
					   doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      throw (typename ElasticConstitutiveLaw<T, Tder>::Err(std::cerr, "doublelinear viscoelastic constitutive law "
			      "is allowed only for rods ad 3D hinges"));
   };

   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return 0;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };

   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
      NO_OP;
   };
};


template<>
class DoubleLinearViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   doublereal dStiffnessPrime;
   doublereal dSecondStiffnessPrime;

 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
					   const doublereal& PStress,
					   doublereal dStiff,
					   doublereal dUpp,
					   doublereal dLow,
					   doublereal dSecondS,
					   doublereal dStiffPrime,
					   doublereal dSecondSPrime)
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff),
     dUpperLimitStrain(dUpp), dLowerLimitStrain(dLow),
     dSecondStiffness(dSecondS),
     dStiffnessPrime(dStiffPrime),
     dSecondStiffnessPrime(dSecondSPrime)
   {
      FDEPrime = dStiffnessPrime;
   };

   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

      typedef DoubleLinearViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress,
                               dStiffness,
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness,
                               dStiffnessPrime,
			       dSecondStiffnessPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear viscoelastic, "
	<< dStiffness << ", "
	<< dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness << ", "
	<< dStiffnessPrime << ", "
	"second damping, " << dSecondStiffnessPrime << ", ";
      return Restart_int(out);
   };

   virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
      Epsilon = Eps;
      EpsilonPrime = EpsPrime;

      doublereal dPreStrain = Get();
      doublereal dCurrStrain = Epsilon-dPreStrain;
      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE = dStiffness;
	 FDEPrime = dStiffnessPrime;
	 F = PreStress+dStiffness*dCurrStrain
	   +dStiffnessPrime*EpsilonPrime;
      } else {
	 FDE = dSecondStiffness;
	 FDEPrime = dSecondStiffnessPrime;

	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress+dStiffness*dUpperLimitStrain
	      +dSecondStiffness*(dCurrStrain-dUpperLimitStrain)
		+dSecondStiffnessPrime*EpsilonPrime;
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {
	    F = PreStress+dStiffness*dLowerLimitStrain
	      +dSecondStiffness*(dCurrStrain-dLowerLimitStrain)
		+dSecondStiffnessPrime*EpsilonPrime;
	 }
      }
   };
};


template<>
class DoubleLinearViscoElasticConstitutiveLaw<Vec3, Mat3x3>
: public ElasticConstitutiveLaw3D {
 private:
   doublereal dStiffness;        /* Isotropa: Eye*dStiffness */
   doublereal dUpperLimitStrain;
   doublereal dLowerLimitStrain;
   doublereal dSecondStiffness;
   doublereal dStiffnessPrime;
   doublereal dSecondStiffnessPrime;

 public:
   DoubleLinearViscoElasticConstitutiveLaw(const TplDriveCaller<Vec3>* pDC,
					   const Vec3& PStress,
					   doublereal dStiff,
					   doublereal dUppLimStrain,
					   doublereal dLowLimStrain,
					   doublereal dSecondStiff,
					   doublereal dStiffPrime,
					   doublereal dSecondStiffPrime)
     : ElasticConstitutiveLaw3D(pDC, PStress),
     dStiffness(dStiff),
     dUpperLimitStrain(dUppLimStrain),
     dLowerLimitStrain(dLowLimStrain),
     dSecondStiffness(dSecondStiff),
     dStiffnessPrime(dStiffPrime),
     dSecondStiffnessPrime(dSecondStiffPrime)
   {
      Mat3x3DEye.Manipulate(FDE, dStiffness);
      Mat3x3DEye.Manipulate(FDEPrime, dStiffnessPrime);
   };

   virtual ~DoubleLinearViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
      ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

      typedef DoubleLinearViscoElasticConstitutiveLaw<Vec3, Mat3x3> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress,
                               dStiffness,
                               dUpperLimitStrain,
                               dLowerLimitStrain,
                               dSecondStiffness,
                               dStiffnessPrime,
			       dSecondStiffnessPrime));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "double linear viscoelastic, "
        << dStiffness << ", "
	<< dUpperLimitStrain << ", "
	<< dLowerLimitStrain << ", "
	<< dSecondStiffness << ", "
	<< dStiffnessPrime << ", "
	<< dSecondStiffnessPrime << ", ";
      return Restart_int(out);
   };

   virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = Zero3) {
      Epsilon = Eps;
      EpsilonPrime = EpsPrime;

      Vec3 PreStrain = Get();
      Vec3 CurrStrain = Epsilon-PreStrain;
      doublereal dCurrStrain = CurrStrain.dGet(3);

      if (dCurrStrain <= dUpperLimitStrain && dCurrStrain >= dLowerLimitStrain) {
	 FDE.Put(3, 3, dStiffness);
	 FDEPrime.Put(3, 3, dStiffnessPrime);
	 F = PreStress+CurrStrain*dStiffness+EpsilonPrime*dStiffnessPrime;
      } else {
	 FDE.Put(3, 3, dSecondStiffness);
	 FDEPrime.Put(3, 3, dSecondStiffnessPrime);

	 if (dCurrStrain > dUpperLimitStrain) {
	    F = PreStress
	      +Vec3(CurrStrain(1)*dStiffness + EpsilonPrime(1)*dStiffnessPrime,
		    CurrStrain(2)*dStiffness + EpsilonPrime(2)*dStiffnessPrime,
		    dUpperLimitStrain*dStiffness
		        +(dCurrStrain-dUpperLimitStrain)*dSecondStiffness
	                +EpsilonPrime(3)*dSecondStiffnessPrime);
	 } else /* if (dCurrStrain < dLowerLimitStrain) */ {
	    F = PreStress
	      +Vec3(CurrStrain.dGet(1)*dStiffness + EpsilonPrime(1)*dStiffnessPrime,
		    CurrStrain.dGet(2)*dStiffness + EpsilonPrime(2)*dStiffnessPrime,
		    dLowerLimitStrain*dStiffness
		        +(dCurrStrain-dLowerLimitStrain)*dSecondStiffness
		        +EpsilonPrime(3)*dSecondStiffnessPrime);
	 }
      }
   };
};

/* DoubleLinearViscoElasticConstitutiveLaw - end */


/* TurbulentViscoElasticConstitutiveLaw - begin */

template <class T, class Tder>
class TurbulentViscoElasticConstitutiveLaw
: public ElasticConstitutiveLaw<T, Tder> {
 public:
   TurbulentViscoElasticConstitutiveLaw(const TplDriveCaller<T>* pDC,
					const T& PStress,
					doublereal = 0.,
					doublereal = 0.,
					doublereal = 0.,
					doublereal = 0.)
     : ElasticConstitutiveLaw<T, Tder>(pDC, PStress) {
      throw (typename ElasticConstitutiveLaw<T, Tder>::Err(std::cerr, "Turbulent viscoelastic constitutive law "
			      "is allowed only for rods"));
   };

   virtual ~TurbulentViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
      return 0;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      return out;
   };

   virtual void Update(const T& /* Eps */ , const T& /* EpsPrime */ = mb_zero<T>()) {
      NO_OP;
   };
};


template<>
class TurbulentViscoElasticConstitutiveLaw<doublereal, doublereal>
: public ElasticConstitutiveLaw1D {
 private:
   doublereal dStiffness;
   doublereal dStiffnessPrime;
   doublereal dTreshold;
   doublereal dParabolicStiffness;

 public:
   TurbulentViscoElasticConstitutiveLaw(const TplDriveCaller<doublereal>* pDC,
					const doublereal& PStress,
					doublereal dStiff,
					doublereal dStiffPrime,
					doublereal dTres,
					doublereal dParabStiff)
     : ElasticConstitutiveLaw1D(pDC, PStress),
     dStiffness(dStiff), dStiffnessPrime(dStiffPrime),
     dTreshold(dTres), dParabolicStiffness(dParabStiff) {
      FDE = dStiffness;
   };

   virtual ~TurbulentViscoElasticConstitutiveLaw(void) {
      NO_OP;
   };

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

   virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
      ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

      typedef TurbulentViscoElasticConstitutiveLaw<doublereal, doublereal> cl;
      SAFENEWWITHCONSTRUCTOR(pCL,
                            cl,
                            cl(pGetDriveCaller()->pCopy(),
                               PreStress,
                               dStiffness,
                               dStiffnessPrime,
                               dTreshold,
                               dParabolicStiffness));

      return pCL;
   };

   virtual std::ostream& Restart(std::ostream& out) const {
      out << "turbulent viscoelastic, "
	<< dStiffness << ", "
	<< dStiffnessPrime << ", "
	<< dTreshold << ", "
	<< dParabolicStiffness << ", ";
      return Restart_int(out);
   };

   virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
      Epsilon = Eps;
      EpsilonPrime = EpsPrime;

      doublereal dPreStrain = Get();

      doublereal d = fabs(EpsilonPrime);
      if (d < dTreshold) {
	 FDEPrime = dStiffnessPrime;
	 F = PreStress+dStiffness*(Epsilon-dPreStrain)
	   +dStiffnessPrime*EpsilonPrime;
      } else {
	 FDEPrime = 2.*dParabolicStiffness*d;
	 F = PreStress+dStiffness*(Epsilon-dPreStrain)
	   +dParabolicStiffness*d*EpsilonPrime;
      }
   };
};

/* TurbulentViscoElasticConstitutiveLaw - end */

/* BiStopCLWrapper - begin */

template <class T, class Tder>
class BiStopCLWrapper
: public ConstitutiveLaw<T, Tder> {
public:
	enum Status { INACTIVE, ACTIVE };

private:
	ConstitutiveLaw<T, Tder> *m_pCL;
	enum Status m_status;
	const DriveCaller *m_pActivatingCondition;
	const DriveCaller *m_pDeactivatingCondition;
	T m_EpsRef;

public:
	BiStopCLWrapper(
			ConstitutiveLaw<T, Tder> *pCL,
			enum Status initialStatus,
			const DriveCaller *pA,
			const DriveCaller *pD
	) : 
	m_pCL(pCL), m_status(initialStatus),
	m_pActivatingCondition(pA), m_pDeactivatingCondition(pD),
	m_EpsRef(mb_zero<T>()) {
		ASSERT(m_pActivatingCondition != 0);
		ASSERT(m_pDeactivatingCondition != 0);
	};

	virtual
	~BiStopCLWrapper(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return m_pCL->GetConstLawType();
	};

	virtual
	ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef BiStopCLWrapper<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL,
			cl,
			cl(m_pCL->pCopy(),
				m_status,
				m_pActivatingCondition->pCopy(),
				m_pDeactivatingCondition->pCopy()));

		return pCL;
	};

	virtual std::ostream&
	Restart(std::ostream& out) const {
		out << "bistop, initial status, ";
		if (m_status == INACTIVE) {
			out << "inactive";

		} else {
			out << "active";
		}
		out << ", ",
			m_pActivatingCondition->Restart(out) << ", ",
			m_pDeactivatingCondition->Restart(out) << ", ";
		return m_pCL->Restart(out);
	};

	virtual void
	Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
		bool ChangeJac(false);

		switch (m_status) {
		case INACTIVE:
			if (m_pActivatingCondition->dGet() == 0.) {
				/* remains inactive: nothing to do */
				break;
			}

			/* activates: change data and ask for jacobian rigeneration */
			m_status = ACTIVE;
			m_EpsRef = ConstitutiveLaw<T, Tder>::Epsilon;
			ChangeJac = true;

		case ACTIVE:
			if (m_pDeactivatingCondition->dGet() != 0.) {
				/* disactivates: reset data and ask for jacobian rigeneration */
				m_status = INACTIVE;
				ConstitutiveLaw<T, Tder>::F = ::mb_zero<T>();
				ConstitutiveLaw<T, Tder>::FDE = ::mb_zero<Tder>();
				ConstitutiveLaw<T, Tder>::FDEPrime = ::mb_zero<Tder>();
#if 0 // skip by now
				throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
#endif
			}

			/* change force as well */
			m_pCL->Update(Eps - m_EpsRef, EpsPrime);
			ConstitutiveLaw<T, Tder>::F = m_pCL->GetF();
			ConstitutiveLaw<T, Tder>::FDE = m_pCL->GetFDE();
			ConstitutiveLaw<T, Tder>::FDEPrime = m_pCL->GetFDEPrime();
			if (ChangeJac) {
				/* if activating, ask for jacobian rigeneration */
#if 0 // skip by now
				throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
#endif
			}
			break;
		}
	};

	virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		if (m_status == ACTIVE) {
			m_pCL->AfterConvergence(Eps - m_EpsRef, EpsPrime);
		}
	};

	// FIXME: OutputAppend() ?
};

/* BiStopCLWrapper - end */

/* helpers */
template <class T>
void
GetPreStress(MBDynParser& HP, T& PreStress)
{
	if (HP.IsKeyWord("prestress")) {
		PreStress = HP.Get(PreStress);
	}
}

template <class T>
TplDriveCaller<T>*
GetPreStrain(const DataManager* pDM, MBDynParser& HP)
{
	if (HP.IsKeyWord("prestrain")) {
		return HP.GetTplDriveCaller<T>();
	}

	TplDriveCaller<T> *pTplDC = 0;

	SAFENEW(pTplDC, ZeroTplDriveCaller<T>);

	return pTplDC;
}

#endif // CONSTLTP_IMPL_H
