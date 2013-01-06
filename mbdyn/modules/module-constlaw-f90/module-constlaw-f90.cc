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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"

#include "usrsub.h"

template <class T, class Tder>
class DummyConstitutiveLaw
: public ConstitutiveLaw<T, Tder> {
private:
public:
	DummyConstitutiveLaw(std::vector<doublereal>& v) {
		NO_OP;
	};

	virtual ~DummyConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;
	};

	virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		NO_OP;
	};
};

template <>
class DummyConstitutiveLaw<doublereal, doublereal>
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	integer size, err;
	std::vector<doublereal>	m_v;

public:
	DummyConstitutiveLaw(const std::vector<doublereal>& v)
	: size(v.size()), err(0), m_v(v)
	{
		us1init_(&size, &m_v[0], &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	virtual ~DummyConstitutiveLaw(void) {
		us1dstr_(&size, &m_v[0], &err);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		typedef DummyConstitutiveLaw<doublereal, doublereal> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_v));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "dummyf90";
		for (std::vector<doublereal>::const_iterator i = m_v.begin();
			i != m_v.end();
			i++)
		{
			out << ", " << *i;
		}
		return out;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

		doublereal dF, dFDE, dFDEP;
		us1updt_(&size, &m_v[0], &Eps, &EpsPrime, &dF, &dFDE, &dFDEP,
			&err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ConstitutiveLaw<doublereal, doublereal>::F = dF;
		ConstitutiveLaw<doublereal, doublereal>::FDE = dFDE;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = dFDEP;
	};

	virtual void AfterConvergence(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		us1aftc_(&size, &m_v[0], &Eps, &EpsPrime, &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};
};

template <>
class DummyConstitutiveLaw<Vec3, Mat3x3>
: public ConstitutiveLaw<Vec3, Mat3x3> {
private:
	integer size, err;
	std::vector<doublereal>	m_v;

public:
	DummyConstitutiveLaw(const std::vector<doublereal>& v)
	: size(v.size()), err(0), m_v(v)
	{
		us3init_(&size, &m_v[0], &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	virtual ~DummyConstitutiveLaw(void) {
		us3dstr_(&size, &m_v[0], &err);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		typedef DummyConstitutiveLaw<Vec3, Mat3x3> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_v));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "dummyf90";
		for (std::vector<doublereal>::const_iterator i = m_v.begin();
			i != m_v.end();
			i++)
		{
			out << ", " << *i;
		}
		return out;
	};

	virtual void Update(const Vec3& Eps, const Vec3& EpsPrime = mb_zero<Vec3>()) {
		ConstitutiveLaw<Vec3, Mat3x3>::Epsilon = Eps;

		doublereal dF[3], dFDE[3*3], dFDEP[3*3];
		us3updt_(&size, &m_v[0], Eps.pGetVec(), EpsPrime.pGetVec(),
			dF, dFDE, dFDEP, &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ConstitutiveLaw<Vec3, Mat3x3>::F = dF;
		ConstitutiveLaw<Vec3, Mat3x3>::FDE = Mat3x3(dFDE, 3);
		ConstitutiveLaw<Vec3, Mat3x3>::FDEPrime = Mat3x3(dFDEP, 3);
	};

	virtual void AfterConvergence(const Vec3& Eps, const Vec3& EpsPrime = mb_zero<Vec3>()) {
		us3aftc_(&size, &m_v[0], Eps.pGetVec(), EpsPrime.pGetVec(),
			&err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};
};

template <>
class DummyConstitutiveLaw<Vec6, Mat6x6>
: public ConstitutiveLaw<Vec6, Mat6x6> {
private:
	integer size, err;
	std::vector<doublereal>	m_v;

public:
	DummyConstitutiveLaw(const std::vector<doublereal>& v)
	: size(v.size()), err(0), m_v(v)
	{
		us6init_(&size, &m_v[0], &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	virtual ~DummyConstitutiveLaw(void) {
		us6dstr_(&size, &m_v[0], &err);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy(void) const {
		ConstitutiveLaw<Vec6, Mat6x6>* pCL = 0;

		typedef DummyConstitutiveLaw<Vec6, Mat6x6> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_v));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "dummyf90";
		for (std::vector<doublereal>::const_iterator i = m_v.begin();
			i != m_v.end();
			i++)
		{
			out << ", " << *i;
		}
		return out;
	};

	virtual void Update(const Vec6& Eps, const Vec6& EpsPrime = mb_zero<Vec6>()) {
		ConstitutiveLaw<Vec6, Mat6x6>::Epsilon = Eps;

		doublereal dE[6], dEP[6], dF[6], dFDE[6*6], dFDEP[6*6];
		Eps.PutTo(dE);
		EpsPrime.PutTo(dEP);
		us6updt_(&size, &m_v[0], dE, dEP, dF, dFDE, dFDEP, &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ConstitutiveLaw<Vec6, Mat6x6>::F = dF;
		ConstitutiveLaw<Vec6, Mat6x6>::FDE = Mat6x6(dFDE, 6);
		ConstitutiveLaw<Vec6, Mat6x6>::FDEPrime = Mat6x6(dFDEP, 6);
	};

	virtual void AfterConvergence(const Vec6& Eps, const Vec6& EpsPrime = mb_zero<Vec6>()) {
		doublereal dE[6], dEP[6];
		Eps.PutTo(dE);
		EpsPrime.PutTo(dEP);
		us6aftc_(&size, &m_v[0], dE, dEP, &err);
		if (err != 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct DummyCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::ELASTIC;

		std::vector<doublereal> v;
		int size = HP.GetInt();
		v.resize(size);
		for (int i = 0; i < size; i++) {
			v[i] = HP.GetReal();
		}

		typedef DummyConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(v));

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new DummyCLR<doublereal, doublereal>;
	if (!SetCL1D("dummyf90", rf1D)) {
		delete rf1D;

		silent_cerr("DummyConstitutiveLaw1D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new DummyCLR<Vec3, Mat3x3>;
	if (!SetCL3D("dummyf90", rf3D)) {
		delete rf3D;

		silent_cerr("DummyConstitutiveLaw3D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec6, Mat6x6> *rf6D = new DummyCLR<Vec6, Mat6x6>;
	if (!SetCL6D("dummyf90", rf6D)) {
		delete rf6D;

		silent_cerr("DummyConstitutiveLaw6D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

