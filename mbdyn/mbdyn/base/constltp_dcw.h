/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef CONSTLTP_DCW_H
#define CONSTLTP_DCW_H

#include "constltp_impl.h"

/* DriveCLWrapper - begin */

template <class T, class Tder>
class DriveCLWrapper
: public ConstitutiveLaw<T, Tder> {
private:
	ConstitutiveLaw<T, Tder> *m_pCL;
	DriveCaller *m_pDC;

public:
	DriveCLWrapper(ConstitutiveLaw<T, Tder> *pCL, DriveCaller *pDC)
	: m_pCL(pCL), m_pDC(pDC)
	{ 
		NO_OP;
	};

	virtual
	~DriveCLWrapper(void) {
		SAFEDELETE(m_pDC);
		SAFEDELETE(m_pCL);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return m_pCL->GetConstLawType();
	};

	virtual
	ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		typedef DriveCLWrapper<T, Tder> DCLW;
		SAFENEWWITHCONSTRUCTOR(pCL,
			DCLW,
			DCLW(m_pCL->pCopy(), m_pDC->pCopy()));

		return pCL;
	};

	virtual std::ostream&
	Restart(std::ostream& out) const {
		out << "drive caller wrap, ", m_pDC->Restart(out) << ", ";
		return m_pCL->Restart(out);
	};

	virtual void
	Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::EpsilonPrime = EpsPrime;

		m_pCL->Update(Eps, EpsPrime);

		doublereal d = m_pDC->dGet();
		ConstitutiveLaw<T, Tder>::F = m_pCL->GetF()*d;

		if (GetConstLawType() & ConstLawType::ELASTIC) {
			ConstitutiveLaw<T, Tder>::FDE = m_pCL->GetFDE()*d;
		}

		if (GetConstLawType() & ConstLawType::VISCOUS) {
			ConstitutiveLaw<T, Tder>::FDEPrime = m_pCL->GetFDEPrime()*d;
		}
	};

	virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		m_pCL->AfterConvergence(Eps, EpsPrime);
	};
};

struct DriveCL1R : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		DriveCaller *pDC = HP.GetDriveCaller();
		ConstitutiveLaw<doublereal, doublereal>* pCL2 = HP.GetConstLaw1D(CLType);
		typedef DriveCLWrapper<doublereal, doublereal> DCLW;
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;
		SAFENEWWITHCONSTRUCTOR(pCL, DCLW, DCLW(pCL2, pDC));
		return pCL;
	};
};

struct DriveCL3R : public ConstitutiveLawRead<Vec3, Mat3x3> {
	virtual ConstitutiveLaw<Vec3, Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		DriveCaller *pDC = HP.GetDriveCaller();
		ConstitutiveLaw<Vec3, Mat3x3>* pCL2 = HP.GetConstLaw3D(CLType);
		typedef DriveCLWrapper<Vec3, Mat3x3> DCLW;
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;
		SAFENEWWITHCONSTRUCTOR(pCL, DCLW, DCLW(pCL2, pDC));
		return pCL;
	};
};

struct DriveCL6R : public ConstitutiveLawRead<Vec6, Mat6x6> {
	virtual ConstitutiveLaw<Vec6, Mat6x6> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		DriveCaller *pDC = HP.GetDriveCaller();
		ConstitutiveLaw<Vec6, Mat6x6>* pCL2 = HP.GetConstLaw6D(CLType);
		typedef DriveCLWrapper<Vec6, Mat6x6> DCLW;
		ConstitutiveLaw<Vec6, Mat6x6>* pCL = 0;
		SAFENEWWITHCONSTRUCTOR(pCL, DCLW, DCLW(pCL2, pDC));
		return pCL;
	};
};

/* DriveCLWrapper - end */

#endif // CONSTLTP_DCW_H
