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

#ifndef CONSTLTP_AXW_H
#define CONSTLTP_AXW_H

#include "constltp_impl.h"

/* AxialCLWrapper - begin */

class AxialCLWrapper
: public ConstitutiveLaw<Vec3, Mat3x3> {
private:
	ConstitutiveLaw<doublereal, doublereal> *m_pCL;
	const Vec3 m_v;
	Vec3 m_F_0;
	Mat3x3 m_FDE_0;

public:
	AxialCLWrapper(ConstitutiveLaw<doublereal, doublereal> *pCL, const Vec3& v)
	: m_pCL(pCL), m_v(v), m_F_0(v), m_FDE_0(v.Tens())
	{ 
		NO_OP;
	};

	virtual
	~AxialCLWrapper(void) {
		SAFEDELETE(m_pCL);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return m_pCL->GetConstLawType();
	};

	virtual
	ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		SAFENEWWITHCONSTRUCTOR(pCL,
			AxialCLWrapper,
			AxialCLWrapper(m_pCL->pCopy(), m_v));

		return pCL;
	};

	virtual std::ostream&
	Restart(std::ostream& out) const {
		out << "axial wrap, " << m_v << ", ";
		return m_pCL->Restart(out);
	};

	virtual void
	Update(const Vec3& Eps, const Vec3& EpsPrime = mb_zero<Vec3>()) {
		ConstitutiveLaw<Vec3, Mat3x3>::Epsilon = Eps;
		ConstitutiveLaw<Vec3, Mat3x3>::EpsilonPrime = EpsPrime;

		m_pCL->Update(m_v*Eps, m_v*EpsPrime);

		ConstitutiveLaw<Vec3, Mat3x3>::F = m_F_0*m_pCL->GetF();

		if (GetConstLawType() & ConstLawType::ELASTIC) {
			ConstitutiveLaw<Vec3, Mat3x3>::FDE = m_FDE_0*m_pCL->GetFDE();
		}

		if (GetConstLawType() & ConstLawType::VISCOUS) {
			ConstitutiveLaw<Vec3, Mat3x3>::FDEPrime = m_FDE_0*m_pCL->GetFDEPrime();
		}
	};

	virtual void AfterConvergence(const Vec3& Eps, const Vec3& EpsPrime = mb_zero<Vec3>()) {
		m_pCL->AfterConvergence(m_v*Eps, m_v*EpsPrime);
	};
};

struct AxialCLR : public ConstitutiveLawRead<Vec3, Mat3x3> {
	virtual ConstitutiveLaw<Vec3, Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		Vec3 v(HP.GetVec3());
		doublereal d = v.Norm();
		v /= d;

		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;
		SAFENEWWITHCONSTRUCTOR(pCL,
			AxialCLWrapper,
			AxialCLWrapper(HP.GetConstLaw1D(CLType), v));

		return pCL;
	};
};

/* AxialCLWrapper - end */

#endif // CONSTLTP_AXW_H
