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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"

class HydraulicDamperCL
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	doublereal m_sigma;
	doublereal m_dotz_L;
	doublereal m_d;

public:
	HydraulicDamperCL(doublereal sigma, doublereal dotz_L, doublereal d)
	: m_sigma(sigma), m_dotz_L(dotz_L), m_d(d) {
		ConstitutiveLaw<doublereal, doublereal>::FDE = 0;
	};

	virtual ~HydraulicDamperCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		typedef HydraulicDamperCL cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_sigma, m_dotz_L, m_d));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "hydraulic damper, " << m_sigma << ", " << m_dotz_L << ", " << m_d;
	};

	virtual void Update(const doublereal& /* Eps */ , const doublereal& dotz) {
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = dotz;

		doublereal abs_dotz = std::abs(dotz);
		if (abs_dotz < m_dotz_L) {
			ConstitutiveLaw<doublereal, doublereal>::F = m_sigma*abs_dotz*dotz;
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 2*m_sigma*abs_dotz;

		} else {
			ConstitutiveLaw<doublereal, doublereal>::F =
				copysign(m_d*(abs_dotz - m_dotz_L) + m_sigma*m_dotz_L*m_dotz_L, dotz);
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_d;
		}
	};
};

/* specific functional object(s) */
struct HydraulicDamperCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		if (HP.IsKeyWord("help")) {
			silent_cerr("HydraulicDamperCL:\n"
				"        hydraulic damper , <sigma> , <dotz_L> , <d>\n"
				"# according to the formula\n"
				"#     f = sigma*|dotz|*dotz                                  for |dotz| <= dotz_L\n"
				"#     f = sign(dotz)*(d*(|dotz| - dotz_L) + sigma*dotz_L^2)  for |dotz| > dotz_L\n"
				<< std::endl);

			if (!HP.IsArg()) {
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}
		}

		CLType = ConstLawType::VISCOUS;

		doublereal sigma = HP.GetReal();
		if (sigma <= 0.) {
			silent_cerr("warning, null or negative parabolic slope "
				"at line " << HP.GetLineData() << std::endl);
		}

		doublereal dotz_L = HP.GetReal();
		if (dotz_L <= 0.) {
			silent_cerr("warning, null or negative limit velocity "
				"at line " << HP.GetLineData() << std::endl);
		}

		doublereal d = HP.GetReal();
		if (d < 0.) {
			silent_cerr("warning, negative slope "
				"at line " << HP.GetLineData() << std::endl);
		}

		typedef HydraulicDamperCL L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(sigma, dotz_L, d));

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
		= new HydraulicDamperCLR;
	if (!SetCL1D("hydraulic" "damper", rf1D)) {
		delete rf1D;

		silent_cerr("HydraulicDamperCL1D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

