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
/*
 * Created by: Pierangelo Masarati <masarati@aero.polimi.it>
 * Modified by: Matteo Fancello <matteo.fancello@gmail.com>
 */

/*
 * Version modified to allow contact without impacts, introducing a threshold for the InitialEpsPrime in order not to have the 
 * dissipative term of the reaction force govern sign, with the effect of imposing adhesive forces. 
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp_impl.h"

class ContactLawCL
: public ElasticConstitutiveLaw<doublereal, doublereal> {
protected:
	const doublereal m_dSign;
	
	const doublereal m_dK;
	const doublereal m_dExp;
	const doublereal m_dC;
	const doublereal m_d_max;

	mutable bool m_bActive;			// is contact ongoing?
	mutable bool m_bToggling;		// toggle m_bActive

public:
	ContactLawCL(const TplDriveCaller<doublereal> *pTplDC,
		doublereal dSign,
		doublereal dK, doublereal dExp, doublereal dC, doublereal d_max)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pTplDC, 0.),
	m_dSign(dSign),
	m_dK(dK), m_dExp(dExp), m_dC(dC), m_d_max(d_max),
	m_bActive(false), m_bToggling(false)
	{
		NO_OP;
	};

	virtual ~ContactLawCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, ContactLawCL,
			ContactLawCL(pGetDriveCaller()->pCopy(),
				m_dSign, m_dK, m_dExp, m_dC, m_d_max));

		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "contact law"
			<< ", sign, " << m_dSign
			<< ", kappa, " << m_dK
			<< ", exp, " << m_dExp
			<< ", damping, " << m_dC
			<< ", d_max, " << m_d_max
			<< ", ", ElasticConstitutiveLaw<doublereal, doublereal>::Restart_int(out);
		return out;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal x = m_dSign*ConstitutiveLaw<doublereal, doublereal>::Epsilon;
		if (x < 0.) {
			if (m_bActive) {
				if (!m_bToggling) {
					m_bToggling = true;
				}
			}

			ConstitutiveLaw<doublereal, doublereal>::F = 0.;
			ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

		} else {
			doublereal xp = m_dSign*ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime;

			if (!m_bActive && !m_bToggling) {
				m_bToggling = true;
			}

			doublereal xn   = std::pow(x, m_dExp);
			doublereal xnm1 = std::pow(x, m_dExp - 1.);
			// doublereal xn   = copysign(std::pow(x, m_dExp), x);
			// doublereal xnm1 = copysign(std::pow(x, m_dExp - 1.), x);

			doublereal df;
			doublereal dfp;
			doublereal dx = x/m_d_max;
			if (dx < 0.) {
				df = 0.;
				dfp = 0.;

			} else if (dx > 1.) {
				df = 1.;
				dfp = 0.;

			} else {
				df = -2.*pow(dx, 3) + 3.*pow(dx, 2);
				dfp = 6.*dx*(1. - dx);
			}

			doublereal force = std::max(0., (xn*m_dK + m_dC*xp*df));
			ConstitutiveLaw<doublereal, doublereal>::F = m_dSign*force;
			if (force == 0.) {
				ConstitutiveLaw<doublereal, doublereal>::FDE = 0.;
				ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;

			} else {
				ConstitutiveLaw<doublereal, doublereal>::FDE = m_dSign*(xnm1*m_dExp*m_dK + m_dC*xp*dfp);
				ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign*m_dC*df;
			}
		}
	};

#if 0
	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return out;
	};
#endif
	

	virtual void AfterConvergence(const doublereal& Eps, const doublereal& EpsPrime = 0.) {

#if 0
		if (m_bActive) {
		doublereal eps2show = m_dSign*ConstitutiveLaw<doublereal, doublereal>::Epsilon ;
		doublereal F2show = m_dSign*ConstitutiveLaw<doublereal, doublereal>::F;
		silent_cout( "	"	 << eps2show  << "	" 	<< F2show << std::endl);
		}	
#endif		

		if (m_bToggling) {
		
#if 0
			silent_cout(">> ContactLawCL::AfterConvergence() "
				"m_bToggling=" << (m_bToggling ? "true" : "false") << " "
				"m_bActive=" << (m_bActive ? "true" : "false") << " "
				"m_dInitialEpsPrime=" << m_dInitialEpsPrime << std::endl);
#endif

			if (m_bActive) {
				m_bActive = false;

			} else {
				m_bActive = true;
			}
			m_bToggling = false;
#if 0
			silent_cout("<< ContactLawCL::AfterConvergence() "
				"m_bToggling=" << (m_bToggling ? "true" : "false") << " "
				"m_bActive=" << (m_bActive ? "true" : "false") << " "
				"m_dInitialEpsPrime=" << m_dInitialEpsPrime << std::endl);
#endif
		}
	};
};

/* specific functional object(s) */
struct ContactLawCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cerr("ContactLawCLR:\n"
				"        contact law ,\n"
				"                [ , sign , { negative | positive | <sign> } , ]\n"
				"                kappa , <stiffness> , (> 0)\n"
				"                exp , <exponent> , (>= 1)\n"
				"                damping , <damping> , (> 0)\n"
				"                d_max , <damping_saturation_displacement> , (> 0)\n"
				"                [ , prestrain , (DriveCaller) <prestrain> ]\n"
				<< std::endl);

			if (!HP.IsArg()) {
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal dSign = -1.;
		if (HP.IsKeyWord("sign")) {
			if (HP.IsKeyWord("positive")) {
				dSign = 1.;
			} else if (HP.IsKeyWord("negative")) {
				dSign = -1.;
			} else {
				doublereal d = HP.GetReal();
				if (d == 0.) {
					silent_cerr("ContactLawCLR: invalid sign " << d
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				dSign = copysign(1., d);
			}
		}
		
		if (!HP.IsKeyWord("kappa")) {
			silent_cerr("ContactLawCLR: \"kappa\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dK = HP.GetReal();
		if (dK <= 0.) {
			silent_cerr("ContactLawCLR: invalid \"kappa\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("exp")) {
			silent_cerr("ContactLawCLR: \"exp\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dExp = HP.GetReal();
		if (dExp < 1.) {
			silent_cerr("ContactLawCLR: invalid \"exp\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		if (!HP.IsKeyWord("damping")) {
			silent_cerr("ContactLawCLR: \"damping\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dC = HP.GetReal();
		if (dC <= 0.) {
			silent_cerr("ContactLawCLR: invalid \"damping\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("damping" "saturation" "displacement")) {
			silent_cerr("ContactLawCLR: \"damping saturation displacement\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal d_max = HP.GetReal();
		if (d_max <= 0.) {
			silent_cerr("ContactLawCLR: invalid \"damping saturation displacement\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* Prestrain */
		TplDriveCaller<doublereal> *pTplDC = GetPreStrain<doublereal>(pDM, HP);

		SAFENEWWITHCONSTRUCTOR(pCL, ContactLawCL,
			ContactLawCL(pTplDC, dSign, dK, dExp, dC, d_max));

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

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new ContactLawCLR;
	if (!SetCL1D("contact" "law", rf1D)) {
		delete rf1D;

		silent_cerr("ContactLawCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}
	
	return 0;
}

