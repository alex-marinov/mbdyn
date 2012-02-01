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
/*
 * Created by: Pierangelo Masarati <masarati@aero.polimi.it>
 * Modified by: Matteo Fancello <matteo.fancello@gmail.com>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp_impl.h"

class ContContactCL
: public ElasticConstitutiveLaw<doublereal, doublereal> {
public:
	enum Type {
		CC_FLORES_ET_AL,
		CC_HUNT_CROSSLEY,
		CC_LANKARANI_NIKRAVESH
	};

protected:
	const doublereal m_dSign;
	
	const doublereal m_dRest;	// coefficient of restitution
			
	ContContactCL::Type m_type; 	// Formulation of continuous contact law: 
					// flores [default]
					// huntcrossley
					// lankaraninikravesh

	const doublereal m_dK;
	const doublereal m_dExp;

	mutable bool m_bActive;			// is contact ongoing?
	mutable bool m_bToggling;		// toggle m_bActive
	mutable doublereal m_dInitialEpsPrime;	// initial contact velocity
	mutable doublereal m_dDissCoef;		// actual dissipation coefficient

public:
	ContContactCL(const TplDriveCaller<doublereal> *pTplDC,
		doublereal dSign, ContContactCL::Type type, doublereal dRest, doublereal dK, doublereal dExp)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pTplDC, 0.),
	m_dSign(dSign), m_dRest(dRest), m_type(type), m_dK(dK), m_dExp(dExp),
	m_bActive(false), m_bToggling(false), m_dInitialEpsPrime(0.)
	{
		NO_OP;
	};

	virtual ~ContContactCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, ContContactCL,
			ContContactCL(pGetDriveCaller()->pCopy(),
				m_dSign, m_type, m_dRest, m_dK, m_dExp));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "continuous contact"
			<< ", sign, " << m_dSign
			<< ", formulation, ";

		switch (m_type) {
		case CC_FLORES_ET_AL:
			out << "flores";
			break;

		case CC_HUNT_CROSSLEY:
			out << "hunt crossley";
			break;

		case CC_LANKARANI_NIKRAVESH:
			out << "lankarani nikravesh";
			break;
		}

		out
			<< ", restitution, " << m_dRest
			<< ", kappa, " << m_dK
			<< ", exp, " << m_dExp
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

			if (!m_bActive) {
				if (!m_bToggling) {
					m_bToggling = true;
					m_dInitialEpsPrime = (xp > 1e-6 ? xp : 1e-6); // FIXME: threshold

					switch (m_type) {
					case CC_FLORES_ET_AL:
						m_dDissCoef = 8./5. * m_dK * (1. - m_dRest) / m_dRest / m_dInitialEpsPrime;
						break;

					case CC_HUNT_CROSSLEY:
						m_dDissCoef = 3./2. * m_dK * (1. - m_dRest) / m_dInitialEpsPrime;
						break;

					case CC_LANKARANI_NIKRAVESH:
						m_dDissCoef =  m_dK * 3./4. * (1. - std::pow(m_dRest, 2)) / m_dInitialEpsPrime;
						break;
					}
				}
			}

			doublereal xn = std::pow(x, m_dExp);
			doublereal xnm1 = std::pow(x, m_dExp - 1.);

			ConstitutiveLaw<doublereal, doublereal>::F = m_dSign*(m_dK*xn + m_dDissCoef *xn*xp);
			ConstitutiveLaw<doublereal, doublereal>::FDE = m_dSign*(m_dExp*m_dK*xnm1 + m_dDissCoef *m_dExp * xnm1*xp);
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign*(m_dDissCoef * xn);

		}
	};

	virtual void AfterConvergence(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		if (m_bToggling) {
#if 0
			silent_cout(">> ContContactCL::AfterConvergence() "
				"m_bToggling=" << (m_bToggling ? "true" : "false") << " "
				"m_bActive=" << (m_bActive ? "true" : "false") << " "
				"m_dInitialEpsPrime=" << m_dInitialEpsPrime << std::endl);
#endif
			if (m_bActive) {
				m_bActive = false;
				m_dInitialEpsPrime = 0.;

			} else {
				m_bActive = true;
			}
			m_bToggling = false;
#if 0
			silent_cout("<< ContContactCL::AfterConvergence() "
				"m_bToggling=" << (m_bToggling ? "true" : "false") << " "
				"m_bActive=" << (m_bActive ? "true" : "false") << " "
				"m_dInitialEpsPrime=" << m_dInitialEpsPrime << std::endl);
#endif
		}
	};
};

/* specific functional object(s) */
struct ContContactCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cerr("ContContactCL:\n"
				"        cont contact,\n"
				"                [ , sign, { negative | positive | <sign> } , ]\n"
				"                [ , formulation, { flores | lankarani nikravesh | hunt crossley } , ]\n"
				"                restitution, <restitution_coefficient>, (0 -> 1)\n"
				"                kappa, <stiffness>, (> 0)\n"
				"                exp, <exponent>, # (>= 1)\n"
				"                [ , prestrain, (DriveCaller)<prestrain> ]\n"
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
					silent_cerr("ContContactCLR: invalid sign " << d
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				dSign = copysign(1., d);
			}
		}
		
		ContContactCL::Type type(ContContactCL::CC_FLORES_ET_AL);
		if (HP.IsKeyWord("formulation")) {
			if (HP.IsKeyWord("flores")) {
				type = ContContactCL::CC_FLORES_ET_AL;

			} else if (HP.IsKeyWord("hunt" "crossley")) {
				type = ContContactCL::CC_HUNT_CROSSLEY;

			} else if (HP.IsKeyWord("lankarani" "nikravesh")) {
				type = ContContactCL::CC_LANKARANI_NIKRAVESH;

			} else {
				silent_cerr("ContContactCLR: invalid formulation"
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}	

		if (!HP.IsKeyWord("restitution")) {
			silent_cerr("ContContactCLR: \"restitution\" expected at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dRest = HP.GetReal();
		if (dRest < 0.) {
			silent_cerr("ContContactCLR: invalid \"restitution\" at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if ((type == ContContactCL::CC_FLORES_ET_AL) && (dRest <= 0.)) {
			silent_cerr("ContContactCLR: null \"restitution\" incompatible with \"Flores\" at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("kappa")) {
			silent_cerr("ContContactCLR: \"kappa\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dK = HP.GetReal();
		if (dK <= 0.) {
			silent_cerr("ContContactCLR: invalid \"kappa\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("exp")) {
			silent_cerr("ContContactCLR: \"exp\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dExp = HP.GetReal();
		if (dExp < 1.) {
			silent_cerr("ContContactCLR: invalid \"exp\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* Prestrain */
		TplDriveCaller<doublereal> *pTplDC = GetPreStrain<doublereal>(pDM, HP);

		SAFENEWWITHCONSTRUCTOR(pCL, ContContactCL,
			ContContactCL(pTplDC, dSign, type, dRest, dK, dExp));

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

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new ContContactCLR;
	if (!SetCL1D("continuous" "contact", rf1D)) {
		delete rf1D;

		silent_cerr("ContContactCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

