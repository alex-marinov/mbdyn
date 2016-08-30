/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
	const doublereal m_dInitialEpsPrimeTol;
	
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
		doublereal dSign, ContContactCL::Type type, doublereal dRest,
		doublereal dK, doublereal dExp, doublereal dInitialEpsPrimeTol)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pTplDC, 0.),
	m_dSign(dSign), m_dInitialEpsPrimeTol(dInitialEpsPrimeTol),
	m_dRest(dRest), m_type(type), m_dK(dK), m_dExp(dExp),
	m_bActive(false), m_bToggling(false),
	m_dInitialEpsPrime(0.)
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
				m_dSign, m_type, m_dRest, m_dK, m_dExp, m_dInitialEpsPrimeTol));

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
			<< ", EpsPrimeTol, " << m_dInitialEpsPrimeTol
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
				m_dInitialEpsPrime = std::max(xp, m_dInitialEpsPrimeTol);

				switch (m_type) {
				case CC_FLORES_ET_AL:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 8./5. * m_dK * (1. - m_dRest) / m_dRest / m_dInitialEpsPrime;
					break;

				case CC_HUNT_CROSSLEY:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 3./2. * m_dK * (1. - m_dRest) / m_dInitialEpsPrime;
					break;

				case CC_LANKARANI_NIKRAVESH:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 3./4. * m_dK * (1. - std::pow(m_dRest, 2)) / m_dInitialEpsPrime;
					break;
				}
			}

//			doublereal xn   = std::pow(x, m_dExp);
//			doublereal xnm1 = std::pow(x, m_dExp - 1.);	

//			doublereal xp_t = ( ( (m_dDissCoef*xp < 0.)&&(abs(m_dDissCoef*xp) > m_dK)) ? abs(m_dK/m_dDissCoef/xp) : 1. );	// Threshold

//			ConstitutiveLaw<doublereal, doublereal>::F = m_dSign*(m_dK*xn + m_dDissCoef *xn* xp * xp_t); 	// FIXME
//			ConstitutiveLaw<doublereal, doublereal>::FDE = m_dSign*(m_dExp*m_dK*xnm1 + m_dDissCoef *m_dExp*xnm1* xp_t);	//FIXME
//			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign*(m_dDissCoef * xn) * xp_t ;
			

			doublereal xn   = std::pow(x, m_dExp);
			doublereal xnm1 = std::pow(x, m_dExp - 1.);			
			
			ConstitutiveLaw<doublereal, doublereal>::F = m_dSign*xn*(m_dK + m_dDissCoef*xp);	// FIXME
			ConstitutiveLaw<doublereal, doublereal>::FDE = m_dSign*xnm1*(m_dExp*m_dK + m_dDissCoef*m_dExp);	
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign*xn*m_dDissCoef;
			
//			silent_cout("hello "  << std::endl);
		}
	};
	
	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return out << " " << m_dInitialEpsPrime << " " << m_dDissCoef;
	};
	

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
			silent_cerr("ContContactCLR:\n"
				"        continuous contact ,\n"
				"                [ , sign , { negative | positive | <sign> } , ]\n"
				"                [ , formulation, { flores | lankarani nikravesh | hunt crossley } , ]\n"
				"                restitution , <restitution_coefficient> , (0 -> 1)\n"
				"                kappa , <stiffness> , (> 0)\n"
				"                exp , <exponent> , (>= 1)\n"
				"                [ , EpsPrimeTol , <EpsPrimeTol> ]\n"
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
		
		doublereal dInitialEpsPrimeTol = 1.e-6;
		if (HP.IsKeyWord("EpsPrimeTol")) {
			dInitialEpsPrimeTol = HP.GetReal();
			if (dInitialEpsPrimeTol < 0.) {
				silent_cerr("ContContactCLR: invalid \"EpsPrimeTol\" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		/* Prestrain */
		TplDriveCaller<doublereal> *pTplDC = GetPreStrain<doublereal>(pDM, HP);

		SAFENEWWITHCONSTRUCTOR(pCL, ContContactCL,
			ContContactCL(pTplDC, dSign, type, dRest, dK, dExp, dInitialEpsPrimeTol));

		return pCL;
	};
};

//---------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------

class ContContact3DCL
: public ElasticConstitutiveLaw<Vec3, Mat3x3> {
public:
	enum Type {
		CC_FLORES_ET_AL = ContContactCL::CC_FLORES_ET_AL,
		CC_HUNT_CROSSLEY = ContContactCL::CC_HUNT_CROSSLEY,
		CC_LANKARANI_NIKRAVESH = ContContactCL::CC_LANKARANI_NIKRAVESH
	};

protected:
	const doublereal m_dSign;
	const doublereal m_dInitialEpsPrimeTol;
	
	const doublereal m_dRest;	// coefficient of restitution
			
	ContContact3DCL::Type m_type; 	// Formulation of continuous contact law: 
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
	ContContact3DCL(const TplDriveCaller<Vec3> *pTplDC,
		doublereal dSign, ContContact3DCL::Type type, doublereal dRest,
		doublereal dK, doublereal dExp, doublereal dInitialEpsPrimeTol)
	: ElasticConstitutiveLaw<Vec3, Mat3x3>(pTplDC, Zero3),
	m_dSign(dSign), m_dInitialEpsPrimeTol(dInitialEpsPrimeTol),
	m_dRest(dRest), m_type(type), m_dK(dK), m_dExp(dExp),
	m_bActive(false), m_bToggling(false), m_dInitialEpsPrime(0.)
	{
		NO_OP;
	};

	virtual ~ContContact3DCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<Vec3, Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, ContContact3DCL,
			ContContact3DCL(pGetDriveCaller()->pCopy(),
				m_dSign, m_type, m_dRest, m_dK, m_dExp, m_dInitialEpsPrimeTol));

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
			<< ", EpsPrimeTol, " << m_dInitialEpsPrimeTol
			<< ", ", ElasticConstitutiveLaw<Vec3, Mat3x3>::Restart_int(out);
		return out;
	};

	virtual void Update(const Vec3& Eps, const Vec3& EpsPrime) {
		ConstitutiveLaw<Vec3, Mat3x3>::Epsilon = Eps - ElasticConstitutiveLaw<Vec3, Mat3x3>::Get();
		ConstitutiveLaw<Vec3, Mat3x3>::EpsilonPrime = EpsPrime;

		doublereal x = m_dSign*ConstitutiveLaw<Vec3, Mat3x3>::Epsilon(3);
		if (x < 0.) {
			if (m_bActive) {
				if (!m_bToggling) {
					m_bToggling = true;
				}
			}

			ConstitutiveLaw<Vec3, Mat3x3>::F = Zero3;
			ConstitutiveLaw<Vec3, Mat3x3>::FDE = Zero3x3;
			ConstitutiveLaw<Vec3, Mat3x3>::FDEPrime = Zero3x3;

		} else {
			doublereal xp = m_dSign*ConstitutiveLaw<Vec3, Mat3x3>::EpsilonPrime(3);

			if (!m_bActive && !m_bToggling) {
				m_bToggling = true;
				m_dInitialEpsPrime = std::max(xp, m_dInitialEpsPrimeTol);
					
				switch (m_type) {
				case CC_FLORES_ET_AL:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 8./5. * m_dK * (1. - m_dRest) / m_dRest / m_dInitialEpsPrime;
					break;

				case CC_HUNT_CROSSLEY:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 3./2. * m_dK * (1. - m_dRest) / m_dInitialEpsPrime;
					break;

				case CC_LANKARANI_NIKRAVESH:
					m_dDissCoef = (std::abs(xp) > m_dInitialEpsPrimeTol) * 3./4. * m_dK * (1. - std::pow(m_dRest, 2)) / m_dInitialEpsPrime;
					break;
				}
			}

//			doublereal xn   = std::pow(x, m_dExp);
//			doublereal xnm1 = std::pow(x, m_dExp - 1.);	

//			doublereal xp_t = ( ( (m_dDissCoef*xp < 0.)&&(abs(m_dDissCoef*xp) > m_dK)) ? abs(m_dK/m_dDissCoef/xp) : 1. );	// Threshold

//			ConstitutiveLaw<doublereal, doublereal>::F = m_dSign*(m_dK*xn + m_dDissCoef *xn* xp * xp_t); 	// FIXME
//			ConstitutiveLaw<doublereal, doublereal>::FDE = m_dSign*(m_dExp*m_dK*xnm1 + m_dDissCoef *m_dExp*xnm1* xp_t);	//FIXME
//			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign*(m_dDissCoef * xn) * xp_t ;
			
			doublereal xn   = std::pow(x, m_dExp);
			doublereal xnm1 = std::pow(x, m_dExp - 1.);			
		
//----------------------FIXME CHECK CORRECTNESS
			ConstitutiveLaw<Vec3, Mat3x3>::F(3) = m_dSign*xn*(m_dK + m_dDissCoef*xp);	// FIXME
			ConstitutiveLaw<Vec3, Mat3x3>::FDE(3, 3) = m_dSign*xnm1*(m_dExp*m_dK + m_dDissCoef*m_dExp);	
			ConstitutiveLaw<Vec3, Mat3x3>::FDEPrime(3, 3) = m_dSign*xn*m_dDissCoef;
			

		}
	};
	
	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return out << " " << m_dInitialEpsPrime << " " << m_dDissCoef;
	};

	virtual void AfterConvergence(const Vec3& Eps, const Vec3& EpsPrime = Zero3) {
		
#if 0
		if (m_bActive) {
			doublereal eps2show = m_dSign*ConstitutiveLaw<Vec3, Mat3x3>::Epsilon(3) ;
			doublereal F2show = m_dSign*ConstitutiveLaw<Vec3, Mat3x3>::F(3);
			silent_cout( "	"	 << eps2show  << "	" 	<< F2show << std::endl);
		}	
#endif		

		if (m_bToggling) {
		
#if 0
			silent_cout(">> ContContact3DCL::AfterConvergence() "
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
			silent_cout("<< ContContact3DCL::AfterConvergence() "
				"m_bToggling=" << (m_bToggling ? "true" : "false") << " "
				"m_bActive=" << (m_bActive ? "true" : "false") << " "
				"m_dInitialEpsPrime=" << m_dInitialEpsPrime << std::endl);
#endif
		}
	};
};

/* specific functional object(s) */
struct ContContact3DCLR : public ConstitutiveLawRead<Vec3, Mat3x3> {
	virtual ConstitutiveLaw<Vec3, Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<Vec3, Mat3x3>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cerr("ContContact3DCL:\n"
				"        continuous contact ,\n"
				"                [ , sign , { negative | positive | <sign> } , ]\n"
				"                [ , formulation , { flores | lankarani nikravesh | hunt crossley } , ]\n"
				"                restitution , <restitution_coefficient>, (0 -> 1)\n"
				"                kappa , <stiffness> , (> 0)\n"
				"                exp , <exponent> , (>= 1)\n"
				"                [ , EpsPrimeTol , <EpsPrimeTol> ]\n"
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
					silent_cerr("ContContact3DCLR: invalid sign " << d
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				dSign = copysign(1., d);
			}
		}
		
		ContContact3DCL::Type type(ContContact3DCL::CC_FLORES_ET_AL);
		if (HP.IsKeyWord("formulation")) {
			if (HP.IsKeyWord("flores")) {
				type = ContContact3DCL::CC_FLORES_ET_AL;

			} else if (HP.IsKeyWord("hunt" "crossley")) {
				type = ContContact3DCL::CC_HUNT_CROSSLEY;

			} else if (HP.IsKeyWord("lankarani" "nikravesh")) {
				type = ContContact3DCL::CC_LANKARANI_NIKRAVESH;

			} else {
				silent_cerr("ContContact3DCLR: invalid formulation"
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}	

		if (!HP.IsKeyWord("restitution")) {
			silent_cerr("ContContact3DCLR: \"restitution\" expected at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dRest = HP.GetReal();
		if (dRest < 0.) {
			silent_cerr("ContContact3DCLR: invalid \"restitution\" at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if ((type == ContContact3DCL::CC_FLORES_ET_AL) && (dRest <= 0.)) {
			silent_cerr("ContContact3DCLR: null \"restitution\" incompatible with \"Flores\" at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("kappa")) {
			silent_cerr("ContContact3DCLR: \"kappa\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dK = HP.GetReal();
		if (dK <= 0.) {
			silent_cerr("ContContact3DCLR: invalid \"kappa\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("exp")) {
			silent_cerr("ContContact3DCLR: \"exp\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dExp = HP.GetReal();
		if (dExp < 1.) {
			silent_cerr("ContContact3DCLR: invalid \"exp\" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		doublereal dInitialEpsPrimeTol = 1.e-6;
		if (HP.IsKeyWord("EpsPrimeTol")) {
			dInitialEpsPrimeTol = HP.GetReal();
			if (dInitialEpsPrimeTol < 0.) {
				silent_cerr("ContContact3DCLR: invalid \"EpsPrimeTol\" at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		/* Prestrain */
		TplDriveCaller<Vec3> *pTplDC = GetPreStrain<Vec3>(pDM, HP);

		SAFENEWWITHCONSTRUCTOR(pCL, ContContact3DCL,
			ContContact3DCL(pTplDC, dSign, type, dRest, dK, dExp, dInitialEpsPrimeTol));

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
	
	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new ContContact3DCLR;
	if (!SetCL3D("continuous" "contact", rf3D)) {
		delete rf3D;

		silent_cerr("ContContact3DCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	} 

	return 0;
}



