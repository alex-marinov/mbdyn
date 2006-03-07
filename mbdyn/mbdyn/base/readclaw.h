/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifndef READCLAW_H
#define READCLAW_H

#include <drive_.h>		/* per TimeDriveCaller */
#include <tpldrive.h>
#include <tpldrive_.h>
#include <constltp_impl.h>
#include <symcltp.h>
#include <ginaccltp.h>
#ifdef USE_GRAALLDAMPER
#include <damper.h>
#endif /* USE_GRAALLDAMPER */
#include <shockabsorber.h>


template <class T>
void GetPreStress(MBDynParser& HP, T& PreStress)
{
	if (HP.IsKeyWord("prestress")) {
		PreStress = HP.Get(PreStress);
	}
}

template <class T>
TplDriveCaller<T>* GetPreStrain(DataManager* pDM,
		MBDynParser& HP,
		T& PreStrain)
{
	if (HP.IsKeyWord("prestrain")) {
		return ReadTplDrive(pDM, HP, PreStrain);
	} else {
		DriveCaller* pDC = NULL;
		SAFENEW(pDC, NullDriveCaller);
		T t(0.);
		TplDriveCaller<T>* pTplDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pTplDC,
   				SingleTplDriveCaller<T>,
   				SingleTplDriveCaller<T>(pDC, t));
		return pTplDC;
	}
}

template <class T, class Tder>
ConstitutiveLaw<T, Tder>* ReadConstLaw(DataManager* pDM,
		MBDynParser& HP,
		ConstLawType::Type& CLType,
		ConstitutiveLaw<T, Tder>*)
{
	DEBUGCOUT("Entering ReadConstLaw" << std::endl);

	const char* sKeyWords[] = {
		"linear" "elastic",
		"linear" "elastic" "isotropic",
		"linear" "elastic" "generic",
		"linear" "elastic" "generic" "axial" "torsion" "coupling",
		"linear" "elastic" "bistop",
		"cubic" "elastic" "generic",
		"log" "elastic",
		"double" "linear" "elastic",
		"isotropic" "hardening" "elastic",
		"contact" "elastic",
		"symbolic" "elastic" "isotropic",

		"linear" "viscous",
		"linear" "viscous" "isotropic",
		"linear" "viscous" "generic",
		"symbolic" "viscous" "isotropic",

		"linear" "viscoelastic",
		"linear" "viscoelastic" "isotropic",
		"linear" "viscoelastic" "generic",
		"linear" "viscoelastic" "generic" "axial" "torsion" "coupling",
		"cubic" "viscoelastic" "generic",
		"doublelinear" "viscoelastic",
		"turbulent" "viscoelastic",
		"linear" "viscoelastic" "bistop",
		"graall" "damper",
		"shock" "absorber",
		"symbolic" "viscoelastic" "isotropic",

		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		LINEARELASTIC = 0,
		LINEARELASTICISOTROPIC,
		LINEARELASTICGENERIC,
		LINEARELASTICGENERICAXIALTORSIONCOUPLING,
		LINEARELASTICBISTOP,
		CUBICELASTICGENERIC,
		LOGELASTIC,
		DOUBLELINEARELASTIC,
		ISOTROPICHARDENINGELASTIC,
		CONTACTELASTIC,
		SYMBOLICELASTICISOTROPIC,

		LINEARVISCOUS,
		LINEARVISCOUSISOTROPIC,
		LINEARVISCOUSGENERIC,
		SYMBOLICVISCOUSISOTROPIC,

		LINEARVISCOELASTIC,
		LINEARVISCOELASTICISOTROPIC,
		LINEARVISCOELASTICGENERIC,
		LINEARVISCOELASTICGENERICAXIALTORSIONCOUPLING,
		CUBICVISCOELASTICGENERIC,
		DOUBLELINEARVISCOELASTIC,
		TURBULENTVISCOELASTIC,
		LINEARVISCOELASTICBISTOP,
		GRAALLDAMPER,
		SHOCKABSORBER,
		SYMBOLICVISCOELASTICISOTROPIC,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	ConstitutiveLaw<T, Tder>* pCL = NULL;
	int CurrKW = HP.GetWord();
	switch (CurrKW) {
	case LINEARELASTIC:
	case LINEARELASTICISOTROPIC: {
		CLType = ConstLawType::ELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("Linear Elastic Isotropic Constitutive Law, stiffness = "
				<< dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearElasticIsotropicConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS));

		break;
	}

	case LINEARELASTICGENERIC: {
		CLType = ConstLawType::ELASTIC;

		DEBUGCOUT("Linear Elastic Generic Constitutive Law" << std::endl);
		Tder S(0.);
		S = HP.Get(S);

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearElasticGenericConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S));

		break;
	}

	case LINEARELASTICGENERICAXIALTORSIONCOUPLING: {
		CLType = ConstLawType::ELASTIC;

		DEBUGCOUT("Linear Elastic Generic Constitutive Law with Axial-Torsion Coupling" << std::endl);
		Tder S(0.);
		S = HP.Get(S);

		/* coefficiente di accoppiamento */
		doublereal dCoupl = HP.GetReal();
		DEBUGCOUT("coupling coefficient: " << dCoupl << std::endl);

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S, dCoupl));

		break;
	}

	case LOGELASTIC: {
		CLType = ConstLawType::ELASTIC;

		DEBUGCOUT("Logaritmic Elastic Constitutive Law" << std::endl);
		doublereal dS = HP.GetReal();
		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LogConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS));

		break;
	}

	case DOUBLELINEARELASTIC: {
		CLType = ConstLawType::ELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("stiffness = " << dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dUpp = HP.GetReal();
		if (dUpp <= 0.) {
			silent_cerr("warning, null or negative upper limit strain at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dLow = HP.GetReal();
		if (dLow >= 0.) {
			silent_cerr("warning, null or positive lower limit strain at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSecondS = HP.GetReal();
		if (dSecondS <= 0.) {
			silent_cerr("warning, null or negative second stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef DoubleLinearElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL,
				L,
				L(pTplDC, PreStress, dS, dUpp, dLow, dSecondS));

		break;
	}

	case ISOTROPICHARDENINGELASTIC: {
		CLType = ConstLawType::ELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("Stiffness = " << dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dE = HP.GetReal();
		DEBUGCOUT("Reference strain = " << dE << std::endl);

		if (dE <= 0.) {
			silent_cerr("error, null or negative reference strain at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef IsotropicHardeningConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dE));

		break;
	}

	case CONTACTELASTIC: {
		CLType = ConstLawType::ELASTIC;

		doublereal dK = HP.GetReal();
		DEBUGCOUT("Stiffness = " << dK << std::endl);

		if (dK <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dGamma = HP.GetReal();
		DEBUGCOUT("Exponent = " << dGamma << std::endl);

		if (dGamma < 1.) {
			silent_cerr("error, exponent < 1. at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef ContactConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dK, dGamma));

		break;
	}

	case SYMBOLICELASTICISOTROPIC:
	case SYMBOLICVISCOUSISOTROPIC:
	case SYMBOLICVISCOELASTICISOTROPIC: {
		CLType = ConstLawType::ELASTIC;

		const char *epsilon = 0;
		const char *epsilonPrime = 0;
		if (HP.IsKeyWord("epsilon")) {
			const char *tmp = HP.GetStringWithDelims();

			if (tmp == 0) {
				silent_cerr("unable to get \"epsilon\" symbol at line " << HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric();
			}
			SAFESTRDUP(epsilon, tmp);
		}

		if (HP.IsKeyWord("epsilon" "prime")) {
			const char *tmp = HP.GetStringWithDelims();

			if (tmp == 0) {
				silent_cerr("unable to get \"epsilonPrime\" symbol at line " << HP.GetLineData() << std::endl);
				throw DataManager::ErrGeneric();
			}
			SAFESTRDUP(epsilonPrime, tmp);

			if (epsilon == 0) {
				CLType = ConstLawType::VISCOUS;
			} else {
				CLType = ConstLawType::VISCOELASTIC;
			}
		}

		/* default - deprecated ... */
		if (epsilon == 0 && epsilonPrime == 0) {
			silent_cerr("need \"epsilon\" or \"epsilon prime\" "
					"at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric();
		}

		if (!HP.IsKeyWord("expression")) {
			silent_cerr("keyword \"expression\" expected at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}

		const char *tmp = HP.GetStringWithDelims();
		if (tmp == 0) {
			silent_cerr("unable to get \"expression\" at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric();
		}
		const char *expression = 0;
		SAFESTRDUP(expression, tmp);

		const char **symbols = 0;
		if (HP.IsKeyWord("symbols")) {
			/* FIXME! */
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		switch (CLType) {
		case ConstLawType::ELASTIC: {
			typedef SymbolicElasticIsotropicConstitutiveLaw<T, Tder> L;
			SAFENEWWITHCONSTRUCTOR(pCL, L,
					L(pTplDC, PreStress,
						epsilon,
						expression, symbols));
			break;
		}

		case ConstLawType::VISCOUS: {
			typedef SymbolicViscousIsotropicConstitutiveLaw<T, Tder> L;
			SAFENEWWITHCONSTRUCTOR(pCL, L,
					L(pTplDC, PreStress,
						epsilonPrime,
						expression, symbols));
			break;
		}

		case ConstLawType::VISCOELASTIC: {
			typedef SymbolicViscoElasticIsotropicConstitutiveLaw<T, Tder> L;
			SAFENEWWITHCONSTRUCTOR(pCL, L,
					L(pTplDC, PreStress,
						epsilon, epsilonPrime,
						expression, symbols));
			break;
		}

		default:
			ASSERT(0);
			throw ErrGeneric();
		}

		break;
	}

	case LINEARVISCOUS:
	case LINEARVISCOUSISOTROPIC: {
		CLType = ConstLawType::VISCOUS;

		doublereal dSP = HP.GetReal();
		DEBUGCOUT("stiffness prime = " << dSP << std::endl);

		if (dSP <= 0.) {
			silent_cerr("warning, null or negative stiffness prime at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress (no prestrain) */
		T PreStress(0.);
		GetPreStress(HP, PreStress);

		typedef LinearViscousIsotropicConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(NULL, PreStress, dSP));

		break;
	}

	case LINEARVISCOUSGENERIC: {
		CLType = ConstLawType::VISCOUS;

		Tder SP(0.);
		SP = HP.Get(SP);

		/* Prestress (no prestrain) */
		T PreStress(0.);
		GetPreStress(HP, PreStress);

		typedef LinearViscousGenericConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(NULL, PreStress, SP));

		break;
	}

	case LINEARVISCOELASTIC:
	case LINEARVISCOELASTICISOTROPIC: {
		CLType = ConstLawType::VISCOELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("Stiffness = " << dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSP = 0.;
		if (HP.IsKeyWord("proportional")) {
			doublereal k = HP.GetReal();
			dSP = k*dS;
		} else {
			dSP = HP.GetReal();
		}
		DEBUGCOUT("stiffness prime = " << dSP << std::endl);

		if (dSP <= 0.) {
			silent_cerr("warning, null or negative stiffness prime at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearViscoElasticIsotropicConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dSP));

		break;
	}

	case LINEARVISCOELASTICGENERIC: {
		CLType = ConstLawType::VISCOELASTIC;

		Tder S(0.);
		S = HP.Get(S);

		Tder SP(0.);
		if (HP.IsKeyWord("proportional")) {
			doublereal k = HP.GetReal();
			SP = S*k;
		} else {
			SP = HP.Get(SP);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearViscoElasticGenericConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S, SP));

		break;
	}

	case LINEARVISCOELASTICGENERICAXIALTORSIONCOUPLING: {
		CLType = ConstLawType::VISCOELASTIC;

		Tder S(0.);
		S = HP.Get(S);

		Tder SP(0.);
		if (HP.IsKeyWord("proportional")) {
			doublereal k = HP.GetReal();
			SP = S*k;
		} else {
			SP = HP.Get(SP);
		}

		/* coefficiente di accoppiamento */
		doublereal dCoupl = HP.GetReal();
		DEBUGCOUT("coupling coefficient: " << dCoupl << std::endl);

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef LinearViscoElasticGenericAxialTorsionCouplingConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S, SP, dCoupl));

		break;
	}

	case CUBICELASTICGENERIC: {
		CLType = ConstLawType::ELASTIC;

		T S1(0.);
		S1 = HP.Get(S1);

		T S2(0.);
		S2 = HP.Get(S2);

		T S3(0.);
		S3 = HP.Get(S3);

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef CubicElasticGenericConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S1, S2, S3));

		break;
	}

	case CUBICVISCOELASTICGENERIC: {
		CLType = ConstLawType::VISCOELASTIC;

		T S1(0.);
		S1 = HP.Get(S1);

		T S2(0.);
		S2 = HP.Get(S2);

		T S3(0.);
		S3 = HP.Get(S3);

		Tder SP(0.);
		SP = HP.Get(SP);

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef CubicViscoElasticGenericConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S1, S2, S3, SP));

		break;
	}

	case DOUBLELINEARVISCOELASTIC: {
		CLType = ConstLawType::VISCOELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("stiffness = " << dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dUpp = HP.GetReal();
		if (dUpp <= 0.) {
			silent_cerr("warning, null or negative upper limit strain at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dLow = HP.GetReal();
		if (dLow >= 0.) {
			silent_cerr("warning, null or positive lower limit strain at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSecondS = HP.GetReal();
		if (dSecondS <= 0.) {
			silent_cerr("warning, null or negative second stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSP = HP.GetReal();
		DEBUGCOUT("stiffness prime = " << dSP << std::endl);

		if (dSP <= 0.) {
			silent_cerr("warning, null or negative stiffness prime at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSecondSP = HP.GetReal();
		DEBUGCOUT("second stiffness prime = " << dSecondSP << std::endl);

		if (dSecondSP <= 0.) {
			silent_cerr("warning, null or negative second stiffness prime at line "
				<< HP.GetLineData() << std::endl);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef DoubleLinearViscoElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL,
				L,
				L(pTplDC, PreStress,
					dS, dUpp, dLow, dSecondS, dSP, dSecondSP));

		break;
	}

	case TURBULENTVISCOELASTIC:	{
		CLType = ConstLawType::VISCOELASTIC;

		doublereal dS = HP.GetReal();
		DEBUGCOUT("Visco-Elastic Turbulent Rod Joint, stiffness = "
				<< dS << std::endl);

		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dParabStiff = HP.GetReal();
		DEBUGCOUT("stiffness prime = " << dParabStiff << std::endl);

		if (dParabStiff <= 0.) {
			silent_cerr("warning, null or negative derivative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dTreshold = 0.;
		if (HP.IsArg()) {
			dTreshold = HP.GetReal(dTreshold);

			/*
			 * Il legame costitutivo ha la forma seguente:
			 *	F = Kp*e + Kd*(de/dt)
			 * con Kp costante e Kd dato dalla seguente legge:
			 *	Kd = cost2                per fabs(de/dt) < Treshold
			 *	Kd = 2*cost1*fabs(de/dt)  per fabs(de/dt) > Treshold
			 * se non viene inserito il valore di treshold, lo si
			 * assume = 0. e quindi il legame e' sempre del secondo
			 * tipo. Altrimenti, se non viene inserita la seconda
			 * costante cost2, si assume che vi sia raccordo tra
			 * i due tipi di legge, ovvero cost2 = cost1*Treshold
			 * altrimenti e' possibile avere un comportamento,
			 * che in prima approssimazione e' valido
			 * per numerosi fluidi, in cui vi e' un salto tra i due
			 * tipi di legge costitutiva. */
		}

		doublereal dSP = dTreshold*dParabStiff;
		if (HP.IsArg()) {
			dSP = HP.GetReal(dSP);
		}

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef TurbulentViscoElasticConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL,
				L,
				L(pTplDC, PreStress,
					dS, dSP, dTreshold, dParabStiff));

		break;
	}

	case LINEARELASTICBISTOP:
	case LINEARVISCOELASTICBISTOP: {
		typedef LinearViscoElasticBiStopConstitutiveLaw<T, Tder> L;
		CLType = ConstLawType::VISCOELASTIC;

		DEBUGCOUT("Linear Viscoelastic Bi Stop Constitutive Law" << std::endl);
		doublereal dS = HP.GetReal();
		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness at line "
				<< HP.GetLineData() << std::endl);
		}

		doublereal dSp = 0.;
		if (CurrKW == LINEARVISCOELASTICBISTOP) {
			dSp = HP.GetReal();
			if (dSp <= 0.) {
				silent_cerr("warning, null or negative stiffness prime at line "
					<< HP.GetLineData() << std::endl);
			}
		}

		typedef typename LinearViscoElasticBiStopConstitutiveLaw<T, Tder>::Status LS;
		LS s = L::INACTIVE;
		if (HP.IsKeyWord("initial" "status")) {
			if (HP.IsKeyWord("active")) {
				s = L::ACTIVE;
			} else if (HP.IsKeyWord("inactive")) {
				s = L::INACTIVE;
			} else {
				silent_cerr("unknown initial status at line "
					<< HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}
		}

		const DriveCaller *pA = HP.GetDriveCaller();
		const DriveCaller *pD = HP.GetDriveCaller();

		/* Prestress and prestrain */
		T PreStress(0.);
		GetPreStress(HP, PreStress);
		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dSp, s, pA, pD));

		break;
	}

	case GRAALLDAMPER: {
#ifdef USE_GRAALLDAMPER
		CLType = ConstLawType::VISCOELASTIC;

		const char* filename = HP.GetFileName();
		DEBUGCOUT("Graall damper input file: \""
				<< filename << "\"" << std::endl);

		doublereal rla = HP.GetReal();
		DEBUGCOUT("Reference length: " << rla << std::endl);

		DriveCaller* pDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pDC,
				TimeDriveCaller,
				TimeDriveCaller(pDM->pGetDrvHdl()));

		T t(1.);
		TplDriveCaller<T>* pTplDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pTplDC,
				SingleTplDriveCaller<T>,
				SingleTplDriveCaller<T>(pDC, t));

		typedef GRAALLDamperConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL,
				L,
				L(pTplDC, rla, filename));

		break;
#else /* USE_GRAALLDAMPER */
		silent_cerr("can't use GRAALL Damper" << std::endl);
		throw ErrGeneric();
#endif /* USE_GRAALLDAMPER */
	}

	/*
	 * Shock absorber per Stefy:
	 *
	 * ``Riprogettazione dell'ammortizzatore del carrello anteriore
	 * di un velivolo di aviazione generale'',
	 * S. Carlucci e S. Gualdi,
	 * A.A. 1997-98
	 */
	case SHOCKABSORBER: {
		CLType = ConstLawType::VISCOELASTIC;

		T PreStrain(0.);
		TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, PreStrain);

		typedef ShockAbsorberConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(pDM, pTplDC, HP));

		break;
	}

    /* Aggiungere altri constitutive laws */

	default:
		silent_cerr("Unknown constitutive law type at line "
			<< HP.GetLineData() << std::endl);

		throw ErrGeneric();
	}

	ASSERT(pCL != NULL);
	return pCL;
} /* End of ReadConstLaw */

#endif /* READCLAW_H */
