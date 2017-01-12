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
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <limits>

#include "dataman.h"
#include "hfluid.h"
#include "hfluid_.h"

/* HydraulicFluid - begin */

HydraulicFluid::HydraulicFluid(unsigned int Label,
	const doublereal dPres0,
	const doublereal dTemp0)
: WithLabel(Label), dPres0(dPres0), dTemp0(dTemp0)
{
	NO_OP;
}

HydraulicFluid::HydraulicFluid(const HydraulicFluid& HF)
: WithLabel(HF.GetLabel()), dPres0(HF.dPres0), dTemp0(HF.dTemp0)
{
	NO_OP;
}

HydraulicFluid::~HydraulicFluid(void)
{
      NO_OP;
}

doublereal
HydraulicFluid::dGetPres0(void) const
{
	ASSERT(dPres0 != -1.);
	return dPres0;
}

doublereal
HydraulicFluid::dGetTemp0(void) const
{
	ASSERT(dTemp0 != -1.);
	return dTemp0;
}

/* HydraulicFluid - end */


HydraulicFluid *
ReadHydraulicFluid(MBDynParser& HP, unsigned int uLabel)
{
	/* tipi di fluidi */
	const char* sKeyWords[] = {
		"uncompressible",		/* deprecated (typo) */
		"incompressible",
		"linear" "compressible",
		"linear" "thermal" "compressible",
		"super",
		"exponential",

		NULL
	};

	enum KeyWords {
		UNKNOWN = -1,

		UNCOMPRESSIBLE = 0,	/* deprecated (typo) */
		INCOMPRESSIBLE,
		LINEARCOMPRESSIBLE,
		LINEARTHERMALCOMPRESSIBLE,
		SUPER,
		EXPONENTIAL,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di drive */
	KeyWords CurrKeyWord = KeyWords(HP.IsKeyWord());
	if (CurrKeyWord == UNKNOWN) {
		CurrKeyWord = INCOMPRESSIBLE;
	}

	HydraulicFluid* pHF = 0;

	switch (CurrKeyWord) {
	case UNCOMPRESSIBLE:	/* deprecated (typo) */
	case INCOMPRESSIBLE: {
		doublereal dDensity(0.);
		if (HP.IsKeyWord("density")) {
			dDensity = HP.GetReal();
		}

		if (dDensity < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("line " << HP.GetLineData()
				<< ": illegal density " << dDensity
				<< " for hydraulic fluid " << uLabel << std::endl);
	  		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal dViscosity(0.);
		if (HP.IsKeyWord("viscosity")) {
			dViscosity = HP.GetReal();
		}

		if (dViscosity < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("line " << HP.GetLineData()
				<< ": illegal viscosity " << dViscosity
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal dPres0(-1.);
		if (HP.IsKeyWord("pressure")) {
			dPres0 = HP.GetReal();
			if (dPres0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference pressure " << dPres0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}
		}

		doublereal dTemp0(-1.);
		if (HP.IsKeyWord("temperature")) {
			dTemp0 = HP.GetReal();
			if (dTemp0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference temperature " << dTemp0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}
		}

		SAFENEWWITHCONSTRUCTOR(pHF,
			IncompressibleHydraulicFluid,
			IncompressibleHydraulicFluid(uLabel,
				dDensity,
				dViscosity,
				dPres0,
				dTemp0));
		} break;

	case LINEARCOMPRESSIBLE:
	case SUPER:
	case EXPONENTIAL: {
		doublereal dDensity(0.);
		doublereal dBeta(0.);
		doublereal dPres0(0.);
		if (HP.IsKeyWord("density")) {
			dDensity = HP.GetReal();

			if (HP.IsKeyWord("sound" "celerity")) {
				doublereal sound = HP.GetReal();
				if (sound < std::numeric_limits<doublereal>::epsilon()) {
					silent_cerr("line " << HP.GetLineData()
						<< ": illegal sound celerity " << sound
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				dBeta = sound*sound*dDensity;
			} else {
				dBeta = HP.GetReal();
				if (std::abs(dBeta) < std::numeric_limits<doublereal>::epsilon()) {
					silent_cerr("line " << HP.GetLineData()
						<< ": illegal bulk modulus " << dBeta
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			if (dDensity < std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal density " << dDensity
					<< " for hydraulic fluid " << uLabel << std::endl);
			}

			dPres0 = HP.GetReal();
			if (dPres0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference pressure " << dPres0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}
		}

		doublereal dViscosity(0.);
		if (HP.IsKeyWord("viscosity")) {
			dViscosity = HP.GetReal();
		}

		if (dViscosity < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("line " << HP.GetLineData()
				<< ": illegal viscosity " << dViscosity
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal dTemp0(-1.);
		if (HP.IsKeyWord("temperature")) {
			dTemp0 = HP.GetReal();
			if (dTemp0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference temperature " << dTemp0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}
		}

		switch (CurrKeyWord) {
		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		case LINEARCOMPRESSIBLE:
			SAFENEWWITHCONSTRUCTOR(pHF,
				LinearCompressibleHydraulicFluid,
				LinearCompressibleHydraulicFluid(uLabel,
					dDensity,
					dBeta,
					dPres0,
					dViscosity,
					dTemp0));
			break;

		case SUPER:
			SAFENEWWITHCONSTRUCTOR(pHF,
				SuperHydraulicFluid,
				SuperHydraulicFluid(uLabel,
					dDensity,
					dBeta,
					dPres0,
					dViscosity,
					dTemp0));
			break;

		case EXPONENTIAL:
			doublereal dPsat = HP.GetReal();

			SAFENEWWITHCONSTRUCTOR(pHF,
				ExpHydraulicFluid,
				ExpHydraulicFluid(uLabel,
					dDensity,
					dBeta,
					dPres0,
					dPsat,
					dViscosity,
					dTemp0));
			break;
		}

		} break;

	case LINEARTHERMALCOMPRESSIBLE: {
		doublereal dDensity(0.);
		doublereal dBeta(0.);
		doublereal dPres0(0.);
		doublereal dAlpha(0.);
		doublereal dTemp0(0.);
		if (HP.IsKeyWord("density")) {
			dDensity = HP.GetReal();

			if (HP.IsKeyWord("sound" "celerity")) {
				doublereal sound = HP.GetReal();
				if (sound < std::numeric_limits<doublereal>::epsilon()) {
					silent_cerr("line " << HP.GetLineData()
						<< ": illegal sound celerity " << sound
						<< " for hydraulic fluid " << uLabel << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				dBeta = sound*sound*dDensity;

			} else {
				dBeta = HP.GetReal();
				if (std::abs(dBeta) < std::numeric_limits<doublereal>::epsilon()) {
					silent_cerr("line " << HP.GetLineData()
						<< ": illegal bulk modulus " << dBeta
						<< " for hydraulic fluid " << uLabel << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

			if (dDensity < std::numeric_limits<doublereal>::epsilon()) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal density " << dDensity
					<< " for hydraulic fluid " << uLabel << std::endl);
			}

			dPres0 = HP.GetReal();
			if (dPres0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference pressure " << dPres0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}

			dAlpha = HP.GetReal();

			dTemp0 = HP.GetReal();
			if (dTemp0 < 0.) {
				silent_cerr("line " << HP.GetLineData()
					<< ": illegal reference temperature " << dTemp0
					<< " for hydraulic fluid " << uLabel << std::endl);
			}
		}

		doublereal dViscosity(0.);
		if (HP.IsKeyWord("viscosity")) {
			dViscosity = HP.GetReal();
		}

		if (dViscosity < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("line " << HP.GetLineData()
				<< ": illegal viscosity " << dViscosity
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		SAFENEWWITHCONSTRUCTOR(pHF,
			LinearCompressibleTHydraulicFluid,
			LinearCompressibleTHydraulicFluid(uLabel,
				dDensity,
				dBeta,
				dPres0,
				dAlpha,
				dTemp0,
				dViscosity));
		} break;

	default:
		silent_cerr("line " << HP.GetLineData()
			<< ": unknown hydraulic fluid type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pHF;
}
