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

/* gestore dell'output */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <sstream>
#include <list>

//#include "dataman.h"
#include "mbpar.h"
#include "units.h"
#include "dataman.h"

/* MBUnits - begin */

const std::unordered_map<const MBUnits::Dimensions, const std::string> DimensionNames ({
	{ MBUnits::Dimensions::Dimensionless , std::string("Dimensionless") },
	{ MBUnits::Dimensions::Boolean , std::string("Boolean") },
	{ MBUnits::Dimensions::Length , std::string("Length") },
	{ MBUnits::Dimensions::Mass , std::string("Mass") },
	{ MBUnits::Dimensions::Time , std::string("Time") },
	{ MBUnits::Dimensions::Current , std::string("Current") },
	{ MBUnits::Dimensions::Temperature , std::string("Temperature") },
	{ MBUnits::Dimensions::Angle , std::string("Angle") },
	{ MBUnits::Dimensions::Area , std::string("Area") },
	{ MBUnits::Dimensions::Force , std::string("Force") },
	{ MBUnits::Dimensions::Velocity , std::string("Velocity") },
	{ MBUnits::Dimensions::Acceleration , std::string("Acceleration") },
	{ MBUnits::Dimensions::AngularVelocity , std::string("Angular velocity") },
	{ MBUnits::Dimensions::AngularAcceleration , std::string("Angular acceleration") },

	{ MBUnits::Dimensions::Momentum , std::string("Momentum") },
	{ MBUnits::Dimensions::MomentaMoment , std::string("Momenta moment") },
	{ MBUnits::Dimensions::MomentumDerivative , std::string("Momentum derivative") },
	{ MBUnits::Dimensions::MomentaMomentDerivative , std::string("Momenta moment derivative") },

	{ MBUnits::Dimensions::LinearStrain , std::string("Linear strain") },
	{ MBUnits::Dimensions::AngularStrain , std::string("Angular strain") },
	{ MBUnits::Dimensions::LinearStrainRate , std::string("Linear strain rate") },
	{ MBUnits::Dimensions::AngularStrainRate , std::string("Angular strain rate") },

	{ MBUnits::Dimensions::StaticMoment , std::string("Static moment") },
	{ MBUnits::Dimensions::MomentOfInertia , std::string("Moment of inertia") },

	{ MBUnits::Dimensions::ForceUnitSpan , std::string("Force per unit span") },

	{ MBUnits::Dimensions::Work , std::string("Work") },
	{ MBUnits::Dimensions::Power , std::string("Power") },
	{ MBUnits::Dimensions::Pressure , std::string("Pressure") },
	{ MBUnits::Dimensions::Moment , std::string("Moment") },
	{ MBUnits::Dimensions::Voltage , std::string("Voltage") },
	{ MBUnits::Dimensions::Charge , std::string("Charge") },
	{ MBUnits::Dimensions::Frequency , std::string("Frequency") },
	{ MBUnits::Dimensions::deg , std::string("deg") },
	{ MBUnits::Dimensions::rad , std::string("rad") },

	/* added later for GetEquationDimension method of DofOwnerOwner class */

	{ MBUnits::Dimensions::MassFlow, std::string("Mass flow")},
	{ MBUnits::Dimensions::Jerk , std::string("Jerk") },
	{ MBUnits::Dimensions::VoltageDerivative , std::string("Voltage derivative") },
	{ MBUnits::Dimensions::UnknownDimension , std::string("Unknown dimension") }
});

/* Costruttore senza inizializzazione */

void MBUnits::ReadOutputUnits(std::ostream& Log, MBDynParser& HP) {
	if (HP.IsKeyWord("MKS")) {
		SetMKSUnits(Units);
	} else if (HP.IsKeyWord("CGS")) {
		SetCGSUnits(Units);
	} else if (HP.IsKeyWord("MMTMS")) {
		SetMMTMSUnits(Units);
	} else if (HP.IsKeyWord("MMKGMS")) {
		SetMMKGMSUnits(Units);
	} else if (HP.IsKeyWord("Custom")) {
		const std::list<Dimensions> BaseUnits ({
			Dimensions::Length,
			Dimensions::Mass,
			Dimensions::Time,
			Dimensions::Current,
			Dimensions::Temperature		 
		});
		for (auto i = BaseUnits.begin(); i != BaseUnits.end(); i++) {
			if (HP.IsKeyWord(DimensionNames.find(*i)->second.c_str())) {
				Units[*i] = HP.GetStringWithDelims();
			} else {
				silent_cerr("Error while reading Custom unit system  at line"
						<< HP.GetLineData()
						<< "\nExpecting the definition of "
						<< DimensionNames.find(*i)->second
						<< " units."
						<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		SetDerivedUnits(Units);
		for (auto i = DimensionNames.begin(); i != DimensionNames.end(); i++) {
			Log << "Unit for " << i->second << ": " << Units[i->first] << std::endl;
		}
	} else {
		silent_cerr("Error while reading the model Units at line"
						<< HP.GetLineData()
						<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void MBUnits::SetDerivedUnits(std::unordered_map<Dimensions, std::string>& Units ) {
	Units[Dimensions::Angle] = "rad";
	Units[Dimensions::Area] = Units[Dimensions::Length] + "^2";
	Units[Dimensions::Force] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Length] + " " +
		Units[Dimensions::Time] + "^-2";
	Units[Dimensions::Velocity] = Units[Dimensions::Length] + " " +
		Units[Dimensions::Time] + "^-1";
	Units[Dimensions::Acceleration] = Units[Dimensions::Length] + " " +
		Units[Dimensions::Time] + "^-2";
	Units[Dimensions::AngularVelocity] = Units[Dimensions::Angle] + " " +
		Units[Dimensions::Time] + "^-1";
	Units[Dimensions::AngularAcceleration] = Units[Dimensions::Angle] + " " +
		Units[Dimensions::Time] + "^-2";

	Units[Dimensions::Momentum] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Velocity];
	Units[Dimensions::MomentaMoment] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Length] + "^2 " + 
		Units[Dimensions::Time] + "^-1";
	Units[Dimensions::MomentumDerivative] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Acceleration];
	Units[Dimensions::MomentaMomentDerivative] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Length] + "^2 " + 
		Units[Dimensions::Time] + "^-2";

	Units[Dimensions::LinearStrain] = Units[Dimensions::Dimensionless];
	Units[Dimensions::AngularStrain] = Units[Dimensions::Angle] + " " +
		Units[Dimensions::Length] + "^-1";
	Units[Dimensions::LinearStrainRate] = Units[Dimensions::Time] + "^-1";
	Units[Dimensions::AngularStrainRate] = Units[Dimensions::Angle] + " " +
		Units[Dimensions::Length] + "^-1 " + 
		Units[Dimensions::Time] + "^-1";

	Units[Dimensions::StaticMoment] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Length];
	Units[Dimensions::MomentOfInertia] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Length] + "^2";

	Units[Dimensions::ForceUnitSpan] = Units[Dimensions::Mass] + " " +
		Units[Dimensions::Time] + "^-2";

	Units[Dimensions::Work] = Units[Dimensions::Force] + " " +
		Units[Dimensions::Length];
	Units[Dimensions::Power] = Units[Dimensions::Force] + " " +
		Units[Dimensions::Velocity];
	Units[Dimensions::Pressure] = Units[Dimensions::Force] + " " +
		Units[Dimensions::Length] + "^-2";
	Units[Dimensions::Moment] = Units[Dimensions::Force] + " " +
		Units[Dimensions::Length];
	Units[Dimensions::Voltage] = Units[Dimensions::Length] + "^2 " +
		Units[Dimensions::Mass] + " " +
		Units[Dimensions::Time] + "^-3 " +
		Units[Dimensions::Current] + "^-1";
	Units[Dimensions::Frequency] = Units[Dimensions::Time] + "^-1";
	Units[Dimensions::Charge] = Units[Dimensions::Time] + " " +
		Units[Dimensions::Current];
	Units[Dimensions::deg] = "deg";
	Units[Dimensions::rad] = "rad";
	Units[Dimensions::MassFlow] = Units[Dimensions::Mass] + " " + 
		Units[Dimensions::Time] + "^-1";
	Units[Dimensions::Jerk] = Units[Dimensions::Mass] + " " +
        Units[Dimensions::Time] + "^-3";
	Units[Dimensions::VoltageDerivative] = Units[Dimensions::Voltage] + " " +
        Units[Dimensions::Time] + "^-1";
	Units[Dimensions::UnknownDimension] = "UnknownDimension";
};

void MBUnits::SetUnspecifiedUnits() {
	for (auto i = DimensionNames.begin(); i != DimensionNames.end(); i++) {
		Units[i->first] = i->second;
	}
}

void MBUnits::SetMKSUnits(std::unordered_map<Dimensions, std::string>& Units) {
	Units[Dimensions::Length] = "m";
	Units[Dimensions::Mass] = "kg";
	Units[Dimensions::Time] = "s";
	Units[Dimensions::Current] = "A";
	Units[Dimensions::Temperature] = "K";
	SetDerivedUnits(Units);
	Units[Dimensions::Force] = "N";
	Units[Dimensions::Moment] = "N m";
	Units[Dimensions::Work] = "J";
	Units[Dimensions::Power] = "W";
	Units[Dimensions::Pressure] = "Pa";
	Units[Dimensions::Voltage] = "V";
	Units[Dimensions::Charge] = "C";
	Units[Dimensions::Frequency] = "Hz";
};

void MBUnits::SetCGSUnits(std::unordered_map<Dimensions, std::string>& Units) {
	Units[Dimensions::Length] = "cm";
	Units[Dimensions::Mass] = "kg";
	Units[Dimensions::Time] = "s";
	Units[Dimensions::Current] = "A";
	Units[Dimensions::Temperature] = "K";
	SetDerivedUnits(Units);
	Units[Dimensions::Force] = "dyn";
	Units[Dimensions::Pressure] = "dyn cm^-2";
	Units[Dimensions::Moment] = "dyn cm";
	Units[Dimensions::Work] = "erg";
	Units[Dimensions::Power] = "erg s^-1";
	Units[Dimensions::Frequency] = "Hz";
	Units[Dimensions::Charge] = "C";
}

void MBUnits::SetMMTMSUnits(std::unordered_map<Dimensions, std::string>& Units) {
	Units[Dimensions::Length] = "mm";
	Units[Dimensions::Mass] = "ton";
	Units[Dimensions::Time] = "ms";
	Units[Dimensions::Current] = "A";
	Units[Dimensions::Temperature] = "K";
	SetDerivedUnits(Units);
	Units[Dimensions::Force] = "N";
	Units[Dimensions::Moment] = "N mm";
	Units[Dimensions::Work] = "N mm";
	Units[Dimensions::Power] = "N mm s^-1";
	Units[Dimensions::Pressure] = "MPa";
	Units[Dimensions::Frequency] = "kHz";
	Units[Dimensions::Charge] = "mC";
}

void MBUnits::SetMMKGMSUnits(std::unordered_map<Dimensions, std::string>& Units) {
	Units[Dimensions::Length] = "mm";
	Units[Dimensions::Mass] = "kg";
	Units[Dimensions::Time] = "ms";
	Units[Dimensions::Current] = "A";
	Units[Dimensions::Temperature] = "K";
	SetDerivedUnits(Units);
	Units[Dimensions::Force] = "kN";
	Units[Dimensions::Moment] = "N m";
	Units[Dimensions::Work] = "N m";
	Units[Dimensions::Power] = "N m ms^-1";
	Units[Dimensions::Pressure] = "GPa";
	Units[Dimensions::Work] = "kN mm";
	Units[Dimensions::Frequency] = "kHz";
	Units[Dimensions::Charge] = "mC";
}


