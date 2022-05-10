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
 * ACKNOWLEDGEMENTS:
 * Support for output with NetCDF is based on a contribution
 * by Patrick Rix <patrick.rix@online.de>
 */

/* gestore dell'output */

#ifndef UNITS_H
#define UNITS_H

/* se #define DEBUG_COUT l'output avviene su cout anziche' nei files */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <unordered_map>


#include "myassert.h"
#include "except.h"
#include "solman.h"
#include "filename.h"

class MBDynParser;

/* MBUnits - begin */

class MBUnits {
public:
	enum struct Dimensions {
		Dimensionless,
		Boolean,
		Length,
		Mass,
		Time,
		Current,
		Temperature,
		Angle,
		Area,
		Force,
		Velocity,
		Acceleration,
		AngularVelocity,
		AngularAcceleration,

		Momentum,
		MomentaMoment,
		MomentumDerivative,
		MomentaMomentDerivative,

		StaticMoment,
		MomentOfInertia,

		LinearStrain,
		AngularStrain,
		LinearStrainRate,
		AngularStrainRate,
		
		ForceUnitSpan,

		Work,
		Power,
		Pressure,
		Moment,
		Voltage,
		Charge,
		Frequency,
		deg,
		rad,

		/* added for GetEquationDimension method of DofOwnerOwner class */
		MassFlow,
		Jerk,
		VoltageDerivative,
		UnknownDimension
	};

public:
	
	void ReadOutputUnits(std::ostream& Log, MBDynParser& HP);
	MBUnits() {
		SetUnspecifiedUnits();
	}

/* Unit system related stuff */
private:
	std::unordered_map<Dimensions, std::string> Units;
	void SetDerivedUnits(std::unordered_map<Dimensions, std::string>& Units );
	void SetUnspecifiedUnits();
	void SetMKSUnits(std::unordered_map<Dimensions, std::string>& Units);
	void SetCGSUnits(std::unordered_map<Dimensions, std::string>& Units);
	void SetMMTMSUnits(std::unordered_map<Dimensions, std::string>& Units);
	void SetMMKGMSUnits(std::unordered_map<Dimensions, std::string>& Units);
}; /* End class MBUnits */

#endif /* UNITS_H */
