/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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
 * Wind profiles
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "windprof.h"

/* PowerLawWindProfile - begin */

PowerLawWindProfile::PowerLawWindProfile(
	const doublereal dZRef,
	const doublereal dVRef,
	const doublereal dPower)
: dZRef(dZRef),
dVRef(dVRef),
dPower(dPower)
{
	ASSERT(dZRef > 0.);
	ASSERT(dVRef > 0.);
	// NOTE: should be < 1?
	ASSERT(dPower > 0.);
}

PowerLawWindProfile::~PowerLawWindProfile(void)
{
	NO_OP;
}

bool
PowerLawWindProfile::GetVelocity(const Vec3& X, Vec3& V) const
{
	doublereal dZ = X(3);
	if (dZ <= 0.) {
		return false;
	}

	V(1) = dVRef*std::pow(dZ/dZRef, dPower);
	V(2) = 0.;
	V(3) = 0.;

	return true;
}

std::ostream&
PowerLawWindProfile::Restart(std::ostream& out) const
{
	return out << "power law"
		<< ", reference elevation, " << dZRef
		<< ", reference velocity, " << dVRef
		<< ", exponent, " << dPower;
}

PowerLawGR::~PowerLawGR(void)
{
	NO_OP;
}

Gust *
PowerLawGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	doublereal dZRef = -1.;
	doublereal dVRef = -1.;
	doublereal dPower = -1.;

	while (HP.IsArg()) {
		if (HP.IsKeyWord("reference" "elevation")) {
			if (dZRef != -1.) {
				silent_cerr("PowerLawWindProfile: "
					"reference elevation provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dZRef = HP.GetReal();

		} else if (HP.IsKeyWord("reference" "velocity")) {
			if (dVRef != -1.) {
				silent_cerr("PowerLawWindProfile: "
					"reference velocity provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dVRef = HP.GetReal();

		} else if (HP.IsKeyWord("exponent")) {
			if (dPower != -1.) {
				silent_cerr("PowerLawWindProfile: "
					"exponent provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dPower = HP.GetReal();

		} else {
			break;
		}
	}

	if (dZRef == -1.) {
		silent_cerr("PowerLawWindProfile: "
			"reference elevation missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dVRef == -1.) {
		silent_cerr("PowerLawWindProfile: "
			"reference velocity missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dPower == -1.) {
		silent_cerr("PowerLawWindProfile: "
			"exponent missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Gust *pG = 0;
	SAFENEWWITHCONSTRUCTOR(pG, PowerLawWindProfile,
		PowerLawWindProfile(dZRef, dVRef, dPower));

	return pG;
}

/* PowerLawWindProfile - end */

/* LogarithmicWindProfile - begin */

LogarithmicWindProfile::LogarithmicWindProfile(
	const doublereal dZRef,
	const doublereal dVRef,
	const doublereal dSurfaceRoughnessLength)
: dZRef(dZRef),
dVRef(dVRef),
dSurfaceRoughnessLength(dSurfaceRoughnessLength),
logZRefZ0(std::log(dZRef/dSurfaceRoughnessLength))
{
	ASSERT(dZRef > 0.);
	ASSERT(dVRef > 0.);
	ASSERT(dSurfaceRoughnessLength > 0.);
}

LogarithmicWindProfile::~LogarithmicWindProfile(void)
{
	NO_OP;
}

bool
LogarithmicWindProfile::GetVelocity(const Vec3& X, Vec3& V) const
{
	doublereal dZ = X(3);
	if (dZ <= 0.) {
		return false;
	}

	V(1) = dVRef*(std::log(dZ/dSurfaceRoughnessLength) - 0.)/(logZRefZ0 - 0.);
	V(2) = 0.;
	V(3) = 0.;

	return true;
}

std::ostream&
LogarithmicWindProfile::Restart(std::ostream& out) const
{
	return out << "logarithmic"
		<< ", reference elevation, " << dZRef
		<< ", reference velocity, " << dVRef
		<< ", surface roughness length, " << dSurfaceRoughnessLength;
}

LogarithmicGR::~LogarithmicGR(void)
{
	NO_OP;
}

Gust *
LogarithmicGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	doublereal dZRef = -1.;
	doublereal dVRef = -1.;
	doublereal dSurfaceRoughnessLength = -1.;

	while (HP.IsArg()) {
		if (HP.IsKeyWord("reference" "elevation")) {
			if (dZRef != -1.) {
				silent_cerr("LogarithmicWindProfile: "
					"reference elevation provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dZRef = HP.GetReal();

		} else if (HP.IsKeyWord("reference" "velocity")) {
			if (dVRef != -1.) {
				silent_cerr("LogarithmicWindProfile: "
					"reference velocity provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dVRef = HP.GetReal();

		} else if (HP.IsKeyWord("surface" "roughness" "length")) {
			if (dSurfaceRoughnessLength != -1.) {
				silent_cerr("LogarithmicWindProfile: "
					"surface roughness length provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dSurfaceRoughnessLength = HP.GetReal();

		} else {
			break;
		}
	}

	if (dZRef == -1.) {
		silent_cerr("LogarithmicWindProfile: "
			"reference elevation missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dVRef == -1.) {
		silent_cerr("LogarithmicWindProfile: "
			"reference velocity missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dSurfaceRoughnessLength == -1.) {
		silent_cerr("LogarithmicWindProfile: "
			"surface roughness length missing "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Gust *pG = 0;
	SAFENEWWITHCONSTRUCTOR(pG, LogarithmicWindProfile,
		LogarithmicWindProfile(dZRef, dVRef, dSurfaceRoughnessLength));

	return pG;
}

/* LogarithmicWindProfile - end */

