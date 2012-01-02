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
 * Wind profiles
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "mbpar.h"

#include "windprof.h"

/* WindProfile - begin */

WindProfile::WindProfile(
	const Vec3& X0,
	const Mat3x3& R0)
: X0(X0), R0(R0)
{
	NO_OP;
}

WindProfile::~WindProfile(void)
{
	NO_OP;
}

/* WindProfile - end */

/* ScalarFuncWindProfile - begin */

ScalarFuncWindProfile::ScalarFuncWindProfile(
	const Vec3& X0,
	const Mat3x3& R0,
	const BasicScalarFunction *sf)
: WindProfile(X0, R0),
sf(sf)
{
	ASSERT(sf != 0);
}

ScalarFuncWindProfile::~ScalarFuncWindProfile(void)
{
	NO_OP;
}

bool
ScalarFuncWindProfile::GetVelocity(const Vec3& X, Vec3& V) const
{
	// dZ is the component of (X - X0) along axis 3 of R0
	doublereal dZ = R0.GetVec(3)*(X - X0);

	// V is projected along axis 1 of R0
	V = R0.GetVec(1)*(*sf)(dZ);

	return true;
}

std::ostream&
ScalarFuncWindProfile::Restart(std::ostream& out) const
{
	return out << "scalar function, "
		<< ", reference position, " << X0
		<< ", reference orientation, " << R0
		<< " # not implemented yet";
}

ScalarFuncGR::~ScalarFuncGR(void)
{
	NO_OP;
}

Gust *
ScalarFuncGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	const BasicScalarFunction *sf = 0;
	Vec3 X0(Zero3);
	bool bGotX0 = false;
	Mat3x3 R0(Eye3);
	bool bGotR0 = false;

	while (HP.IsArg()) {
		if (HP.IsKeyWord("reference" "position")) {
			if (bGotX0) {
				silent_cerr("ScalarFuncWindProfile: "
					"reference position provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			X0 = HP.GetVecAbs(AbsRefFrame);

		} else if (HP.IsKeyWord("reference" "orientation")) {
			if (bGotR0) {
				silent_cerr("ScalarFuncWindProfile: "
					"reference orientation provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			R0 = HP.GetRotAbs(AbsRefFrame);

		} else {
			break;
		}
	}

	if (!HP.IsArg()) {
		silent_cerr("ScalarFuncWindProfile: "
			"function expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	sf = ParseScalarFunction(HP, const_cast<DataManager *const>(pDM));
	if (sf == 0) {
		silent_cerr("ScalarFuncWindProfile: "
			"unable to parse scalar function "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Gust *pG = 0;
	SAFENEWWITHCONSTRUCTOR(pG, ScalarFuncWindProfile,
		ScalarFuncWindProfile(X0, R0, sf));

	return pG;
}

/* ScalarFuncWindProfile - end */

/* PowerLawWindProfile - begin */

PowerLawWindProfile::PowerLawWindProfile(
	const Vec3& X0,
	const Mat3x3& R0,
	const doublereal dZRef,
	const DriveCaller *pVRef,
	const doublereal dPower)
: WindProfile(X0, R0),
dZRef(dZRef),
VRef(pVRef),
dPower(dPower)
{
	ASSERT(dZRef > 0.);
	ASSERT(pVRef != 0);
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
	// dZ is the component of (X - X0) along axis 3 of R0
	doublereal dZ = R0.GetVec(3)*(X - X0);
	if (dZ <= 0.) {
		return false;
	}

	// V is projected along axis 1 of R0
	V = R0.GetVec(1)*(VRef.dGet()*std::pow(dZ/dZRef, dPower));

	return true;
}

std::ostream&
PowerLawWindProfile::Restart(std::ostream& out) const
{
	out << "power law"
		<< ", reference position, " << X0
		<< ", reference orientation, " << R0
		<< ", reference elevation, " << dZRef
		<< ", reference velocity, ", VRef.pGetDriveCaller()->Restart(out)
		<< ", exponent, " << dPower;
	return out;
}

PowerLawGR::~PowerLawGR(void)
{
	NO_OP;
}

Gust *
PowerLawGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	doublereal dZRef = -1.;
	DriveCaller *pVRef = 0;
	doublereal dPower = -1.;
	Vec3 X0(Zero3);
	bool bGotX0 = false;
	Mat3x3 R0(Eye3);
	bool bGotR0 = false;

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
			if (pVRef != 0) {
				silent_cerr("PowerLawWindProfile: "
					"reference velocity provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			pVRef = HP.GetDriveCaller();

		} else if (HP.IsKeyWord("exponent")) {
			if (dPower != -1.) {
				silent_cerr("PowerLawWindProfile: "
					"exponent provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dPower = HP.GetReal();

		} else if (HP.IsKeyWord("reference" "position")) {
			if (bGotX0) {
				silent_cerr("PowerLawWindProfile: "
					"reference position provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			X0 = HP.GetVecAbs(AbsRefFrame);

		} else if (HP.IsKeyWord("reference" "orientation")) {
			if (bGotR0) {
				silent_cerr("PowerLawWindProfile: "
					"reference orientation provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			R0 = HP.GetRotAbs(AbsRefFrame);

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

	if (pVRef == 0) {
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
		PowerLawWindProfile(X0, R0, dZRef, pVRef, dPower));

	return pG;
}

/* PowerLawWindProfile - end */

/* LogarithmicWindProfile - begin */

LogarithmicWindProfile::LogarithmicWindProfile(
	const Vec3& X0,
	const Mat3x3& R0,
	const doublereal dZRef,
	const DriveCaller *pVRef,
	const doublereal dSurfaceRoughnessLength)
: WindProfile(X0, R0),
dZRef(dZRef),
VRef(pVRef),
dSurfaceRoughnessLength(dSurfaceRoughnessLength),
logZRefZ0(std::log(dZRef/dSurfaceRoughnessLength))
{
	ASSERT(dZRef > 0.);
	ASSERT(pVRef != 0);
	ASSERT(dSurfaceRoughnessLength > 0.);
}

LogarithmicWindProfile::~LogarithmicWindProfile(void)
{
	NO_OP;
}

bool
LogarithmicWindProfile::GetVelocity(const Vec3& X, Vec3& V) const
{
	// dZ is the component of (X - X0) along axis 3 of R0
	doublereal dZ = R0.GetVec(3)*(X - X0);
	if (dZ <= 0.) {
		return false;
	}

	// V is projected along axis 1 of R0
	V = R0.GetVec(1)*(VRef.dGet()*(std::log(dZ/dSurfaceRoughnessLength) - 0.)/(logZRefZ0 - 0.));

	return true;
}

std::ostream&
LogarithmicWindProfile::Restart(std::ostream& out) const
{
	out << "logarithmic"
		<< ", reference position, " << X0
		<< ", reference orientation, " << R0
		<< ", reference elevation, " << dZRef
		<< ", reference velocity, ", VRef.pGetDriveCaller()->Restart(out)
		<< ", surface roughness length, " << dSurfaceRoughnessLength;
	return out;
}

LogarithmicGR::~LogarithmicGR(void)
{
	NO_OP;
}

Gust *
LogarithmicGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	doublereal dZRef = -1.;
	DriveCaller *pVRef = 0;
	doublereal dSurfaceRoughnessLength = -1.;
	Vec3 X0(Zero3);
	bool bGotX0 = false;
	Mat3x3 R0(Eye3);
	bool bGotR0 = false;

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
			if (pVRef != 0) {
				silent_cerr("LogarithmicWindProfile: "
					"reference velocity provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			pVRef = HP.GetDriveCaller();

		} else if (HP.IsKeyWord("surface" "roughness" "length")) {
			if (dSurfaceRoughnessLength != -1.) {
				silent_cerr("LogarithmicWindProfile: "
					"surface roughness length provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			dSurfaceRoughnessLength = HP.GetReal();

		} else if (HP.IsKeyWord("reference" "position")) {
			if (bGotX0) {
				silent_cerr("PowerLawWindProfile: "
					"reference position provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			X0 = HP.GetVecAbs(AbsRefFrame);

		} else if (HP.IsKeyWord("reference" "orientation")) {
			if (bGotR0) {
				silent_cerr("PowerLawWindProfile: "
					"reference orientation provided twice "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			R0 = HP.GetRotAbs(AbsRefFrame);

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

	if (pVRef == 0) {
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
		LogarithmicWindProfile(X0, R0,
			dZRef, pVRef, dSurfaceRoughnessLength));

	return pG;
}

/* LogarithmicWindProfile - end */

