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
 * Wind profiles from
 * TurbSim User's Guide, Revised February 1, 2007 for Version 1.21
 * B.J. Jonkman and M.L. Buhl, Jr.
 * Technical Report NREL/TP-500-41136 April 2007
 */

#ifndef WINDPROF_H
#define WINDPROF_H

#include "aerodyn.h"
#include "ScalarFunctions.h"

/* WindProfile - begin */

class WindProfile : public Gust {
protected:
	const Vec3 X0;
	const Mat3x3 R0;
	
public:
	WindProfile(const Vec3& X0, const Mat3x3& R0);
	virtual ~WindProfile(void);
};

/* WindProfile - end */

/* ScalarFuncWindProfile - begin */

class ScalarFuncWindProfile : public WindProfile {
protected:
	const BasicScalarFunction *sf;

public:
	ScalarFuncWindProfile(const Vec3& X0, const Mat3x3& R0,
		const BasicScalarFunction *sf);
	virtual ~ScalarFuncWindProfile(void);
	virtual bool GetVelocity(const Vec3& X, Vec3& V) const;
	virtual std::ostream& Restart(std::ostream& out) const;
};

struct ScalarFuncGR : public GustRead {
public:
	virtual ~ScalarFuncGR(void);
	virtual Gust *
	Read(const DataManager* pDM, MBDynParser& HP);
};

/* ScalarFuncWindProfile - end */

/* PowerLawWindProfile - begin */

class PowerLawWindProfile : public WindProfile {
protected:
	const doublereal dZRef;
	DriveOwner VRef;
	const doublereal dPower;

public:
	PowerLawWindProfile(const Vec3& X0, const Mat3x3& R0,
		const doublereal dZRef, const DriveCaller *pVRef,
		const doublereal dPower);
	virtual ~PowerLawWindProfile(void);
	virtual bool GetVelocity(const Vec3& X, Vec3& V) const;
	virtual std::ostream& Restart(std::ostream& out) const;
};

struct PowerLawGR : public GustRead {
public:
	virtual ~PowerLawGR(void);
	virtual Gust *
	Read(const DataManager* pDM, MBDynParser& HP);
};

/* PowerLawWindProfile - end */

/* LogarithmicWindProfile - begin */

class LogarithmicWindProfile : public WindProfile {
protected:
	const doublereal dZRef;
	DriveOwner VRef;
	const doublereal dSurfaceRoughnessLength;

	const doublereal logZRefZ0;
	// TODO: stability function
	
public:
	LogarithmicWindProfile(const Vec3& X0, const Mat3x3& R0,
		const doublereal dZRef, const DriveCaller *pVRef,
		const doublereal dSurfaceRoughnessLength);
	virtual ~LogarithmicWindProfile(void);
	virtual bool GetVelocity(const Vec3& X, Vec3& V) const;
	virtual std::ostream& Restart(std::ostream& out) const;
};

struct LogarithmicGR : public GustRead {
public:
	virtual ~LogarithmicGR(void);
	virtual Gust *
	Read(const DataManager* pDM, MBDynParser& HP);
};

/* LogarithmicWindProfile - end */

#endif // WINDPROF_H

