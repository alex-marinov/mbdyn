/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "gust.h"
#include "dataman.h"
#include "drive_.h"

#include "windprof.h"

/* Gust - begin */

Gust::~Gust(void)
{
	NO_OP;
}

void
Gust::SetAirProperties(const AirProperties *pap)
{
	ASSERT(pap != 0);
	pAP = pap;
}

Vec3
Gust::GetVelocity(const Vec3& X) const
{
	Vec3 V(Zero3);
	GetVelocity(X, V);
	return V;
}

GustRead::~GustRead(void)
{
	NO_OP;
}

/* Gust - end */

/* Gust1D - begin */

class Gust1D : public Gust {
private:
	Vec3 FrontDir;
	Vec3 GustDir;
	doublereal dVRef;
	DriveOwner Time;
	DriveOwner GustProfile;

public:
	Gust1D(const Vec3& f, const Vec3& g, const doublereal& v,
			DriveCaller *pT, DriveCaller *pG);
	~Gust1D(void);

	virtual std::ostream& Restart(std::ostream& out) const;
	bool GetVelocity(const Vec3& X, Vec3& V) const;
};

Gust1D::Gust1D(const Vec3& f, const Vec3& g, const doublereal& v,
	DriveCaller *pT, DriveCaller *pG)
: FrontDir(f),
GustDir(g),
dVRef(v),
Time(pT),
GustProfile(pG)
{
	ASSERT(pT != NULL);
	ASSERT(pG != NULL);
}

Gust1D::~Gust1D(void)
{
	NO_OP;
}

std::ostream&
Gust1D::Restart(std::ostream& out) const
{
	out << "front 1D, ",
		FrontDir.Write(out, ", ")
		<< ", ", GustDir.Write(out, ", ")
		<< ", " << dVRef
		<< ", ", GustProfile.pGetDriveCaller()->Restart(out);
	return out;
}

bool
Gust1D::GetVelocity(const Vec3& X, Vec3& V) const
{
	doublereal x = FrontDir*X + dVRef*Time.dGet();
	doublereal v = GustProfile.dGet(x);
	V = GustDir*v;
	return true;
}

/* reads a front 1D gust */
struct Gust1DGR : public GustRead {
public:
	virtual ~Gust1DGR(void);
	virtual Gust *
	Read(const DataManager* pDM, MBDynParser& HP);
};

Gust1DGR::~Gust1DGR(void)
{
	NO_OP;
}

Gust *
Gust1DGR::Read(const DataManager* pDM, MBDynParser& HP)
{
	/* front direction */
	Vec3 f = HP.GetVecAbs(::AbsRefFrame);

	/* gust velocity direction */
	Vec3 g = HP.GetVecAbs(::AbsRefFrame);

	/* reference velocity */
	doublereal v = HP.GetReal();

	/* time drive caller
	 * FIXME: not needed if v = 0 */
	DriveCaller *pT = NULL;
	SAFENEWWITHCONSTRUCTOR(pT, TimeDriveCaller,
			TimeDriveCaller(pDM->pGetDrvHdl()));

	/* gust profile drive caller */
	DriveCaller *pP = HP.GetDriveCaller();

	/* gust */
	Gust *pG = 0;
	SAFENEWWITHCONSTRUCTOR(pG, Gust1D, Gust1D(f, g, v, pT, pP));

	return pG;
}

/* Gust1D - end */

/* bag that contains functions to parse gusts */

typedef std::map<std::string, GustRead *, ltstrcase> GustFuncMapType;
static GustFuncMapType GustFuncMap;

struct GustWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return GustFuncMap.find(s) != GustFuncMap.end();
	};
};

static GustWordSetType GustWordSet;

bool
SetGustData(const char *name, GustRead *rf)
{
	pedantic_cout("registering gust \"" << name << "\"" << std::endl);
	return GustFuncMap.insert(GustFuncMapType::value_type(name, rf)).second;
}

// read gust data
Gust *
ReadGustData(const DataManager* pDM, MBDynParser& HP)
{
	DEBUGCOUTFNAME("ReadGustData()");

	const char *s = HP.IsWord(GustWordSet);
	if (s == 0) {
		return 0;
	}

	GustFuncMapType::iterator func = GustFuncMap.find(std::string(s));
	if (func == GustFuncMap.end()) {
		silent_cerr("unknown gust type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

static unsigned done;

void
InitGustData(void)
{
	if (::done++ > 0) {
		return;
	}

	SetGustData("front" "1D", new Gust1DGR);
	SetGustData("scalar" "function", new ScalarFuncGR);
	SetGustData("power" "law", new PowerLawGR);
	SetGustData("logarithmic", new LogarithmicGR);

	/* NOTE: add here initialization of new built-in drive callers;
	 * alternative ways to register new custom drive callers are:
	 * - call SetDriveData() from anywhere in the code
	 * - write a module that calls SetDriveData() from inside a function
	 *   called module_init(), and run-time load it using "module load"
	 *   in the input file.
	 */
}

void
DestroyGustData(void)
{
	if (::done == 0) {
		silent_cerr("DestroyGustData() called once too many" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--::done > 0) {
		return;
	}

	/* free stuff */
	for (GustFuncMapType::iterator i = GustFuncMap.begin(); i != GustFuncMap.end(); ++i) {
		delete i->second;
	}
	GustFuncMap.clear();
}

