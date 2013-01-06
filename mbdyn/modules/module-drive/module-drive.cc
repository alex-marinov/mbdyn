/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "drive.h"

class DummyDriveCaller : public DriveCaller {
private:
	doublereal dConst;

public:
	DummyDriveCaller(const doublereal &d);
	virtual ~DummyDriveCaller(void);
 
	/* Copia */
	virtual DriveCaller* pCopy(void) const;
 
	/* Scrive il contributo del DriveCaller al file di restart */   
	virtual std::ostream& Restart(std::ostream& out) const;
 
	inline doublereal dGet(const doublereal& /* dVar */ ) const;
	inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
};

DummyDriveCaller::DummyDriveCaller(const doublereal &d)
: DriveCaller(0),
dConst(d)
{
	NO_OP;
}

DummyDriveCaller::~DummyDriveCaller(void)
{
	NO_OP;
}

DriveCaller *
DummyDriveCaller::pCopy(void) const
{
	return new DummyDriveCaller(dConst);
}

std::ostream&
DummyDriveCaller::Restart(std::ostream& out) const
{
	return out << "dummy, " << dConst;
}

inline doublereal 
DummyDriveCaller::dGet(const doublereal& /* dVar */ ) const 
{
	return dConst; 
}

inline doublereal
DummyDriveCaller::dGet(void) const
{
	return dConst;
}

inline bool
DummyDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

inline doublereal 
DummyDriveCaller::dGetP(const doublereal& /* dVar */ ) const
{
	return 0.;
}

inline doublereal 
DummyDriveCaller::dGetP(void) const
{
	return 0.;
}

/* prototype of the functional object: reads a drive caller */
struct DummyDCR : public DriveCallerRead {
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) {
		doublereal d = HP.GetReal();
		return new DummyDriveCaller(d);
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	DriveCallerRead	*rf = new DummyDCR;

	if (!SetDriveData("dummy", rf)) {
		delete rf;

		silent_cerr("DummyDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

