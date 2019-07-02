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
 AUTHOR: Reinhard Resch <r.resch@a1.net>
        Copyright (C) 2015(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"
#endif

#include <vector>
#include <drive.h>

#include <myassert.h>
#include <except.h>
#include <mynewmem.h>
#include <dataman.h>

#include "module-switch_drive.h"

class SwitchDriveCaller : public DriveCaller
{
public:
	enum Type
	{
		SDC_FLOOR,
		SDC_CEIL,
		SDC_NEAREST
	};
	explicit SwitchDriveCaller(const DriveHandler* pDH,  
                               const DriveOwner& oSwitch, 
                               Type eType,
                               const std::vector<DriveOwner>& drives);
	virtual ~SwitchDriveCaller();
	bool bIsDifferentiable(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	doublereal dGet(const doublereal& dVar) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	DriveCaller* pCopy(void) const;

private:
	const DriveCaller* pGetDrive() const;

    DriveOwner oSwitch;
	Type eType;
	typedef std::vector<DriveOwner>::const_iterator iterator;
	const std::vector<DriveOwner> rgDrives;
};

struct SwitchDriveDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

SwitchDriveCaller::SwitchDriveCaller(const DriveHandler* pDH,
                                     const DriveOwner& oSwitch,
                                     Type eType, 
                                     const std::vector<DriveOwner>& drives)
: DriveCaller(pDH),
  oSwitch(oSwitch),
  eType(eType),
  rgDrives(drives)
{
	NO_OP;
}

SwitchDriveCaller::~SwitchDriveCaller()
{
	NO_OP;
}

const DriveCaller* SwitchDriveCaller::pGetDrive() const
{
    const doublereal dSwitch = oSwitch.dGet();

	integer iDrive;

	switch (eType) {
	case SDC_FLOOR:
		iDrive = floor(dSwitch);
		break;

	case SDC_CEIL:
		iDrive = ceil(dSwitch);
		break;

	case SDC_NEAREST:
		iDrive = round(dSwitch);
		break;

	default:
		ASSERT(false);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer iMaxIndex = rgDrives.size() - 1;

	if (iDrive < 0 || iDrive > iMaxIndex) {
		silent_cerr("drive caller(" << GetLabel()
				<< "): argument " << iDrive
                << " from drive caller " << oSwitch.pGetDriveCaller()->GetLabel()
				<< " is not in range [0:" << iMaxIndex << "]"
				<< std::endl);

		throw ErrIndexOutOfRange(iDrive, 0, iMaxIndex, MBDYN_EXCEPT_ARGS);
	}

	return rgDrives[iDrive].pGetDriveCaller();
}

doublereal
SwitchDriveCaller::dGet(const doublereal& dVar) const
{
	return pGetDrive()->dGet(dVar);
}

doublereal SwitchDriveCaller::dGetP(const doublereal& dVar) const
{
	return pGetDrive()->dGetP(dVar);
}

bool SwitchDriveCaller::bIsDifferentiable(void) const
{
	for (iterator i = rgDrives.begin(); i != rgDrives.end(); ++i) {
		if (!i->bIsDifferentiable()) {
			return false;
		}
	}

	return true;
}

/* Restart */
std::ostream&
SwitchDriveCaller::Restart(std::ostream& out) const
{
	out << "switch, " << (eType == SDC_FLOOR ? "floor, " : "ceil, ")
		<< rgDrives.size() << ", ";

	for (iterator i = rgDrives.begin(); i != rgDrives.end(); ++i) {
		i->pGetDriveCaller()->Restart(out);

		if (rgDrives.end() - i > 1) {
			out << ", ";
		}
	}

	return out;
}


DriveCaller*
SwitchDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
			SwitchDriveCaller,
			SwitchDriveCaller(pGetDrvHdl(), oSwitch, eType, rgDrives));

	return pDC;
}

DriveCaller *
SwitchDriveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "switch");

	const DriveHandler* pDrvHdl = 0;

	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	/* driver legato ad un grado di liberta' nodale */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no DOF dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pDC = 0;

    const DriveOwner oSwitch(HP.GetDriveCaller());

	SwitchDriveCaller::Type eType;

	if (HP.IsKeyWord("floor")) {
		eType = SwitchDriveCaller::SDC_FLOOR;
	} else if (HP.IsKeyWord("ceil")) {
		eType = SwitchDriveCaller::SDC_CEIL;
	} else if (HP.IsKeyWord("nearest")) {
		eType = SwitchDriveCaller::SDC_NEAREST;
	} else {
		silent_cerr("Keyword \"floor\" or \"ceil\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer iNumDrives = HP.GetInt();

	if (iNumDrives < 1) {
		silent_cerr("At least one drive caller expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<DriveOwner> rgDrives;

	rgDrives.reserve(iNumDrives);

	for (int i = 0; i < iNumDrives; ++i) {
		rgDrives.push_back(HP.GetDriveCaller());
	}

	SAFENEWWITHCONSTRUCTOR(pDC,
		SwitchDriveCaller,
		SwitchDriveCaller(pDrvHdl, oSwitch, eType, rgDrives));

	return pDC;
}

bool
switch_drive_set()
{
	DriveCallerRead	*rf = new SwitchDriveDCR;

	if (!SetDriveCallerData("switch", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!switch_drive_set()) {
		silent_cerr("switch_drive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // ! STATIC_MODULES
