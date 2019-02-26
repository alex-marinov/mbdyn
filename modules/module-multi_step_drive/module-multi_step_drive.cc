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
 AUTHOR: Reinhard Resch <r.resch@secop.com>
        Copyright (C) 2015(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"
#endif

#include <limits>
#include <vector>
#include <drive.h>

#include <myassert.h>
#include <except.h>
#include <dataman.h>
#include "module-multi_step_drive.h"

class MultiStepDrive : public DriveCaller
{
public:
	struct StepRecord
	{
		StepRecord(doublereal x=0., doublereal y=0.)
		:x(x), y(y) { }
		doublereal x;
		doublereal y;
	};

	explicit MultiStepDrive(const DriveHandler* pDH, const std::vector<StepRecord>& drives);
	virtual ~MultiStepDrive();
	bool bIsDifferentiable(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	doublereal dGet(const doublereal& dVar) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	DriveCaller* pCopy(void) const;

private:
	typedef std::vector<StepRecord>::const_iterator iterator;
	const std::vector<StepRecord> rgSteps;
};

struct MultiStepDriveDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

MultiStepDrive::MultiStepDrive(const DriveHandler* pDH, const std::vector<StepRecord>& rgSteps)
: DriveCaller(pDH),
  rgSteps(rgSteps)
{
	NO_OP;
}

MultiStepDrive::~MultiStepDrive()
{
	NO_OP;
}

doublereal
MultiStepDrive::dGet(const doublereal& dVar) const
{
	ASSERT(rgSteps.size() >= 1u);

	if (dVar < rgSteps.front().x) {
		return rgSteps.front().y;
	}

	for (iterator i = rgSteps.begin(); i != rgSteps.end() - 1; ++i) {
		if (dVar >= i->x && dVar < (i + 1)->x) {
			return i->y;
		}
	}

	return rgSteps.back().y;
}

doublereal MultiStepDrive::dGetP(const doublereal&) const
{
	return 0.;
}

bool MultiStepDrive::bIsDifferentiable(void) const
{
	return false;
}

/* Restart */
std::ostream&
MultiStepDrive::Restart(std::ostream& out) const
{
	out << "multi step, " << rgSteps.size() << ", ";

	for (iterator i = rgSteps.begin(); i != rgSteps.end(); ++i) {
		out << i->x << ", " << i->y;

		if (rgSteps.end() - i > 1) {
			out << ", ";
		}
	}

	return out;
}


DriveCaller*
MultiStepDrive::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
			MultiStepDrive,
			MultiStepDrive(pGetDrvHdl(), rgSteps));

	return pDC;
}

DriveCaller *
MultiStepDriveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "multi step");

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

	const integer iNumSteps = HP.GetInt();

	if (iNumSteps < 1) {
		silent_cerr("At least one step expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<MultiStepDrive::StepRecord> rgSteps;

	rgSteps.reserve(iNumSteps);

	doublereal xPrev = -std::numeric_limits<doublereal>::max();

	for (int i = 0; i < iNumSteps; ++i) {
		const doublereal x = HP.GetReal();

		if (x <= xPrev) {
			silent_cerr("X-values must be in ascending order at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal y = HP.GetReal();

		rgSteps.push_back(MultiStepDrive::StepRecord(x, y));

		xPrev = x;
	}

	SAFENEWWITHCONSTRUCTOR(pDC,
		MultiStepDrive,
		MultiStepDrive(pDrvHdl, rgSteps));

	return pDC;
}

bool
multi_step_drive_set()
{
	DriveCallerRead	*rf = new MultiStepDriveDCR;

	if (!SetDriveCallerData("multi" "step", rf)) {
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
		silent_cerr("multi_step_drive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // ! STATIC_MODULES
