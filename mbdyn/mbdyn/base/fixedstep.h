/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* fixed step file driver */

#ifndef FIXEDSTEP_H
#define FIXEDSTEP_H

#include <drive.h>

/* FixedStepFileDrive - begin */

class FixedStepFileDrive : public FileDrive {
protected:
	doublereal dT0;
	doublereal dDT;
	integer iNumSteps;
	bool bLinear;
	bool bPadZeroes;
	Bailout boWhen;

	doublereal* pd;
	doublereal** pvd;

public:
	FixedStepFileDrive(unsigned int uL, const DriveHandler* pDH,
			const char* const sFileName, integer is, integer id,
			doublereal t0, doublereal dt,
			bool bl, bool pz, Drive::Bailout bo);
	virtual ~FixedStepFileDrive(void);

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

/* FixedStepFileDrive - end */

class DataManager;
class MBDynParser;

struct FixedStepDR : public DriveRead {
public:
	virtual Drive *
	Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP);
};

#endif /* FIXEDSTEP_H */

