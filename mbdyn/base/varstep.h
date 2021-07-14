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

/* variable step file driver */

#ifndef VARSTEP_H
#define VARSTEP_H

#include "drive.h"

/* VariableStepFileDrive - begin */

class VariableStepFileDrive : public FileDrive {
protected:
	integer iNumSteps;
	integer iCurrStep;
	bool bLinear;
	bool bPadZeroes;
	Bailout boWhen;

	doublereal* pd;
	doublereal** pvd;

public:
	VariableStepFileDrive(unsigned int uL, const DriveHandler* pDH,
			const char* const sFileName, integer id,
			bool bl, bool pz, Drive::Bailout bo);
	virtual ~VariableStepFileDrive(void);

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

/* VariableStepFileDrive - end */

class DataManager;
class MBDynParser;

struct VariableStepDR : public DriveRead {
public:
	virtual Drive *
	Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP);
};

#endif // VARSTEP_H

