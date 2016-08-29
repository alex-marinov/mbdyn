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

 /*
 * With the contribution of Ankit Aggarwal <ankit.ankit.aggarwal@gmail.com>
 * during Google Summer of Code 2016
 */
 

#ifndef TIMESTEPCONTROL_H
#define TIMESTEPCONTROL_H

#include "mbconfig.h"

#include <map>

class TimeStepControl;
#include "solver.h"
#include "drive_.h"

extern void InitTimeStepData(void);
extern void DestroyTimeStepData(void);


class TimeStepControl {
protected:
	doublereal dCurrTimeStep;

public:
	TimeStepControl(void) { NO_OP; };
	virtual ~TimeStepControl(void) { NO_OP; };
	virtual doublereal dGetNewStepTime(StepIntegrator::StepChange currStep, doublereal iPerformedIters) = 0;
	virtual void SetDriveHandler(const DriveHandler* driveHandler) = 0;
	virtual void Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep) = 0;
};

extern TimeStepControl *ReadTimeStepData(Solver *s, MBDynParser& HP);

class TimeStepRead {
public:
	virtual ~TimeStepRead(void) { NO_OP; };
	virtual TimeStepControl *Read(Solver *s, MBDynParser& HP) = 0;
};

class  NoChange : public TimeStepControl {
private:
	Solver *s;

public:
	NoChange(Solver * s);
	~NoChange(void) { NO_OP; };
	doublereal dGetNewStepTime(StepIntegrator::StepChange currStep, doublereal iPerformedIters);
	void SetDriveHandler(const DriveHandler* driveHandler) { NO_OP; };
 	void Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep);
};

class ChangeStep : public TimeStepControl {
private:
	Solver *s;
	DriveCaller* pStrategyChangeDrive;
	doublereal dMinTimeStep;
	DriveOwner MaxTimeStep;

public:
	ChangeStep(Solver *s, DriveCaller* pStrategyChangeDrive);
	~ChangeStep(void) {
		SAFEDELETE(pStrategyChangeDrive);
	};

	doublereal dGetNewStepTime(StepIntegrator::StepChange currStep, doublereal iPerformedIters);
	void SetDriveHandler(const DriveHandler* driveHandler);
	void Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep);
};

class Factor : public TimeStepControl {
private:
	Solver *s;
	doublereal dReductionFactor;
	doublereal iStepsBeforeReduction;
	doublereal dRaiseFactor;
	doublereal iStepsBeforeRaise;
	doublereal iMinIters;
	doublereal iMaxIters;
	bool bLastChance;
	doublereal iStepsAfterReduction;
	doublereal iStepsAfterRaise;
	doublereal iWeightedPerformedIters;
	doublereal dMinTimeStep;
	DriveOwner MaxTimeStep;

public:
	Factor(Solver *s,
		doublereal dReductionFactor,
		doublereal iStepsBeforeReduction,
		doublereal dRaiseFactor,
		doublereal iStepsBeforeRaise,
		doublereal iMinIters,
		doublereal iMaxIters);
	~Factor(void) { NO_OP; };
	doublereal dGetNewStepTime(StepIntegrator::StepChange Why, doublereal iPerformedIters);
	void SetDriveHandler(const DriveHandler* driveHandler) { NO_OP; };
	void Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep);
};

#endif // TIMESTEPCONTROL_H
