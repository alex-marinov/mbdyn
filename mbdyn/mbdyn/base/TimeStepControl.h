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
 

#ifndef TIME_STEP_H
#define TIME_STEP_H

class TimeStepControl;
#include "mbconfig.h"
#include "output.h"
#include "solman.h"
#include <map>
#include "solver.h"
#include "drive_.h"

extern void InitTimeStepData(void);
extern void DestroyTimeStepData(void);
class TimeStepControl {

public:
	TimeStepControl(){};
	virtual ~TimeStepControl(void){};
	virtual doublereal dGetNewStepTime(StepIntegrator::StepChange currStep , doublereal iPerformedIters) = 0;
	virtual void SetDriveHandler(const DriveHandler* driveHandler) = 0;
	virtual void Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep) = 0;
};

extern TimeStepControl * ReadTimeStepData(Solver * s, MBDynParser& HP);

class TimeStepRead {

public:
	virtual ~TimeStepRead(void){};
	virtual TimeStepControl *
	Read( Solver * s, MBDynParser& HP) = 0;
};

class  NoChange : public TimeStepControl {
private:
	Solver * s;
	doublereal dCurrTimeStep;
public:
	NoChange(Solver * s);
	doublereal dGetNewStepTime(StepIntegrator::StepChange currStep , doublereal iPerformedIters);
	void SetDriveHandler(const DriveHandler* driveHandler){}
 	void Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep);
	~NoChange(){};
};


class ChangeStep : public TimeStepControl {
private:
	Solver * s;
	DriveCaller* pStrategyChangeDrive;
	DriveOwner MaxTimeStep;
	doublereal dCurrTimeStep;
public:
	ChangeStep(Solver * s,DriveCaller* pStrategyChangeDrive);
	doublereal dGetNewStepTime(StepIntegrator::StepChange currStep , doublereal iPerformedIters);
	
	void SetDriveHandler(const DriveHandler* driveHandler);

	void Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep);
	~ChangeStep(){
		SAFEDELETE(pStrategyChangeDrive);
	};
};

class Factor : public TimeStepControl {
private:
	Solver * s;
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
	DriveOwner MaxTimeStep;

	doublereal dCurrTimeStep;

public:
	Factor( Solver * s ,doublereal dReductionFactor,
			doublereal iStepsBeforeReduction,
			doublereal dRaiseFactor,
			doublereal iStepsBeforeRaise,
			doublereal iMinIters,
			doublereal iMaxIters );
	doublereal dGetNewStepTime(StepIntegrator::StepChange Why , doublereal iPerformedIters);
	void SetDriveHandler(const DriveHandler* driveHandler){}
	void Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep);
	~Factor(){};
};

/********** READER FOR FACTOR **************/

/***************SET AND CLEAR FOR MAP ************/

#endif
