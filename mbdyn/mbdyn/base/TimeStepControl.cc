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

#include "mbconfig.h"

#include "TimeStepControl.h"
#include "solver_impl.h"

typedef std::map<std::string, TimeStepRead *> TimeStepFuncMapType; 

TimeStepFuncMapType TimeStepReadMap;

struct TimeStepWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return TimeStepReadMap.find(s) != TimeStepReadMap.end();
	};
};

TimeStepWordSetType TimeStepWordSet;

bool 
SetTimeStepData(const char *name, TimeStepRead * sr) {
	return TimeStepReadMap.insert(std::make_pair(name, sr)).second;
}

NoChange::NoChange(Solver *s): s(s) {
	NO_OP;
}

doublereal
NoChange::dGetNewStepTime(StepIntegrator::StepChange currStep, doublereal iPerformedIters)
{
	return dCurrTimeStep;
}

void
NoChange::Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep)
{
	dCurrTimeStep = dInitialTimeStep;

	doublereal dInitialMaxTimeStep;
	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();

	} else{
		dInitialMaxTimeStep =  MaxTimeStep.dGet();
	}

	if (dMinTimeStep > dInitialMaxTimeStep) {
		silent_cerr("warning: minimum time step " << dMinTimeStep << " greater than maximum (initial) time step " << dInitialMaxTimeStep << " (ignored)" << std::endl);
	}

	if (dInitialTimeStep < dMinTimeStep) {
		silent_cerr("warning: (constant) time step " << dInitialTimeStep << " less than minimum time step " << dMinTimeStep << " (ignored)" << std::endl);
	}

	if (dInitialTimeStep > dInitialMaxTimeStep) {
		silent_cerr("warning: (constant) time step " << dInitialTimeStep << " greater than maximum (initial) time step " << dInitialMaxTimeStep << " (ignored)" << std::endl);
	}
}

struct NoChangeTSR : public TimeStepRead {
public:
	TimeStepControl *Read(Solver *s, MBDynParser& HP);
};

TimeStepControl *
NoChangeTSR::Read(Solver *s, MBDynParser& HP)
{
	// FIXME: the (constant) time step could be parsed here
	return new NoChange(s);
}

ChangeStep::ChangeStep(Solver *s, DriveCaller* pStrategyChangeDrive)
: s(s),
pStrategyChangeDrive(pStrategyChangeDrive),
dMinTimeStep(-1.)
{
	NO_OP;
}

doublereal
ChangeStep::dGetNewStepTime(StepIntegrator::StepChange currStep, doublereal iPerformedIters)
{
	doublereal dNewStep = pStrategyChangeDrive->dGet();
	doublereal dMaxTimeStep = MaxTimeStep.dGet();;

	// are we sure we intend to change the time step if repeat is requested?
	if (currStep == StepIntegrator::REPEATSTEP && dNewStep == dCurrTimeStep) {
		return (dCurrTimeStep = dMinTimeStep/2.);
	}

	return (dCurrTimeStep = std::max(std::min(dNewStep, dMaxTimeStep), dMinTimeStep));
}

void
ChangeStep::SetDriveHandler(const DriveHandler* driveHandler)
{
	pStrategyChangeDrive->SetDrvHdl(driveHandler);
}

void
ChangeStep::Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep)
{
	this->dCurrTimeStep = dInitialTimeStep;
	this->MaxTimeStep = MaxTimeStep;
	this->dMinTimeStep = dMinTimeStep;

	doublereal dInitialMaxTimeStep ;
	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();

	} else{
		dInitialMaxTimeStep =  MaxTimeStep.dGet();
	}

	if (dMinTimeStep > dInitialMaxTimeStep) {
		silent_cerr("error: minimum time step " << dMinTimeStep << " greater than maximum (initial) time step " << dInitialMaxTimeStep << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dInitialTimeStep < dMinTimeStep) {
		silent_cerr("error: initial time step " << dInitialTimeStep << " less than minimum time step " << dMinTimeStep << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dInitialTimeStep > dInitialMaxTimeStep) {
		silent_cerr("error: initial time step " << dInitialTimeStep << " greater than maximum (initial) time step " << dInitialMaxTimeStep << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

class ChangeStepTSR : public TimeStepRead {
public:	
	TimeStepControl *Read(Solver *s, MBDynParser& HP);	
};

TimeStepControl *
ChangeStepTSR::Read(Solver *s, MBDynParser& HP)
{
	return new ChangeStep(s , HP.GetDriveCaller(true));
}

Factor::Factor(Solver *s,
	doublereal dReductionFactor,
	doublereal iStepsBeforeReduction,
	doublereal dRaiseFactor,
	doublereal iStepsBeforeRaise,
	doublereal iMinIters,
	doublereal iMaxIters)
: s(s),
dReductionFactor(dReductionFactor),
iStepsBeforeReduction(iStepsBeforeReduction),
dRaiseFactor(dRaiseFactor),
iStepsBeforeRaise(iStepsBeforeRaise),
iMinIters(iMinIters),
iMaxIters(iMaxIters),
bLastChance(false),
iStepsAfterReduction(0),
iStepsAfterRaise(0),
iWeightedPerformedIters(0),
dMinTimeStep(::dDefaultMinTimeStep)
{
	NO_OP;
}

doublereal
Factor::dGetNewStepTime(StepIntegrator::StepChange Why, doublereal iPerformedIters)
{
	doublereal dMaxTimeStep = MaxTimeStep.dGet();
	
	switch (Why) {
	case StepIntegrator::REPEATSTEP:
		if (dCurrTimeStep*dReductionFactor > dMinTimeStep) {
			if (bLastChance == true) {
				bLastChance = false;
			}
			iStepsAfterReduction = 0;
			return (dCurrTimeStep = std::min(dCurrTimeStep*dReductionFactor, dMaxTimeStep));

		}

		if (bLastChance == false) {
			bLastChance = true;
			return (dCurrTimeStep = dMinTimeStep);

		}

		return (dCurrTimeStep = std::min(dCurrTimeStep*dReductionFactor, dMaxTimeStep));

	case StepIntegrator::NEWSTEP:
		iStepsAfterReduction++;
		iStepsAfterRaise++;

		iWeightedPerformedIters = (10*iPerformedIters + 9*iWeightedPerformedIters)/10;

		if (iPerformedIters > iMaxIters) {
			iStepsAfterReduction = 0;
			bLastChance = false;
			return (dCurrTimeStep = std::min(std::max(dCurrTimeStep*dReductionFactor, dMinTimeStep), dMaxTimeStep));

		}

		if (iPerformedIters <= iMinIters
				&& iStepsAfterReduction >iStepsBeforeReduction
				&& iStepsAfterRaise > iStepsBeforeRaise
				&& dCurrTimeStep < dMaxTimeStep)
		{
			iStepsAfterRaise = 0;
			iWeightedPerformedIters = 0;
			return (dCurrTimeStep = std::min(dCurrTimeStep*dRaiseFactor, dMaxTimeStep));
		}

		return dCurrTimeStep;

	default:
		// Should Not Reach Over here
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}


void
Factor::Init(integer iMaxIterations, doublereal dMinTimeStep, const DriveOwner& MaxTimeStep, doublereal dInitialTimeStep)
{
	this->dCurrTimeStep = dInitialTimeStep;
	this->MaxTimeStep = MaxTimeStep;
	this->dMinTimeStep = dMinTimeStep;

	if (iMaxIters <= iMinIters) {
		silent_cerr("error: maximum number of iterations " << iMaxIters << " less than or equal to minimum " << iMinIters << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dInitialMaxTimeStep;
	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();
	} else{

		dInitialMaxTimeStep = MaxTimeStep.dGet();
	}

	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(ConstDriveCaller)
			&& MaxTimeStep.dGet() == std::numeric_limits<doublereal>::max())
	{
		silent_cerr("warning: maximum time step not set; the initial time step value " << dInitialTimeStep << " will be used" << std::endl);
	}

	if (dMinTimeStep == ::dDefaultMinTimeStep) {
		silent_cerr("warning: minimum time step not set; the initial time step value " << dInitialTimeStep << " will be used" << std::endl);
		dMinTimeStep = dInitialTimeStep;
	}

	if (dMinTimeStep == dInitialMaxTimeStep) {
		silent_cerr("warning: minimum time step " << dMinTimeStep << "is equal to (initial) maximum time step " << dInitialMaxTimeStep << "; no time step adaptation will take place" << std::endl);

	} else if (dMinTimeStep == dInitialMaxTimeStep) {
		silent_cerr("error: minimum time step " << dMinTimeStep << "is greater than (initial) maximum time step " << dInitialMaxTimeStep << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

class FactorTSR : public TimeStepRead {
public:
	TimeStepControl *Read(Solver *s, MBDynParser& HP);	
};

TimeStepControl *
FactorTSR::Read(Solver *s, MBDynParser& HP)
{
	integer iMinIters = 1;
	integer iMaxIters = 0;
	doublereal dReductionFactor = 1.;
	try {
		dReductionFactor = HP.GetReal(1., HighParser::range_gt_le<doublereal>(0., 1.));

	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("error: invalid reduction factor " << e.Get() << " (must be positive and less than (or equal to) one) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	integer iStepsBeforeReduction = 1;
	try {
		iStepsBeforeReduction = HP.GetInt(1, HighParser::range_ge<integer>(1));

	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("error: invalid number of steps before reduction " << e.Get() << " (must be greater than zero) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	doublereal dRaiseFactor = 1.;
	try {
		dRaiseFactor = HP.GetReal(1., HighParser::range_ge<doublereal>(1.));

	} catch (HighParser::ErrValueOutOfRange<doublereal> e) {
		silent_cerr("error: invalid increase factor " << e.Get() << " (must be greater than (or equal to) one) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	integer iStepsBeforeRaise = 1;
	try {
		iStepsBeforeRaise = HP.GetInt(1, HighParser::range_ge<integer>(1));

	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("error: invalid number of steps before increase " << e.Get() << " (must be greater than zero) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	try {
		iMinIters = HP.GetInt(1, HighParser::range_ge<integer>(1));

	} catch (HighParser::ErrValueOutOfRange<integer> e) {
		silent_cerr("error: invalid number of minimum iterations " << e.Get() << " (must be greater than zero) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
		throw e;
	}

	iMaxIters = 0;
	if (HP.IsArg()) {
		try {
			iMaxIters = HP.GetInt(1, HighParser::range_ge<integer>(1));

		} catch (HighParser::ErrValueOutOfRange<integer> e) {
			silent_cerr("error: invalid number of maximum iterations " << e.Get() << " (must be greater than zero) [" << e.what() << "] at line " << HP.GetLineData() << std::endl);
			throw e;
		}
	}

	// RETURN THE new FACTOR
	return new Factor(s, dReductionFactor,
		iStepsBeforeReduction,
		dRaiseFactor,
		iStepsBeforeRaise,
		iMinIters,
		iMaxIters);
}

void
InitTimeStepData(void)
{
	SetTimeStepData("no" "change", new NoChangeTSR);
	SetTimeStepData("change", new ChangeStepTSR);
	SetTimeStepData("factor", new FactorTSR);
}

void
DestroyTimeStepData(void)
{
	for (TimeStepFuncMapType::iterator it = TimeStepReadMap.begin(); it != TimeStepReadMap.end(); ++it) {
		delete it->second;
	}
	TimeStepReadMap.clear();
}

TimeStepControl * 
ReadTimeStepData(Solver *s, MBDynParser& HP)
{
	const char *str = HP.IsWord(TimeStepWordSet);
	TimeStepFuncMapType::iterator func = TimeStepReadMap.find(std::string(str));
	if (func == TimeStepReadMap.end()) {
		silent_cerr("unknown TimeStepControl type " << s << "; using \"no change\"" << std::endl);
		func = TimeStepReadMap.find(std::string("no" "change"));		
	}
	ASSERT(func != TimeStepReadMap.end());
	return func->second->Read(s, HP);
}

