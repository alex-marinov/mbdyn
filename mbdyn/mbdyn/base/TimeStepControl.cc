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

typedef std::map<std::string,TimeStepRead * > TimeStepFuncMapType; 

TimeStepFuncMapType TimeStepReadMap;

struct TimeStepWordSetType : public HighParser::WordSet{

	bool IsWord(const std::string& s) const{
		return TimeStepReadMap.find(s) != TimeStepReadMap.end();
	};
};

TimeStepWordSetType TimeStepWordSet;

bool 
SetTimeStepData(const char *name , TimeStepRead * sr){
	return TimeStepReadMap.insert(std::make_pair(name,sr)).second;
}


NoChange::NoChange(Solver * s): s(s){
}

doublereal
NoChange::dGetNewStepTime(StepIntegrator::StepChange currStep , doublereal iPerformedIters) {
	
	return dCurrTimeStep;
}

void
NoChange::Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep){
	
	this->dCurrTimeStep = dCurrTimeStep;
	doublereal dMinTimeStep = s->dGetdMinTimeStep();
	doublereal dDefaultMinTimeStep = -1.;
	if (dMinTimeStep != dDefaultMinTimeStep) {
			silent_cerr("\"min time step\" only allowed with variable time step (ignored)." <<std::endl);
	}
	dMinTimeStep = dInitialTimeStep;

	doublereal dInitialMaxTimeStep ;
	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();
	}
	else{
		dInitialMaxTimeStep =  MaxTimeStep.dGet();
	}

	if (dInitialMaxTimeStep != std::numeric_limits<doublereal>::max()) {
		silent_cerr("\"max time step\" only allowed with variable time step (ignored)." <<std::endl);
	}
	s->dSetMinTimeStep(dMinTimeStep);

}

struct NoChangeSTR : public TimeStepRead {
	TimeStepControl *
	Read(Solver * s, MBDynParser& HP);
};

TimeStepControl *
NoChangeSTR::Read(Solver * s, MBDynParser& HP){
	return new NoChange(s);
}

ChangeStep::ChangeStep(Solver *s , DriveCaller* pStrategyChangeDrive):
s(s),
pStrategyChangeDrive(pStrategyChangeDrive)
{}

doublereal
ChangeStep::dGetNewStepTime(StepIntegrator::StepChange currStep , doublereal iPerformedIters){
	doublereal dNewStep = pStrategyChangeDrive->dGet();
	doublereal maxTimeStepval = MaxTimeStep.dGet();;
	
	doublereal dMinTimeStep = s->dGetdMinTimeStep();

	if (currStep == StepIntegrator::REPEATSTEP
		&& dNewStep == dCurrTimeStep){
		return dMinTimeStep/2.;
	}
	return std::max(std::min(dNewStep, maxTimeStepval), dMinTimeStep);
}

void
ChangeStep::SetDriveHandler(const DriveHandler* driveHandler){
	pStrategyChangeDrive->SetDrvHdl(driveHandler);
}


void
ChangeStep::Init(doublereal iMaxIterations,doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep){
	this->dCurrTimeStep = dCurrTimeStep;
	this->MaxTimeStep = MaxTimeStep;
	doublereal dMinTimeStep = s->dGetdMinTimeStep();
	doublereal dInitialMaxTimeStep ;
	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();
	}
	else{
		dInitialMaxTimeStep =  MaxTimeStep.dGet();
	}
	if (dMinTimeStep > dInitialMaxTimeStep) {
		silent_cerr("inconsistent min/max time step"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

class ChangeStepSTR : public TimeStepRead{

public:	
	TimeStepControl *
	Read(Solver * s, MBDynParser& HP);	
};

TimeStepControl *
ChangeStepSTR::Read(Solver * s, MBDynParser& HP){
	return new ChangeStep(s , HP.GetDriveCaller(true));
}

Factor::Factor(Solver * s ,doublereal dReductionFactor,
			doublereal iStepsBeforeReduction,
			doublereal dRaiseFactor,
			doublereal iStepsBeforeRaise,
			doublereal iMinIters,
			doublereal iMaxIters ):
s(s),
dReductionFactor(dReductionFactor),
iStepsBeforeReduction(iStepsBeforeReduction),
dRaiseFactor(dRaiseFactor),
iStepsBeforeRaise(iStepsBeforeRaise),
iMinIters(iMinIters),
iMaxIters(iMaxIters),
bLastChance(false),
iStepsAfterReduction(0),
iStepsAfterRaise(0),
iWeightedPerformedIters(0)
{}

doublereal
Factor::dGetNewStepTime(StepIntegrator::StepChange Why , doublereal iPerformedIters){
	doublereal dMinTimeStep = s->dGetdMinTimeStep();
	doublereal maxTimeStepval = MaxTimeStep.dGet();
	
	if (Why == StepIntegrator::REPEATSTEP) {
		if (dCurrTimeStep*dReductionFactor > dMinTimeStep){
			if (bLastChance == true) {
				bLastChance = false;
			}
			iStepsAfterReduction = 0;
			return std::min(dCurrTimeStep*dReductionFactor, maxTimeStepval);

		} else {
			if (bLastChance == false) {
				bLastChance = true;
				return dMinTimeStep;
			} else {
					return std::min(dCurrTimeStep*dReductionFactor, maxTimeStepval);
				}
		}
	}

	if (Why == StepIntegrator::NEWSTEP) {
		iStepsAfterReduction++;
		iStepsAfterRaise++;

		iWeightedPerformedIters = (10*iPerformedIters + 9*iWeightedPerformedIters)/10;

		if (iPerformedIters > iMaxIters) {
			iStepsAfterReduction = 0;
			bLastChance = false;
			return std::min(std::max(dCurrTimeStep*dReductionFactor, dMinTimeStep), maxTimeStepval);

		} else if (iPerformedIters <= iMinIters
				&& iStepsAfterReduction >iStepsBeforeReduction
				&& iStepsAfterRaise > iStepsBeforeRaise
				&& dCurrTimeStep < maxTimeStepval)
			{
				iStepsAfterRaise = 0;
				iWeightedPerformedIters = 0;
				return std::min(dCurrTimeStep*dRaiseFactor, maxTimeStepval);
			}
		return dCurrTimeStep;
	}
	//Should Not Reach Over here
	return 0.;
}


void
Factor::Init(doublereal iMaxIterations , doublereal dCurrTimeStep ,DriveOwner & MaxTimeStep , doublereal dInitialTimeStep){

	this->dCurrTimeStep = dCurrTimeStep;
	this->MaxTimeStep = MaxTimeStep;
	doublereal dDefaultMinTimeStep = -1.;
	doublereal dInitialMaxTimeStep ;
	doublereal dMinTimeStep = s->dGetdMinTimeStep();

	if(iMaxIters <= iMinIters) {
		silent_cerr("warning, "
			<< "strategy maximum number "
			<< "of iterations "
			<< "is <= minimum: "
			<< iMaxIters << " <= "
			<< iMinIters << "; "
			<< "the maximum global iteration value "
			<< iMaxIterations << " "
			<< "will be used"
			<< std::endl);
		iMaxIters = iMaxIterations;
	}

	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(PostponedDriveCaller)) {
		dInitialMaxTimeStep = std::numeric_limits<doublereal>::max();
	}
	else{
		dInitialMaxTimeStep =  MaxTimeStep.dGet();
	}

	if (typeid(*MaxTimeStep.pGetDriveCaller()) == typeid(ConstDriveCaller)
			&& MaxTimeStep.dGet() == std::numeric_limits<doublereal>::max()) {
			silent_cerr("warning, "
				<< "maximum time step not set and strategy "
				<< "factor selected:\n"
				<< "the initial time step value will be used as "
				<< "maximum time step: "
				<< "max time step = "
				<< dInitialTimeStep
				<< std::endl);
	}
	if (dMinTimeStep == dDefaultMinTimeStep) {
		silent_cerr("warning, "
			<< "minimum time step not set and strategy "
			<< "factor selected:\n"
			<< "the initial time step value will be used as "
			<< "minimum time step: "
			<< "min time step = "
			<< dInitialTimeStep
			<< std::endl);
		dMinTimeStep = dInitialTimeStep;
	}
	if (dMinTimeStep == dInitialMaxTimeStep) {
		silent_cerr("error, "
			<< "minimum and maximum time step are equal, but "
			<< "strategy factor has been selected:\n"
			<< "this is almost meaningless"
			<< std::endl);
	}
	if (dMinTimeStep > dInitialMaxTimeStep) {
		silent_cerr("error, "
			<< "minimum and maximum time step are equal, but "
			<< "strategy factor has been selected:\n"
			<< "this is meaningless - bailing out"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	s->dSetMinTimeStep(dMinTimeStep);
}

class FactorSTR: public TimeStepRead {

public:
	TimeStepControl *
	Read(Solver * s, MBDynParser& HP);	
};

TimeStepControl *
FactorSTR::Read( Solver * s, MBDynParser& HP){

	integer iMinIters = 1;
	integer iMaxIters = 0;
	doublereal dReductionFactor = HP.GetReal();
	if (dReductionFactor >= 1.) {
		silent_cerr("warning, "
			"illegal reduction factor "
			"at line " << HP.GetLineData()
			<< "; default value 1. "
			"(no reduction) will be used"
			<< std::endl);
		dReductionFactor = 1.;
	}

	doublereal iStepsBeforeReduction = HP.GetInt();
	if (iStepsBeforeReduction <= 0) {
		silent_cerr("warning, "
			"illegal number of steps "
			"before reduction at line "
			<< HP.GetLineData()
			<< "; default value 1 will be "
			"used (it may be dangerous)"
			<< std::endl);
		iStepsBeforeReduction = 1;
	}

	doublereal dRaiseFactor = HP.GetReal();
	if (dRaiseFactor <= 1.) {
		silent_cerr("warning, "
			"illegal raise factor at line "
			<< HP.GetLineData()
			<< "; default value 1. "
			"(no raise) will be used"
			<< std::endl);
		dRaiseFactor = 1.;
	}

	doublereal iStepsBeforeRaise = HP.GetInt();
	if (iStepsBeforeRaise <= 0) {
		silent_cerr("warning, "
			"illegal number of steps "
			"before raise at line "
			<< HP.GetLineData()
			<< "; default value 1 will be "
			"used (it may be dangerous)"
			<< std::endl);
		iStepsBeforeRaise = 1;
	}

	iMinIters = HP.GetInt();
	if (iMinIters <= 0) {
		silent_cerr("warning, "
			"illegal minimum number "
			"of iterations at line "
			<< HP.GetLineData()
			<< "; default value 0 will be "
			"used (never raise)"
			<< std::endl);
		iMinIters = 1;
	}
	iMaxIters = 0;
	if (HP.IsArg()) {
		iMaxIters = HP.GetInt();
		if (iMaxIters <= 0) {
			silent_cerr("warning, "
				"illegal mmaximim number "
				"of iterations at line "
				<< HP.GetLineData()
				<< "; default value will be "
				"used"
				<< std::endl);
		}
	}

	// RETURN THE new FACTOR
	return new Factor(s,dReductionFactor,
						iStepsBeforeReduction,
						dRaiseFactor,
						iStepsBeforeRaise,
						iMinIters,
						iMaxIters);
	
}

void
InitTimeStepData(void){
	SetTimeStepData("nochange", new NoChangeSTR);
	SetTimeStepData("change", new ChangeStepSTR);
	SetTimeStepData("factor", new FactorSTR);
}

void
DestroyTimeStepData(void){
	for(TimeStepFuncMapType::iterator it = TimeStepReadMap.begin(); it != TimeStepReadMap.end(); it++){
		delete it->second;
	}
	TimeStepReadMap.clear();
}


TimeStepControl * 
ReadTimeStepData( Solver * s, MBDynParser& HP){
	const char *str = HP.IsWord(TimeStepWordSet);
	TimeStepFuncMapType::iterator func = TimeStepReadMap.find(std::string(str));
	if(func == TimeStepReadMap.end()){
		silent_cerr("unknown Time step type " << s <<std::endl);
		func = TimeStepReadMap.find(std::string("nochange"));		
	}
	return func->second->Read(s , HP );
}

