/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
	Copyright (C) 2016(-2016) all rights reserved.
	The copyright of this patch is transferred
	to Pierangelo Masarati and Paolo Mantegazza
	for use in the software MBDyn as described 
	in the GNU Public License version 2.1

*/
#include "mbconfig.h" 
#include "mbdynFMI.h"
#include "stepsol.h"

#include <sstream>
#include <iostream>

#include <fmilib.h>

#include <FMI2/fmi2_import.h>
//#include <FMI2/fmi2_xml_model_description.h>

#include <FMI/fmi_version.h>
#include <FMI2/fmi2_types.h>
#include <FMI2/fmi2_functions.h>
#include <FMI2/fmi2_enums.h>
//#include <FMI2/fmi2_capi.h>
#include <FMI/fmi_util.h>
#include "fmi2_import_impl.h"

//#include <FMI2/fmi2_capi_impl.h>
//#include <FMI2/fmi2_capi.h>

#include <FMI1/fmi1_import.h>
//#include <FMI1/fmi1_xml_model_description.h>

//#include "../FMI/fmi_import_context_impl.h"
//#include <FMI1/fmi1_capi_impl.h>

#include <FMI/fmi_version.h>
#include <FMI1/fmi1_types.h>
#include <FMI1/fmi1_functions.h>
#include <FMI1/fmi1_enums.h>
//#include <FMI1/fmi1_capi.h>
#include <FMI/fmi_util.h>
#include "fmi1_import_impl.h"

#define STATUSCHECK(a)\
	if (a!=0) {\
		printf("Error in %s, %d\n", __func__, __LINE__); \
		throw ErrGeneric(MBDYN_EXCEPT_ARGS); \
	}

void do_event_iteration(fmi2_import_t *fmu, fmi2_event_info_t *eventInfo)
{
    eventInfo->newDiscreteStatesNeeded = fmi2_true;
    eventInfo->terminateSimulation     = fmi2_false;
    while (eventInfo->newDiscreteStatesNeeded && !eventInfo->terminateSimulation) {
        fmi2_import_new_discrete_states(fmu, eventInfo);
    }
}

static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status) {
	silent_cout("stepFinished is called with fmiStatus = "<< fmi2_status_to_string(status)<<std::endl);
}


int annotation_start_handle(void *context, const char *parentName, void *parent, const char *elm, const char **attr) {
        int i = 0;
        printf("Annotation element %s start (tool: %s, parent:%s)\n", elm, parentName,
                parent?fmi2_import_get_variable_name((fmi2_import_variable_t*)parent):"model");
        while(attr[i]) {
                printf("Attribute %s = %s\n", attr[i], attr[i+1]);
                i+=2;
        }
        printf("Annotation data:\n");
        return 0;
}

int annotation_data_handle(void* context, const char *s, int len) {
        int i;
        for(i = 0; i < len; i++){
                printf("%c", s[i]);
	}
        return 0;
}

int annotation_end_handle(void *context, const char *elm) {
        printf("\nAnnotation element %s end\n", elm);
        return 0;
}


fmi2_xml_callbacks_t annotation_callbacks = {
        annotation_start_handle,
        annotation_data_handle,
        annotation_end_handle, NULL};



void fmu2::parseXML(fmi_import_context_t* context, const char* dirPath){

	fmu = fmi2_import_parse_xml(context, dirPath, &annotation_callbacks);
	if (!fmu){
		silent_cout("Error parsing XML 2.0\n");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	} else { 
		silent_cout("XML Parsed successfully\n");
	}
}

void fmu2::setCallBackFunction(){
	callBackFunctions.logger = fmi2_log_forwarding;
	callBackFunctions.allocateMemory = calloc;
	callBackFunctions.freeMemory = free;
	callBackFunctions.stepFinished = stepFinished;
	callBackFunctions.componentEnvironment = 0;
}

void fmu2::ImportCreateDLL(){

        if(simType != fmi2_import_get_fmu_kind(fmu)-1){
                silent_cout("This FMU does not support the specified simulation type");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }


	if(simType == IMPORT){
		jmstatus = fmi2_import_create_dllfmu(fmu, fmi2_fmu_kind_me, &callBackFunctions);
	}
	else if (simType == COSIM){
		jmstatus = fmi2_import_create_dllfmu(fmu, fmi2_fmu_kind_cs, &callBackFunctions);
	}
	else{
		std::cout<<"Unknown Type";
	}

	STATUSCHECK(jmstatus);
}

int fmu2::GetNumOfContinousStates(void){
	numOfContStates = fmi2_import_get_number_of_continuous_states(fmu);
	
	return numOfContStates;
}

int fmu2::GetNumOfEventIndicators(void){
	return fmi2_import_get_number_of_event_indicators(fmu);
}

int fmu2::GetNumOfVar(void){
	fmi2_import_variable_list_t* vl;
	vl = fmi2_import_get_variable_list(fmu, 0);
	return fmi2_import_get_variable_list_size(vl);
}

void fmu2::Initialize(double dTol, double time, double rTol){
	jmstatus = fmi2_import_instantiate(fmu, "Test ME model instance",fmi2_model_exchange,0,0);
	STATUSCHECK(jmstatus);

	if (rTol ==0){
		rTol = 0.001;
	}

	fmistatus = fmi2_import_setup_experiment(fmu, dTol, rTol, time, fmi2_false, 0.0);
	STATUSCHECK(fmistatus);
	fmistatus = fmi2_import_enter_initialization_mode(fmu);
	STATUSCHECK(fmistatus);
	fmistatus = fmi2_import_exit_initialization_mode(fmu);
	STATUSCHECK(fmistatus);
	
	do_event_iteration(fmu, &eventInfo);
	fmi2_import_enter_continuous_time_mode(fmu);
	fmistatus = fmi2_import_set_time(fmu, time);
	STATUSCHECK(fmistatus);

	eventInfo.newDiscreteStatesNeeded           = fmi2_false;
        eventInfo.terminateSimulation               = fmi2_false;
        eventInfo.nominalsOfContinuousStatesChanged = fmi2_false;
        eventInfo.valuesOfContinuousStatesChanged   = fmi2_true;
        eventInfo.nextEventTimeDefined              = fmi2_false;
        eventInfo.nextEventTime                     = -0.0;

}

bool fmu2::SupportsDirectionalDerivatives(){
	unsigned int capability = 0;

	if(simType == COSIM){
		fmi2_capabilities_enu_t id = fmi2_cs_providesDirectionalDerivatives;
		capability = fmi2_xml_get_capability(fmu->md, id);
	} else if (simType == IMPORT){
		fmi2_capabilities_enu_t id = fmi2_me_providesDirectionalDerivatives;
		capability = fmi2_xml_get_capability(fmu->md, id);
	}

	if (!capability){
		silent_cout("This FMU does not support Directional Derivatives\n");
		return false;
	}
	
	silent_cout("This FMU supports Directional Derivatives\n");
	return true;


}

void fmu2::GetDirectionalDerivatives(FullMatrixHandler* jacobian, int* inputStatesVRef, int inputLength, double* seedVector){
//Seed Vector (the state value)

	fmi2_import_variable_list_t* derivativeList = fmi2_import_get_derivatives_list(fmu);
	size_t derivativeListSize = fmi2_import_get_variable_list_size(derivativeList);

	fmi2_import_variable_t *var;
	fmi2_import_real_variable_t *varCont;
	fmi2_value_reference_t* derList = new fmi2_value_reference_t[derivativeListSize];
	fmi2_value_reference_t* stateList = new fmi2_value_reference_t[derivativeListSize + inputLength];

	for (unsigned int i=0; i<derivativeListSize; i++){
		var          = fmi2_import_get_variable(derivativeList, i);
		derList[i]   = fmi2_import_get_variable_vr(var);
		varCont      = fmi2_import_get_real_variable_derivative_of((fmi2_import_real_variable_t*)var);
		stateList[i] = fmi2_import_get_variable_vr((fmi2_import_variable_t*)varCont);
	}

	for (int i=numOfContStates; i<inputLength+numOfContStates; i++){
		stateList[i] = inputStatesVRef[i-numOfContStates];
	}

	fmi2_real_t output[numOfContStates + inputLength];	

	for (int i=0; i<numOfContStates; i++){
		fmistatus = fmi2_import_get_directional_derivative(fmu, stateList, numOfContStates + inputLength, &derList[i], 1, output , seedVector);
		STATUSCHECK(fmistatus);

		for (int j=0; j<numOfContStates + inputLength; j++){
			jacobian->PutCoef(i+1, j+1, output[j]);
		}
	}
	delete[] derList;
	delete[] stateList;
}



void fmu2::EventIndicatorInit(void){
	nEventIndicators = fmi2_import_get_number_of_event_indicators(fmu);

	eventIndicators = new fmi2_real_t [nEventIndicators];
	eventIndicatorsPrev = new fmi2_real_t [nEventIndicators];

	fmistatus = fmi2_import_get_event_indicators(fmu, eventIndicatorsPrev, nEventIndicators);

	if(fmistatus){
		silent_cout("This FMU does triggers events. Warning: Event triggers not supported by MBDyn\n");
	}	
}

void fmu2::SetRelativeTol(double dTol){
	relativeTolerance = (fmi1_real_t)dTol;
}

void fmu2::SetTime(double time){
        currTime = time;
        fmistatus = fmi2_import_set_time(fmu, currTime);
	STATUSCHECK(fmistatus);	
}

void fmu2::SetStates(double* states){
	callEventUpdate = fmi2_false;
	terminateSimulation = fmi2_false;
	fmistatus = fmi2_import_set_continuous_states(fmu, states, numOfContStates);
	STATUSCHECK(fmistatus);

	fmistatus = fmi2_import_completed_integrator_step(fmu, fmi2_true, &callEventUpdate, &terminateSimulation);
	STATUSCHECK(fmistatus);

}

bool fmu2::CheckInput(const std::string s){
	v = fmi2_import_get_variable_by_name(fmu, s.c_str());
	if(v==NULL){
		silent_cout(s.c_str()<<" is not defined in XML, hence ");
		return false;
	}
	return ((fmi2_import_get_causality(v))==2);
}

int fmu2::GetRefValueFromString(const char* s){
	v = fmi2_import_get_variable_by_name(fmu, s);
        return  static_cast<int>(fmi2_import_get_variable_vr(v));

}

double fmu2::GetStateFromRefValue(unsigned int i){
	double realValue[1];
        fmi1_value_reference_t valueReference = {i};
        fmistatus = fmi2_import_get_real(fmu, &valueReference, 1, realValue);
	STATUSCHECK(fmistatus);

        return realValue[0];

}

fmu2::~fmu2(void){
        
	fmi_import_free_context(context);

	if (simType == IMPORT ){
		delete[] eventIndicators;
	        delete[] eventIndicatorsPrev;

		fmi2_import_free_instance(fmu);
		fmi2_import_destroy_dllfmu(fmu);
		fmi2_import_free(fmu);

	} else {

		fmistatus = fmi2_import_terminate(fmu);
		STATUSCHECK(fmistatus);
		fmi2_import_free_instance(fmu);
		fmi2_import_destroy_dllfmu(fmu);
		fmi2_import_free(fmu);

	}

}


void fmu2::GetStateDerivatives(double* derivatives){
    
        fmistatus = fmi2_import_get_derivatives(fmu, derivatives, numOfContStates);

	STATUSCHECK(fmistatus);

}

void fmu2::GetStates(double* states){
        fmistatus = fmi2_import_get_continuous_states(fmu, states, numOfContStates);
	STATUSCHECK(fmistatus);
}

bool fmu2::CheckInterrupts(double currTime, double* states){
        int zeroCrossningEvent = 0;

        fmistatus = fmi2_import_get_event_indicators(fmu, eventIndicators, nEventIndicators);
	STATUSCHECK(fmistatus);
        for (unsigned int k = 0; k < nEventIndicators; k++) {
              if (eventIndicators[k]*eventIndicatorsPrev[k] < 0) {
                      zeroCrossningEvent = 1;
                      break;
              }
        }

	if ( zeroCrossningEvent || (eventInfo.nextEventTimeDefined && currTime == eventInfo.nextEventTime)) {
		fmistatus = fmi2_import_enter_event_mode(fmu);
		do_event_iteration(fmu, &eventInfo);
		fmistatus = fmi2_import_enter_continuous_time_mode(fmu);
		fmistatus = fmi2_import_get_continuous_states(fmu, states, numOfContStates);
		fmistatus = fmi2_import_get_event_indicators(fmu, eventIndicatorsPrev, nEventIndicators);


		return true;
        } else {
                return false;
        }
}

void fmu2::SetValuesByVariable(const std::string s, double value){
        fmi2_value_reference_t ref;
        fmi2_variable_alias_kind_enu_t aliasKind;
        fmi2_base_type_enu_t baseType;

        fmi2_import_variable_t* v = fmi2_import_get_variable_by_name(fmu, s.c_str());

        ref = fmi2_import_get_variable_vr(v);
        aliasKind = fmi2_import_get_variable_alias_kind(v);
        baseType = fmi2_import_get_variable_base_type(v);

        if (aliasKind == fmi2_variable_is_not_alias){
                if(baseType != fmi2_base_type_bool){
                        value = -value;
		}
                else{
                        value = (bool)!value;
		}
        }

        if (baseType == fmi2_base_type_real){
                fmistatus = fmi2_import_set_real(fmu, &ref, 1, &value);
	}
        else if (baseType == fmi2_base_type_int){
                fmistatus = fmi2_import_set_integer(fmu, &ref, 1, (const int*)&value);
	}
        else if (baseType == fmi2_base_type_bool){
                fmistatus = fmi2_import_set_boolean(fmu, &ref, 1, (const fmi2_boolean_t*)&value);
	}
        else{
                silent_cerr("Input type not Supported");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        STATUSCHECK(fmistatus);
}


void fmu2::InitializeAsSlave(const char* location, double tstart, double tend){

	std::stringstream resourceLocation;
	std::string strLocation = location;

	resourceLocation << "file://" << strLocation;
	fmi2_boolean_t visible = fmi2_false;
	jmstatus = fmi2_import_instantiate(fmu, "Model for CS", fmi2_cosimulation, resourceLocation.str().c_str(), visible);
	STATUSCHECK(jmstatus);

	if (tend == std::numeric_limits<double>::max()){
		fmistatus = fmi2_import_setup_experiment(fmu, fmi2_true, relativeTolerance, tstart, fmi2_false, 0);
	} else {
		fmistatus = fmi2_import_setup_experiment(fmu, fmi2_true, relativeTolerance, tstart, fmi2_true, tend);
	}

	STATUSCHECK(fmistatus);

        fmistatus = fmi2_import_enter_initialization_mode(fmu);
	STATUSCHECK(fmistatus);

        fmistatus = fmi2_import_exit_initialization_mode(fmu);
	STATUSCHECK(fmistatus);
	
}

void fmu2::CSPropogate(double tcur, double dt){
	fmistatus = fmi2_import_do_step(fmu, tcur, dt, fmi2_true);
	STATUSCHECK(fmistatus);
}


fmu::~fmu(void){
	NO_OP;
}

fmu1::~fmu1(void){
	fmi_xml_free_context(context);

	if(simType == IMPORT){

		fmi1_capi_free_dll(fmu->capi);

		delete[] eventIndicators;
		delete[] eventIndicatorsPrev;
		//	delete[] currStates;

		delete[] deriv;
		delete[] vrs;

		fmi1_import_free_model_instance(fmu);
		fmi1_import_destroy_dllfmu(fmu);
		fmi1_import_free(fmu);

	} else {

		fmistatus = fmi1_import_terminate_slave(fmu);
		STATUSCHECK(fmistatus);
		fmi1_import_free_slave_instance(fmu);

	}
}


void fmu1::parseXML(fmi_import_context_t* context, const char* dirPath){
	fmu = fmi1_import_parse_xml(context, dirPath);
        if (!fmu){
                silent_cout("Error parsing XML\n");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        } else {
                silent_cout("XML Parsed successfully\n");
        }
}

void fmu1::setCallBackFunction(){
	callBackFunctions.logger = fmi1_log_forwarding;
	callBackFunctions.allocateMemory = calloc;
	callBackFunctions.freeMemory = free;
	silent_cout("Callback Functions Set.\n");
}

void fmu1::ImportCreateDLL(void){
	jmstatus = fmi1_import_create_dllfmu(fmu, callBackFunctions, 1);
	STATUSCHECK(jmstatus);
	
	silent_cout("Successfully created the DLL loading mechanism.\n");
}

int fmu1::GetNumOfContinousStates(void){
	numOfContStates = fmi1_import_get_number_of_continuous_states(fmu);
//	currStates      = new double [numOfContStates];
//	currStatesDer   = new double [numOfContStates];	
	vrs             = new fmi1_value_reference_t [numOfContStates];
	deriv           = new fmi1_real_t [numOfContStates];
	return numOfContStates;
}

void fmu1::Initialize(double dTol, double time, double rTol){

	fmi1_boolean_t tolControlled = fmi1_true;

	if (dTol ==0 ){
		dTol = fmi1_import_get_default_experiment_tolerance(fmu);
		tolControlled = fmi1_false;
	}
	currTime  = time;

	STATUSCHECK(fmi1_status_ok);

	jmstatus  = fmi1_import_instantiate_model(fmu, "Test ME model instance");
	STATUSCHECK(jmstatus);

	fmistatus = fmi1_import_initialize(fmu, tolControlled, dTol, &eventInfo);
	STATUSCHECK(fmistatus);

	fmistatus = fmi1_import_set_time(fmu, currTime);
	STATUSCHECK(fmistatus);

// Necessary to initialize model as capi
	fmi1_capi_load_dll(fmu->capi);

	jmstatus = fmi1_capi_load_fcn(fmu->capi);	
	STATUSCHECK(jmstatus);

#ifdef DEBUG
	fmi1_boolean_t loggingOn = fmi1_true;
#else
	fmi1_boolean_t loggingOn = fmi1_false;
#endif

	fmi1_capi_instantiate_model(fmu->capi, "noName", fmi1_import_get_GUID(fmu), loggingOn); 

	intermediateResults = fmi1_false;
	silent_cout("FMU Initialized\n");
}

int fmu1::GetNumOfEventIndicators(void){
	return fmi1_import_get_number_of_event_indicators(fmu);
}

bool fmu1::CheckInterrupts(double currTime, double* states){
	int zeroCrossningEvent = 0;

	fmistatus = fmi1_import_get_event_indicators(fmu, eventIndicators, nEventIndicators);
	STATUSCHECK(fmistatus);
	for (unsigned int k = 0; k < nEventIndicators; k++) {
  	      if (eventIndicators[k]*eventIndicatorsPrev[k] < 0) {
        	      zeroCrossningEvent = 1;
                      break;
              }
       	}

	if ( zeroCrossningEvent || (eventInfo.upcomingTimeEvent && currTime == eventInfo.nextEventTime)) {
		fmistatus = fmi1_import_eventUpdate(fmu, intermediateResults, &eventInfo);
		fmistatus = fmi1_import_get_continuous_states(fmu, states, numOfContStates);
		fmistatus = fmi1_import_get_event_indicators(fmu, eventIndicatorsPrev, nEventIndicators);
		return true;
	} else {
		return false;
	}	
}

bool fmu1::CheckInput(const std::string s){
	v = fmi1_import_get_variable_by_name(fmu, s.c_str());
	if(v==NULL){
		silent_cout(s.c_str()<<" is not defined in XML, hence ");
		return false;
	}

	return ((fmi1_import_get_causality(v))==0);
}



void fmu1::SetValuesByVariable(const std::string s, double value){
	fmi1_value_reference_t ref;
	fmi1_variable_alias_kind_enu_t aliasKind;
        fmi1_base_type_enu_t baseType;

	fmi1_import_variable_t* v = fmi1_import_get_variable_by_name(fmu, s.c_str());

	ref = fmi1_import_get_variable_vr(v);
	aliasKind = fmi1_import_get_variable_alias_kind(v);
	baseType = fmi1_import_get_variable_base_type(v);

        if (aliasKind == fmi1_variable_is_negated_alias){
		if(baseType != fmi1_base_type_bool){
			value = -value;
		}
		else {
			value = (bool)!value; 
		}
	}

	if (baseType == fmi1_base_type_real){
		fmistatus = fmi1_import_set_real(fmu, &ref, 1, &value);
	}
	else if (baseType == fmi1_base_type_int){
		fmistatus = fmi1_import_set_integer(fmu, &ref, 1, (const int*)&value);
	}
	else if (baseType == fmi1_base_type_bool){
		fmistatus = fmi1_import_set_boolean(fmu, &ref, 1, (const fmi1_boolean_t*)&value);
	}
	else{
		silent_cerr("Input type not Supported");
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	STATUSCHECK(fmistatus);
}

int fmu1::GetRefValueFromString(const char* s){
	v = fmi1_import_get_variable_by_name(fmu, s);
	return	static_cast<int>(fmi1_import_get_variable_vr(v));
}

void fmu1::GetStateDerivatives(double* derivatives){

	fmistatus = fmi1_import_get_derivatives(fmu, derivatives, numOfContStates);
	STATUSCHECK(fmistatus);
}

void fmu1::SetTime(double time){
	currTime = time;
	fmistatus = fmi1_import_set_time(fmu, currTime);
	STATUSCHECK(fmistatus);
}

void fmu1::SetStates(double* states){

	fmistatus = fmi1_capi_get_state_value_references(fmu->capi, vrs, numOfContStates);
	STATUSCHECK(fmistatus);
	fmistatus = fmi1_import_set_real(fmu, vrs, numOfContStates, states);
	STATUSCHECK(fmistatus);

}

double fmu1::GetStateFromRefValue(unsigned int i){
	double realValue[1];
	fmi1_value_reference_t valueReference = {i};
	fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, realValue);	
	STATUSCHECK(fmistatus);
	
	return realValue[0];
} 

bool fmu1::SupportsDirectionalDerivatives(){
	return false;
}

void fmu1::GetDirectionalDerivatives(FullMatrixHandler* jacobian, int* vector, int vectorLength, double* seedVector){
	NO_OP;
}

int fmu1::GetNumOfVar(){
	fmi1_import_variable_list_t* vl;
	vl = fmi1_import_get_variable_list(fmu);
	return fmi1_import_get_variable_list_size(vl);
}

void fmu1::SetRelativeTol(double dTol){
	relativeTolerance = (fmi1_real_t)dTol;
}

void fmu1::EventIndicatorInit(void){
	nEventIndicators = fmi1_import_get_number_of_event_indicators(fmu);

	eventIndicators = new fmi1_real_t [nEventIndicators];
        eventIndicatorsPrev = new fmi1_real_t [nEventIndicators];

	fmistatus = fmi1_import_get_event_indicators(fmu, eventIndicatorsPrev, nEventIndicators);	
	if(fmistatus){	
		silent_cout("This FMU does triggers events. Warning: Event triggers not supported by MBDyn\n");
	}
}

void fmu1::CSPropogate(double tcur, double dt){
	fmistatus = fmi1_import_do_step(fmu, tcur, dt, fmi1_true);
	STATUSCHECK(fmistatus);
}


void fmu1::GetStates(double* states){

	fmistatus = fmi1_import_get_continuous_states(fmu, states, numOfContStates);
	if(!fmistatus){	
		fmistatus = fmi1_import_get_nominal_continuous_states(fmu, states, numOfContStates);
	}
	if(!fmistatus){
		for(int i=0; i<numOfContStates; i++){
			states[i] = 0;
		}
	}
	
//	STATUSCHECK(fmistatus);
	
}

void fmu1::InitializeAsSlave(const char* location, double tstart, double tend){
	std::stringstream resourceLocation;
	std::string strLocation = location;

	resourceLocation << "file://" << strLocation;
	jmstatus = fmi1_import_instantiate_slave(fmu, "Test CS model instance", resourceLocation.str().c_str(), "", 0.0, fmi1_false, fmi1_false);
	STATUSCHECK(jmstatus);
// INITIALIZE PRINTS SOME WARNINGS.
	fmistatus = fmi1_import_initialize_slave(fmu, 1.0, fmi1_false, 10.0);
	STATUSCHECK(fmistatus);

}
