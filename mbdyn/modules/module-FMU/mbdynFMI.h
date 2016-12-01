/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
        Copyright (C) 2016(-2016) all rights reserved.
        The copyright of this patch is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described 
        in the GNU Public License version 2.1

*/


#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <errno.h>
#include <fmilib.h>

#include "fmi_zip_unzip.h"

//#include <JM/jm_types.h>
//#include <JM/jm_callbacks.h>
//#include <JM/jm_portability.h>

#include <FMI2/fmi2_import.h>
#include <FMI/fmi_import_context.h>

#include "mbconfig.h"
#include "stepsol.h"
 

//#include <fmi1_import_impl.h>
//#include <fmi1_import_variable_list_impl.h>

class fmu {
    public:
        enum SimulationTypes{
                IMPORT,
                COSIM
        };	

	int simType; 
	fmi_import_context_t* context;

	fmu(fmi_import_context_t* text, int type){
		context = text;
		simType = type;
	}

        virtual void parseXML(fmi_import_context_t* context, const char* dirPath) = 0;
        virtual void setCallBackFunction() = 0;

	virtual void ImportCreateDLL(void) = 0;

	virtual bool CheckInput(const std::string) = 0;

	virtual int GetNumOfContinousStates(void) = 0;
	virtual int GetNumOfEventIndicators(void) = 0;
	virtual int GetNumOfVar(void) = 0;

	virtual void Initialize(double dTol, double time, double rTol) = 0;
	virtual void EventIndicatorInit(void) = 0;

	virtual void SetRelativeTol(double dTol) = 0; //Not Req
	virtual void SetTime(double time) = 0;
	virtual void SetStates(double* states) = 0;
	virtual void SetValuesByVariable(const std::string, double) = 0;

	virtual int GetRefValueFromString(const char* s) = 0;
	virtual double GetStateFromRefValue(unsigned int i) = 0;
	virtual void GetStateDerivatives(double*) = 0;
	virtual void GetStates(double*) = 0;
	virtual void GetDirectionalDerivatives(FullMatrixHandler *,int*, int, double *) = 0;
	virtual bool SupportsDirectionalDerivatives() = 0;

	virtual bool CheckInterrupts(double, double*) = 0;

	virtual void InitializeAsSlave(const char*, double, double) = 0;
	virtual void CSPropogate(double tcur, double dt) = 0;


	virtual ~fmu(void);
};

class fmu1 :public fmu{
    private:
        fmi1_callback_functions_t callBackFunctions;
        fmi1_import_t* fmu;

	fmi1_status_t fmistatus;	
	jm_status_enu_t jmstatus;	

	fmi1_real_t currTime;
	fmi1_real_t relativeTolerance;
	fmi1_event_info_t eventInfo;

	size_t nEventIndicators;
	fmi1_real_t* eventIndicators;
        fmi1_real_t* eventIndicatorsPrev;
	fmi1_import_variable_t* v;

	int numOfContStates;
	fmi1_real_t* deriv;
        fmi1_value_reference_t *vrs;

	fmi1_boolean_t intermediateResults;


    public:
        void parseXML(fmi_import_context_t* context, const char* dirPath);
        void setCallBackFunction();

	fmu1(fmi_import_context_t* text, int type):fmu::fmu(text, type){
		NO_OP;
	}

	void ImportCreateDLL(void);

	bool CheckInput(const std::string);

	int GetNumOfContinousStates(void);
	int GetNumOfEventIndicators(void);
	int GetNumOfVar(void);

	void Initialize(double dTol, double time, double rTol);
	void EventIndicatorInit(void);

	void SetRelativeTol(double dTol);
	void SetTime(double time);
	void SetStates(double* states);
	void SetValuesByVariable(const std::string, double);

	int GetRefValueFromString(const char* s);
	double GetStateFromRefValue(unsigned int i);
	void GetStateDerivatives(double*);
	void GetStates(double *);
	void GetDirectionalDerivatives(FullMatrixHandler *, int*, int, double*);
	bool SupportsDirectionalDerivatives();

	bool CheckInterrupts(double, double*);

	void InitializeAsSlave(const char* location, double tstart, double tend);
	void CSPropogate(double tcur, double dt);

	virtual ~fmu1(void);

};

class fmu2 : public fmu{
    private:
        fmi2_callback_functions_t callBackFunctions;
        fmi2_import_t* fmu;

	fmi2_status_t fmistatus;	
	jm_status_enu_t jmstatus;

	fmi2_real_t currTime;
	fmi2_real_t relativeTolerance;
	fmi2_event_info_t eventInfo;

	size_t nEventIndicators;
	fmi2_real_t* eventIndicators;
        fmi2_real_t* eventIndicatorsPrev;
	fmi2_import_variable_t* v;
	fmi2_boolean_t callEventUpdate;
	fmi2_boolean_t terminateSimulation;

	int numOfContStates;
        fmi2_value_reference_t *vrs;

	fmi2_boolean_t intermediateResults;


    public:
        void parseXML(fmi_import_context_t* context, const char* dirPath);
        void setCallBackFunction();

	fmu2(fmi_import_context_t* text, int type):fmu::fmu(text, type){
		NO_OP;
	}

	void ImportCreateDLL(void);

	bool CheckInput(const std::string);

	int GetNumOfContinousStates(void);
	int GetNumOfEventIndicators(void);
	int GetNumOfVar(void);

	void Initialize(double dTol, double time, double rTol);
	void EventIndicatorInit(void);

	void SetRelativeTol(double dTol); //Not Req
	void SetTime(double time);
	void SetStates(double* states);
	void SetValuesByVariable(const std::string, double);

	int GetRefValueFromString(const char* s);
	double GetStateFromRefValue(unsigned int i);
	void GetStateDerivatives(double*);
	void GetStates(double *);
	void GetDirectionalDerivatives(FullMatrixHandler *, int*, int, double*);
	bool SupportsDirectionalDerivatives();

	bool CheckInterrupts(double, double*);

	void InitializeAsSlave(const char* , double, double);
	void CSPropogate(double tcur, double dt);

	virtual ~fmu2(void);
};

