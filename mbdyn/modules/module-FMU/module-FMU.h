/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
        Copyright (C) 2016(-2017) all rights reserved.
        The copyright of this patch is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described 
        in the GNU Public License version 2.1

*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "drive.h"
#include "privdrive.h"

#include "mbdynFMI_config.h"



class ModuleFMU
: virtual public Elem, public UserDefinedElem {
private:
       	fmu* model;
        char FMUlocation[1000];
        const char* simType;
	double currTime;
	double initialTime;
	double timeStep;
	double endTime;

        fmi_import_context_t* context;
        fmi_version_enu_t version;

        jm_callbacks callbacks;
        jm_status_enu_t status;
	
        typedef std::map<std::string, const DriveCaller *> strDriveCon;
        typedef std::map<int, const PrivDriveCaller*> intDriveCon;
        strDriveCon drivesContainer;
	intDriveCon privDrivesIndex;
	
	int numOfContinousStates;
	int numOfEventIndicators;
	DataManager* pDM;
	
	int* statesOrder;
	double* currState;
	double* stateDerivatives;
	int* jacobianInputVector;
	bool* privDriveArray;
	int privDriveLength;

	bool directionalFlag;	
	double *seedVector;	

public:
        ModuleFMU(unsigned uLabel, const DofOwner *pDO,
                DataManager* pDM, MBDynParser& HP);
        virtual ~ModuleFMU(void);

        virtual void Output(OutputHandler& OH) const;
        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr);
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr);
        unsigned int iGetNumPrivData(void) const;
        unsigned int iGetPrivDataIdx(const char *s) const; //
        doublereal dGetPrivData(unsigned int i) const; //
        int iGetNumConnectedNodes(void) const;
        void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
        void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                SimulationEntity::Hints *ph);
        std::ostream& Restart(std::ostream& out) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr);
        SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType(unsigned int i) const;

	fmu::SimulationTypes SIMTYPE;

};

