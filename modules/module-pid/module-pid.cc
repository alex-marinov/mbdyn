// Author: Matteo Daniele <matteo.daniele@polimi.it>
// pid module for mbdyn
#include "mbconfig.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "module-pid.h"

Pid::Pid(unsigned int uLabel, const DofOwner *pDO,
            DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout("\nModule: PID\n"
        "Author: Matteo Daniele <matteo.daniele@polimi.it>\n"
        "Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
        "		Politecnico di Milano					\n"
        "		http://www.aero.polimi.it				\n"
        "									\n"
        "	All rights reserved						\n"
        "Implementation of a PID controller with anti windup\n"
        "Input:\n"
        "   Reference value\n"
        "   Measured value\n"
        "   Kp: proportional constant\n"
        "   Ki: integral constant\n"
        "   Kd: derivative constant\n"
        "   Kn: derivative washout constant\n"
        "   Ii0: initial value of the integrator (integrator)\n"
        "   Id0: initial value of the integrator (derivative)\n"
        "   Ks: anti windup scale factor\n"
        "   Ymin: saturation min value\n"
        "   Ymax: saturation max value\n"
        "\n"
        "Output:\n"
        "element label\n"
        "input errror\n"
        "output after saturation\n"
        "output before saturation\n"
        "proportional output\n"
        "integral output\n"
        "derivative output\n"
        "Integrator intergral\n"
        "Derivative integral\n"
        "\n"
        "Private data:\n"
        "output after saturation\n"
        "output before saturation\n"
        << std::endl);

        if (!HP.IsArg()){
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }

    // constructor begins:
    // get timestep
    DoTime.Set(new TimeStepDriveCaller(pDM->pGetDrvHdl()));

    // gather input parameters
    DEBUGCOUT( "reading PID(" << uLabel << ") element..." << std::endl);
    //
    //bSetPoint = false;
    if (HP.IsKeyWord("setpoint")){
        pSetPoint = HP.GetDriveCaller();
        bSetPoint = true;
    }
    else {
        //SetPoint = 0.0;
        DEBUGCOUT( "PID(" << uLabel << "): setpoint not provided, using " << SetPoint << " as setpoint value" << std::endl);
    }
    //
    if (HP.IsKeyWord("measure")){
        pMeasure = HP.GetDriveCaller();
    }
    if (HP.IsKeyWord("Kp")){
        Kp = HP.GetReal();
    }
    if (HP.IsKeyWord("Ki")){
        Ki = HP.GetReal();
    }
    if (HP.IsKeyWord("Kd")){
        Kd = HP.GetReal();
    }
    //
    if (HP.IsKeyWord("Kn")){
        Kn = HP.GetReal();
    }
    else{
        //Kn = 100.0;
        DEBUGCOUT( "PID(" << uLabel << "): Kn not provided, using " << Kn << " as default value" << std::endl);
    }
    //
    if (HP.IsKeyWord("Ii0")){
        Ii0 = HP.GetReal();
    }
    else{
        //Ii0 = 0.0;
        DEBUGCOUT( "PID(" << uLabel << "): Ii0 for the integrator not provided, using " << Ii0 << " as initial integrator value" << std::endl);
    }
    //
    if (HP.IsKeyWord("Id0")){
        Id0 = HP.GetReal();
    }
    else {
        //Id0 = 0.0;
        DEBUGCOUT( "PID(" << uLabel << "): Id0 for the derivator not provided, using " << Id0 << " as initial derivator value" << std::endl);
    }
    //
    if (HP.IsKeyWord("Ks")){
        Ks = HP.GetReal();
    }
    else {
        //Ks = 1.0;
        DEBUGCOUT( "PID(" << uLabel << "): anti windup scale factor not provided, using " << Ks << " as value" << std::endl);
    }
    //
    //bGotYmin = false;
    //bGotYmax = false;

    if (HP.IsKeyWord("Ymin")){
        pYmin = HP.GetDriveCaller();
        bGotYmin = true;
    }
    else {
        //Ymin = -2.e32; //-std::numeric_limits<doublereal>::max();
        DEBUGCOUT( "PID(" << uLabel << "): lower limit saturation not provided, using " << Ymin << " as minimum value" << std::endl);
    }

    if (HP.IsKeyWord("Ymax")){
        pYmax = HP.GetDriveCaller();
        bGotYmax = true;
    }
    else {
        //Ymax = 2.e32; //std::numeric_limits<doublereal>::max();
        DEBUGCOUT( "PID(" << uLabel << "): upper limit saturation not provided, using " << Ymax << " as maximum value" << std::endl);
    }

    // INITIALIZATION
    if (bSetPoint){
        SetPoint = pSetPoint->dGet();
    }

    Measure  = pMeasure->dGet();

    // initialize integrals
    Ii = Ii0;
    Id = Id0;

    // initialize outputs
    YpOut = 0.0;
    YiOut = Ii0;
    YdOut = -Kn*Id0;
    YbOut = YpOut + YiOut + YdOut;

    if (bGotYmin){
        Ymin  = pYmin->dGet();
    }
    if (bGotYmax){
        Ymax  = pYmax->dGet();
    }

    YOut  = std::max(std::min(YbOut, Ymax), Ymin);


    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    out << "PID: " << uLabel
        << " " << Kp
        << " " << Ki
        << " " << Kd
        << " " << Kn
        << " " << Ii0
        << " " << Id0
        << " " << Ks
        << std::endl;
}

Pid::~Pid()
{
    NO_OP;
}

void Pid::Output(OutputHandler& OH) const
{
    if (bToBeOutput())
    {
        std::ostream& out = OH.Loadable();
        out << std::setw(8) << GetLabel() // 1: label
            << " " << InputError          // 2: error
            << " " << YOut                // 3: output
            << " " << YpOut               // 4: proportional output
            << " " << YiOut               // 5: integral output
            << " " << YdOut               // 6: derivative output
            << " " << YbOut               // 7: output before saturation
            << " " << Ii                  // 8: integral (integrator)
            << " " << Id                  // 9: integral (derivator)
            << std::endl;
    }
}

void Pid::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

void Pid::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

//contributo al file di restart
std::ostream& Pid::Restart(std::ostream& out) const
{
    return out << "# not implemented yet" << std::endl;
}

VariableSubMatrixHandler&
Pid::AssJac(VariableSubMatrixHandler& Workmat,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
    Workmat.SetNullMatrix();
    return Workmat;
}

unsigned int Pid::iGetNumPrivData(void) const
{
    // how many private data
    // output
    // output before saturation
    return 2;
}

unsigned int Pid::iGetPrivDataIdx(const char* s) const
{
    ASSERT(s != NULL);

    struct
    {
        const char* s;
        int i;
    } sPrivData[] = {
        {"yout", YOUT},
        {"ybout", YBOUT},
        {0}
    };

    for (int i = 0; sPrivData[i].s != 0; i++)
    {
        if (strcasecmp(s,sPrivData[i].s) == 0)
        {
            return sPrivData[i].i;
        }
    }

    return 0;

}

doublereal Pid::dGetPrivData(unsigned int i) const
{
    if (i<=0 || i >= LASTPRIVDATA)
    {
        silent_cerr("PID("<<GetLabel()<<"): "
        "private data "<< i << "not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    else
    {
        switch (i)
        {
            case YOUT   : return YOut;
            case YBOUT  : return YbOut;
        }
    }

    return 0.;
}

int Pid::iGetNumConnectedNodes(void) const
{
    return 0;
}

void Pid::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    connectedNodes.resize(0);
}

void Pid::SetValue(DataManager* pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    NO_OP;
}

unsigned Pid::iGetInitialNumDof(void) const
{
    return 0;
}

SubVectorHandler&
Pid::AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr)
{
    WorkVec.ResizeReset(0);
    // get setpoint
    if (bSetPoint) {
        SetPoint = pSetPoint->dGet();
    }
    // get measure
    Measure = pMeasure->dGet();
    // saturations
    if (bGotYmin) {
        Ymin = pYmin->dGet();
    }
    if (bGotYmax) {
        Ymax = pYmax->dGet();
    }

    // retrieve time step
    dt = DoTime.dGet();
    // compute input error
    InputError = SetPoint - Measure;
    // error derivative
    dInputError = (InputError-PrevInputErr)/dt;
    PrevInputErr = InputError;

    DEBUGCOUT( "InputError: " << InputError << std::endl);
    DEBUGCOUT( "dt: " << dt << std::endl);
    DEBUGCOUT( "dInputError/dt: " << dInputError << std::endl);

    // proportional control
    YpOut = Kp*InputError;
    // integral control
    // anti windup
    e2 = InputError - Ks*(YbOut-YOut);
    //  update output
    YiOut = Ii + Ki*e2*dt;
    // update integral
    Ii = YiOut;
    // derivative control
    // compute internal error
    e1 = Kd*InputError - Id;
    // update output
    YdOut = Kn*e1;
    // update integral
    Id = Id + YdOut*dt;
    // compute output before saturation
    YbOut = YpOut + YiOut + YdOut;
    // compute output after saturation
    YOut = std::max(std::min(YbOut, Ymax), Ymin);

    return WorkVec;
}

VariableSubMatrixHandler&
Pid::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                    const VectorHandler& XCurr)
{
    WorkMat.SetNullMatrix();
    return WorkMat;
}

SubVectorHandler&
Pid::InitialAssRes(SubVectorHandler& WorkVec,
                    const VectorHandler& XCurr)
{
    WorkVec.Resize(0);
    return WorkVec;
}

////////////////////////////////////////////////////////////////////////


bool PIDSet(void)
{
    UserDefinedElemRead *rf = new UDERead<Pid>;

    if (!SetUDE("pid", rf))
    {
        delete rf;
        return false;
    }

    return true;
}

//#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
    if (!PIDSet())
    {
        silent_cerr("pid: "
                    "module_init("<< module_name << ") "
                    "failed" << std::endl);
        return -1;
    }

    return 0;

}
//#endif // ! STATIC_MODULES
