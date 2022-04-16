// Author: Matteo Daniele <matteo.daniele@polimi.it>
// pid module for mbdyn
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "userelem.h"
#include "dataman.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "grapes-integrator.h"

Integrator::Integrator(unsigned int uLabel, const DofOwner *pDO,
            DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout("\nModule: Integrator\n"
        "Author: Matteo Daniele <matteo.daniele@polimi.it>\n"
        "Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
        "		Politecnico di Milano					\n"
        "		http://www.aero.polimi.it				\n"
        "									\n"
        "	All rights reserved						\n"
        "Implementation of a Integrator\n"
        "Input:\n"
        "   Input value\n"
        "   I0: initial integral value\n"
        "\n"
        "Output:\n"
        "element label\n"
        "integrated value\n"
        "\n"
        "Private data:\n"
        << std::endl);

        if (!HP.IsArg()){
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }

    // constructor begins:
    // get timestep
    DoTime.Set(new TimeStepDriveCaller(pDM->pGetDrvHdl()));

    // gather input parameters
    DEBUGCOUT( "reading INTEGRATOR(" << uLabel << ") element..." << std::endl);
    //
    //bSetPoint = false;
    if (HP.IsKeyWord("input")){
        pInput = HP.GetDriveCaller();
        bInput = true;
    }
    else {
        DEBUGCERR( "INTEGRATOR(" << uLabel << "): input not provided"<< std::endl);
    }
    //
    if (HP.IsKeyWord("I0")){
        Ii0 = HP.GetReal();
    }
    else{
        DEBUGCOUT( "INTEGRATOR(" << uLabel << "): Ii0 for the integrator not provided, using " << Ii0 << " as initial integrator value" << std::endl);
    }
    if (HP.IsKeyWord("rule")){
        if (HP.IsKeyWord("trapezoidal")){
            bTrapezoidal = true;
            bLinear = false;
        }
    }

    // INITIALIZATION
    if (bInput){
        inputValue = pInput->dGet();
    }

    prevInputValue = Ii0;
    // initialize
    prevIntegral = Ii0;
    YOut = Ii0;
    


    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    out << "INTEGRATOR: " << uLabel
        << " " << Ii0
        << " " << inputValue
        << std::endl;
}

Integrator::~Integrator()
{
    NO_OP;
}

void Integrator::Output(OutputHandler& OH) const
{
    if (bToBeOutput())
    {
        std::ostream& out = OH.Loadable();
        out << std::setw(4) << GetLabel() // 1: label
            << " " << prevIntegral        // 2: proportional output
            << " " << YOut                // 3: integral output            
            << " " << inputValue          // 4: input
            << std::endl;
    }
}

void Integrator::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

void Integrator::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

//contributo al file di restart
std::ostream& Integrator::Restart(std::ostream& out) const
{
    return out << "# not implemented yet" << std::endl;
}

VariableSubMatrixHandler&
Integrator::AssJac(VariableSubMatrixHandler& Workmat,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
    Workmat.SetNullMatrix();
    return Workmat;
}

unsigned int Integrator::iGetNumPrivData(void) const
{
    // how many private data
    // output
    // output before saturation
    return 2;
}

unsigned int Integrator::iGetPrivDataIdx(const char* s) const
{
    ASSERT(s != NULL);

    struct
    {
        const char* s;
        int i;
    } sPrivData[] = {
        {"yout", YOUT},
        {"yprev", YPREV},
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

doublereal Integrator::dGetPrivData(unsigned int i) const
{
    if (i<=0 || i >= LASTPRIVDATA)
    {
        silent_cerr("INTEGRATOR("<<GetLabel()<<"): "
        "private data "<< i << "not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    else
    {
        switch (i)
        {
            case YOUT   : return YOut;
            case YPREV  : return prevIntegral;
        }
    }

    return 0.;
}

int Integrator::iGetNumConnectedNodes(void) const
{
    return 0;
}

void Integrator::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    connectedNodes.resize(0);
}

void Integrator::SetValue(DataManager* pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    NO_OP;
}

unsigned Integrator::iGetInitialNumDof(void) const
{
    return 0;
}

SubVectorHandler&
Integrator::AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr)
{
    WorkVec.ResizeReset(0);
    // get input
    prevInputValue = inputValue;
    if (bInput) {
        inputValue = pInput->dGet();
    }
    // retrieve time step
    dt = DoTime.dGet();

    DEBUGCOUT( "InputValue: " << inputValue << std::endl);
    DEBUGCOUT( "dt: " << dt << std::endl);
    
    prevIntegral = YOut;
    bTrapezoidal ? YOut = Ii0+(prevInputValue+inputValue)*dt/2 : YOut = Ii0+inputValue*dt;
    
    return WorkVec;
}

VariableSubMatrixHandler&
Integrator::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                    const VectorHandler& XCurr)
{
    WorkMat.SetNullMatrix();
    return WorkMat;
}

SubVectorHandler&
Integrator::InitialAssRes(SubVectorHandler& WorkVec,
                    const VectorHandler& XCurr)
{
    WorkVec.Resize(0);
    return WorkVec;
}

////////////////////////////////////////////////////////////////////////


bool IntegratorSet(void)
{
    UserDefinedElemRead *rf = new UDERead<Integrator>;

    if (!SetUDE("Integrator", rf))
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
    if (!IntegratorSet())
    {
        silent_cerr("Integrator: "
                    "module_init("<< module_name << ") "
                    "failed" << std::endl);
        return -1;
    }

    return 0;

}
//#endif // ! STATIC_MODULES
