// Author: Matteo Daniele <matteo.daniele@polimi.it>
// pid module for mbdyn
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "userelem.h"
#include "dataman.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "grapes-adder.h"

AdderNode::AdderNode(unsigned int uLabel, const DofOwner *pDO,
            DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout("\nModule: AdderNode\n"
        "Author: Matteo Daniele <matteo.daniele@polimi.it>\n"
        "Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
        "		Politecnico di Milano					\n"
        "		http://www.aero.polimi.it				\n"
        "									\n"
        "	All rights reserved						\n"
        "Implementation of a AdderNode\n"
        "Input:\n"
        "   number of inputs\n"
        "   list of the inputs preceeded by the addition/subtraction term\n"
        "\n"
        "Output:\n"
        "adder node output\n"
        "\n"
        "Private data:\n"
        "adder node output\n"
        << std::endl);

        if (!HP.IsArg()){
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }

    // constructor begins:
    // get timestep
    DoTime.Set(new TimeStepDriveCaller(pDM->pGetDrvHdl()));

    // gather input parameters
    DEBUGCOUT( "reading ADDERNODE(" << uLabel << ") element..." << std::endl);
    //
    //bSetPoint = false;
    if (HP.IsKeyWord("ninputs")){
        nInputs = HP.GetInt();
        pInputs.reserve(nInputs);
        inputValues.reserve(nInputs);
    }
    //
    if (HP.IsKeyWord("input")){
        if (HP.IsKeyWord("+")){
            pInputs.push_back(HP.GetDriveCaller());
            addOrSub.push_back(1.0);
        }
        else if (HP.IsKeyWord("-")){
            pInputs.push_back(HP.GetDriveCaller());
            addOrSub.push_back(-1.0);
        }
    }

    // INITIALIZATION
    for (int k=0;k<nInputs;k++){
        inputValues[k] = pInputs[k]->dGet()*addOrSub[k];
        YOut += inputValues[k];
    }

    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    out << "ADDERNODE: " << uLabel
        << " " << nInputs
        << " " << YOut
        << std::endl;
}

AdderNode::~AdderNode()
{
    NO_OP;
}

void AdderNode::Output(OutputHandler& OH) const
{
    if (bToBeOutput())
    {
        std::ostream& out = OH.Loadable();
        out << std::setw(4) << GetLabel() // 1: label
            << " " << nInputs        // 2: number of inputs
            << " " << YOut           // 3: adder output           
            << std::endl;
    }
}

void AdderNode::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

void AdderNode::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

//contributo al file di restart
std::ostream& AdderNode::Restart(std::ostream& out) const
{
    return out << "# not implemented yet" << std::endl;
}

VariableSubMatrixHandler&
AdderNode::AssJac(VariableSubMatrixHandler& Workmat,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
    Workmat.SetNullMatrix();
    return Workmat;
}

unsigned int AdderNode::iGetNumPrivData(void) const
{
    // how many private data
    // output
    // output before saturation
    return 2;
}

unsigned int AdderNode::iGetPrivDataIdx(const char* s) const
{
    ASSERT(s != NULL);

    struct
    {
        const char* s;
        int i;
    } sPrivData[] = {
        {"yout", YOUT},
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

doublereal AdderNode::dGetPrivData(unsigned int i) const
{
    if (i<=0 || i >= LASTPRIVDATA)
    {
        silent_cerr("ADDERNODE("<<GetLabel()<<"): "
        "private data "<< i << "not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    else
    {
        switch (i)
        {
            case YOUT   : return YOut;
        }
    }

    return 0.;
}

int AdderNode::iGetNumConnectedNodes(void) const
{
    return 0;
}

void AdderNode::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    connectedNodes.resize(0);
}

void AdderNode::SetValue(DataManager* pDM,
    VectorHandler& X, VectorHandler& XP,
    SimulationEntity::Hints *ph)
{
    NO_OP;
}

unsigned AdderNode::iGetInitialNumDof(void) const
{
    return 0;
}

SubVectorHandler&
AdderNode::AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr)
{
    WorkVec.ResizeReset(0);
    // get input
    for (int k=0;k<nInputs;k++){
        inputValues[k] = pInputs[k]->dGet()*addOrSub[k];
        YOut += inputValues[k];
    }

    // retrieve time step
    // dt = DoTime.dGet();

    //DEBUGCOUT( "InputValue: " << inputValue << std::endl);
    //DEBUGCOUT( "dt: " << dt << std::endl);
    
    return WorkVec;
}

VariableSubMatrixHandler&
AdderNode::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                    const VectorHandler& XCurr)
{
    WorkMat.SetNullMatrix();
    return WorkMat;
}

SubVectorHandler&
AdderNode::InitialAssRes(SubVectorHandler& WorkVec,
                    const VectorHandler& XCurr)
{
    WorkVec.Resize(0);
    return WorkVec;
}

////////////////////////////////////////////////////////////////////////


bool AdderNodeSet(void)
{
    UserDefinedElemRead *rf = new UDERead<AdderNode>;

    if (!SetUDE("AdderNode", rf))
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
    if (!AdderNodeSet())
    {
        silent_cerr("AdderNode: "
                    "module_init("<< module_name << ") "
                    "failed" << std::endl);
        return -1;
    }

    return 0;

}
//#endif // ! STATIC_MODULES
