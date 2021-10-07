// Author: Matteo Daniele <matteo.daniele@polimi.it>
// display module for mbdyn
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <iostream>
#include <cmath>
#include <fstream>
#include "display_mainwindow.h"
#include "module-display.h"

Display::Display(unsigned int uLabel, const DofOwner *pDO,
                    DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout("\nModule: DISPLAY\n"
        "Author: Matteo Daniele <matteo.daniele@polimi.it>\n"
        "Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
        "		Politecnico di Milano					\n"
        "		http://www.aero.polimi.it				\n"
        "									\n"
        "	All rights reserved						\n"
        "Implementation of a Display for visualization of MBDyn data during simulation\n"
        << std::endl);
    }

    if (!HP.IsArg()){
        throw NoErr(MBDYN_EXCEPT_ARGS);
    }

    // set time step drive caller
    DoTimeStep.Set(new TimeStepDriveCaller(pDM->pGetDrvHdl()));
    // set time elapsed drive caller
    DoTimeElapsed.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));

    // prepare window input
    if (HP.IsKeyWord("number" "of" "inputs"))
    {
        uNSignals = HP.GetInt();
    }

    // resize array of pointers to drive callers and
    // the value of signals to the number of signals that will be displayed
    pSignals.resize(uNSignals);
    // resize array of signals values to numebr
    dSignals.resize(uNSignals);

    // prepare the drive caller pointers for each signal
    for (int i = 0; i < uNSignals; i++)
    {
        if (HP.IsKeyWord("signal")){
            pSignals[i] = HP.GetDriveCaller();
        }
    }

    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    out << "Display: " << uLabel;
    for (int i = 0; i < uNSignals; i++){
        out << " " << pSignals[i]->GetLabel();
    }
    out << std::endl;

    t0 = DoTimeElapsed.dGet();

    // inizializzazione finestra con tutti i fronzoli,
    // i salamelecchi, gli inchini, i ricchi premi e i cotillon
    MainWindow appWindow;
    appWindow.show();
    //applicationDisplay.exec(); 


}

Display::~Display()
{
    //Display::~applicationDisplay();
    NO_OP;
}

void Display::Output(OutputHandler& OH) const
{
    if (bToBeOutput())
    {
        std::ostream& out = OH.Loadable();
        out << std::setw(8) << GetLabel();   // 1: label
        for(int i = 0; i < uNSignals; i++){
            out << " " << pSignals[i]->GetLabel();  // label of the drive callers
        }
        out << std::endl;
    }
}

void Display::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

void Display::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 0;
    *piNumCols = 0;
}

// restart file
std::ostream& Display::Restart(std::ostream& out) const
{
    return out << "# not implemented yet" << std::endl;
}

VariableSubMatrixHandler& Display::AssJac(VariableSubMatrixHandler& Workmat,
                                            doublereal dCoef,
                                            const VectorHandler& XCurr,
                                            const VectorHandler& XPrimeCurr)
{
    Workmat.SetNullMatrix();
    return Workmat;
}

unsigned int Display::iGetNumPrivData(void) const
{
    // how many private data
    // label of the display 
    return 1;
}

unsigned int Display::iGetPrivDataIdx(const char* s) const
{
    ASSERT(s != NULL);

    struct
    {
        const char* s;
        int i;
    } sPrivData[] = {
        {"label", DISPLAYLABEL},
        {0}
    };

    for (int i = 0; sPrivData[i].s != 0; i++)
    {
        if (strcasecmp(s, sPrivData[i].s) == 0)
        {
            return sPrivData[i].i;
        }
    }

    return 0;
}

doublereal Display::dGetPrivData(unsigned int i) const
{
    if (i<=0 || i >= LASTPRIVDATA)
    {
        silent_cerr("Display("<<GetLabel()<<"): "
        "private data " << i << "not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    else
    {
        switch (i)
        {
        case DISPLAYLABEL : return GetLabel();
        // TO BE IMPLEMENTED IF NEEDED
        }
    }
    return 0.;
}

int Display::iGetNumConnectedNodes(void) const
{
    return 0;
}

void Display::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
    connectedNodes.resize(0);
}

void Display::SetValue(DataManager* pDM,
                        VectorHandler& X, VectorHandler& XP,
                        SimulationEntity::Hints *ph)
{
    NO_OP;
}

unsigned Display::iGetInitialNumDof(void) const
{
    return 0;
}

SubVectorHandler& Display::AssRes(SubVectorHandler& WorkVec,
                                    doublereal dCoef,
                                    const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr)
{
    WorkVec.ResizeReset(0);

    // retrieve time step
    dt = DoTimeStep.dGet();

    // store the value of the attached drive callers
    for (int i = 0; i < uNSignals; i++)
    {
        dSignals[i] = pSignals[i]->dGet();
        // update window graph goes here
    }

    return WorkVec;
}

VariableSubMatrixHandler& 
Display::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                        const VectorHandler& XCurr)
{
    WorkMat.SetNullMatrix();
    return WorkMat;
}

SubVectorHandler&
Display::InitialAssRes(SubVectorHandler& WorkVec,
                        const VectorHandler& XCurr)
{
    WorkVec.Resize(0);
    return WorkVec;
}       


bool DisplaySet(void)
{
    UserDefinedElemRead *rf = new UDERead<Display>;

    if (!SetUDE("display", rf))
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
    if (!DisplaySet())
    {
        silent_cerr("display: "
                    "module_init("<< module_name << ") "
                    "failed" << std::endl);
        return -1;
    }

    return 0;

}
//#endif // ! STATIC_MODULES
