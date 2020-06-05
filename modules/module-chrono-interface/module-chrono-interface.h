
#ifndef MODULE_CHRONO_INTERFACE_H
#define MODULE_CHRONO_INTERFACE_H

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"

/*start class: buffer for interfaces*/
// send(): data: MBDyn->buffer->C::E
// recv(): data: C::E->buffer->MBDyn



/* start class ChronoInterfaceBaseELem: a base class*/
// functions: refer to "/base/extforce.h" class ExtForce;
// functions that introduces methods to handle the simulation ("bsae/simentity.h")

class ChronoInterfaceBaseElem
    : virtual public Elem,
      public UserDefinedElem
{
protected:
    Converged c;
    // exchange after predict ?
    bool bSendAfterPredict;
    
    int iCoupling;

    // some parameters to record the coupling states
    // interation counter
    mutable int iCouplingCounter;
    mutable bool bFirstSend; 
    mutable bool bFirstRecv; 
    enum SendWhen
    {
        SEND_FIRST_TIME,
        SEND_REGULAR,
        SEND_AFTER_CONVERGENCE
    };
    enum RecvWhen
    {
        RECV_FIRST_TIME,
        RECV_REGULAR,
        RECV_AFTER_CONVERGENCE
    };

    // functions for data exchanging
    void Send(ChronoInterfaceBaseElem::SendWhen sendwhen);// data from MBDyn to C::E;
    void Recv(ChronoInterfaceBaseElem::RecvWhen recvwhen);// data from C::E to MBDyn;

public:
    // 0: loose interface
	// 1: tight coupling
	// >1: exchange every iCoupling iterations // TO DO 
    enum
    {
        LOOSE = 0,
        TIGHT = 1,
        MULTIRATE = 2 // TO DO 
    };

    // constructor
    ChronoInterfaceBaseElem(unsigned uLabel,     // Label
                        const DofOwner *pDO, // ？
                        DataManager *pDM,    // a lot of information for solvers (nodes, elements, solver...)
                        MBDynParser &HP);    // for obtain data from mbdyn script);
    virtual ~ChronoInterfaceBaseElem(void);



    /* functions that introduces methods to handle the simulation: start*/
    virtual void SetValue(DataManager *pDM,
                          VectorHandler &X, VectorHandler &XP,
                          SimulationEntity::Hints *h = 0);
    virtual void Update(const VectorHandler &XCurr,
                        const VectorHandler &XprimeCurr);
    virtual void AfterConvergence(const VectorHandler &X,
                                  const VectorHandler &XP);
    virtual void AfterPredict(VectorHandler &X,
                              VectorHandler &XP);
    unsigned int iGetNumPrivData(void) const;
    /* functions that introduces methods to handle the simulation: start*/



    /* functions for element, set the Jac and Res: start*/
    // Initial Assemble
    virtual unsigned int iGetInitialNumDof(void) const;
    virtual void
    InitialWorkSpaceDim(integer *piNumRows, integer *piNumCols) const;
    VariableSubMatrixHandler &
    InitialAssJac(VariableSubMatrixHandler &WorkMat,
                  const VectorHandler &XCurr);
    SubVectorHandler &
    InitialAssRes(SubVectorHandler &WorkVec, const VectorHandler &XCurr);
    // Assemble
    virtual void WorkSpaceDim(integer *piNumRows,
                              integer *piNumCols) const;
    VariableSubMatrixHandler &
    AssJac(VariableSubMatrixHandler &WorkMat,
           doublereal dCoef,
           const VectorHandler &XCurr,
           const VectorHandler &XPrimeCurr);
    SubVectorHandler &
    AssRes(SubVectorHandler &WorkVec,
           doublereal dCoef,
           const VectorHandler &XCurr,
           const VectorHandler &XPrimeCurr);
    /* functions for element, set the Jac and Res: end*/
    


    /*Miscellaneous refers to the module-template2.cc and force.h: start*/
    virtual void Output(OutputHandler &OH) const; // for output;    
    int iGetNumConnectedNodes(void) const; 
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
    std::ostream &Restart(std::ostream &out) const; // module information
    /*Miscellaneous refers to the module-template2.cc and force.h: end*/


    /* data information of MBDyn side, refer to strext.h: start */
public:
    //struc includes all the information of a point
    struct PointData{      
        unsigned uLabel;
        const StructNode *pNode;
        Vec3 Offset;
        Vec3 F;
        Vec3 M;
    };

protected:
    std::vector<PointData> m_Points; // now only support one node


    /* data information of MBDyn side, refer to strext.h: end */
};

// Read struct and read function; can use the UDERead <> instead;
/*
struct ChronoInterfaceElemRead : public UserDefinedElemRead{
    virtual ~ChronoInterfaceElemRead(void) { NO_OP; };
    virtual ChronoInterfaceBaseElem *
    Read(unsigned uLabel, const DofOwner* pDO,
         DataManager* const pDM, MBDynParser& HP) const; // why const ?
};

ChronoInterfaceBaseElem *
ChronoInterfaceElemRead::Read(unsigned uLabel, const DofOwner *pDO,
                              DataManager *const pDM, MBDynParser &HP) const
{
    return new ChronoInterfaceBaseElem(uLabel, pDO, pDM, HP);
}*/
#endif
