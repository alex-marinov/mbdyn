/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
  * With the contribution of Runsen Zhang <runsen.zhang@polimi.it>
  * during Google Summer of Code 2020
  */
/* name rules for MBDyn-chrono::engine:
functions/variables: MBDyn_CE_XxxYyy (eg. MBDyn_CE_NodesNum)
pointer: pMBDyn_CE_XxxYyy
structures/class: MBDYN_CE_XXXYYY (eg. MBDYN_CE_POINTDATA)
temporary variables: mbdynce_xxx (eg. mbdynce_point)
*/

#ifndef MODULE_CHRONO_INTERFACE_H
#define MODULE_CHRONO_INTERFACE_H

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>
#include <vector>

#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"
#include "converged.h"

#include "mbdyn_ce.h" // interfaces functions to C::E;




/* start class ChronoInterfaceBaseELem: a base class*/
// functions: refer to "/base/extforce.h" class ExtForce;
// functions that introduces methods to handle the simulation ("bsae/simentity.h")

class ChronoInterfaceBaseElem
    : virtual public Elem,
      public UserDefinedElem
{
private:
    DataManager *m_pDM;
    // A function that transfer doublereal Vec3 to double Vec3 [data between MBDyn and C::E should use the same type]
    void MBDyn_CE_Vec3D(const Vec3& mbdynce_Vec3, double *mbdynce_temp, double MBDyn_CE_CELengthScale);
    // A function that transfer doublereal Mat3x3 to double Mat3x3 [data between MBDyn and C::E should use the same type]
    void MBDyn_CE_Mat3x3D(const Mat3x3& mbdynce_Mat3x3, double *mbdynce_temp);
    double MBDyn_CE_CalculateError(); // calculate the error of coupling force.

protected:
    std::vector<double> MBDyn_CE_CouplingKinematic; // for coupling motion
    std::vector<double> MBDyn_CE_CouplingDynamic; // for coupling forces
    std::vector<double> MBDyn_CE_CouplingDynamic_pre; // for coupling forces in last iterations.
    double *pMBDyn_CE_CouplingKinematic_x = NULL; // consistent with the external struc force element
    double *pMBDyn_CE_CouplingKinematic_R = NULL;
    double *pMBDyn_CE_CouplingKinematic_xp = NULL;
    double *pMBDyn_CE_CouplingKinematic_omega = NULL;
    double *pMBDyn_CE_CouplingKinematic_xpp = NULL;
    double *pMBDyn_CE_CouplingKinematic_omegap = NULL;
    double *pMBDyn_CE_CEFrame = NULL; // the position [3] and the orietation [9] of chrono ground coordinate
    double *pMBDyn_CE_CouplingDynamic_f = NULL;
    double *pMBDyn_CE_CouplingDynamic_m = NULL;
    double *pMBDyn_CE_CouplingDynamic_f_pre = NULL;
    double *pMBDyn_CE_CouplingDynamic_m_pre = NULL;
    struct{
        unsigned Size_Kinematic;
        unsigned Size_Dynamic;
    } MBDyn_CE_CouplingSize;
    // some parameters about the convergence
    unsigned MBDyn_CE_CouplingIter_Max;
    unsigned MBDyn_CE_CouplingIter_Count;
    double MBDyn_CE_Coupling_Tol; 

protected:
    Converged MBDyn_CE_CEModel_Converged; // denote whether the coupling variables are converged
    bool bMBDyn_CE_CEModel_DoStepDynamics;         // detect whether CEModel is needed to be simulated and sends back data
    bool bMBDyn_CE_FirstSend;      // whether the current residual is the first or not..
    
public:
    pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel = NULL;
    std::vector<double> MBDyn_CE_CEModel_Data; // for reload C::E data in the tight coupling scheme
    // Coupling nodes information
    struct MBDYN_CE_POINTDATA {
        unsigned MBDyn_CE_uLabel;
        unsigned MBDyn_CE_CEBody_Label; //coupling bodies in C::E model
        const StructNode *pMBDyn_CE_Node;
        Vec3 MBDyn_CE_Offset;
		Vec3 MBDyn_CE_F;
		Vec3 MBDyn_CE_M;
    };

    std::vector<MBDYN_CE_CEMODELDATA> MBDyn_CE_CEModel_Label; //IDs of coupling bodies and motors in C::E model
    double MBDyn_CE_CEScale[4]; // the Unit used in Chrono::Engine. 1 Unit(m) in MBDyn = MBDyn_CE_CEScale * Unit() in Chrono::Engine;

protected:
    double time_step;
    std::vector<MBDYN_CE_POINTDATA> MBDyn_CE_Nodes;
    unsigned MBDyn_CE_NodesNum; // for now, only the case of one node
public:
    // 0: loose interface
    // 1: tight coupling
    // >1: exchange every iCoupling iterations // TO DO
    enum MBDyn_CE_COUPLING
    {
        COUPLING_NONE = -2,
        COUPLING_STSTAGGERED = -1, // meanless, only for consistent with strext
        COUPLING_LOOSE = 0,
        COUPLING_TIGHT = 1,
        //COUPLING_MULTIRATE >1 // TO DO
    };
    int MBDyn_CE_CouplingType;

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
    //unsigned int iGetNumPrivData(void) const;
    /* functions that introduces methods to handle the simulation: start*/



    /* functions for element, set the Jac and Res: start*/
    // Initial Assemble
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

    /* functions for coupling variables: start*/
    void MBDyn_CE_SendDataToBuf(); // write the kinematic variables to the vector.
    void MBDyn_CE_RecvDataFromBuf(); // read the dynamic variables from the Buf.
    /* functions for coupling variables: start*/

    /*Miscellaneous refers to the module-template2.cc and force.h: start*/
    
    virtual void Output(OutputHandler &OH) const;                                              // for output;
    std::ostream &Restart(std::ostream &out) const; // module information
    virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	}; // force style
    /*Miscellaneous refers to the module-template2.cc and force.h: end*/
};

/*// Read struct and read function; can use the UDERead <> instead;

struct ChronoInterfaceElemRead : public UserDefinedElemRead{
    virtual ~ChronoInterfaceElemRead(void) { NO_OP; };
    virtual ChronoInterfaceBaseElem *
    Read(unsigned uLabel, const DofOwner* pDO,
         DataManager* const pDM, MBDynParser& HP) const; 
};

ChronoInterfaceBaseElem *
ChronoInterfaceElemRead::Read(unsigned uLabel, const DofOwner *pDO,
                              DataManager *const pDM, MBDynParser &HP) const
{
    // Read element: obtain information from MBDyn script
    return new ChronoInterfaceBaseElem(uLabel, pDO, pDM, HP);
}*/
#endif
