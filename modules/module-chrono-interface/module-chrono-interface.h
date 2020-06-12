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


#ifndef MODULE_CHRONO_INTERFACE_H
#define MODULE_CHRONO_INTERFACE_H

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "elem.h"
#include "strnode.h"
#include "dataman.h"
#include "userelem.h"

#include "mbdyn_ce.h" // interfaces functions to C::E;


/* start class ChronoInterfaceBaseELem: a base class*/
// functions: refer to "/base/extforce.h" class ExtForce;
// functions that introduces methods to handle the simulation ("bsae/simentity.h")

class ChronoInterfaceBaseElem
    : virtual public Elem,
      public UserDefinedElem
{
public:
    pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel = NULL;
    // 0: loose interface
    // 1: tight coupling
    // >1: exchange every iCoupling iterations // TO DO
    enum 
    {
        COUPLING_NONE = -2,
        COUPLING_STSTAGGERED = -1,
        COUPLING_LOOSE = 0,
        COUPLING_TIGHT = 1,
        COUPLING_MULTIRATE = 2 // TO DO
    } iCoupling;

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
