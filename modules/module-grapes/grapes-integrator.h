/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
/* module-pid
 * Author: Matteo Daniele
 *
 * Copyright (C) 2008-2021
 *
 * Matteo Daniele <matteo.daniele@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */
// Author: Matteo Daniele <matteo.daniele@polimi.it>
// pid module for mbdyn
#ifndef GRAPES_INTEGRATOR_H
#define GRAPES_INTEGRATOR_H

extern bool IntegratorSet(void);
#include "elem.h"
#include "dataman.h"
#include "userelem.h"
#include "drive_.h"

class Integrator : virtual public Elem, public UserDefinedElem
{
    private:

        // time drive owner
        DriveOwner DoTime;
        doublereal dt;

        // measured drive caller
        const DriveCaller* pInput;
        bool bInput;
        doublereal inputValue;
        doublereal prevInputValue;
        // previous integral value
        doublereal prevIntegral = 0.0;
        // output
        doublereal YOut;
        // initial value of the integrals
        doublereal Ii0 = 0.0;
        // misc
        bool bTrapezoidal = false;
        bool bLinear = true;
        ////////////////////////////////////////////////////////////////////////


    public:

        Integrator( unsigned int uL, const DofOwner* pDO,
                DataManager* pDM, MBDynParser& HP);
        virtual ~Integrator();

        ////////////////////////////////////////////////////////////////////////
        enum PrivData {
            YOUT = 1,
            YPREV,
            LASTPRIVDATA
        };

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        //contributo al file di restart
        std::ostream& Restart(std::ostream& out) const;
        // assemblaggio yacobiano
        VariableSubMatrixHandler& AssJac(   VariableSubMatrixHandler& WorkMat,
                                            doublereal dCoef,
                                            const VectorHandler& /* XCurr */ ,
                                            const VectorHandler& /* XPrimeCurr */ );
        // assemblaggio residuo
        SubVectorHandler& AssRes(   SubVectorHandler& WorkVec,
                                    doublereal /* dCoef */ ,
                                    const VectorHandler& /* XCurr */ ,
                                    const VectorHandler& /* XPrimeCurr */ );
        unsigned int iGetNumPrivData(void) const;
        unsigned int iGetPrivDataIdx(const char* s) const;
        doublereal dGetPrivData(unsigned int i) const;
        int iGetNumConnectedNodes(void) const;
        void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
        void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
            SimulationEntity::Hints *ph);
        virtual unsigned int iGetInitialNumDof(void) const;
        // output
        virtual void Output(OutputHandler& OH) const;
        /* Contributo allo jacobiano durante l'assemblaggio iniziale */
        VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                                const VectorHandler& /* XCurr */ );
        /* Contributo al residuo durante l'assemblaggio iniziale */
        SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
                                        const VectorHandler& /* XCurr */ );

};

#endif // ! GRAPES_INTEGRATOR_H
