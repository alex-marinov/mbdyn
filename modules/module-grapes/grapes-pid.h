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
#ifndef GRAPES_PID_H
#define GRAPES_PID_H

extern bool PIDSet(void);
#include "mbconfig.h"
#include "elem.h"
#include "dataman.h"
#include "userelem.h"
#include "drive_.h"

class Pid : virtual public Elem, public UserDefinedElem
{
    private:

        // time drive owner
        DriveOwner DoTime;
        doublereal dt;

        // reference drive caller
        const DriveCaller* pSetPoint;
        // measured drive caller
        const DriveCaller* pMeasure;
        // reference and measured value
        doublereal SetPoint = 0.0;
        bool bSetPoint = false;
        doublereal Measure;
        // error
        doublereal PrevInputErr = 0.0;
        doublereal InputError;
        doublereal dInputError;
        doublereal e1;  // derivator
        doublereal e2;  // anti windup
        // output
        doublereal YOut;
        // PID constants
        doublereal Kp; // proportional
        doublereal Ki; // integral
        doublereal Kd; // derivative
        doublereal Kn = 100.0; // derivative washout constant
        // initial value of the integrals
        doublereal Ii0 = 0.0;
        doublereal Id0 = 0.0;
        // anti wind-up
        doublereal Ks = 1.0; // wind-up scale factor
        // wind-up saturation
        const DriveCaller* pYmin;
        const DriveCaller* pYmax;
        doublereal Ymin = std::numeric_limits<doublereal>::min();
        doublereal Ymax = std::numeric_limits<doublereal>::max();
        bool bGotYmin = false;
        bool bGotYmax = false;

        // output
        doublereal YpOut; // proportional output
        doublereal YiOut; // integral output
        doublereal YdOut; // derivative output
        doublereal YbOut; // output before saturation

        // integral value
        doublereal Ii;
        doublereal Id;

        ////////////////////////////////////////////////////////////////////////


    public:

        Pid( unsigned int uL, const DofOwner* pDO,
                DataManager* pDM, MBDynParser& HP);
        virtual ~Pid();



        ////////////////////////////////////////////////////////////////////////
        enum PrivData {
            YOUT = 1,
            YBOUT,
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

#endif // ! GRAPES_PID_H
