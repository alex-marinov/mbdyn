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
/* module-display
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
// simulation data display for mbdyn
#ifndef MODULE_DISPLAY_H
#define MODULE_DISPLAY_H

extern bool DisplaySet(void);
#include "mbconfig.h"
#include "elem.h"
#include "dataman.h"
#include "userelem.h"
#include "drive_.h"
#include "display_mainwindow.h"

class Display : virtual public Elem, public UserDefinedElem
{
    private:
        // time step drive owner
        DriveOwner DoTimeStep;
        // time drive owner
        DriveOwner DoTimeElapsed;
        // time step
        doublereal dt;
        // initial time of simulation
        doublereal t0;
        // drive callers are added (how many drive callers we want)
        // number of signals connected to display (updated at 
        // the beginning of the simulation)
        unsigned int uNSignals;
        // array of pointers to drive callers
        std::vector<const DriveCaller*> pSignals;
        // array of values read by drive callers
        std::vector<doublereal> dSignals;

        QApplication applicationDisplay(int argc, char* argv[]);

    public:
        Display(unsigned int uL, const DofOwner* pDo, DataManager* pDM, MBDynParser& HP);
        virtual ~Display();

        enum PrivData {

            DISPLAYLABEL = 1,
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

#endif // ! MODULE_DISPLAY_H