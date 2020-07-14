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
functions/variables for C::E Model:MBDyn_CE_CEModel_XxxYyy
*/

#ifndef MBDYN_CE_H
#define MBDYN_CE_H

#include"chrono/ChConfig.h"
#include <vector>

extern "C" {
// IDs of coupling bodies and motors in C::E
struct MBDYN_CE_CEMODELDATA{
        unsigned MBDyn_CE_CEBody_Label;
        unsigned MBDyn_CE_CEMotor_Label;
    };
//opaque pointer to C::E system
typedef  void* pMBDyn_CE_CEModel_t;

pMBDyn_CE_CEModel_t MBDyn_CE_CEModel_Init(std::vector<double> & MBDyn_CE_CEModel_DataSave, // space for data save
                                          const double* pMBDyn_CE_CEFrame, const double* MBDyn_CE_CEScale, // coupling information: ground ref and scale for units
                                          std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label, // coupling information: coupling bodies
                                          const int & MBDyn_CE_CouplingType); // coupling information: coupling type

// destroy
void
MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel);

// C::E models receive coupling motion from the buffer
void MBDyn_CE_CEModel_RecvFromBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel,
                                  const std::vector<double> &MBDyn_CE_CouplingKinematic,
                                  const unsigned &MBDyn_CE_NodesNum,
                                  const std::vector<MBDYN_CE_CEMODELDATA> &MBDyn_CE_CEModel_Label,
                                  double time_step);

// C::E models send coupling forces to the buffer
void MBDyn_CE_CEModel_SendToBuf(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, std::vector<double> &MBDyn_CE_CouplingDynamic, 
                                double* pMBDyn_CE_CEFrame,
                                const unsigned& MBDyn_CE_NodesNum,
                                const double* MBDyn_CE_CEScale,
                                const std::vector<MBDYN_CE_CEMODELDATA> & MBDyn_CE_CEModel_Label);

// update CEModel, and do time integration.
void
MBDyn_CE_CEModel_DoStepDynamics(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, double time_step);

// save CEModel at current step for reloading it in the tight coupling scheme
// (before advance())
int
MBDyn_CE_CEModel_DataSave(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data);

// reload data in the tight coupling scheme at each iteration
int
MBDyn_CE_CEModel_DataReload(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, 
                        std::vector<double> & MBDyn_CE_CEModel_Data);
// destroy
//extern void
//MBDyn_CE_destroy(MBDyn_CE_t *);

// add if needed
//extern void
//MBDyn_CE_AfterPredict(MBDyn_CE_t *);

// add arguments as needed
//extern void
//MBDyn_CE_Exchange(MBDyn_CE_t *, double *x, double *R, double *f, double *m);

// add if needed
//extern void
//MBDyn_CE_AfterConvergence(MBDyn_CE_t *);
}
#endif // MBDYN_CE_H
