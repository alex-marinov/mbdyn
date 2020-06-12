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


#ifndef MBDYN_CE_H
#define MBDYN_CE_H

#include"chrono/ChConfig.h"

extern "C" {
//opaque pointer to C::E system
typedef  void* pMBDyn_CE_CEModel_t;

pMBDyn_CE_CEModel_t MBDyn_CE_CEModel_Init();

// destroy
extern void
MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel_t);

// update CEModel
extern void
MBDyn_CE_CEModel_Update(pMBDyn_CE_CEModel_t, double time_step);

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

/*class MBDyn_CE_CEModel
{
private:
  
public:
  MBDyn_CE_CEModel();
  ~MBDyn_CE_CEModel();
  // opaque struct to C::E system
  struct MBDYN_CE_CEMODEL_DATA{
    double F[3]; //forces from CEModel;
    double M[3]; //torques from CEModel;
  }MBDyn_CE_CEModel_Data;

  enum MBDYN_CE_CEMODEL_SIMCMD
  {
    CEMODEL_SAVE_DATA=1,
    CEMODEL_FIRST_SEND=2,
    CEMODEL_REGULAR_STEP=3
  }MBDyn_CE_CEModel_SimCmd;

  void getMBDyn_CE_CEModel_Data();
};

MBDyn_CE_CEModel::MBDyn_CE_CEModel()
{
}

MBDyn_CE_CEModel::~MBDyn_CE_CEModel()
{
}*/

#endif // MBDYN_CE_H
