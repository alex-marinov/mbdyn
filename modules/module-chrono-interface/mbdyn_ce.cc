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

#include "mbdyn_ce.h"
#include "chrono/ChConfig.h"
#include "chrono_parallel/physics/ChSystemParallel.h"

using namespace chrono;
using namespace chrono::collision;

extern "C" void
MBDyn_CE_CEModel_Create(ChSystemParallelNSC *pMBDyn_CE_CEModel);

extern "C" pMBDyn_CE_CEModel_t MBDyn_CE_CEModel_Init()
{
	std::cout << "Initial MBDyn_CE_CEModel pointer:\n";
	ChSystemParallelNSC *pMBDyn_CE_CEModel = new ChSystemParallelNSC;
	if(pMBDyn_CE_CEModel==NULL)
	{
		std::cout << "Initial MBDyn_CE_CEModel pointer fails\n";
	}
	MBDyn_CE_CEModel_Create(pMBDyn_CE_CEModel);
	return pMBDyn_CE_CEModel;
}

extern "C" void
MBDyn_CE_CEModel_Destroy(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel)
{
	std::cout << "destroy the CE_model...\n";
	if(pMBDyn_CE_CEModel!=NULL)
	{
		// must convert to the correct type
		// delete  (int *) pMBDyn_CE_CEModel;
		delete  (ChSystemParallelNSC *) pMBDyn_CE_CEModel;
		pMBDyn_CE_CEModel = NULL;
	}
}

extern "C" void
MBDyn_CE_CEModel_Update(pMBDyn_CE_CEModel_t pMBDyn_CE_CEModel, double time_step)
{
	std::cout << "update CE_models:\n";
	if (pMBDyn_CE_CEModel!=NULL)
	{   // it's not a good idea to do this convert
		// it's dangerous, how to avoid?
		ChSystemParallelNSC* temp=(ChSystemParallelNSC *) pMBDyn_CE_CEModel;
		temp->DoStepDynamics(time_step);
		std::cout << temp->GetChTime()<<"\n";
	}
}