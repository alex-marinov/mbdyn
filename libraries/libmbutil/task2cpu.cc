/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "myassert.h"
#include "task2cpu.h"

Task2CPU Task2CPU::oGlobalState;

Task2CPU::Task2CPU()
{
#ifdef USE_PTHREAD_SETAFFINITY_NP
     CPU_ZERO(&oCPUSet);
#endif
}

Task2CPU::~Task2CPU()
{
}

int Task2CPU::iGetCPU(int iCPU) const
{
     ASSERT(iCPU >= 0);
     ASSERT(iCPU < iGetMaxSize());
     
#ifdef USE_PTHREAD_SETAFFINITY_NP
     return CPU_ISSET(iCPU, &oCPUSet);
#else
     return 0;
#endif
}

void Task2CPU::SetCPU(int iCPU)
{
     ASSERT(iCPU >= 0);
     ASSERT(iCPU < iGetMaxSize());
     
#ifdef USE_PTHREAD_SETAFFINITY_NP
     CPU_SET(iCPU, &oCPUSet);
#endif
}

void Task2CPU::ClearCPU(int iCPU)
{
#ifdef USE_PTHREAD_SETAFFINITY_NP
     CPU_CLR(iCPU, &oCPUSet);
#endif
}

int Task2CPU::iGetFirstCPU() const
{
     int iCPU;
     
     for (iCPU = 0; iCPU < iGetMaxSize(); ++iCPU) {
	  if (iGetCPU(iCPU)) {
	       break;
	  }
     }

     ASSERT(iCPU >= 0);
     ASSERT(iCPU <= iGetMaxSize());

     return iCPU;
}

int Task2CPU::iGetNextCPU(int iCPU) const
{
     ASSERT(iCPU >= 0);
     ASSERT(iCPU < iGetMaxSize());
     
     for ( ++iCPU; iCPU < iGetMaxSize(); ++iCPU) {
	  if (iGetCPU(iCPU)) {
	       break;
	  }
     }

     ASSERT(iCPU >= 0);
     ASSERT(iCPU <= iGetMaxSize());
     
     return iCPU;
}

int Task2CPU::iGetCount() const
{
#ifdef USE_PTHREAD_SETAFFINITY_NP
     return CPU_COUNT(&oCPUSet);
#else
     return 0;
#endif
}

int Task2CPU::iGetMaxSize()
{
#ifdef USE_PTHREAD_SETAFFINITY_NP
     return CPU_SETSIZE;
#else
     return 0;
#endif
}

bool Task2CPU::bGetAffinity()
{
#ifdef USE_PTHREAD_SETAFFINITY_NP     
     return 0 == pthread_getaffinity_np(pthread_self(), sizeof(oCPUSet), &oCPUSet);
#else
     return false;
#endif
}

bool Task2CPU::bSetAffinity() const
{
#ifdef USE_PTHREAD_SETAFFINITY_NP     
     return 0 == pthread_setaffinity_np(pthread_self(), sizeof(oCPUSet), &oCPUSet);
#else
     return false;
#endif
}

const Task2CPU& Task2CPU::GetGlobalState()
{
     return oGlobalState;
}

void Task2CPU::SetGlobalState(const Task2CPU& oState)
{
     oGlobalState = oState;
}


