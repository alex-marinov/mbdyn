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

#ifndef TASK2CPU_H
#define TASK2CPU_H

#include "mbconfig.h"
#include "ac/pthread.h"

class Task2CPU {
public:     
     Task2CPU();
     ~Task2CPU();

     void SetCPU(int iCPU);
     void ClearCPU(int iCPU);
     int iGetCPU(int iCPU) const;
     int iGetCount() const;
     int iGetFirstCPU() const;
     int iGetNextCPU(int iCPU) const;
     static int iGetMaxSize();
     bool bSetAffinity() const;
     bool bGetAffinity();
     static const Task2CPU& GetGlobalState();
     static void SetGlobalState(const Task2CPU& oState);
     
private:     
#ifdef USE_PTHREAD_SETAFFINITY_NP
     cpu_set_t oCPUSet;
#endif
     static Task2CPU oGlobalState;
};

#endif /* TASK2CPU_H */
