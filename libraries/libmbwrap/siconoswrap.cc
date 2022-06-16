/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#include "mbconfig.h"

#ifdef USE_SICONOS
#include "siconoswrap.h"

SiconosSolutionManager::SiconosSolutionManager(NumericsMatrix_types eMatType, integer iDim, integer iNumNz)
     :oRowMap(iDim), x(iDim, nullptr, &oRowMap), A(eMatType, iDim, iDim, iNumNz, &oRowMap)
{
}

SiconosSolutionManager::~SiconosSolutionManager()
{
}

#ifdef DEBUG
void SiconosSolutionManager::IsValid() const
{
     oRowMap.IsValid();
     A.IsValid();
     x.IsValid();
}
#endif

SiconosMatrixHandler* SiconosSolutionManager::pMatHdl() const
{
     return &A;
}

SiconosVectorHandler* SiconosSolutionManager::pResHdl() const
{
     return &x;
}

SiconosVectorHandler* SiconosSolutionManager::pSolHdl() const
{
     return &x;
}

bool SiconosSolutionManager::bGetConditionNumber(doublereal& dCond) const
{
     return false;
}

void SiconosSolutionManager::MatrReset()
{
     NM_set_LU_factorized(A.pGetMatrix(), false);
}

void SiconosSolutionManager::MatrInitialize()
{
     MatrReset();
}

void SiconosSolutionManager::Solve()
{
#ifdef DEBUG
     MyVectorHandler f(x);
#endif

     integer ierr = NM_gesv(A.pGetMatrix(), x.pdGetVec(), true);

     if (0 != ierr) {
          throw LinearSolver::ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }

#ifdef DEBUG
     A.MatVecDecMul(f, x);
     DEBUGCERR("||A*x-b||=" << f.Norm() << "\n");
#endif
}

SiconosDenseSolutionManager::SiconosDenseSolutionManager(integer iDim)
     :SiconosSolutionManager(NM_DENSE, iDim, -1)
{
}

SiconosDenseSolutionManager::~SiconosDenseSolutionManager()
{
}

SiconosSparseSolutionManager::SiconosSparseSolutionManager(integer iDim, integer iNumNz)
     :SiconosSolutionManager(NM_SPARSE, iDim, iNumNz)
{
}

SiconosSparseSolutionManager::~SiconosSparseSolutionManager()
{
}
#endif
