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

#ifdef USE_TRILINOS
#undef HAVE_BLAS
#undef HAVE_BOOL
#include "aztecoowrap.h"

AztecOOSolutionManager::AztecOOSolutionManager(integer Dim, integer iMaxIter, doublereal dTol, integer iVerbose)
     :x(Dim, Comm),
      b(Dim, Comm),
      A(Dim, Dim, 1, Comm),
      problem(A, x, b),
      solver(problem),
      iMaxIter(iMaxIter),
      dTol(dTol)
{
}

AztecOOSolutionManager::~AztecOOSolutionManager(void)
{
}

#ifdef DEBUG
void AztecOOSolutionManager::IsValid(void) const
{
     A.IsValid();
     x.IsValid();
     b.IsValid();
}
#endif

void AztecOOSolutionManager::MatrReset(void)
{
     A.Reset();
}

void AztecOOSolutionManager::Solve(void)
{
     A.FillComplete();
     
     integer ierr = solver.Iterate(iMaxIter, dTol);

     if (ierr != 0) {
          silent_cerr("AztecOO warning: linear solver did not converge");
     }
}

MatrixHandler* AztecOOSolutionManager::pMatHdl(void) const
{
     return &A;
}

VectorHandler* AztecOOSolutionManager::pResHdl(void) const
{
     return &b;
}

VectorHandler* AztecOOSolutionManager::pSolHdl(void) const
{
     return &x;
}

bool AztecOOSolutionManager::bGetConditionNumber(doublereal& dCond) const
{
     return false;
}
#endif
