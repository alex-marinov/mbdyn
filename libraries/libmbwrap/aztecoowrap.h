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

#ifndef ___AZTEC_OO_SOLUTION_MANAGER_H__INCLUDED__
#define ___AZTEC_OO_SOLUTION_MANAGER_H__INCLUDED__

#ifdef USE_TRILINOS
#include "solman.h"
#undef HAVE_BLAS // FIXME: mbconfig.h should not be included from solman.h
#undef HAVE_BOOL
#include "epetravh.h"
#include "epetraspmh.h"
#include <Epetra_SerialComm.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

class AztecOOSolutionManager: public SolutionManager {
public:
     AztecOOSolutionManager(integer Dim, integer iMaxIter, doublereal dTol, integer iVerbose);

     virtual ~AztecOOSolutionManager(void);

#ifdef DEBUG
     virtual void IsValid(void) const override;
#endif

     virtual void MatrReset(void) override;

     virtual void Solve(void) override;

     virtual MatrixHandler* pMatHdl(void) const override;

     virtual VectorHandler* pResHdl(void) const override;

     virtual VectorHandler* pSolHdl(void) const override;

     virtual bool bGetConditionNumber(doublereal& dCond) const override;
     
private:
     Epetra_SerialComm Comm;
     mutable EpetraVectorHandler x, b;
     mutable EpetraSparseMatrixHandler A;
     Epetra_LinearProblem problem;
     AztecOO solver;
     const integer iMaxIter;
     const doublereal dTol;
};

#endif
#endif
