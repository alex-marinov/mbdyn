/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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
  Copyright (C) 2020(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifndef STRUMPACK_SOLUTION_MANAGER_HH
#define STRUMPACK_SOLUTION_MANAGER_HH

#ifdef USE_STRUMPACK
#include <iostream>
#include <vector>

#include <StrumpackSparseSolver.hpp>

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"

#include "spmapmh.h"
#include "sp_gradient_spmh.h"

class StrumpackSolver: public LinearSolver {
public:
     explicit StrumpackSolver(SolutionManager* pSolMan,
			      integer iNumThreads,
			      integer iNumIter,
			      unsigned uSolverFlags = 0u,
			      integer iVerbose = 0);
     ~StrumpackSolver(void);

     virtual void Reset(void) override;
     virtual void Solve(void) const override;

     virtual void MakeCompactForm(class SparseMatrixHandler& mh,
				  std::vector<doublereal>& Ax,
				  std::vector<integer>& Ar,
				  std::vector<integer>& Ac,
				  std::vector<integer>& Ap) const override;

     void ResetSymbolic();
     
private:
     mutable strumpack::StrumpackSparseSolver<doublereal, integer> oSolver;
     mutable integer iNumNonZeros;
     mutable bool bDoOrderingStep;
};

template <typename MatrixHandlerType>
class StrumpackSolutionManager: public SolutionManager {
public:
     explicit StrumpackSolutionManager(integer iDim,
				       integer iNumThreads,
				       integer iNumIter,
				       unsigned uSolverFlags = 0u,
				       integer iVerbose = 0);
     virtual ~StrumpackSolutionManager(void);
#ifdef DEBUG
     virtual void IsValid(void) const;
#endif /* DEBUG */
     virtual void MatrInitialize() override;
     virtual void MatrReset(void) override;
     virtual void Solve(void) override;
     virtual MatrixHandler* pMatHdl(void) const override;
     virtual MyVectorHandler* pResHdl(void) const override;
     virtual MyVectorHandler* pSolHdl(void) const override;

private:
     inline StrumpackSolver* pGetSolver() const;

     mutable MatrixHandlerType A;
     mutable MyVectorHandler x, b;
     std::vector<integer> Ai, Ap;
     std::vector<doublereal> Ax;
};
#endif
#endif 
