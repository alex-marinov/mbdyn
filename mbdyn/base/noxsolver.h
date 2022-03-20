/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
  */

#ifndef __NOX_SOLVER_H__INCLUDED__
#define __NOX_SOLVER_H__INCLUDED__

#ifdef USE_TRILINOS
#include "nonlin.h"

struct NoxSolverParameters: public CommonNonlinearSolverParam {
     NoxSolverParameters();

     enum SolverFlags {
          JACOBIAN_NEWTON_KRYLOV = 0x10,
          JACOBIAN_NEWTON        = 0x20,
          JACOBIAN_OPERATOR_MASK = JACOBIAN_NEWTON_KRYLOV | JACOBIAN_NEWTON,
          ALGORITHM_LINESEARCH_BASED   = 0x40,
          ALGORITHM_TRUST_REGION_BASED = 0x80,
          ALGORITHM_INEXACT_TRUST_REGION_BASED = 0x100,
          ALGORITHM_TENSOR_BASED               = 0x200,
          ALGORITHM_MASK = ALGORITHM_LINESEARCH_BASED |
          ALGORITHM_TRUST_REGION_BASED |
          ALGORITHM_INEXACT_TRUST_REGION_BASED |
          ALGORITHM_TENSOR_BASED
     };
     
     doublereal dNewtonKrylovPerturbation;
     doublereal dWrmsRelTol;
     doublereal dWrmsAbsTol;
     doublereal dTolLinSol;
     integer iMaxIterLinSol;
};

NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam);

#endif
#endif
