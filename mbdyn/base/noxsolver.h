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
          SOLVER_LINESEARCH_BASED   = 0x40,
          SOLVER_TRUST_REGION_BASED = 0x80,
          SOLVER_INEXACT_TRUST_REGION_BASED = 0x100,
          SOLVER_TENSOR_BASED               = 0x200,
          SOLVER_MASK = SOLVER_LINESEARCH_BASED |
          SOLVER_TRUST_REGION_BASED |
          SOLVER_INEXACT_TRUST_REGION_BASED |
          SOLVER_TENSOR_BASED,
          DIRECTION_NEWTON = 0x400,
          DIRECTION_STEEPEST_DESCENT = 0x800,
          DIRECTION_NONLINEAR_CG = 0x1000,
          DIRECTION_BROYDEN = 0x2000,
          DIRECTION_MASK = DIRECTION_NEWTON |
          DIRECTION_STEEPEST_DESCENT |
          DIRECTION_BROYDEN,
          FORCING_TERM_CONSTANT =  0x4000,
          FORCING_TERM_TYPE1    =  0x8000,
          FORCING_TERM_TYPE2    = 0x10000,
          FORCING_TERM_MASK = FORCING_TERM_CONSTANT |
          FORCING_TERM_TYPE1 |
          FORCING_TERM_TYPE2,
          LINESEARCH_BACKTRACK    = 0x20000,
          LINESEARCH_POLYNOMIAL   = 0x40000,
          LINESEARCH_MORE_THUENTE = 0x80000,
          LINESEARCH_MASK = LINESEARCH_BACKTRACK |
          LINESEARCH_POLYNOMIAL |
          LINESEARCH_MORE_THUENTE,
          LINEAR_SOLVER_GMRES    = 0x100000,
          LINEAR_SOLVER_CG       = 0x200000,
          LINEAR_SOLVER_CGS      = 0x400000,
          LINEAR_SOLVER_TFQMR    = 0x800000,
          LINEAR_SOLVER_BICGSTAB = 0x1000000,
          LINEAR_SOLVER_MASK = LINEAR_SOLVER_GMRES |
          LINEAR_SOLVER_CG |
          LINEAR_SOLVER_CGS |
          LINEAR_SOLVER_TFQMR |
          LINEAR_SOLVER_BICGSTAB,
          RECOVERY_STEP_TYPE_CONST     = 0x2000000,
          RECOVERY_STEP_TYPE_LAST_STEP = 0x4000000,
          RECOVERY_STEP_TYPE_MASK = RECOVERY_STEP_TYPE_CONST |
          RECOVERY_STEP_TYPE_LAST_STEP,
          USE_PRECOND_AS_SOLVER = 0x8000000,
          SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN = 0x10000000,
          SUFFICIENT_DEC_COND_ARED_PRED = 0x20000000,
          SUFFICIENT_DEC_COND_MASK = SUFFICIENT_DEC_COND_ARMIJO_GOLDSTEIN |
                                      SUFFICIENT_DEC_COND_ARED_PRED
     };

     doublereal dWrmsRelTol;
     doublereal dWrmsAbsTol;
     doublereal dTolLinSol;
     doublereal dMinStep;
     doublereal dRecoveryStep;
     doublereal dForcingTermMinTol;
     doublereal dForcingTermMaxTol;
     doublereal dForcingTermAlpha;
     doublereal dForcingTermGamma;
     integer iMaxIterLinSol;
     integer iKrylovSubSpaceSize;
     integer iMaxIterLineSearch;
     integer iInnerIterBeforeAssembly;
};

NonlinearSolver*
pAllocateNoxNonlinearSolver(const NonlinearSolverTestOptions& oSolverOpt,
                            const NoxSolverParameters& oParam);

#endif
#endif
