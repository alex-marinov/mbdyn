/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
  *
  * Copyright (C) 2003
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

#ifndef NONLIN_H
#define NONLIN_H

#include<external.h>
#include<nonlinpb.h>
#include<solman.h>  
#include <ac/float.h>
#include<vector>

/*
 * Directory tree rationale:
 * 
 * everything is based on NonlinearSolver; direct methods,
 * e.g. NewtonRaphsonSolver, inherit directly from NonlinearSolver
 * matrix-free methods are based on MatrixFreeSolver, which
 * inherits from NonlinearSolver, and may require Preconditioner.
 * The available matrix-free methods are BiCGStab and Gmres,
 * which requires UpHessMatrix.
 *
 * - nonlin.h 			NonlinearSolver
 *   |--- nr.h			NewtonRaphsonSolver
 *   |- precond.h		Preconditioner
 *   |  |- precond_.h 		FullJacobianPr
 *   |- mfree.h 		MatrixFreeSolver
 *      |--- bicg.h		BiCGStab
 *      |--- gmres.h		Gmres, UpHessMatrix
 */

class SolverDiagnostics {
protected:
 	unsigned OutputFlags;

	enum {
		OUTPUT_NONE		= 0x0000,

		OUTPUT_ITERS		= 0x0001,
		OUTPUT_RES		= 0x0002,
		OUTPUT_SOL		= 0x0004,
		OUTPUT_JAC		= 0x0008,
		OUTPUT_MSG		= 0x0010,

		OUTPUT_DEFAULT		= OUTPUT_MSG,

		OUTPUT_MASK		= 0x00FF
	};
public:

	SolverDiagnostics(unsigned OF = OUTPUT_DEFAULT);
	virtual ~SolverDiagnostics(void);
	
	void SetOutputFlags(unsigned OF);
	void AddOutputFlags(unsigned OF);
	void DelOutputFlags(unsigned OF);
		
	inline bool outputIters(void) const {
		return (OutputFlags & OUTPUT_ITERS);
	};
 
	inline bool outputRes(void) const {
		return (OutputFlags & OUTPUT_RES);
	};
 
	inline bool outputSol(void) const {
		return (OutputFlags & OUTPUT_SOL);
	};
 
	inline bool outputJac(void) const {
		return (OutputFlags & OUTPUT_JAC);
	};

        /*
	 * all messages not protected behind any other condition
	 * must be protected by a "if (outputMsg())" condition
	 */
	inline bool outputMsg(void) const {
		return (OutputFlags & OUTPUT_MSG);
	};
};

 
class NonlinearSolver : public SolverDiagnostics
{
public:
 	class ErrSimulationDiverged {};
 	class NoConvergence {};
	class ConvergenceOnSolution {};
	class ErrGeneric {};

protected:
	integer Size;
	integer TotJac;	
#ifdef USE_EXTERNAL	
	External::ExtMessage ExtStepType;
#endif /* USE_EXTERNAL */

#ifdef __HACK_SCALE_RES__
	VectorHandler* pScale; 
#endif /* __HACK_SCALE_RES__ */

	virtual doublereal MakeTest(const VectorHandler& Vec) = 0;

public:
	NonlinearSolver(void);
		
#ifdef __HACK_SCALE_RES__
	virtual void SetScale(const VectorHandler* pScl);
#endif /* __HACK_SCALE_RES__ */

	virtual ~NonlinearSolver(void);

	virtual void Solve(const NonlinearProblem* NLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr
#ifdef MBDYN_X_CONVSOL
			, const doublereal& SolTol,
			doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
			) = 0;

	virtual integer TotalAssembledJacobian(void);
	
#ifdef USE_EXTERNAL
	void SetExternal(const External::ExtMessage Ty);
	
protected:
	void SendExternal(void);
#endif /* USE_EXTERNAL */
};

#endif /* NONLIN_H */

