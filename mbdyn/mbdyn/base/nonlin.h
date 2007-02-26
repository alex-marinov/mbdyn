/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
  * Copyright (C) 2003-2007
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

#ifndef NONLIN_H
#define NONLIN_H

#include <solverdiagnostics.h>
#include <external.h>
#include <nonlinpb.h>
#include <solman.h>  
#include <ac/float.h>
#include <vector>

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

/* Needed for callback declaration; defined in <mbdyn/base/solver.h> */
class Solver;
 
class NonlinearSolverTest {
public:
	enum Type {
		NONE,
		NORM,
		MINMAX,

		LASTNONLINEARSOLVERTEST
	};

	virtual ~NonlinearSolverTest(void);

	/* loops over the vector Vec */
	virtual doublereal MakeTest(Solver *pS, const integer& Size,
			const VectorHandler& Vec, bool bResidual = false);

	/* tests a single value, and returns the measure accordingly */
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const = 0;

	/* merges results of multiple tests */
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const = 0;

	/* post-processes the test */
	virtual doublereal TestPost(const doublereal& dRes) const;

	/* scales a single value */
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestNone : public NonlinearSolverTest {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal MakeTest(Solver *pS, integer Size,
			const VectorHandler& Vec, bool bResidual = false);
};

class NonlinearSolverTestNorm : public NonlinearSolverTest {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal TestPost(const doublereal& dRes) const;
};

class NonlinearSolverTestMinMax : public NonlinearSolverTest {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
};

class NonlinearSolverTestScale : public NonlinearSolverTest {
protected:
	const VectorHandler* pScale; 
	
public:
	NonlinearSolverTestScale(const VectorHandler* pScl = 0);
	virtual ~NonlinearSolverTestScale(void);
	virtual void SetScale(const VectorHandler* pScl);
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleNorm : public NonlinearSolverTestScale,
	public NonlinearSolverTestNorm {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleMinMax : public NonlinearSolverTestScale,
	public NonlinearSolverTestMinMax {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolver : public SolverDiagnostics
{
public:
 	class ErrSimulationDiverged {};
 	class NoConvergence {};
	class ConvergenceOnSolution {};
	class ErrGeneric {};

	enum Type {
		UNKNOWN = -1,

		NEWTONRAPHSON,
		MATRIXFREE,

		DEFAULT = NEWTONRAPHSON,

		LASTSOLVERTYPE
	};

protected:
	integer Size;
	integer TotJac;	
	bool bHonorJacRequest;
#ifdef USE_MPI
	bool bParallel;
#endif /* USE_MPI */
	NonlinearSolverTest *pResTest;
	NonlinearSolverTest *pSolTest;
#ifdef USE_EXTERNAL	
	External::ExtMessage ExtStepType;
#endif /* USE_EXTERNAL */

	virtual doublereal MakeResTest(Solver* pS, const VectorHandler& Vec);
	virtual doublereal MakeSolTest(Solver* pS, const VectorHandler& Vec);

public:
	NonlinearSolver(bool JacReq = false);

	virtual void SetTest(NonlinearSolverTest *pr, NonlinearSolverTest *ps);
		
	virtual ~NonlinearSolver(void);

	virtual void Solve(const NonlinearProblem *pNLP,
			Solver *pS,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr,
			const doublereal& SolTol,
			doublereal& dSolErr) = 0;

	virtual integer TotalAssembledJacobian(void);
	
#ifdef USE_EXTERNAL
	void SetExternal(const External::ExtMessage Ty);
	
protected:
	void SendExternal(void);
#endif /* USE_EXTERNAL */
};

#endif /* NONLIN_H */

