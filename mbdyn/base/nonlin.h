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
 
 /* 
  *
  * Copyright (C) 2003-2017
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

#ifndef NONLIN_H
#define NONLIN_H

#include "solverdiagnostics.h"
#include "external.h"
#include "nonlinpb.h"
#include "solman.h"

#include <iomanip>
#include <chrono>
#include <cfloat>
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
class InverseSolver;
 
class NonlinearSolverTest {
public:
	enum Type {
		NONE,
		NORM,
		MINMAX,
		RELNORM,
		SEPNORM,

		LASTNONLINEARSOLVERTEST
	};

	virtual ~NonlinearSolverTest(void);

        virtual Type GetType() const=0;

	/* loops over the vector Vec */
	virtual doublereal MakeTest(Solver *pS, const integer& Size,
			const VectorHandler& Vec, bool bResidual = false,
			doublereal dScaleAlgEqu = 1.,
			doublereal* pTestDiff=0);

	/* tests a single value, and returns the measure accordingly */
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const = 0;

	/* merges results of multiple tests */
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const = 0;

	/* post-processes the test */
	virtual doublereal TestPost(const doublereal& dRes) const;

	/* scales a single value */
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;

	/* returns pointer to global vector for absolute residual*/
	virtual VectorHandler* GetAbsRes() { return 0;}

	/* returns pointer to map of dimension and corresponding equations */
	virtual std::map<OutputHandler::Dimensions, std::set<integer>>* GetDimMap() { return 0; };
};

class NonlinearSolverTestNone : virtual public NonlinearSolverTest {
public:
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal MakeTest(Solver *pS, integer Size,
			const VectorHandler& Vec, bool bResidual = false,
			doublereal* pTestDiff=0);
};

class NonlinearSolverTestNorm : virtual public NonlinearSolverTest {
public:
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal TestPost(const doublereal& dRes) const;
};

class NonlinearSolverTestRelNorm : virtual public NonlinearSolverTest {
public:
	MyVectorHandler AbsRes;
	virtual VectorHandler* GetAbsRes();
		virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal TestPost(const doublereal& dRes) const;

	/* loops over Res and AbsRes */
	virtual doublereal MakeTest(Solver *pS, const integer& Size,
			const VectorHandler& Vec, bool bResidual = false,
			doublereal dScaleAlgEqu = 1., doublereal* pTestDiff=0);
};

class NonlinearSolverTestSepNorm : virtual public NonlinearSolverTest {
public:
	/* Indices for corresponding dimensions */
	std::map<OutputHandler::Dimensions, std::set<integer>> MapOfDimensionIndices;
	virtual std::map<OutputHandler::Dimensions, std::set<integer>>* GetDimMap();

	/* Vector of the absolute values */
	MyVectorHandler AbsRes;
	virtual VectorHandler* GetAbsRes();

		virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual doublereal TestPost(const doublereal& dRes) const;

	/* loops over Res components seperately */
	virtual doublereal MakeTest(Solver *pS, const integer& Size,
			const VectorHandler& Vec, bool bResidual = false,
			doublereal dScaleAlgEqu = 1., doublereal* pTestDiff=0);
};


class NonlinearSolverTestMinMax : virtual public NonlinearSolverTest {
public:
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
};

class NonlinearSolverTestScale : virtual public NonlinearSolverTest {
protected:
	const VectorHandler* pScale; 
	
public:
	NonlinearSolverTestScale(const VectorHandler* pScl = 0);
	virtual ~NonlinearSolverTestScale(void);
	virtual void SetScale(const VectorHandler* pScl);
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleNorm : virtual public NonlinearSolverTestScale,
	virtual public NonlinearSolverTestNorm {
public:
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleRelNorm : virtual public NonlinearSolverTestScale,
	virtual public NonlinearSolverTestNorm {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleSepNorm : virtual public NonlinearSolverTestScale,
	virtual public NonlinearSolverTestSepNorm {
public:
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestScaleMinMax : virtual public NonlinearSolverTestScale,
	virtual public NonlinearSolverTestMinMax {
public:
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;
	virtual void TestMerge(doublereal& dResCurr,
			const doublereal& dResNew) const;
	virtual const doublereal& dScaleCoef(const integer& iIndex) const;
};

class NonlinearSolverTestRange : public NonlinearSolverTest {
protected:
	integer m_iFirstIndex;
	integer m_iLastIndex;
	NonlinearSolverTest *m_pTest;

	bool bIsValid(const integer& iIndex) const;

public:
	NonlinearSolverTestRange(NonlinearSolverTest *pTest, integer iFirstIndex = -1, integer iLastIndex = -1);
	virtual ~NonlinearSolverTestRange(void);

#if 0
	virtual doublereal MakeTest(Solver *pS, const integer& Size,
		const VectorHandler& Vec, bool bResidual = false);
#endif
        virtual Type GetType() const;
	virtual void TestOne(doublereal& dRes, const VectorHandler& Vec,
			const integer& iIndex, doublereal dCoef) const;

	virtual void TestMerge(doublereal& dResCurr, const doublereal& dResNew) const;

	virtual doublereal TestPost(const doublereal& dRes) const;

	virtual const doublereal& dScaleCoef(const integer& iIndex) const;

	void SetRange(integer iFirstIndex, integer iLastIndex);
};

struct NonlinearSolverOptions
{
	bool bHonorJacRequest;

	enum ScaleFlags {
		SCALE_ALGEBRAIC_EQUATIONS_NO = 0,
		SCALE_ALGEBRAIC_EQUATIONS_YES = 1
	} eScaleFlags;

	doublereal dScaleAlgebraic;

	NonlinearSolverOptions(bool bHonorJacRequest = false,
		enum ScaleFlags eScaleFlags = SCALE_ALGEBRAIC_EQUATIONS_NO,
		doublereal dScaleAlgebraic = 1.);
};

class NonlinearSolver : public SolverDiagnostics, protected NonlinearSolverOptions
{
public:
 	class ErrSimulationDiverged : MBDynErrBase {
	public:
		ErrSimulationDiverged(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
 	class NoConvergence : public MBDynErrBase {
	public:
		NoConvergence(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
 	class TimeStepLimitExceeded: public NoConvergence {
 	public:
 		TimeStepLimitExceeded(MBDYN_EXCEPT_ARGS_DECL): NoConvergence(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
 	};
 	class MaxResidualExceeded: public NoConvergence {
 	public:
 		MaxResidualExceeded(MBDYN_EXCEPT_ARGS_DECL): NoConvergence(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
 	};
	class ConvergenceOnSolution : public MBDynErrBase {
	public:
		ConvergenceOnSolution(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

	enum Type {
		UNKNOWN = -1,

		NEWTONRAPHSON,
		MATRIXFREE,
		LINESEARCH,
                BFGS,
		DEFAULT = NEWTONRAPHSON,

		LASTSOLVERTYPE
	};

protected:
	integer Size;
	integer TotJac;	

#ifdef USE_MPI
	bool bParallel;
#endif /* USE_MPI */
	NonlinearSolverTest *pResTest;
	NonlinearSolverTest *pSolTest;
#ifdef USE_EXTERNAL	
	External::ExtMessage ExtStepType;
#endif /* USE_EXTERNAL */

	virtual bool MakeSolTest(Solver* pS,
		const VectorHandler& Vec,
		const doublereal& dTol,
		doublereal& dTest);

	enum CPUTimeType {
		CPU_RESIDUAL,
		CPU_JACOBIAN,
		CPU_LINEAR_SOLVER,
		CPU_LAST_TYPE
	};

	doublereal dGetCondMax()const { return dMaxCond; }
	doublereal dGetCondMin()const { return dMinCond; }
	doublereal dGetCondAvg()const { return dSumCond / iNumCond; }
        inline std::chrono::nanoseconds dGetTimeCPU(CPUTimeType eType) const;
	inline void AddCond(doublereal dCond);
        inline void ResetCond();
        inline void AddTimeCPU(std::chrono::nanoseconds dTime, CPUTimeType eType);

    friend class CPUStopWatch;
                
    class CPUStopWatch  {
    public:
        explicit CPUStopWatch(NonlinearSolver& oSolver, CPUTimeType eType)
            :oSolver(oSolver),
             eType(eType),
	     eStatus(SWST_INACTIVE),	     
             dStartTimeCPU(std::chrono::nanoseconds(0)),
             dElapsedCPU(std::chrono::nanoseconds(0)) {
        }

        ~CPUStopWatch() {
            Toc();
        }
        
	std::chrono::time_point<std::chrono::high_resolution_clock> Tic() {
            using namespace std::chrono;
	    
	    dStartTimeCPU = high_resolution_clock::now();
	    eStatus = SWST_ACTIVE;

            return dStartTimeCPU;
        }

	 std::chrono::time_point<std::chrono::high_resolution_clock> Tic(CPUStopWatch& oOther) {
            ASSERT(&oSolver == &oOther.oSolver);
	    
            dStartTimeCPU = oOther.Toc();
	    eStatus = SWST_ACTIVE;
            
            return dStartTimeCPU;
        }
        
        std::chrono::time_point<std::chrono::high_resolution_clock> Toc() {
            using namespace std::chrono;
	    
            time_point<high_resolution_clock> dEndTimeCPU = high_resolution_clock::now();;
	    
	    if (eStatus == SWST_ACTIVE) {
		 dElapsedCPU = dEndTimeCPU - dStartTimeCPU;
		 oSolver.AddTimeCPU(dElapsedCPU, eType);
            }

	    eStatus = SWST_INACTIVE;

            return dEndTimeCPU;
        }

        std::chrono::nanoseconds dGetElapsedCPU() const {
            using namespace std::chrono;
	    
            if (eStatus == SWST_ACTIVE) {
		 return duration_cast<nanoseconds>(high_resolution_clock::now() - dStartTimeCPU);
            } else {
                return dElapsedCPU;
            }
        }
                   
        friend inline std::ostream& operator<<(std::ostream& os, const CPUStopWatch& oWatch) {
            return oWatch.Print(os);
        }
	 
    private:
        std::ostream& Print(std::ostream& os) const {
	     using namespace std::chrono;
	     typedef duration<float, std::ratio<1, 1> > FloatSec;
	     
	     auto flags = os.flags();
	     auto prec = os.precision();
	     os.setf(std::ios::scientific);
	     os.precision(2);
	     os << FloatSec(dGetElapsedCPU()).count() << "s/" << FloatSec(oSolver.dGetTimeCPU(eType)).count() << "s";	     
	     os.flags(flags);
	     os.precision(prec);
	     
	     return os;
        }
                   
        NonlinearSolver& oSolver;
        const CPUTimeType eType;
	enum StatusType {
	      SWST_ACTIVE,
	      SWST_INACTIVE
	} eStatus;
	 
        std::chrono::time_point<std::chrono::high_resolution_clock> dStartTimeCPU;
        std::chrono::nanoseconds dElapsedCPU;
    };
    
private:
	integer iNumCond;
	doublereal dMaxCond;
	doublereal dMinCond;
	doublereal dSumCond;
        std::chrono::nanoseconds dTimeCPU[CPU_LAST_TYPE];

public:
	explicit NonlinearSolver(const NonlinearSolverOptions& options);

	virtual void SetTest(NonlinearSolverTest *pr, NonlinearSolverTest *ps);
		
	virtual ~NonlinearSolver(void);

	virtual bool MakeResTest(Solver* pS,
		const NonlinearProblem *pNLP,
		const VectorHandler& Vec,
		const doublereal& dTol,
		doublereal& dTest,
		doublereal& dTestDiff);
	
	virtual void Solve(const NonlinearProblem *pNLP,
			Solver *pS,
			const integer iMaxIter,
			const doublereal& Tol,
			integer& iIterCnt,
			doublereal& dErr,
			const doublereal& SolTol,
			doublereal& dSolErr) = 0;

	virtual integer TotalAssembledJacobian(void);

	virtual NonlinearSolverTest* pGetResTest(void)	{
		return pResTest;
	}
	
	virtual NonlinearSolverTest* pGetSolTest(void)	{
		return pSolTest;
	}

        enum NonlinearSolverHintReal {
                LINESEARCH_LAMBDA_MAX,
                LINESEARCH_LAMBDA_CURR,
                NONLINEAR_SOLVER_LAST_HINT_REAL
        };

        enum NonlinearSolverHintInteger {
                LINESEARCH_ITERATION_CURR,
                NONLINEAR_SOLVER_LAST_HINT_INTEGER
        };


    void SetNonlinearSolverHint(NonlinearSolverHintReal eType, doublereal dHint) {
        ASSERT(eType >= 0);
        ASSERT(eType < NONLINEAR_SOLVER_LAST_HINT_REAL);

        oSolverHints.rgRealVal[eType] = dHint;
    }
    
    doublereal GetNonlinearSolverHint(NonlinearSolverHintReal eType) const {
        ASSERT(eType >= 0);
        ASSERT(eType < NONLINEAR_SOLVER_LAST_HINT_REAL);
        
        return oSolverHints.rgRealVal[eType];
    }

    void SetNonlinearSolverHint(NonlinearSolverHintInteger eType, integer iHint) {
        ASSERT(eType >= 0);
        ASSERT(eType < NONLINEAR_SOLVER_LAST_HINT_INTEGER);

        oSolverHints.rgIntegerVal[eType] = iHint;
    }
    
    integer GetNonlinearSolverHint(NonlinearSolverHintInteger eType) const {
        ASSERT(eType >= 0);
        ASSERT(eType < NONLINEAR_SOLVER_LAST_HINT_INTEGER);
        
        return oSolverHints.rgIntegerVal[eType];
    }

     std::ostream& PrintSolverTime(std::ostream& os) const;
     
#ifdef USE_EXTERNAL
	void SetExternal(const External::ExtMessage Ty);
	
protected:
	void SendExternal(void);
#endif /* USE_EXTERNAL */

private:
        struct SolverHints {
           doublereal rgRealVal[NONLINEAR_SOLVER_LAST_HINT_REAL];
           integer rgIntegerVal[NONLINEAR_SOLVER_LAST_HINT_INTEGER];
        } oSolverHints;
};

inline void
NonlinearSolver::AddCond(doublereal dCond)
{
	iNumCond++;
	dSumCond += dCond;

	if (dCond > dMaxCond) {
		dMaxCond = dCond;
	}

	if (dCond < dMinCond) {
		dMinCond = dCond;
	}
}

inline void
NonlinearSolver::ResetCond()
{
    iNumCond = 0;
    dSumCond = 0;
    dMaxCond = 0;
    dMinCond = std::numeric_limits<doublereal>::max();
}
        
inline std::chrono::nanoseconds
NonlinearSolver::dGetTimeCPU(CPUTimeType eType) const
{
	ASSERT(eType >= 0);
	ASSERT(eType < CPU_LAST_TYPE);

	return dTimeCPU[eType];
}

inline void
NonlinearSolver::AddTimeCPU(std::chrono::nanoseconds dTime, CPUTimeType eType)
{
	ASSERT(eType >= 0);
	ASSERT(eType < CPU_LAST_TYPE);

	dTimeCPU[eType] += dTime;
}

#endif /* NONLIN_H */

