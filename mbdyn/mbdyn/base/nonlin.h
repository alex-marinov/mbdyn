/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

class NonlinearSolver
{

public:
 	class ErrSimulationDiverged{};
 	class NoConvergence{};
	class ConvergenceOnSolution{};
	class ErrGeneric{};
protected:
	integer Size;
	
	bool foutIters, foutRes, foutJac, foutSol;	
#ifdef USE_EXTERNAL	
	External::ExtMessage ExtStepType;
#endif /* USE_EXTERNAL */

#ifdef __HACK_SCALE_RES__
	VectorHandler* pScale; 
#endif /* __HACK_SCALE_RES__ */



public:
	NonlinearSolver(void): 
		Size(0),
		foutIters(false),
		foutRes(false),
		foutJac(false),
		foutSol(false)
#ifdef USE_EXTERNAL
		, ExtStepType(External::ERROR)  
#endif /* USE_EXTERNAL */
#ifdef __HACK_SCALE_RES__
		, VectorHandler* pScale(NULL) 
#endif /* __HACK_SCALE_RES__ */
		{ };
		
#ifdef __HACK_SCALE_RES__
	virtual void SetScale(const VectorHandler* pScl) {
		pScale = pScl;
		return;
	}  
#endif /* __HACK_SCALE_RES__ */

	virtual void SetOutputFlag(bool fIt, bool fRes, bool fJac, bool fSol) {
		foutIters = fIt;
		foutRes = fRes;
		foutJac = fJac;
		foutSol = fSol;
		return;
	};
		
	virtual ~NonlinearSolver(void) {};
	
	virtual void Solve(const NonlinearProblem* NLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal Toll,
			const doublereal SolToll,
			integer& iIterCnt,
			doublereal& dErr
#ifdef MBDYN_X_CONVSOL
			, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
			) = 0;

	
#ifdef USE_EXTERNAL
	void SetExternal(const External::ExtMessage Ty)
	{
		ExtStepType = Ty;
		return;
	
	};

protected:
	void SendExternal(void)
	{
		switch (ExtStepType) {
		
		case (External::EMPTY):
			External::SendNull();
			break;

		case (External::REGULAR):
			External::SendRegular();
			break;
			
		case (External::CLOSE):
			External::SendClose();
			break;
	
		case (External::ERROR):
		default:
			External::SendError();
		}
		return;
	};
#endif /* USE_EXTERNAL */			
};
 
class NewtonRaphsonSolver : public  NonlinearSolver
{
	SolutionManager* pSM;
	VectorHandler* 	pRes;
	VectorHandler* 	pSol;
	MatrixHandler*  pJac;
	flag fTrueNewtonRaphson;
	integer IterationBeforeAssembly;
	

public:

	NewtonRaphsonSolver(const flag fTNR, 
			const integer IterBfAss);
	
	~NewtonRaphsonSolver(void) { };
	
	void Solve(const NonlinearProblem* NLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal Toll,
			const doublereal SolToll,
			integer& iIterCnt,
			doublereal& dErr
#ifdef MBDYN_X_CONVSOL
			, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
			);
			
private:
	doublereal MakeTest(const VectorHandler& Vec);
};

class Preconditioner
{
public: 

	enum PrecondType{
		FULLJACOBIAN,
		UNKNOWN
	};
	
	virtual ~Preconditioner(void) { };
	
	virtual void Precond(VectorHandler& v, 
			SolutionManager* pSM) const = 0;
};


class FullJacobianPr : public Preconditioner
{	

public:
	FullJacobianPr(void) {};
	
	~FullJacobianPr(void) {};
	
	void Precond(VectorHandler& v, 
			SolutionManager* pSM) const {
		pSM->ChangeResPoint(v.pdGetVec());
		pSM->Solve();
	   	if (pSM->pSolHdl() != &v) {
			v = *(pSM->pSolHdl());
		}
	};
	 
};

const doublereal defaultTau = 1.e-7;
const doublereal defaultGamma = 0.9;
const doublereal defaultEtaMax = 0.9;

class BiCGMatrixFreeSolver : public NonlinearSolver
{
	SolutionManager* pSM;
	Preconditioner* pPM;
	VectorHandler* 	pRes;
	doublereal IterToll;
	integer MaxLinIt;
	doublereal Tau;
	doublereal gamma;
	doublereal etaMax; 
	integer PrecondIter; 
	bool fBuildMat;
	NonlinearProblem* pPrevNLP;
	
public:
	BiCGMatrixFreeSolver(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt); 

	~BiCGMatrixFreeSolver(void) { };
	
	void Solve(const NonlinearProblem* NLP,
			SolutionManager* pSolMan,
			const integer iMaxIter,
			const doublereal Toll,
			const doublereal SolToll,
			integer& iIterCnt,
			doublereal& dErr
#ifdef MBDYN_X_CONVSOL
			, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
			);

private:
	doublereal MakeTest(const VectorHandler& Vec);
};




#endif /* NONOLIN_H */
