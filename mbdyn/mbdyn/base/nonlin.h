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
#include <ac/float.h>
#include<vector>

class NonlinearSolver
{

public:
 	class ErrSimulationDiverged{};
 	class NoConvergence{};
	class ConvergenceOnSolution{};
	class ErrGeneric{};
protected:
	integer Size;
	integer TotJac;	
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
		TotJac(0),
		foutIters(false),
		foutRes(false),
		foutJac(false),
		foutSol(false)
#ifdef USE_EXTERNAL
		, ExtStepType(External::ERROR)  
#endif /* USE_EXTERNAL */
#ifdef __HACK_SCALE_RES__
		, pScale(NULL) 
#endif /* __HACK_SCALE_RES__ */
		{ };
		
#ifdef __HACK_SCALE_RES__
	virtual void SetScale(const VectorHandler* pScl) {
		pScale = (VectorHandler *)pScl;
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

	virtual integer TotalAssembledJacobian(void) {
		return TotJac;
	};
	
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
	
	virtual void Precond(VectorHandler& b,
			VectorHandler& x, 
			SolutionManager* pSM) const = 0;
};


class FullJacobianPr : public Preconditioner
{	

public:
	FullJacobianPr(void) {};
	
	~FullJacobianPr(void) {};
	
	void Precond(VectorHandler& b,
			VectorHandler& x, 
			SolutionManager* pSM) const {
		
		if (pSM->pSolHdl() != pSM->pResHdl()) {	
			pSM->ChangeResPoint(b.pdGetVec());
			pSM->ChangeSolPoint(x.pdGetVec());
			pSM->Solve();
		} else {
			x = b;
			pSM->ChangeResPoint(x.pdGetVec());
			pSM->Solve();
		}
	};
	 
};

const doublereal defaultTau = sqrt(DBL_EPSILON);
const doublereal defaultGamma = 0.9;

class MatrixFreeSolver : public NonlinearSolver
{

public:
	enum SolverType {
		BICGSTAB,
		GMRES
	};	

protected:
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
	const NonlinearProblem* pPrevNLP;
	
public:
	MatrixFreeSolver(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt,
			doublereal etaMx); 

	~MatrixFreeSolver(void) { };
	
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

protected:
	doublereal MakeTest(const VectorHandler& Vec);
};


class BiCGStab : public MatrixFreeSolver
{

public:

	BiCGStab(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt,
			doublereal etaMx) 
	: MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx) {};
	
	~BiCGStab(void) { };
	
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
			);
};


class UpHessMatrix 
{
	
	std::vector<doublereal> M;
	integer Size;

public:	
	UpHessMatrix(integer n) : M(n*n+1), Size(n) {};
	void Reset(doublereal d = 0.) {
		for (unsigned int i=0; i < M.size(); i++) M[i] = 0;
	};
	doublereal& operator() (const integer i, const integer j) {
		return M[i*Size+j];
	};
	
	doublereal operator() (const integer i, const integer j) const {
		return M[i*Size+j];
	};

};


class Gmres : public MatrixFreeSolver
{

public:

	Gmres(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt,
			doublereal etaMx) 
	:MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx) {};
	
	~Gmres(void) { };
	
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
			);
			
private:

	void GeneratePlaneRotation(const doublereal &dx, const doublereal &dy, 
			doublereal &cs, doublereal &sn) const 
	{
		if (fabs(dy) < DBL_EPSILON) {
    			cs = 1.0;
    			sn = 0.0;
  		} else if (fabs(dy) > fabs(dx)) {
    			doublereal temp = dx / dy; 
    			sn = 1.0 / sqrt( 1.0 + temp*temp );
    			cs = temp * sn;
  		} else {
    			doublereal temp = dy / dx; 
    			cs = 1.0 / sqrt( 1.0 + temp*temp );
    			sn = temp * cs;
  		}
	};
			
	void ApplyPlaneRotation(doublereal &dx, doublereal &dy, 
			const doublereal &cs, const doublereal &sn) const 
	{ 
  		doublereal temp = cs * dx + sn * dy; 
  		dy = -sn * dx + cs * dy;
		dx  = temp;
		return;  
	};
			
	
	void Backsolve(VectorHandler& x, integer sz,  UpHessMatrix& H, 
			VectorHandler& s, MyVectorHandler* v) 
	{ 
 
		for (int i = sz; i >= 0; i--) {
    			s.fPutCoef(i+1, s.dGetCoef(i+1) / H(i,i));
    			for (int j = i - 1; j >= 0; j--)
      				s.fDecCoef(j+1, H(j,i) * s.dGetCoef(i+1));
  		}

  		for (int j = 0; j <= sz; j++) {
    			x.ScalarAddMul(v[j], s.dGetCoef(j+1));
		}
		return;
	};

					
};






#endif /* NONOLIN_H */
