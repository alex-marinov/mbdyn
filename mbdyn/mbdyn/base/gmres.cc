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
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
 
#include <gmres.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <umfpackwrap.h>
#include <unistd.h>
#include <output.h>

UpHessMatrix::UpHessMatrix(integer n)
: M(n*n+1), Size(n)
{
	NO_OP;
}

UpHessMatrix::~UpHessMatrix(void)
{
	NO_OP;
}
	
void
UpHessMatrix::Reset(doublereal d)
{
	for (unsigned int i = 0; i < M.size(); i++) {
		M[i] = 0;
	}
}

doublereal&
UpHessMatrix::operator() (const integer i, const integer j)
{
	return M[i*Size+j];
}
	
doublereal
UpHessMatrix::operator() (const integer i, const integer j) const
{
	return M[i*Size+j];
}

Gmres::Gmres(const Preconditioner::PrecondType PType, 
		const integer iPStep,
		doublereal ITol,
		integer MaxIt,
		doublereal etaMx) 
: MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx)
{
	NO_OP;
}
	
Gmres::~Gmres(void)
{
	NO_OP;
}
	
void
Gmres::GeneratePlaneRotation(const doublereal &dx, const doublereal &dy, 
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
}

void
Gmres::ApplyPlaneRotation(doublereal &dx, doublereal &dy, 
		const doublereal &cs, const doublereal &sn) const 
{ 
	doublereal temp = cs * dx + sn * dy; 
	dy = -sn * dx + cs * dy;
	dx = temp;
}

void
Gmres::Backsolve(VectorHandler& x, integer sz,  UpHessMatrix& H, 
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
}

void
Gmres::Solve(const NonlinearProblem* pNLP,
		SolutionManager* pSolMan,
		const integer iMaxIter,
		const doublereal Tol,
		const doublereal SolTol,
		integer& iIterCnt,
		doublereal& dErr
#ifdef MBDYN_X_CONVSOL
		, doublereal& dSolErr
#endif /* MBDYN_X_CONVSOL  */	
	    )
{
	ASSERT(pNLP != NULL);
	ASSERT(pSM != NULL);
	
	iIterCnt = 0;

#ifdef MBDYN_X_CONVSOL
	dSolErr = 0.;
#endif /* MBDYN_X_CONVSOL  */	

	/* external nonlinear iteration */	
	
	/* riassembla sempre lo jacobiano se l'integratore e' nuovo */
	if (pNLP != pPrevNLP) {
		fBuildMat = true;
	}
	
	pPrevNLP = pNLP;
	
	pSM  = pSolMan;
        pRes = pSM->pResHdl();
	Size = pRes->iGetSize();

	doublereal eta = etaMax;
	doublereal rateo = 0.;
	doublereal Fnorm = 1.;
        doublereal resid;
	doublereal norm1, norm2;
	VectorHandler* pr;
	MyVectorHandler s(MaxLinIt+1), cs(MaxLinIt+1), sn(MaxLinIt+1), w(Size);
	MyVectorHandler* v = NULL;
	SAFENEWARR(v, MyVectorHandler, MaxLinIt+1); 
	MyVectorHandler vHat(Size); 
	MyVectorHandler dx(Size); 
	UpHessMatrix H(MaxLinIt+1);
	integer TotalIter = 0;

#ifdef DEBUG_ITERATIVE
	std::cout << " Gmres New Step " <<std::endl; 			
#endif /* DEBUG_ITERATIVE */

	while (true) {

#ifdef 	USE_EXTERNAL 	
		SendExternal();
#endif /* USE_EXTERNAL */
		
		pRes->Reset(0.);
      		pNLP->Residual(pRes);
		
      		if (foutRes) {
	 		std::cout << "Residual:" << std::endl;
	 		std::cout << iIterCnt <<std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< pRes->dGetCoef(iTmpCnt) << std::endl;
			}
      		}

		dErr = MakeTest(*pRes);		

#ifdef DEBUG_ITERATIVE
      		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (dErr < Tol) {
	 		return;
      		}
      		if (!isfinite(dErr)) {
			THROW(ErrSimulationDiverged());
		}
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
		}
          	rateo = dErr*dErr/Fnorm;

#ifdef DEBUG_ITERATIVE
		std::cerr << "rateo " << rateo << std::endl;
#endif /* DEBUG_ITERATIVE */

		Fnorm = dErr*dErr;
		
      		iIterCnt++;

		/* inner iteration to solve the linear system */	
	
		/* Gmres(m) Iterative solver */
		DEBUGCOUT("Using Gmres(m) iterative solver" << std::endl);

        	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */
        	/* N.B. *pRes = -F(0) */ 
		
		pr = pRes;
        	doublereal LocTol = eta * dErr;
 
#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocTol " << LocTol << std::endl;
#endif /* DEBUG_ITERATIVE */
	
		resid = dErr;
		dx.Reset(0.);
		
		if (fBuildMat) {
			pSM->MatrInit(0.);
			pNLP->Jacobian(pSM->pMatHdl());
			fBuildMat = false;
			TotalIter = 0;
			TotJac++;
			
#ifdef DEBUG_ITERATIVE
			std::cerr << "Jacobian " << std::endl;
#endif /* DEBUG_ITERATIVE */

		}

		int i = 0;
        	v[0].Resize(Size);
		v[0].ScalarMul(*pr, 1./resid);
		s.Reset(0.);
		cs.Reset(0.);
		sn.Reset(0.);
		s.fPutCoef(1, resid);
		while ((i < MaxLinIt)) {

#ifdef DEBUG_ITERATIVE
			std::cout << "v[i]:" << i << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< v[i].dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pPM->Precond(v[i], vHat, pSM); 

#ifdef DEBUG_ITERATIVE
			std::cout << "vHat:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< vHat.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			pNLP->EvalProd(Tau, *pr, vHat, w);

#ifdef DEBUG_ITERATIVE
			std::cout << "w:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< w.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			norm1 = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm1: " << norm1 << std::endl; 
#endif /* DEBUG_ITERATIVE */

			for (int k = 0; k <= i; k++) {
        			H(k, i) = w.InnerProd(v[k]);

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k,i): " << k << " " << i << " "<< H(k,i) << std::endl; 
#endif /* DEBUG_ITERATIVE */

				w.ScalarAddMul(v[k], -H(k, i));
      			}
			H(i+1, i) = w.Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "w.Norm(): " << w.Norm() << std::endl;
			std::cerr << "H(i+1, i): " << i+1 << " " << i << " " << H(i+1, i) << std::endl;
#endif /* DEBUG_ITERATIVE */

			norm2 = H(i+1, i); 

#ifdef DEBUG_ITERATIVE
			std::cerr << "norm2: " << norm2 << std::endl;
#endif /* DEBUG_ITERATIVE */

			v[i+1].Resize(Size);
			
			/*  Reorthogonalize? */
			if  (.001*norm2/norm1 < DBL_EPSILON) {
    				for (int k = 0;  k <= i; k++) {       
					doublereal hr = v[k].InnerProd(w);
        				H(k,i) = H(k,i)+hr;
        				w.ScalarAddMul(v[k], -hr);

#ifdef DEBUG_ITERATIVE
					std::cerr << "Reorthogonalize: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				}
				H(i+1, i) = w.Norm();
			}
			if (fabs(H(i+1, i)) > DBL_EPSILON) { 
				v[i+1].ScalarMul(w, 1./H(i+1, i));

#ifdef DEBUG_ITERATIVE
				std::cout << "v[i+1]:" << i+1 << std::endl;
#endif /* DEBUG_ITERATIVE */

	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< v[i+1].dGetCoef(iTmpCnt) << std::endl;
				}
			} else {
				/* happy breakdown !! */

#ifdef DEBUG_ITERATIVE
				std::cerr << "happy breakdown: "  << std::endl;
#endif /* DEBUG_ITERATIVE */

				v[i+1].Reset(0.);
			} 
			for (int k = 0; k < i; k++) {
        			ApplyPlaneRotation(H(k,i), H(k+1,i), cs.dGetCoef(k+1), sn.dGetCoef(k+1));

#ifdef DEBUG_ITERATIVE
				std::cerr << "H(k, i): " << k << " " << i << " " << H(k, i) << std::endl;
				std::cerr << "H(k+1, i): " << k+1 << " " << i << " " << H(k+1, i) << std::endl;
#endif /* DEBUG_ITERATIVE */

			}
							
			GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i+1), sn(i+1));

#ifdef DEBUG_ITERATIVE
			std::cerr << "cs(i): " << cs.dGetCoef(i+1) << std::endl;
			std::cerr << "sn(i): " << sn.dGetCoef(i+1) << std::endl;
#endif /* DEBUG_ITERATIVE */

      			ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i+1), sn(i+1));
      			ApplyPlaneRotation(s(i+1), s(i+2), cs(i+1), sn(i+1));
				
			if ((resid = fabs(s.dGetCoef(i+2))) < LocTol) {

#ifdef DEBUG_ITERATIVE
				std::cerr << "resid 1: " << resid  << std::endl;
#endif /* DEBUG_ITERATIVE */

				Backsolve(dx, i, H, s, v);
				pPM->Precond(dx, dx, pSM); 

#ifdef DEBUG_ITERATIVE
				std::cout << "dx:" << std::endl;
	 			for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    				std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
						<< dx.dGetCoef(iTmpCnt) << std::endl;
				}
#endif /* DEBUG_ITERATIVE */

        			break;
      			}

#ifdef DEBUG_ITERATIVE
			std::cerr << "resid 2 : " << resid << std::endl;
#endif /* DEBUG_ITERATIVE */

			TotalIter++;
			i++;
		}
		if (i == MaxLinIt) {
			Backsolve(dx, MaxLinIt, H, s, v);
			pPM->Precond(dx, dx, pSM);				
		        std::cerr << "Iterative inner solver didn't converge."
				<< " Continuing..." << std::endl;
		}

		/* se ha impiegato troppi passi riassembla lo jacobiano */
		
		if (TotalIter >= PrecondIter) {
			fBuildMat = true;
		}
		/* calcola il nuovo eta */
		doublereal etaNew = gamma * rateo;
		doublereal etaBis;
		if ((etaBis = gamma*eta*eta) > .1) {
			etaNew = (etaNew > etaBis) ? etaNew : etaBis;
		}
		eta = (etaNew < etaMax) ? etaNew : etaMax;
		/* prevent oversolving */
		etaBis = .5*Tol/dErr;
		eta = (eta > etaBis) ? eta : etaBis;

#ifdef DEBUG_ITERATIVE
		std::cerr << "eta " << eta << std::endl;
#endif /* DEBUG_ITERATIVE */
		
		if (foutSol) {      
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": "
					<< dx.dGetCoef(iTmpCnt) << std::endl;
			}
		}		
		
		
		if (foutIters) {
#ifdef USE_MPI
			if (dynamic_cast<SchurSolutionManager*> (pSM) && (MBDynComm.Get_rank() == 0)) {
#endif /* USE_MPI */
				std::cerr << "\tIteration " << iIterCnt
					<< " " << dErr << " J"
					<< std::endl;
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		
      		pNLP->Update(&dx);

		
#ifdef MBDYN_X_CONVSOL
		if (SolTol > 0.) {
			dSolErr = MakeTest(dx);
        		if (dSolErr < SolTol) {
				THROW(ConvergenceOnSolution());
			}
      		}
#endif /* MBDYN_X_CONVSOL */
	}

	SAFEDELETEARR(v);
}

