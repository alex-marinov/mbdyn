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

 /*****************************************************************************
 *                                                                           *
 *                    BiCGStab ITERATIVE SOLUTION MANAGER                    *
 *                                                                           *
 *****************************************************************************/

/* Author: Giuseppe Quaranta <quaranta@aero.polimi.it> */

#ifndef BiCGStabSolutionManager_hh
#define BiCGStabSolutionManager_hh

#include <ac/iostream>

#include <myassert.h>
#include <mynewmem.h>
#include <solman.h>
#include <spmapmh.h>
class Preconditioner;

class IterativeSolutionManager: public SolutionManager {
private:
	
	integer Size; 
	SolutionManager* pSM;
	SolutionDataManager* pDM;
	Preconditioner* pP;         /* Preconditioner class */
	MatrixHandler*  pA;
	VectorHandler *pxVH, *pbVH;
	VectorHandler *pXCurr, *pXPrimeCurr;
	
	doublereal Toll;
	int MaxIt;
        doublereal Tau;
        
	void ComputeTau();

public:
	template <class S> IterativeSolutionManager(integer iDim,
	                         doublereal Tollerance,
				 int max,
				 SolutionDataManager* pDM,
				 VectorHandler* pX,
				 VectorHandler* pXP, 			 
				 S* SM,
				 integer WorkSpaceSize = 0,
				 const doublereal dPivot = 1.0);
	
	~IterativeSolutionManager(void);
	
	virtual void IsValid(void) const {
		NO_OP;
	};

	/* Inizializzatore generico */
	virtual void MatrInit(const doublereal& d = 0.);
	
	/* Risolve il sistema */
	virtual void Solve(const doublereal dCoef);

   	/* sposta il puntatore al vettore del residuo */
   	virtual void ChangeResPoint(doublereal* pRes) {
		std::cerr << "IterativeSolutionManager::ChangeResPoint: "
			<< "you should not be here !!"
			<< "Aborting..." << std::endl;
		THROW(ErrGeneric());
	};

  
   	/* sposta il puntatore al vettore del residuo */
   	virtual void ChangeSolPoint(doublereal* pSol) {
		std::cerr << "IterativeSolutionManager::ChangeSolPoint: "
			<< "you should not be here !!"
			<< "Aborting..." << std::endl;
		THROW(ErrGeneric());
	};


	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const {
		return pA;
	};

	/* Rende disponibile l'handler per il termine noto */
	virtual VectorHandler* pResHdl(void) const {
		return pbVH;
	};

	/* Rende disponibile l'handler per la soluzione */
	virtual VectorHandler* pSolHdl(void) const {
		return pxVH;
	};

};

class Preconditioner {
private:
       MatrixHandler* pA;
       SolutionManager* pSM;
       int Size;
public:

      	Preconditioner(integer iDim, MatrixHandler *pa, 
		SolutionManager* psm)
	:Size(iDim),
	pA(pa),
	pSM(psm) 
	{};
       	
	~Preconditioner(void){
		NO_OP;
	};
	
        void solve(VectorHandler &bVH, VectorHandler &xVH) {
		MyVectorHandler TmpVH(Size);
		pSM->ChangeResPoint(bVH.pdGetVec());
		pSM->ChangeSolPoint(xVH.pdGetVec());
		pSM->Solve(0.);
	};
};

#endif /* BiCGStabSolutionManager_hh */
