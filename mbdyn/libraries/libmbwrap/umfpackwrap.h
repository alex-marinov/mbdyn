/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2002
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* December 2001 
 * Modified to add a Sparse matrix in row form and to implement methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 1996-2002
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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

#ifndef UmfpackSparseLUSolutionManager_hh
#define UmfpackSparseLUSolutionManager_hh

#ifdef USE_UMFPACK

#include <ac/iostream>
#include <vector>

extern "C" {
#include <umfpack.h>
}

#include <myassert.h>
#include <mynewmem.h>
#include <solman.h>
#include <spmapmh.h>

#warning "Umfpack3SparseLUSolutionManager => UmfpackSparseLUSolutionManager"
#define Umfpack3SparseLUSolutionManager UmfpackSparseLUSolutionManager

class UmfpackSparseLUSolutionManager: public SolutionManager {
private:
	mutable SpMapMatrixHandler A;
	MyVectorHandler *xVH, *bVH;
	std::vector<double> x;
	std::vector<double> b;
	std::vector<double> Ax;
	std::vector<int> Ai;
	std::vector<int> Ap;
	doublereal* pdRhs;
	doublereal* pdSol;
	

	void * Symbolic;
	double Control[UMFPACK_CONTROL];
	double Info[UMFPACK_INFO];
	void * Numeric;
	
	bool HasBeenReset;

	bool PrepareSymbolic(void);
	
public:
	UmfpackSparseLUSolutionManager(integer Dim, integer /* unused */ = 0, 
			doublereal dPivot = -1.);
	virtual ~UmfpackSparseLUSolutionManager(void);
	virtual void IsValid(void) const {
		NO_OP;
	};

	/* Inizializzatore generico */
	virtual void MatrInit(const doublereal& d = 0.);
	
	/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
	virtual void Solve(const doublereal /* dCoef */ = 0.);

	/* Bacward Substitution */
	void BackSub(doublereal t_iniz = 0.);
   
   	/* sposta il puntatore al vettore del residuo */
   	void ChangeResPoint(doublereal* pRes){
		pdRhs = pRes;
	};
   
   	/* sposta il puntatore al vettore del residuo */
   	void ChangeSolPoint(doublereal* pSol){
		pdSol = pSol;
	};
	
	/* Rende disponibile l'handler per la matrice */
	virtual SpMapMatrixHandler* pMatHdl(void) const;

	/* Rende disponibile l'handler per il termine noto */
	virtual MyVectorHandler* pResHdl(void) const;

	/* Rende disponibile l'handler per la soluzione */
	virtual MyVectorHandler* pSolHdl(void) const;
};

#endif /* USE_UMFPACK */

#endif /* UmfpackSparseLUSolutionManager_hh */

