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

#ifndef GMRES_H
#define GMRES_H

#include <mfree.h>

class UpHessMatrix 
{
	std::vector<doublereal> M;
	integer Size;

public:	
	UpHessMatrix(integer n);
	~UpHessMatrix(void);

	void Reset(doublereal d = 0.);

	doublereal& operator() (const integer i, const integer j);
	doublereal operator() (const integer i, const integer j) const;
};

class Gmres : public MatrixFreeSolver
{
public:
	Gmres(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt,
			doublereal etaMx);
	
	~Gmres(void);
	
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
			doublereal &cs, doublereal &sn) const;
			
	void ApplyPlaneRotation(doublereal &dx, doublereal &dy, 
			const doublereal &cs, const doublereal &sn) const;
			
	void Backsolve(VectorHandler& x, integer sz,  UpHessMatrix& H, 
			VectorHandler& s, MyVectorHandler* v);
};

#endif /* GMRES_H */

