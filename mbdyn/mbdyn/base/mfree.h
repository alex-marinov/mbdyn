/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
  * Copyright (C) 2004
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */

#ifndef MFREE_H
#define MFREE_H

#include <nonlin.h>
#include <precond.h>

class MatrixFreeSolver : public NonlinearSolver
{
public:
	enum SolverType {
		UNKNOWN = -1,

		BICGSTAB,
		GMRES,

		DEFAULT = GMRES,

		LAST
	};	

protected:
	Preconditioner* pPM;
	VectorHandler* 	pRes;
	doublereal IterTol;
	integer MaxLinIt;
	doublereal Tau;
	doublereal gamma;
	doublereal etaMax; 
	integer PrecondIter; 
	bool bBuildMat;
	bool honorJacRequest;
	const NonlinearProblem* pPrevNLP;
	
public:
	MatrixFreeSolver(const Preconditioner::PrecondType PType, 
			const integer iPStep,
			doublereal ITol,
			integer MaxIt,
			doublereal etaMx,
			doublereal T,
			bool JacReq = false); 

	~MatrixFreeSolver(void);
};

#endif /* MFREE_H */

