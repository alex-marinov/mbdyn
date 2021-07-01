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
  * interfaccia dello Step integrator con il nonlinear solver 
  */
  
#ifndef NONLINPB_H
#define NONLINPB_H

#include <solman.h>

class NonlinearSolverTest;

class NonlinearProblem
{
public:
	/* Distruttore */
   	virtual ~NonlinearProblem(void) { };

	virtual void Residual(VectorHandler* pRes, VectorHandler* pAbsRes=0) const = 0;

	virtual void Jacobian(MatrixHandler* pJac) const = 0;
	
	virtual void Update(const VectorHandler* pSol) const = 0;

	/* scale factor for tests */
	virtual doublereal TestScale(const NonlinearSolverTest *pTest,
								 doublereal& dAlgebraicEquations) const = 0;
	
	virtual void EvalProd(doublereal Tau, const VectorHandler& f0,
			const VectorHandler& w, VectorHandler& z) const = 0;
};   

#endif /* NONLINPB_H */

