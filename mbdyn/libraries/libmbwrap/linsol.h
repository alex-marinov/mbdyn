/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* linear solver generico */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

/* Integrator - begin */

class LinSol
{
public:
   	enum SolverType {
		HARWELL_SOLVER = 0,
		MESCHACH_SOLVER,
		Y12_SOLVER,
                UMFPACK_SOLVER,
                UMFPACK_CC_SOLVER,
		EMPTY_SOLVER,

		LAST_SOLVER
	};

protected:
	SolverType CurrSolver;
	integer iWorkSpaceSize;
   	doublereal dPivotFactor;

public:
	static SolverType defaultSolver;

	LinSol(void);
   	virtual ~LinSol(void);
	void Read(HighParser &HP, bool bAllowEmpty = false);
 
	SolverType GetSolver(void) const;
	integer iGetWorkSpaceSize(void) const;
	const doublereal& dGetPivotFactor(void) const;

	SolutionManager *const
	LinSol::GetSolutionManager(integer iNLD, integer iLWS = 0) const;
};

extern const char *psSolverNames[];

/* Integrator - end */

#endif /* LINEARSOLVER_H */

