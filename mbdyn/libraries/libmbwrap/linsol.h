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

/* linear solver generico */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

/* Integrator - begin */

class LinSol {
public:
   	enum SolverType {
		HARWELL_SOLVER = 0,
		MESCHACH_SOLVER,
		Y12_SOLVER,
                UMFPACK_SOLVER,
		SUPERLU_SOLVER,
		LAPACK_SOLVER,
		TAUCS_SOLVER,
		NAIVE_SOLVER,
		EMPTY_SOLVER,

		LAST_SOLVER
	};

	enum {
		SOLVER_FLAGS_NONE = 0x00U,
		SOLVER_FLAGS_ALLOWS_MAP = 0x01U,
		SOLVER_FLAGS_ALLOWS_CC = 0x02U,
		SOLVER_FLAGS_ALLOWS_DIR = 0x04U,
		SOLVER_FLAGS_TYPE_MASK = SOLVER_FLAGS_ALLOWS_MAP|SOLVER_FLAGS_ALLOWS_CC|SOLVER_FLAGS_ALLOWS_DIR,
		SOLVER_FLAGS_ALLOWS_MT = 0x10U,
		SOLVER_FLAGS_ALLOWS_COLAMD = 0x20U
	};
	
protected:
	SolverType CurrSolver;
	unsigned solverFlags;
	/*
	 * number of threads
	 * currently used by:
	 * 	SparseLU
	 */
	unsigned nThreads;
	/*
	 * max workspace size
	 * currently used by:
	 * 	Y12
	 */
	integer iWorkSpaceSize;
	/*
	 * block size
	 * currently used by:
	 * 	Umfpack
	 */
	unsigned blockSize;
	/*
	 * pivot factor
	 * currently used by:
	 * 	SparseLU	0.0 -> 1.0
	 * 	Y12		0.0=no, 1.0=full
	 * 	Umfpack		0.0 -> 1.0
	 */
   	doublereal dPivotFactor;

public:
	static SolverType defaultSolver;

	LinSol(void);
   	virtual ~LinSol(void);
	void Read(HighParser &HP, bool bAllowEmpty = false);
 
	SolverType GetSolver(void) const;
	unsigned GetSolverFlags(void) const;
	unsigned GetSolverFlags(SolverType t) const;
	const char *const GetSolverName(void) const;
	const char *const GetSolverName(SolverType t) const;
	bool SetSolver(SolverType t, unsigned f = SOLVER_FLAGS_NONE);
	bool SetSolverFlags(unsigned f);
	bool AddSolverFlags(unsigned f);
	bool MaskSolverFlags(unsigned f);
	bool SetNumThreads(unsigned nt);
	integer iGetWorkSpaceSize(void) const;
	const doublereal& dGetPivotFactor(void) const;

	SolutionManager *const
	GetSolutionManager(integer iNLD, integer iLWS = 0) const;
};

/* Integrator - end */

#endif /* LINEARSOLVER_H */

