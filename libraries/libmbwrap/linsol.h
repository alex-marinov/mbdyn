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

/* linear solver generico */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "solman.h"
#ifdef USE_MPI
#include "mbcomm.h"
#endif
/* Integrator - begin */

class LinSol {
public:
   	enum SolverType {
		EMPTY_SOLVER = 0,
		HARWELL_SOLVER,
		LAPACK_SOLVER,
		NAIVE_SOLVER,
		SUPERLU_SOLVER,
		TAUCS_SOLVER,
                UMFPACK_SOLVER,
                KLU_SOLVER,
		Y12_SOLVER,
                PARDISO_SOLVER,
                PARDISO_64_SOLVER,
                PASTIX_SOLVER,
                QR_SOLVER,
                SPQR_SOLVER,
		STRUMPACK_SOLVER,
		WATSON_SOLVER,
                AZTECOO_SOLVER,
                AMESOS_SOLVER,
                SICONOS_SPARSE_SOLVER,
                SICONOS_DENSE_SOLVER,
		LAST_SOLVER
	};

     enum SolverFlags: unsigned {
		SOLVER_FLAGS_NONE = 0x00U,
		SOLVER_FLAGS_ALLOWS_MAP = 0x01U,
		SOLVER_FLAGS_ALLOWS_CC = 0x02U,
		SOLVER_FLAGS_ALLOWS_DIR = 0x04U,
		SOLVER_FLAGS_ALLOWS_GRAD = 0x08U,
		SOLVER_FLAGS_TYPE_MASK = SOLVER_FLAGS_ALLOWS_MAP|SOLVER_FLAGS_ALLOWS_CC|SOLVER_FLAGS_ALLOWS_DIR|SOLVER_FLAGS_ALLOWS_GRAD,
		SOLVER_FLAGS_ALLOWS_MT_FCT = 0x10U,
		SOLVER_FLAGS_ALLOWS_MT_ASS = 0x20U,
		//permutations
		SOLVER_FLAGS_ALLOWS_COLAMD = 0x40U,
		SOLVER_FLAGS_ALLOWS_MMDATA = 0x80U,
		SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE = 0x100U,
		SOLVER_FLAGS_ALLOWS_KING = 0x200U,
		SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION = 0x400U,
		SOLVER_FLAGS_ALLOWS_MDAPLUSAT = 0x800U,
		SOLVER_FLAGS_ALLOWS_SLOAN = 0x1000U,
                SOLVER_FLAGS_ALLOWS_AMD = 0x2000U,
                SOLVER_FLAGS_ALLOWS_METIS = 0x4000U,
                SOLVER_FLAGS_ALLOWS_GIVEN = 0x8000U,
		SOLVER_FLAGS_ALLOWS_SCOTCH = 0x10000U,
		SOLVER_FLAGS_PERM_MASK = 
			SOLVER_FLAGS_ALLOWS_COLAMD |
			SOLVER_FLAGS_ALLOWS_MMDATA |
			SOLVER_FLAGS_ALLOWS_MDAPLUSAT |
			SOLVER_FLAGS_ALLOWS_REVERSE_CUTHILL_MC_KEE |
			SOLVER_FLAGS_ALLOWS_KING |
			SOLVER_FLAGS_ALLOWS_SLOAN |
                        SOLVER_FLAGS_ALLOWS_NESTED_DISSECTION |
                        SOLVER_FLAGS_ALLOWS_AMD |
                        SOLVER_FLAGS_ALLOWS_METIS |
                        SOLVER_FLAGS_ALLOWS_GIVEN |
		        SOLVER_FLAGS_ALLOWS_SCOTCH,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_HSS = 0x20000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_BLR = 0x40000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_HODLR = 0x80000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_SVD = 0x100000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_PQRCP =  0x200000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRCP =  0x400000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_TQRCP =  0x800000U,
		SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRRT = 0x1000000U,
		SOLVER_FLAGS_COMPRESSION_MASK =
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_HSS |
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_BLR | 
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_HODLR | 
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_SVD |
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_PQRCP | 
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRCP |
		        SOLVER_FLAGS_ALLOWS_COMPRESSION_TQRCP |
                        SOLVER_FLAGS_ALLOWS_COMPRESSION_RQRRT,
                SOLVER_FLAGS_ALLOWS_PRECOND_LAPACK    = 0x10000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_UMFPACK   = 0x20000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_KLU       = 0x30000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_ILUT      = 0x40000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_SUPERLU   = 0x50000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_MUMPS     = 0x60000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_SCALAPACK = 0x70000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_DSCPACK   = 0x80000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_PARDISO   = 0x90000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_PARAKLETE = 0xA0000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_TAUCS     = 0xB0000000U,
                SOLVER_FLAGS_ALLOWS_PRECOND_CSPARSE   = 0xC0000000U,
                SOLVER_FLAGS_PRECOND_MASK             = 0xF0000000U
	};

	/* solver data */
	struct solver_t {
		const char *const	s_name;
		const char *const	s_alias;
		enum SolverType		s_type;
		unsigned		s_flags;
		unsigned		s_default_flags;
		doublereal		s_pivot_factor;
		doublereal		s_drop_tolerance;
	};

protected:
	SolverType currSolver;
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
	 * 	Naive		1.0 -> 0.0
	 */
   	doublereal dPivotFactor;

	/*
	 * drop tolerance
	 * currently used by:
	 *	Umfpack		>= 0.0 (0.0: use all)
	 */
	doublereal dDropTolerance;

	/*
	 * matrix scaling
	 * currently used by:
	 *	Naive, KLU, Umfpack
	 */
	SolutionManager::ScaleOpt scale;

        doublereal dLowRankCompressTol; /* Tolerance for low-rank kernels (PaStiX) */
        doublereal dLowRankCompressMinRatio; /* Min ratio for rank w.r.t. strict rank (PaStiX) */
	/*
	 * maximum number of iterations for iterative refinement
	 * used only by:
	 *  Umfpack
	 */
	integer iMaxIter;
        doublereal dTolRes; // Used by AztecOO iterative linear solver
        integer iVerbose;
public:
	static SolverType defaultSolver;

	LinSol(void);
   	virtual ~LinSol(void);
 
	SolverType GetSolver(void) const;
	const char *const GetSolverName(void) const;
	unsigned GetSolverFlags(void) const;
	unsigned GetNumThreads(void) const;
	integer iGetWorkSpaceSize(void) const;
	const doublereal& dGetPivotFactor(void) const;
	const doublereal& dGetDropTolerance(void) const;
	unsigned GetBlockSize(void) const;
	const SolutionManager::ScaleOpt& GetScale(void) const{ return scale; }
	integer GetMaxIterations(void) const;

	const char *const GetSolverName(SolverType t) const;
	unsigned GetSolverFlags(SolverType t) const;

	bool SetSolver(SolverType t, unsigned f = SOLVER_FLAGS_NONE);
	bool SetSolverFlags(unsigned f);
        bool AddSolverFlags(unsigned mask, unsigned flag);
	bool SetNumThreads(unsigned nt);
	bool SetWorkSpaceSize(integer);
	bool SetPivotFactor(const doublereal &d);
	bool SetDropTolerance(const doublereal &d);
	bool SetBlockSize(unsigned bs);
	bool SetScale(const SolutionManager::ScaleOpt& scale);
        bool SetLowRankCompressTol(const doublereal& d);
	bool SetLowRankCompressMinRatio(const doublereal& d);
	bool SetMaxIterations(integer iMaxIter);
        bool SetTolerance(doublereal dToleranceRes);
        doublereal dGetTolerance() const { return dTolRes; }
        integer iGetMaxIterations() const { return iMaxIter; }
        bool SetVerbose(integer iVerb);
	SolutionManager *const
	GetSolutionManager(integer iNLD,
#ifdef USE_MPI
                           MPI::Intracomm& oComm,
#endif
                           integer iLWS = 0) const;
};

extern const LinSol::solver_t solver[];

/* LinSol - end */

#endif // LINEARSOLVER_H

