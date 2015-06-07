/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 *                            SOLUTION MANAGER                               *
 *                                                                           *
 *****************************************************************************/

/* Pierangelo Masarati */


#ifndef SOLMAN_H
#define SOLMAN_H

#include <cmath>
#include <iostream>
#include "ac/f2c.h"

/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "vh.h"
#include "mh.h"

#include "matvec3.h"

/* classi virtuali dichiarate in questo file */
class MatrixHandler;    /* gestore matrice */
class SparseMatrixHandler;
class VectorHandler;    /* gestore vettore */
class SolutionManager;  /* gestore della soluzione */

/* classi usate in questo file */
class SubMatrixHandler;
class FullSubMatrixHandler;
class SparseSubMatrixHandler;
class VariableSubMatrixHandler;
class SubVectorHandler;
class Vec3;
class Mat3x3;
class LinearSolver;

/* SolutionDataManager - begin */

class SolutionDataManager {
public:
	virtual ~SolutionDataManager(void);

	struct ChangedEquationStructure: public MBDynErrBase {
	public:
		ChangedEquationStructure(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	
	/* Collega il DataManager ed il DriveHandler ai vettori soluzione */
	virtual void
	LinkToSolution(VectorHandler& XCurr,
		VectorHandler& XPrimeCurr) = 0;

	/* Assembla il residuo */
	virtual void AssRes(VectorHandler& ResHdl, doublereal dCoef) = 0;
};

/* SolutionDataManager - end */


/* SolutionManager - begin */

class SolutionManager {
public:
	// whether matrix scaling should be performed, and when, using dgeequ()
	enum ScaleWhen {
		SCALEW_NEVER = 0,
		SCALEW_ONCE,
		SCALEW_ALWAYS
	};

	enum ScaleAlgorithm {
		SCALEA_NONE,
		SCALEA_UNDEF,
		SCALEA_ROW_MAX,
		SCALEA_ROW_SUM,
		SCALEA_COL_MAX,
		SCALEA_COL_SUM,
		SCALEA_LAPACK,
		SCALEA_ITERATIVE,
		SCALEA_ROW_MAX_COL_MAX
	};

	enum ScaleFlags {
		SCALEF_DEFAULT = 0x0u,
		SCALEF_WARN    = 0x1u,
		SCALEF_VERBOSE = 0x2u,
		SCALEF_COND_NUM_1 = 0x4u,
		SCALEF_COND_NUM_INF = 0x8u,
		SCALEF_COND_NUM = SCALEF_COND_NUM_1 | SCALEF_COND_NUM_INF
	};

	struct ScaleOpt {
		ScaleOpt(ScaleWhen when = SCALEW_NEVER,
				 ScaleAlgorithm alg = SCALEA_UNDEF,
				 integer iMaxIter = 100,
				 doublereal dTol = sqrt(std::numeric_limits<doublereal>::epsilon()),
				 unsigned flags = SCALEF_DEFAULT):
			when(when),
			algorithm(alg),
			iMaxIter(iMaxIter),
			dTol(dTol),
			uFlags(flags)
		{

		}

		ScaleWhen when;
		ScaleAlgorithm algorithm;
		integer iMaxIter;
		doublereal dTol;
		unsigned uFlags;
	};
protected:
	LinearSolver *pLS;

public:
	SolutionManager(void);

	/* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
	virtual ~SolutionManager(void);

	virtual void
	LinkToSolution(VectorHandler& XCurr,
		VectorHandler& XPrimeCurr);

#ifdef DEBUG
	virtual void IsValid(void) const = 0;
#endif /* DEBUG */

	/* Inizializzatore generico */
	virtual void MatrReset(void) = 0;

	/* Inizializzatore "speciale" */
	virtual void MatrInitialize(void);

	/* Risolve il sistema */
	virtual void Solve(void) = 0;
	virtual void SolveT(void);

	/* Rende disponibile l'handler per la matrice */
	virtual MatrixHandler* pMatHdl(void) const = 0;

	/* Rende disponibile l'handler per il termine noto */
	virtual VectorHandler* pResHdl(void) const = 0;

	/* Rende disponibile l'handler per la soluzione (e' lo stesso
	 * del termine noto, ma concettualmente sono separati) */
	virtual VectorHandler* pSolHdl(void) const = 0;

   	/* sposta il puntatore al vettore del residuo */
   	doublereal *pdSetResVec(doublereal* pd);
   
   	/* sposta il puntatore al vettore della soluzione */
   	doublereal *pdSetSolVec(doublereal* pd);

   	/* return true if the condition number is available */
   	bool bGetConditionNumber(doublereal& dCond) const;
};

/* SolutionManager - end */

#endif /* SOLMAN_H */

