/* $Header$ */
/*
 * HmFe (C) is a FEM analysis code.
 *
 * Copyright (C) 1996-2017
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
 * Copyright (C) 2001-2017
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
 * Wsmp is used by permission; please read its Copyright,
 * License and Availability note.
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_WSMP
#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "wsmpwrap.h"

extern "C" {
     extern void wgsmp_(int *N,
			int *IA,
			int *JA,
			double *AVALS,
			double *B,
			int *LDB,
			int *NRHS,
			double *RMISC,
			int *IPARAM,
			double *DPARAM);
     extern void wsetmaxthrds_(int *);
     extern void wsmp_clear_();
}


/* WsmpSolver - begin */

WsmpSolver::WsmpSolver(const integer size, const doublereal dPivot,
		       const unsigned blockSize, const unsigned nt,
		       const SolutionManager::ScaleOpt& oScale, integer iMaxIterRefine)
     : LinearSolver(0),
       iSize(size),
       Axp(0),
       Aip(0),
       App(0),
       ldb(size),
       nrhs(1),
       rmisc(0),
       iNumNzPrev(-1),
       Symbolic(false)
{
     int tnt = nt;
     wsetmaxthrds_(&tnt);
     iparm[0] = 0;
     iparm[1] = 0;
     iparm[2] = 0;

     wgsmp_(&iSize, 0, 0, 0, 0, &ldb, &nrhs, rmisc, iparm, dparm);

     /* CSC format */
     iparm[3] = 1;
     /* 0 index */
     iparm[4] = 0;
     /* pivot */

     iparm[5] = iMaxIterRefine;
     iparm[6] = 1;

     switch (oScale.when) {
     case SolutionManager::SCALEW_NEVER:
	  iparm[7] = 2;
	  iparm[8] = -1;
	  break;
     default:
	  switch (oScale.algorithm) {
	  case SolutionManager::SCALEA_ROW_MAX:
	       iparm[8] = 1;
	       break;
	  case SolutionManager::SCALEA_COL_MAX:
	       iparm[8] = 3;
	       break;
	  case SolutionManager::SCALEA_ROW_MAX_COL_MAX:
	       iparm[8] = 2;
	       break;
	  case SolutionManager::SCALEA_ITERATIVE:
	       iparm[7] = 1;
	       break;
	  default:
	       ;
	  }
     }

     constexpr doublereal dTolRefine = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9);

     dparm[5] = dTolRefine;

     if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
	  /*
	   * 1.0: true partial pivoting
	   * 0.0: treated as 1.0
	   *
	   * default: 0.1
	   */
	  iparm[10] = 1;
	  dparm[10] = dPivot;
     }

     if (blockSize > 0) {
	  iparm[31] = blockSize;
     }
}

WsmpSolver::~WsmpSolver(void)
{
     wsmp_clear_();
}

void
WsmpSolver::Reset(void)
{
     bHasBeenReset = true;
}

void
WsmpSolver::Solve(void) const
{
     if (bHasBeenReset) {
	  const_cast<WsmpSolver *>(this)->Factor();
	  bHasBeenReset = false;
     }

     int status;

     /*
      * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
      */

     iparm[1] = 3;
     iparm[2] = 3;

     int tsize = iSize;

     wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);

     status = iparm[63];

     if (status != 0) {
	  silent_cerr("Wsmp back substitution failed" << std::endl);

	  /* de-allocate memory */
	  wsmp_clear_();

	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if (iparm[5] > 0) {
	  iparm[1] = 4;
	  iparm[2] = 4;

	  wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);

	  status = iparm[63];

	  if (status != 0) {
	       silent_cerr("Wsmp iterative refinement failed" << std::endl);

	       /* de-allocate memory */
	       wsmp_clear_();

	       throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	  }
     }
}

void
WsmpSolver::Factor(void)
{
     int status;

     /*
      * NOTE: Axp, Aip, App should have been set by * MakeCompactForm()
      */

     if (!Symbolic && !bPrepareSymbolic()) {
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     iparm[1] = 2;
     iparm[2] = 2;

     int tsize = iSize;
	  
     wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);

     status = iparm[63];

     if (status != 0) {
	  silent_cerr("Wsmp numeric factorization failed" << std::endl);

	  /* de-allocate memory */
	  wsmp_clear_();
	  
	  throw ErrFactor(-1, MBDYN_EXCEPT_ARGS);
     }
}

void
WsmpSolver::MakeCompactForm(SparseMatrixHandler& mh,
                            std::vector<doublereal>& Ax,
                            std::vector<integer>& Ai,
                            std::vector<integer>& Ac,
                            std::vector<integer>& Ap) const
{
     if (!bHasBeenReset) {
	  return;
     }

     if (iNumNzPrev != mh.Nz()) {
	  Symbolic = false;
     }

     iNumNzPrev = mh.MakeCompressedColumnForm(Ax, Ai, Ap, 0);

     Axp = &(Ax[0]);
     Aip = &(Ai[0]);
     App = &(Ap[0]);
}

bool
WsmpSolver::bPrepareSymbolic(void)
{
     int status;

     iparm[1] = 1;
     iparm[2] = 1;

     int tsize = iSize;

     wgsmp_(&tsize, App, Aip, Axp, LinearSolver::pdRhs, &ldb, &nrhs, rmisc, iparm, dparm);

     status = iparm[63];

     if (status != 0) {
	  silent_cerr("Wsmp symbolic factorization failed" << std::endl);

	  /* de-allocate memory */
	  wsmp_clear_();

	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     Symbolic = true;

     return true;
}

/* WsmpSolver - end */

/* WsmpSparseSolutionManager - begin */
template <typename MatrixHandlerType>
WsmpSparseSolutionManager<MatrixHandlerType>::WsmpSparseSolutionManager(integer Dim,
									doublereal dPivot,
									const unsigned blockSize,
									const unsigned nt,
									const SolutionManager::ScaleOpt& oScale,
									integer iMaxIterRefine)
     : A(Dim, Dim),
       xb(Dim),
       xbVH(Dim, &xb[0])
{
     SAFENEWWITHCONSTRUCTOR(pLS, WsmpSolver,
			    WsmpSolver(Dim, dPivot, blockSize, nt, oScale, iMaxIterRefine));

     (void)pLS->pdSetResVec(&xb[0]);
     (void)pLS->pdSetSolVec(&xb[0]);
     pLS->SetSolutionManager(this);
}

template <typename MatrixHandlerType>
WsmpSparseSolutionManager<MatrixHandlerType>::~WsmpSparseSolutionManager(void)
{
     NO_OP;
}

template <typename MatrixHandlerType>
void
WsmpSparseSolutionManager<MatrixHandlerType>::MatrReset(void)
{
     pLS->Reset();
}

template <typename MatrixHandlerType>
void
WsmpSparseSolutionManager<MatrixHandlerType>::MakeCompressedColumnForm(void)
{
     pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
template <typename MatrixHandlerType>
void
WsmpSparseSolutionManager<MatrixHandlerType>::Solve(void)
{
     MakeCompressedColumnForm();
     pLS->Solve();
}

/* Rende disponibile l'handler per la matrice */
template <typename MatrixHandlerType>
MatrixHandler*
WsmpSparseSolutionManager<MatrixHandlerType>::pMatHdl(void) const
{
     return &A;
}

/* Rende disponibile l'handler per il termine noto */
template <typename MatrixHandlerType>
MyVectorHandler*
WsmpSparseSolutionManager<MatrixHandlerType>::pResHdl(void) const
{
     return &xbVH;
}

/* Rende disponibile l'handler per la soluzione */
template <typename MatrixHandlerType>
MyVectorHandler*
WsmpSparseSolutionManager<MatrixHandlerType>::pSolHdl(void) const
{
     return &xbVH;
}

/* WsmpSparseSolutionManager - end */
template <class CC>
WsmpSparseCCSolutionManager<CC>::WsmpSparseCCSolutionManager(integer Dim,
							     doublereal dPivot,
							     const unsigned& blockSize,
							     const unsigned nt,
							     const SolutionManager::ScaleOpt& oScale,
							     integer iMaxIterRefine)
     : WsmpSparseSolutionManager<SpMapMatrixHandler>(Dim, dPivot, blockSize, nt, oScale, iMaxIterRefine),
       CCReady(false),
       Ac(0)
{
     NO_OP;
}

template <class CC>
WsmpSparseCCSolutionManager<CC>::~WsmpSparseCCSolutionManager(void)
{
     if (Ac) {
	  SAFEDELETE(Ac);
     }
}

template <class CC>
void
WsmpSparseCCSolutionManager<CC>::MatrReset(void)
{
     pLS->Reset();
}

template <class CC>
void
WsmpSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
     if (!CCReady) {
	  pLS->MakeCompactForm(A, Ax, Ai, Adummy, Ap);

	  if (Ac == 0) {
	       SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
	  }

	  CCReady = true;
     }
}

/* Inizializzatore "speciale" */
template <class CC>
void
WsmpSparseCCSolutionManager<CC>::MatrInitialize()
{
     CCReady = false;

     if (Ac) {
	  SAFEDELETE(Ac); // Needed for DirCColMatrixHandler
	  Ac = nullptr;
     }
     
     MatrReset();
}

/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
WsmpSparseCCSolutionManager<CC>::pMatHdl(void) const
{
     if (!CCReady) {
	  return &A;
     }

     ASSERT(Ac != 0);
     return Ac;
}

template class WsmpSparseSolutionManager<SpMapMatrixHandler>;
template class WsmpSparseSolutionManager<SpGradientSparseMatrixHandler>;
template class WsmpSparseCCSolutionManager<CColMatrixHandler<0> >;
template class WsmpSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* WsmpSparseCCSolutionManager - end */

#endif /* USE_WSMP */
