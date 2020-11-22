/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2018(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifndef __PASTIX_SOLUTION_MANAGER_H__INCLUDED__
#define __PASTIX_SOLUTION_MANAGER_H__INCLUDED__

#ifdef USE_PASTIX

#include <iostream>
#include <vector>
#include <complex>

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "spmapmh.h"

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient_spmh.h"
#endif

#define MPI_COMM_WORLD 0

extern "C" {
#include <pastix.h>
}

class PastixSolver: public LinearSolver {
private:
     struct SpMatrix: spmatrix_t {
	  SpMatrix();
	  SpMatrix(const SpMatrix&)=delete;	 
	  ~SpMatrix();

	  SpMatrix& operator=(const SpMatrix&)=delete;

	  bool MakeCompactForm(const SparseMatrixHandler& mh);
	  
	  doublereal* pAx() const { return reinterpret_cast<doublereal*>(values); }
	  pastix_int_t* pAi() const { return rowptr; }
	  pastix_int_t* pAp() const { return colptr; }
	  
     private:
	  template <typename T>
	  static inline T* pAllocate(T* pMem, size_t nSize);
	  void Allocate(size_t iNumNZ, size_t iMatSize);
	  integer iNumNonZeros;
     };
     
     pastix_int_t iparm[IPARM_SIZE];
     doublereal dparm[DPARM_SIZE];
     pastix_data_t* pastix_data;
     mutable SpMatrix spm;
     mutable bool bDoOrdering;
     
public:
     explicit PastixSolver(SolutionManager* pSM,
			   integer iDim,
			   integer iNumIter,
			   integer iNumThreads,
			   integer iVerbose);
     ~PastixSolver();

#ifdef DEBUG
     void IsValid(void) const;
#endif /* DEBUG */

     void Initialize() { bDoOrdering = true; }
     virtual void Solve(void) const;
     SpMatrix& MakeCompactForm(SparseMatrixHandler& mh);
};

template <typename MatrixHandlerType>
class PastixSolutionManager: public SolutionManager {
private:
    std::vector<doublereal> x;
    std::vector<doublereal> b;
    mutable MyVectorHandler xVH, bVH;
    ScaleOpt scale;
    MatrixScaleBase* pMatScale;

protected:
    mutable MatrixHandlerType A;

    PastixSolver* pGetSolver() { return static_cast<PastixSolver*>(pLS); }

    template <typename MH>
    void ScaleMatrixAndRightHandSide(MH &mh);

    template <typename MH>
    MatrixScale<MH>& GetMatrixScale();

    void ScaleSolution(void);

public:
    PastixSolutionManager(integer iDim,
			  integer iNumThreads,
			  integer iNumIter,
			  const ScaleOpt& scale = ScaleOpt(),
			  integer iVerbose = 0);
    virtual ~PastixSolutionManager(void);
#ifdef DEBUG
    virtual void IsValid(void) const;
#endif /* DEBUG */
    virtual void MatrReset(void);
    virtual void MatrInitialize(void);
    virtual void Solve(void);
    virtual void MakeCompressedColumnForm(void);
    virtual MatrixHandler* pMatHdl(void) const;
    virtual VectorHandler* pResHdl(void) const;
    virtual VectorHandler* pSolHdl(void) const;
};

#endif

#endif
