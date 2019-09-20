/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2019
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
  AUTHOR: Reinhard Resch <r.resch@a1.net>
  Copyright (C) 2019(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

/*
  References:
  Numerical recipes in C: the art of scientific computing / William H. Press [et al.]. – 2nd ed.
  ISBN 0-521-43108-5
*/

#ifndef QrSolutionManager_hh
#define QrSolutionManager_hh

#include <iostream>
#include <vector>

#if defined(USE_SUITESPARSE_QR)
#include "SuiteSparseQR.hpp"
#endif

#include "myassert.h"
#include "mynewmem.h"
#include "ls.h"
#include "solman.h"
#include "fullmh.h"

#if defined(USE_SUITESPARSE_QR)
#include "spmapmh.h"
#endif

class QrDenseSolver: public LinearSolver {
public:
        struct Data {
                explicit Data(integer iDim)
                        :r(iDim, iDim),
                         qt(iDim, iDim),
                         b(iDim),
                         c(iDim),
                         d(iDim),
                         t(iDim) {
                }

                FullMatrixHandler r, qt;
                MyVectorHandler b, c, d, t;
        };

private:
        Data* pData;

        void Factor(void);
        void Rotate(integer i, doublereal a, doublereal b);
        void SolveQR(void) const;

public:
        explicit QrDenseSolver(Data* pData);
        ~QrDenseSolver(void);

        void Reset(void);
        void Solve(void) const;
        void SolveR(void) const;

        void InitQR();
        void UpdateQR(VectorHandler& u, VectorHandler& v);
};

class QrDenseSolutionManager: public QrSolutionManager {
private:
        mutable QrDenseSolver::Data oData;
        inline QrDenseSolver* pGetLS() const;

public:
        explicit QrDenseSolutionManager(integer Dim);
        virtual ~QrDenseSolutionManager(void);
#ifdef DEBUG
        virtual void IsValid(void) const {
                NO_OP;
        };
#endif /* DEBUG */
        virtual void MatrReset(void);
        virtual void Solve(void);
        virtual MatrixHandler* pMatHdl(void) const;
        virtual MyVectorHandler* pResHdl(void) const;
        virtual MyVectorHandler* pSolHdl(void) const;
        virtual VectorHandler&
        MatVecOp(MatVecOpType op,
                 VectorHandler& a,
                 const VectorHandler& b) const;        
        virtual void SolveR(void);
        virtual void InitQR(void);
        virtual void UpdateQR(VectorHandler& u, VectorHandler& v);
};

#if defined(USE_SUITESPARSE_QR)

class CholModCommon {
public:
        CholModCommon() {
                cholmod_l_start(&oCommon);
        }

        ~CholModCommon() {
                cholmod_l_finish(&oCommon);
        }

        operator cholmod_common*() {
                return &oCommon;
        }
        
private:
        cholmod_common oCommon;
};

class CholModVectorHandler: public VectorHandler {
public:
        explicit CholModVectorHandler(integer iSize, CholModCommon& oCommon);
        explicit CholModVectorHandler(cholmod_dense* v, CholModCommon& oCommon);
        virtual ~CholModVectorHandler(void);

#ifdef DEBUG
        virtual void IsValid(void) const;
#endif /* DEBUG */

        inline doublereal* pdGetVec(void) const {
                return v ? reinterpret_cast<doublereal*>(v->x) : nullptr;
        }

        CholModVectorHandler& operator=(const CholModVectorHandler& rhs) {
                VectorHandler::operator=(rhs);
                return *this;
        }       

        CholModVectorHandler& operator=(const VectorHandler& rhs) {
                VectorHandler::operator=(rhs);
                return *this;
        }
        
        CholModVectorHandler& operator=(cholmod_dense* pRhs);
        
        virtual integer iGetSize(void) const;

        virtual void Reset(void);

        virtual void Resize(integer iNewSize);

        virtual void ResizeReset(integer);

        virtual void PutCoef(integer iRow, const doublereal& dCoef);

        virtual void IncCoef(integer iRow, const doublereal& dCoef);

        virtual void DecCoef(integer iRow, const doublereal& dCoef);

        virtual const doublereal& dGetCoef(integer iRow) const;

        virtual const doublereal& operator () (integer iRow) const;

        virtual doublereal& operator () (integer iRow);

        operator cholmod_dense*() { return v; }
        operator const cholmod_dense*() const { return v; }

private:
        void Free();
        
        cholmod_dense* v;
        CholModCommon& oCommon;
};

class QrSparseSolver: public LinearSolver {
public:
        class Data {
        public:
                explicit Data(integer iSize, unsigned flags);
                ~Data();
                void Reset(void);
                void Factor(void);
                void SolveR(void);
                void SolveQR(void);               
                void UpdateQR(VectorHandler& u, VectorHandler& v);
                VectorHandler&
                MatVecOp(QrSolutionManager::MatVecOpType op,
                         VectorHandler& a,
                         const VectorHandler& b);
                MatrixHandler* pMatHdl();
                VectorHandler* pResHdl();                
                VectorHandler* pSolHdl();
                
        private:
                void MakeCompressedColumnForm(void);                
                void InitQR(void);
                void FactorClean(void);
                void MatrixClean(void);
                static void UTSolve(const cholmod_sparse* const pR, cholmod_dense* const pX);
                static void Permute(VectorHandler& X, SuiteSparse_long* const perm);
                template <typename T>
                static void MatVecOp(const cholmod_sparse* const pMat, const T& op);
                static void MatTVecIncMul(const cholmod_sparse* pMat, VectorHandler& out, const VectorHandler& in);
                static void MatTVecDecMul(const cholmod_sparse* pMat, VectorHandler& out, const VectorHandler& in);
                static void MatVecDecMul(const cholmod_sparse* pMat, VectorHandler& out, const VectorHandler& in);

                unsigned ordering;
                CholModCommon oCommon; // Must be initialized before any CholModVectorHandler!
                SpMapMatrixHandler A;
                FullMatrixHandler Q, R;
                CholModVectorHandler B, X, V;
                
                enum State {
                        ST_MAT_MAP_HDL,
                        ST_MAT_COMPR_COL,
                        ST_FACT_SPARSE,                        
                        ST_FACT_UPDATE
                } eState;
                
                cholmod_sparse* pA;
                cholmod_sparse* pQ;
                cholmod_sparse* pR;
                SuiteSparse_long* pE;
                std::vector<SuiteSparse_long> Einv;
                SuiteSparse_long iNumNonZeros;
                std::vector<doublereal> w;
        };

private:
        void Factor(void) const;
        
        Data* const pData;
        
public:
        explicit QrSparseSolver(Data* pData);
        ~QrSparseSolver(void);

        void Reset(void);
        void Solve(void) const;
        void SolveR(void) const;

        void InitQR();
        void UpdateQR(VectorHandler& u, VectorHandler& v);
};

class QrSparseSolutionManager: public QrSolutionManager {
private:
        mutable QrSparseSolver::Data oData;
        inline QrSparseSolver* pGetLS() const;

public:
        explicit QrSparseSolutionManager(integer Dim, unsigned flags);
        virtual ~QrSparseSolutionManager(void);
#ifdef DEBUG
        virtual void IsValid(void) const;
#endif
        virtual void MatrReset(void);
        virtual void Solve(void);
        virtual MatrixHandler* pMatHdl(void) const;
        virtual VectorHandler* pResHdl(void) const;
        virtual VectorHandler* pSolHdl(void) const;
        virtual VectorHandler&
        MatVecOp(MatVecOpType op,
                 VectorHandler& a,
                 const VectorHandler& b) const;        
        virtual void SolveR(void);
        virtual void InitQR(void);
        virtual void UpdateQR(VectorHandler& u, VectorHandler& v);

};
#endif

#endif /* QrSolutionManager_hh */
