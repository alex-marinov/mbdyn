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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "ac/lapack.h"
#include "ac/qrupdate.h"
#include "myassert.h"
#include "solman.h"
#include "qrwrap.h"
#include "linsol.h"

QrDenseSolver::QrDenseSolver(Data* pData)
        : LinearSolver(nullptr),
          pData(pData)
{

}

QrDenseSolver::~QrDenseSolver(void)
{
}

void
QrDenseSolver::Reset(void)
{
        bHasBeenReset = true;
}

void
QrDenseSolver::Solve(void) const
{
        if (bHasBeenReset) {
                const_cast<QrDenseSolver *>(this)->Factor();
                bHasBeenReset = false;
        }

        SolveQR();
}

void QrDenseSolver::SolveR() const
{
        const auto& r = pData->r;
        auto& b = pData->b;
        const auto& d = pData->d;
        const integer n = r.iGetNumCols();

        b(n) /= d(n);

        for (integer i = n - 1; i >= 1; i--) {
                doublereal sum = 0.0;

                for (integer j = i + 1; j <= n; j++) {
                        sum += r(i, j) * b(j);
                }

                b(i) = (b(i) - sum) / d(i);
        }
}

void QrDenseSolver::SolveQR() const
{
        const auto& r = pData->r;
        auto& b = pData->b;
        const auto& c = pData->c;
        const integer n = r.iGetNumCols();

        for (integer j = 1; j < n; j++) {
                doublereal sum = 0.0;

                for (integer i = j; i <= n; i++) {
                        sum += r(i, j) * b(i);
                }

                doublereal tau = sum / c(j);

                for (integer i = j; i <= n; i++) {
                        b(i) -= tau * r(i, j);
                }
        }

        SolveR();
}

void
QrDenseSolver::Factor(void)
{
        auto& r = pData->r;
        auto& d = pData->d;
        auto& c = pData->c;

        const integer n = r.iGetNumCols();

        for (integer k = 1; k < n; k++) {
                doublereal scale = 0.0;

                for (integer i = k; i <= n; i++) {
                        scale = std::max(scale, std::fabs(r(i, k)));
                }

                if (scale == 0.0) {
                        throw ErrNoPivot(k, MBDYN_EXCEPT_ARGS);
                } else {
                        for (integer i = k; i <= n; i++) {
                                r(i, k) /= scale;
                        }

                        doublereal sum = 0.0;

                        for (integer i = k; i <= n; i++) {
                                sum += r(i, k) * r(i, k);
                        }

                        doublereal sigma = copysign(sqrt(sum), r(k, k));

                        r(k, k) += sigma;
                        c(k) = sigma * r(k, k);
                        d(k) = -scale * sigma;

                        for (integer j = k + 1; j <= n; j++) {
                                sum = 0.0;
                                for (integer i = k; i <= n; i++) {
                                        sum += r(i, k) * r(i, j);
                                }

                                doublereal tau = sum / c(k);

                                for (integer i = k; i <= n; i++) {
                                        r(i, j) -= tau * r(i, k);
                                }
                        }
                }
        }

        d(n) = r(n, n);

        if (d(n) == 0.0) {
                throw ErrNoPivot(n, MBDYN_EXCEPT_ARGS);
        }
}

// Compute QR decomposition of Q * R + u * v^T = Q * (R + Q^T * u * v^T) = Q * (R + t * v^T)
void QrDenseSolver::UpdateQR(VectorHandler& u, VectorHandler& v)
{
        auto& r = pData->r;
        auto& d = pData->d;
        auto& t = pData->t;
        const auto& qt = pData->qt;
        const integer n = r.iGetNumCols();

        t.Reset();

        qt.MatVecIncMul(t, u); // t = Q^T * u

#ifdef DEBUG
        FullMatrixHandler A(n);

        for (integer i = 1; i <= n; ++i) {
                for (integer j = 1; j <= n; ++j) {
                        A(i, j) = 0;

                        for (integer k = 1; k <= n; ++k) {
                                A(i, j) += pData->qt(k, i) * (pData->r(k, j) + t(k) * v(j));
                        }
                }
        }
#endif
        integer k;

        for (k = n; k >= 1; k--) {
                if (t(k)) {
                        break;
                }
        }

        if (k < 1) {
                k=1;
        }

        for (integer i = k - 1; i >= 1; i--) {
                Rotate(i, t(i), -t(i + 1));

                if (t(i) == 0.0) {
                        t(i) = fabs(t(i + 1));
                } else if (fabs(t(i)) > fabs(t(i + 1))) {
                        t(i) = std::fabs(t(i)) * std::sqrt(1.0 + std::pow(t(i + 1) / t(i), 2));
                } else {
                        t(i) = fabs(t(i + 1)) * sqrt(1.0 + std::pow(t(i) / t(i + 1), 2));
                }
        }

        for (integer j = 1; j <= n; j++) {
                r(1, j) += t(1) * v(j);
        }

        for (integer i = 1; i < k; i++) {
                Rotate(i, r(i, i), -r(i + 1, i));
        }

        for (integer i = 1; i <= n; ++i) {
                if (r(i, i) == 0.) {
                        throw ErrNoPivot(i, MBDYN_EXCEPT_ARGS);
                }

                d(i) = r(i, i);
        }

#ifdef DEBUG
        FullMatrixHandler B(n);

        for (integer i = 1; i <= n; ++i) {
                for (integer j = 1; j <= n; ++j) {
                        B(i, j) = 0;

                        for (integer k = 1; k <= n; ++k) {
                                B(i, j) += pData->qt(k, i) * pData->r(k, j);
                        }
                }
        }

        doublereal normA = A.Norm();
        doublereal normB = B.Norm();
        doublereal dTol = std::max(normA, normB) * sqrt(std::numeric_limits<doublereal>::epsilon());
        doublereal dMaxErr = 0;

        for (integer i = 1; i <= n; ++i) {
                for (integer j = 1; j <= n; ++j) {
                        dMaxErr = std::max(dMaxErr, fabs(A(i, j) - B(i, j)));
                }
        }

        DEBUGCERR("dMaxErr=" << dMaxErr << std::endl);
        DEBUGCERR("dTol=" << dTol << std::endl);
        ASSERT(dMaxErr < dTol);
#endif
}

void QrDenseSolver::Rotate(integer i, doublereal a, doublereal b)
{
        auto& r = pData->r;
        auto& qt = pData->qt;
        const integer n = r.iGetNumCols();

        integer j;
        doublereal c, fact, s, w, y;

        if (a == 0.0) {
                c = 0.0;
                s = (b >= 0.0 ? 1.0 : -1.0);
        } else if (fabs(a) > fabs(b)) {
                fact = b / a;
                c = copysign(1.0 / sqrt(1.0 + (fact * fact)), a);
                s = fact * c;
        } else {
                fact = a / b;
                s = copysign(1.0 / sqrt(1.0 + (fact * fact)), b);
                c = fact * s;
        }

        for (j = i; j <= n; j++) {
                y = r(i, j);
                w = r(i + 1, j);
                r(i, j) = c * y - s * w;
                r(i + 1, j) = s * y + c * w;
        }

        for (j = 1; j <= n; j++) {
                y = qt(i, j);
                w = qt(i + 1, j);
                qt(i, j) = c * y - s * w;
                qt(i + 1, j) = s * y + c*w;
        }
}

void QrDenseSolver::InitQR()
{
        auto& r = pData->r;
        auto& qt = pData->qt;
        const auto& c = pData->c;
        const auto& d = pData->d;
        const integer n = r.iGetNumCols();

#ifdef DEBUG
        FullMatrixHandler A(pData->r);
#endif

        Factor();

        if (qt.iGetNumRows() != n) {
                ASSERT(qt.iGetNumRows() == 0);
                ASSERT(qt.iGetNumCols() == 0);

                qt.Resize(n, n);
        }

        for (integer i = 1; i <= n; i++) { /* Form Q explicitly. */
                for (integer j = 1; j <= n; j++) {
                        qt(i, j) = 0.0;
                }

                qt(i, i) = 1.0;
        }

        for (integer k = 1; k < n; k++) {
                if (c(k)) {
                        for (integer j=1; j <= n; j++) {
                                doublereal sum = 0.0;

                                for (integer i = k; i <= n; i++) {
                                        sum += r(i, k) * qt(i, j);
                                }

                                sum /= c(k);

                                for (integer i = k; i <= n; i++) {
                                        qt(i, j) -= sum * r(i, k);
                                }
                        }
                }
        }

        for (integer i = 1; i <= n; i++) { /* Form R explicitly. */
                r(i, i) = d(i);

                for (integer j = 1; j < i; j++) {
                        r(i, j) = 0.0;
                }
        }

#ifdef DEBUG
        FullMatrixHandler B(n);

        for (integer i = 1; i <= n; ++i) {
                for (integer j = 1; j <= n; ++j) {
                        B(i, j) = 0;

                        for (integer k = 1; k <= n; ++k) {
                                B(i, j) += pData->qt(k, i) * pData->r(k, j);
                        }
                }
        }

        doublereal dTol = std::max(A.Norm(), B.Norm()) * sqrt(std::numeric_limits<doublereal>::epsilon());
        doublereal dMaxErr = 0;

        for (integer i = 1; i <= n; ++i) {
                for (integer j = 1; j <= n; ++j) {
                        dMaxErr = std::max(dMaxErr, fabs(A(i, j) - B(i, j)));
                }
        }

        DEBUGCERR("dMaxErr(Q * R - A)=" << dMaxErr << std::endl);
        DEBUGCERR("dTol(Q * R - A)=" << dTol << std::endl);
        ASSERT(dMaxErr < dTol);
#endif
}

QrDenseSolutionManager::QrDenseSolutionManager(integer iDim)
        : oData(iDim)
{
        SAFENEWWITHCONSTRUCTOR(pLS,
                               QrDenseSolver,
                               QrDenseSolver(&oData));

        pLS->pdSetResVec(oData.b.pdGetVec());
        pLS->pdSetSolVec(oData.b.pdGetVec());

        pLS->SetSolutionManager(this);
}

QrDenseSolutionManager::~QrDenseSolutionManager(void)
{
        NO_OP;
}

QrDenseSolver* QrDenseSolutionManager::pGetLS() const
{
        auto pSolver = dynamic_cast<QrDenseSolver*>(pLS);

        ASSERT(pSolver != nullptr);

        return pSolver;
}

void
QrDenseSolutionManager::MatrReset(void)
{
        pLS->Reset();
}

void
QrDenseSolutionManager::Solve(void)
{
        pLS->Solve();
}

MatrixHandler*
QrDenseSolutionManager::pMatHdl(void) const
{
        return &oData.r;
}

MyVectorHandler*
QrDenseSolutionManager::pResHdl(void) const
{
        return &oData.b;
}

MyVectorHandler*
QrDenseSolutionManager::pSolHdl(void) const
{
        return &oData.b;
}

VectorHandler& QrDenseSolutionManager::MatVecOp(MatVecOpType op,
                                                VectorHandler& a,
                                                const VectorHandler& b) const
{
        switch (op) {
        case OP_A_MINUS_Q_B:
                return oData.qt.MatTVecDecMul(a, b);

        case OP_A_PLUS_QT_B:
                return oData.qt.MatVecIncMul(a, b);

        case OP_A_MINUS_R_B:
                return oData.r.MatVecDecMul(a, b);

        case OP_A_MINUS_RT_B:
                return oData.r.MatTVecDecMul(a, b);

        default:
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

void QrDenseSolutionManager::InitQR(void)
{
        pGetLS()->InitQR();
}

void QrDenseSolutionManager::UpdateQR(VectorHandler& u, VectorHandler& v)
{
        pGetLS()->UpdateQR(u, v);
}

void QrDenseSolutionManager::SolveR()
{
        pGetLS()->SolveR();
}

#if defined(USE_SUITESPARSE_QR)

CholModVectorHandler::CholModVectorHandler(integer iSize, CholModCommon& oCommon)
        :oCommon(oCommon)
{
        static_assert(sizeof(doublereal) == sizeof(double));

        v = cholmod_l_zeros(iSize, 1, CHOLMOD_REAL, oCommon);

        if (!v) {
                ErrMemory(MBDYN_EXCEPT_ARGS);
        }
}

CholModVectorHandler::CholModVectorHandler(cholmod_dense* v, CholModCommon& oCommon)
        :v(v), oCommon(oCommon)
{
}

CholModVectorHandler::~CholModVectorHandler(void)
{
        Free();
}

CholModVectorHandler& CholModVectorHandler::operator=(cholmod_dense* rhs)
{
        if (rhs != v) {
                Free();
                v = rhs;
        }

        return *this;
}

#ifdef DEBUG
void CholModVectorHandler::IsValid(void) const
{
        ASSERT(cholmod_l_check_common(oCommon));
        ASSERT(cholmod_l_check_dense(v, oCommon));
}
#endif /* DEBUG */

integer CholModVectorHandler::iGetSize(void) const
{
        return v ? v->nrow : 0;
}

void CholModVectorHandler::Reset(void)
{
        auto p = CholModVectorHandler::pdGetVec();

        if (p) {
                std::fill(p, p + v->nrow, 0.);
        }
}

void CholModVectorHandler::Resize(integer iNewSize)
{
        integer iCurrSize = iGetSize();
        
        if (iNewSize == iCurrSize) {
                return;
        }

        CholModVectorHandler oTmp(iNewSize, oCommon);

        if (iCurrSize) {
                auto p = CholModVectorHandler::pdGetVec();
                std::copy(p, p + std::min(iCurrSize, iNewSize), oTmp.pdGetVec());
        }
        
        std::swap(v, oTmp.v);
}

void CholModVectorHandler::ResizeReset(integer iNewSize)
{
        if (iNewSize != iGetSize()) {
                CholModVectorHandler oTmp(iNewSize, oCommon);

                std::swap(v, oTmp.v);
        } else {
                Reset();
        }
}

void CholModVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
        (*this)(iRow) = dCoef;
}

void CholModVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
        (*this)(iRow) += dCoef;
}

void CholModVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
        (*this)(iRow) -= dCoef;
}

const doublereal& CholModVectorHandler::dGetCoef(integer iRow) const
{
        return (*this)(iRow);
}

const doublereal& CholModVectorHandler::operator () (integer iRow) const
{
        ASSERT(iRow >= 1);
        ASSERT(iRow <= iGetSize());

        --iRow;

        return CholModVectorHandler::pdGetVec()[iRow];
}

doublereal& CholModVectorHandler::operator () (integer iRow)
{
        ASSERT(iRow >= 1);
        ASSERT(iRow <= iGetSize());

        --iRow;

        return CholModVectorHandler::pdGetVec()[iRow];
}

void CholModVectorHandler::Free()
{
        if (v) {
                cholmod_l_free_dense(&v, oCommon);
                ASSERT(!v);
        }
}

QrSparseSolver::Data::Data(integer iSize, unsigned flags)
        :ordering(0),
         A(iSize),
         Q(iSize), R(iSize),
         B(iSize, oCommon),
         X(iSize, oCommon),
         V(iSize, oCommon),
         eState(ST_MAT_MAP_HDL),
         pA(nullptr),
         pQ(nullptr),
         pR(nullptr),
         pE(nullptr),
         Einv(iSize),
         iNumNonZeros(-1),
         w(2 * iSize)
{
        switch (flags & LinSol::SOLVER_FLAGS_PERM_MASK) {
        case LinSol::SOLVER_FLAGS_ALLOWS_AMD:
                ordering = SPQR_ORDERING_AMD;
                break;
        case LinSol::SOLVER_FLAGS_ALLOWS_METIS:
                ordering = SPQR_ORDERING_METIS;
                break;
        case LinSol::SOLVER_FLAGS_ALLOWS_GIVEN:
                ordering = SPQR_ORDERING_GIVEN;
                break;
        case LinSol::SOLVER_FLAGS_ALLOWS_COLAMD:
                ordering = SPQR_ORDERING_COLAMD;
                break;
        default:
                ordering = SPQR_ORDERING_DEFAULT;
        }
}

QrSparseSolver::Data::~Data()
{
        FactorClean();
        MatrixClean();
}

void QrSparseSolver::Data::Reset(void)
{
        eState = ST_MAT_MAP_HDL;
}

VectorHandler&
QrSparseSolver::Data::MatVecOp(QrSolutionManager::MatVecOpType op,
                               VectorHandler& a,
                               const VectorHandler& b)
{
        if (eState < ST_FACT_SPARSE) {
                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        switch (op) {
        case QrSolutionManager::OP_A_MINUS_Q_B:
        case QrSolutionManager::OP_A_PLUS_QT_B:
                if (eState < ST_FACT_UPDATE) {
                        switch (op) {
                        case QrSolutionManager::OP_A_MINUS_Q_B:
                                MatVecDecMul(pQ, a, b);
                                break;

                        case QrSolutionManager::OP_A_PLUS_QT_B:
                                MatTVecIncMul(pQ, a, b);
                                break;

                        default:
                                ASSERT(0);
                        }
                } else {
                        // dense matrix vector product
                        switch (op) {
                        case QrSolutionManager::OP_A_MINUS_Q_B:
                                Q.MatVecDecMul(a, b);
                                break;

                        case QrSolutionManager::OP_A_PLUS_QT_B:
                                Q.MatTVecIncMul(a, b);
                                break;

                        default:
                                ASSERT(0);
                        }
                }
                break;

        case QrSolutionManager::OP_A_MINUS_R_B:
        case QrSolutionManager::OP_A_MINUS_RT_B:
                const VectorHandler* pb;
                
                if (pE && op == QrSolutionManager::OP_A_MINUS_R_B) {                        
                        for (integer i = 1; i <= b.iGetSize(); ++i) {
                                V(i) = b(pE[i - 1] + 1);
                        }
                        
                        pb = &V;
                } else {
                        pb = &b;
                }
                
                if (eState < ST_FACT_UPDATE) {
                        switch (op) {
                        case QrSolutionManager::OP_A_MINUS_R_B:
                                MatVecDecMul(pR, a, *pb);
                                break;
                        case QrSolutionManager::OP_A_MINUS_RT_B:
                                MatTVecDecMul(pR, a, *pb);                                
                                break;
                        default:
                                ASSERT(0);
                        }
                } else {
                        // dense matrix vector product
                        switch (op) {
                        case QrSolutionManager::OP_A_MINUS_R_B:
                                R.MatVecDecMul(a, *pb);
                                break;

                        case QrSolutionManager::OP_A_MINUS_RT_B:
                                R.MatTVecDecMul(a, *pb);
                                break;

                        default:
                                ASSERT(0);
                        }
                }

                if (pE && op == QrSolutionManager::OP_A_MINUS_RT_B) {
                        Permute(a, pE);
                }
                
                break;
                
        default:
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return a;
}

MatrixHandler* QrSparseSolver::Data::pMatHdl()
{
        return &A;
}

VectorHandler* QrSparseSolver::Data::pResHdl()
{
        return &B;
}

VectorHandler* QrSparseSolver::Data::pSolHdl()
{
        return &X;
}

void QrSparseSolver::Data::MakeCompressedColumnForm(void)
{
        const size_t nr = A.iGetNumRows();
        const size_t nc = A.iGetNumCols();
        const size_t nz = A.Nz();

        if (pA && (pA->nzmax != nz || pA->nrow != nr || pA->ncol != nc)) {
                MatrixClean();
        }

        if (!pA) {
                pA = cholmod_l_allocate_sparse(nr, nc, nz, 1, 1, 0, CHOLMOD_REAL, oCommon);

                if (!pA) {
                        throw ErrMemory(MBDYN_EXCEPT_ARGS);
                }
        }

        ASSERT(nz == pA->nzmax);
        ASSERT(nr == pA->nrow);
        ASSERT(nc == pA->ncol);

        auto const colptr = reinterpret_cast<SuiteSparse_long*>(pA->p);
        auto const rowptr = reinterpret_cast<SuiteSparse_long*>(pA->i);
        auto const values = reinterpret_cast<doublereal*>(pA->x);

        SuiteSparse_long j = -1, k = 0;

        for (const auto& elem: A) {
                ASSERT(elem.iRow >= 0);
                ASSERT(elem.iRow < A.iGetNumRows());
                ASSERT(elem.iCol >= 0);
                ASSERT(elem.iCol < A.iGetNumCols());

                if (j != elem.iCol) {
                        ASSERT(j == -1 || elem.iCol == j + 1);

                        colptr[elem.iCol] = k;
                        j = elem.iCol;
                }

                rowptr[k] = elem.iRow;
                values[k] = elem.dCoef;

#ifdef DEBUG
                if (!std::isfinite(elem.dCoef)) {
                        throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
                }
#endif
                ++k;
        }

        colptr[j + 1] = k;

        eState = ST_MAT_COMPR_COL;
}

void QrSparseSolver::Data::MatrixClean(void)
{
        if (pA) {
                cholmod_l_free_sparse(&pA, oCommon);

                ASSERT(!pA);
        }

        eState = ST_MAT_MAP_HDL;
}

void QrSparseSolver::Data::FactorClean(void)
{
        if (pQ) {
                cholmod_l_free_sparse(&pQ, oCommon);
                ASSERT(!pQ);
        }

        if (pR) {
                cholmod_l_free_sparse(&pR, oCommon);
                ASSERT(!pR);
        }

        if (pE) {
                cholmod_l_free(pA->nrow, sizeof(*pE), pE, oCommon);
                pE = nullptr;
        }
}

void QrSparseSolver::Data::Factor(void)
{
        if (eState < ST_MAT_COMPR_COL) {
                MakeCompressedColumnForm();
        }

        FactorClean();
                
        SuiteSparseQR<doublereal>(ordering,
                                  SPQR_NO_TOL,
                                  pA->ncol,
                                  pA,
                                  &pQ,
                                  &pR,
                                  &pE,
                                  oCommon);

        if (pE) {
                ASSERT(Einv.size() == pA->ncol);

                for (size_t i = 0; i < Einv.size(); ++i) {
                        Einv[pE[i]] = i;
                }
        }
        
        eState = ST_FACT_SPARSE;
}

void QrSparseSolver::Data::UpdateQR(VectorHandler& u, VectorHandler& v)
{
        if (eState < ST_FACT_UPDATE) {
                InitQR();
        }

        ASSERT(Q.iGetNumCols() == R.iGetNumRows());
        ASSERT(u.iGetSize() == Q.iGetNumRows());
        ASSERT(v.iGetSize() == R.iGetNumCols());
        ASSERT(w.size() == 2 * static_cast<size_t>(Q.iGetNumCols()));

        if (pE) {
                Permute(v, &Einv.front());
        }
        
        integer M = Q.iGetNumRows();
        integer N = R.iGetNumCols();
        integer K = Q.iGetNumCols();
        integer LDQ = Q.iGetNumRows();
        integer LDR = R.iGetNumRows();

        __FC_DECL__(dqr1up)(&M,
                            &N,
                            &K,
                            Q.pdGetMat(),
                            &LDQ,
                            R.pdGetMat(),
                            &LDR,
                            u.pdGetVec(),
                            v.pdGetVec(),
                            &w[0]);
}

void QrSparseSolver::Data::UTSolve(const cholmod_sparse* const pR, cholmod_dense* const pX)
{
        // Upper triangular solve: GNU-Octave dSparse.cc
        const SuiteSparse_long nc = pR->ncol;

        auto const cidx = reinterpret_cast<const SuiteSparse_long*>(pR->p);
        auto const ridx = reinterpret_cast<const SuiteSparse_long*>(pR->i);
        auto const data = reinterpret_cast<const doublereal*>(pR->x);
        auto const rhs = reinterpret_cast<doublereal*>(pX->x);

        for (SuiteSparse_long k = nc - 1; k >= 0; k--)
        {
                if (rhs[k] != 0.)
                {
                        if (ridx[cidx[k + 1] - 1] != k || data[cidx[k + 1] - 1] == 0.)
                        {
                                throw ErrFactor(k, MBDYN_EXCEPT_ARGS);
                        }

                        const doublereal tmp = rhs[k] / data[cidx[k + 1] - 1];

                        rhs[k] = tmp;

                        for (SuiteSparse_long i = cidx[k]; i < cidx[k + 1] - 1; i++)
                        {
                                rhs[ridx[i]] -= tmp * data[i];
                        }
                }
        }
}

void QrSparseSolver::Data::Permute(VectorHandler& X, SuiteSparse_long* const perm)
{
        // Permute vector: Pastix bvec_dcompute.c
        auto const x = X.pdGetVec();
        const SuiteSparse_long m = X.iGetSize();

        for(SuiteSparse_long k = 0; k < m; k++) {
                SuiteSparse_long i = k;
                SuiteSparse_long j = perm[i];

                if (j < 0) {
                        continue;
                }

                perm[i] = -j - 1;

                while(j != k) {
                        std::swap(x[j], x[k]);

                        i = j;
                        j = perm[i];
                        perm[i] = -j - 1;

                        ASSERT( (j != i) && (j >= 0) );
                }
        }

        for(SuiteSparse_long k = 0; k < m; k++) {
                ASSERT(perm[k] < 0);
                perm[k] = -perm[k] - 1;
        }
}

template <typename T>
void QrSparseSolver::Data::MatVecOp(const cholmod_sparse* const pMat, const T& op)
{
        const auto cidx = reinterpret_cast<const SuiteSparse_long*>(pMat->p);
        const auto ridx = reinterpret_cast<const SuiteSparse_long*>(pMat->i);
        const auto data = reinterpret_cast<const doublereal*>(pMat->x);
        const SuiteSparse_long n = pMat->nrow;

        for (SuiteSparse_long j = 0; j < n; ++j) {
                for (SuiteSparse_long k = cidx[j]; k < cidx[j + 1]; ++k) {
                        op(ridx[k], j, data[k]);
                }
        }
}

void QrSparseSolver::Data::MatTVecIncMul(const cholmod_sparse* pMat, VectorHandler& out, const VectorHandler& in)
{
        ASSERT(out.iGetSize() == static_cast<ssize_t>(pMat->ncol));
        ASSERT(in.iGetSize() == static_cast<ssize_t>(pMat->nrow));

        doublereal* const pout = out.pdGetVec();
        const doublereal* const pin = in.pdGetVec();

        MatVecOp(pMat,
                 [pout, pin] (SuiteSparse_long i,
                              SuiteSparse_long j,
                              doublereal d) {
                         pout[j] += d * pin[i];
                 });
}

void QrSparseSolver::Data::MatTVecDecMul(const cholmod_sparse* pMat, VectorHandler& out, const VectorHandler& in)
{
        ASSERT(out.iGetSize() == static_cast<ssize_t>(pMat->ncol));
        ASSERT(in.iGetSize() == static_cast<ssize_t>(pMat->nrow));

        doublereal* const pout = out.pdGetVec();
        const doublereal* const pin = in.pdGetVec();

        MatVecOp(pMat,
                 [pout, pin] (SuiteSparse_long i,
                              SuiteSparse_long j,
                              doublereal d) {
                         pout[j] -= d * pin[i];
                 });
}

void QrSparseSolver::Data::MatVecDecMul(const cholmod_sparse* pMat,
                                        VectorHandler& out,
                                        const VectorHandler& in)
{
        ASSERT(out.iGetSize() == static_cast<ssize_t>(pMat->nrow));
        ASSERT(in.iGetSize() == static_cast<ssize_t>(pMat->ncol));

        doublereal* const pout = out.pdGetVec();
        const doublereal* const pin = in.pdGetVec();

        MatVecOp(pMat,
                 [pout, pin] (SuiteSparse_long i,
                              SuiteSparse_long j,
                              doublereal d) {
                         pout[i] -= d * pin[j];
                 });
}

void QrSparseSolver::Data::InitQR(void)
{
        if (eState < ST_FACT_SPARSE) {
                Factor();
        }

        ASSERT(pA);
        ASSERT(pQ);
        ASSERT(pR);

        ASSERT(static_cast<size_t>(Q.iGetNumRows()) == pQ->nrow);
        ASSERT(static_cast<size_t>(Q.iGetNumCols()) == pQ->ncol);

        const auto Qcptr = reinterpret_cast<const SuiteSparse_long*>(pQ->p);
        const auto Qridx = reinterpret_cast<const SuiteSparse_long*>(pQ->i);
        const auto Qval = reinterpret_cast<const doublereal*>(pQ->x);

        Q.Reset();

        for (size_t j = 0; j < pQ->ncol; ++j) {
                for (SuiteSparse_long k = Qcptr[j]; k < Qcptr[j + 1]; ++k) {
                        Q.PutCoef(Qridx[k] + 1, j + 1, Qval[k]);
                }
        }

#ifdef DEBUG
        auto prec = std::cerr.precision();

        std::cerr.precision(16);

        DEBUGCERR("A=[" << std::endl);

        A.Print(std::cerr, MatrixHandler::MAT_PRINT_FULL);

        DEBUGCERR("];\nQ=[" << std::endl);

        Q.Print(std::cerr, MatrixHandler::MAT_PRINT_FULL);
#endif

        const auto Rcptr = reinterpret_cast<const SuiteSparse_long*>(pR->p);
        const auto Rridx = reinterpret_cast<const SuiteSparse_long*>(pR->i);
        const auto Rval = reinterpret_cast<const doublereal*>(pR->x);

        R.Reset();

        for (size_t j = 0; j < pR->ncol; ++j) {
                for (SuiteSparse_long k = Rcptr[j]; k < Rcptr[j + 1]; ++k) {
                        R.PutCoef(Rridx[k] + 1, j + 1, Rval[k]);
                }
        }

#ifdef DEBUG
        {
#if DEBUG > 1
                DEBUGCERR("];\nR=[" << std::endl);

                R.Print(std::cerr, MatrixHandler::MAT_PRINT_FULL);

                DEBUGCERR("];\n");

                std::cerr.precision(prec);
#endif
                FullMatrixHandler Atest(A.iGetNumRows());

                Q.MatMatMul(Atest, R);

                doublereal dErrA = 0;

                for (integer i = 1; i <= Atest.iGetNumRows(); ++i) {
                        for (integer j = 1; j <= Atest.iGetNumCols(); ++j) {
                                dErrA = std::max(dErrA, fabs(Atest.dGetCoef(i, j) - A.dGetCoef(i, j)));
                        }
                }

                dErrA /= A.Norm();

                constexpr doublereal dTolRel = sqrt(std::numeric_limits<doublereal>::epsilon());

                FullMatrixHandler Qtest(A.iGetNumRows());

                Q.MatTMatMul(Qtest, Q);

                for (integer i = 1; i <= Qtest.iGetNumCols(); ++i) {
                        Qtest.DecCoef(i, i, 1);
                }

                doublereal dErrQ = 0;

                for (integer i = 1; i <= Qtest.iGetNumRows(); ++i) {
                        for (integer j = 1; j <= Qtest.iGetNumCols(); ++j) {
                                dErrQ = std::max(dErrQ, fabs(Qtest.dGetCoef(i, j)));
                        }
                }

                DEBUGCERR("dErrA=" << dErrA << std::endl);
                DEBUGCERR("dErrQ=" << dErrQ << std::endl);

                ASSERT(dErrA < dTolRel);
                ASSERT(dErrQ < dTolRel);
        }
#endif

        eState = ST_FACT_UPDATE;
}

void QrSparseSolver::Data::SolveR(void)
{
        if (eState < ST_FACT_SPARSE) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

#ifdef DEBUG
        for (integer i = 1; i <= X.iGetSize(); ++i) {
                if (!std::isfinite(X(i))) {
                        throw ErrFactor(i, MBDYN_EXCEPT_ARGS);
                }
        }
#endif

        if (eState < ST_FACT_UPDATE) {
                // sparse linear solver
                UTSolve(pR, X);
        } else {
                // dense linear solver
                char UPLO = 'U';
                char TRANS = 'N';
                char DIAG = 'N';
                integer N = R.iGetNumCols();
                integer NRHS = 1;
                const doublereal* A = R.pdGetMat();
                integer LDA = R.iGetNumRows();
                doublereal* B = X.pdGetVec();
                integer LDB = X.iGetSize();
                integer INFO;

                __FC_DECL__(dtrtrs)(&UPLO,
                                    &TRANS,
                                    &DIAG,
                                    &N,
                                    &NRHS,
                                    A,
                                    &LDA,
                                    B,
                                    &LDB,
                                    &INFO);

                if (INFO > 0) {
                        throw ErrFactor(INFO, MBDYN_EXCEPT_ARGS);
                } else if (INFO < 0) {
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }

        if (pE) {
                Permute(X, pE);
        }
        
#ifdef DEBUG
        for (integer i = 1; i <= X.iGetSize(); ++i) {
                if (!std::isfinite(X(i))) {
                        throw ErrFactor(i, MBDYN_EXCEPT_ARGS);
                }
        }
#endif
}

void QrSparseSolver::Data::SolveQR(void)
{
        if (eState < ST_FACT_SPARSE) {
                Factor();
        }

        X.Reset();
        
        if (eState < ST_FACT_UPDATE) {
                MatTVecIncMul(pQ, X, B);
        } else {
                Q.MatTVecIncMul(X, B);
        }

        SolveR();
}

QrSparseSolver::QrSparseSolver(Data* pData)
        :pData(pData)
{
}

QrSparseSolver::~QrSparseSolver(void)
{
}

void QrSparseSolver::Reset(void)
{
        pData->Reset();
}

void QrSparseSolver::Factor(void) const
{
        pData->Factor();
}

void QrSparseSolver::Solve(void) const
{
        pData->SolveQR();
}

void QrSparseSolver::SolveR(void) const
{
        pData->SolveR();
}

void QrSparseSolver::InitQR()
{
        Factor();
}

void QrSparseSolver::UpdateQR(VectorHandler& u, VectorHandler& v)
{
        pData->UpdateQR(u, v);
}

QrSparseSolutionManager::QrSparseSolutionManager(integer Dim, unsigned flags)
        :oData(Dim, flags)
{
        SAFENEWWITHCONSTRUCTOR(pLS,
                               QrSparseSolver,
                               QrSparseSolver(&oData));

        pLS->SetSolutionManager(this);
}

QrSparseSolutionManager::~QrSparseSolutionManager(void)
{
}

QrSparseSolver* QrSparseSolutionManager::pGetLS() const
{
        return dynamic_cast<QrSparseSolver*>(pLS);
}

#ifdef DEBUG
void QrSparseSolutionManager::IsValid(void) const
{
                NO_OP;
};
#endif /* DEBUG */

void QrSparseSolutionManager::MatrReset(void)
{
        pLS->Reset();
}

void QrSparseSolutionManager::Solve(void)
{
        pLS->Solve();
}

MatrixHandler* QrSparseSolutionManager::pMatHdl(void) const
{
        return oData.pMatHdl();
}

VectorHandler* QrSparseSolutionManager::pResHdl(void) const
{
        return oData.pResHdl();
}

VectorHandler* QrSparseSolutionManager::pSolHdl(void) const
{
        return oData.pSolHdl();
}

void QrSparseSolutionManager::SolveR(void)
{
        pGetLS()->SolveR();
}

void QrSparseSolutionManager::InitQR(void)
{
        pGetLS()->InitQR();
}

void QrSparseSolutionManager::UpdateQR(VectorHandler& u, VectorHandler& v)
{
        pGetLS()->UpdateQR(u, v);
}

VectorHandler&
QrSparseSolutionManager::MatVecOp(MatVecOpType op,
                                  VectorHandler& a,
                                  const VectorHandler& b) const
{
        return oData.MatVecOp(op, a, b);
}

#endif
