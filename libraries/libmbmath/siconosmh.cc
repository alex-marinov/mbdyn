/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"

#ifdef USE_SICONOS
#include <iomanip>
#include "sp_gradient.h"
#include "myassert.h"
#include "siconosmh.h"

SiconosRowMap::SiconosRowMap(integer iSize)
{
     rgRowMap.reserve(iSize);

     for (integer i = 1; i <= iSize; ++i) {
          rgRowMap.push_back(i);
     }
}

#ifdef DEBUG
void SiconosRowMap::IsValid() const
{
     std::vector<bool> rgFlags(iGetSize(), false);

     for (integer i = 1; i <= iGetSize(); ++i) {
          ASSERT(iGetIndex(i) >= 1);
          ASSERT(iGetIndex(i) <= iGetSize());
          ASSERT(!rgFlags[iGetIndex(i) - 1]);

          rgFlags[iGetIndex(i) - 1] = true;
     }
}
#endif

void SiconosRowMap::SetIndex(std::vector<integer>&& rgIndex)
{
#ifdef DEBUG
     IsValid();
#endif

     rgRowMap = std::move(rgIndex);

#ifdef DEBUG
     IsValid();
#endif
}

SiconosVectorHandler::SiconosVectorHandler(integer iSize, doublereal* pdTmpVec, const SiconosRowMap* pRowMap)
     :MyVectorHandler(iSize, pdTmpVec), pRowMap(pRowMap)
{

}

SiconosVectorHandler::~SiconosVectorHandler()
{
}

void SiconosVectorHandler::Resize(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif

     MyVectorHandler::Resize(iNewSize);

#ifdef DEBUG
     IsValid();
#endif
}

void SiconosVectorHandler::ResizeReset(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif

     MyVectorHandler::ResizeReset(iNewSize);

#ifdef DEBUG
     IsValid();
#endif
}

void SiconosVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
     (*this)(iRow) = dCoef;
}

void SiconosVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
     (*this)(iRow) += dCoef;
}

void SiconosVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
     (*this)(iRow) -= dCoef;
}

const doublereal& SiconosVectorHandler::dGetCoef(integer iRow) const
{
     return (*this)(iRow);
}

const doublereal& SiconosVectorHandler::operator() (integer iRow) const
{
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());
     ASSERT(pRowMap->iGetSize() == iGetSize());
     ASSERT(pRowMap->iGetIndex(iRow) >= 1);
     ASSERT(pRowMap->iGetIndex(iRow) <= iGetSize());

     return MyVectorHandler::operator()(pRowMap->iGetIndex(iRow));
}

doublereal& SiconosVectorHandler::operator() (integer iRow)
{
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());
     ASSERT(pRowMap->iGetSize() == iGetSize());
     ASSERT(pRowMap->iGetIndex(iRow) >= 1);
     ASSERT(pRowMap->iGetIndex(iRow) <= iGetSize());

     return MyVectorHandler::operator()(pRowMap->iGetIndex(iRow));
}

void SiconosVectorHandler::SetRowMap(const SiconosRowMap* pRowMapNew)
{
#ifdef DEBUG
     IsValid();
#endif

     pRowMap = pRowMapNew;

#ifdef DEBUG
     IsValid();
#endif
}

SiconosMatrixHandler::SiconosMatrixHandler(const SiconosMatrixHandler& oMH)
     :iNumNz(oMH.iNumNz),
      pRowMap(oMH.pRowMap)
{
     pMatrix = static_cast<NumericsMatrix*>(malloc(sizeof(*pMatrix)));

     if (!pMatrix) {
          throw ErrMemory(MBDYN_EXCEPT_ARGS);
     }

     NM_null(pMatrix);

     NM_copy(oMH.pMatrix, pMatrix);
}

SiconosMatrixHandler::SiconosMatrixHandler(NM_types eStorageType, integer iNumRows, integer iNumCols, integer iNumNz, const SiconosRowMap* pRowMap)
     :iNumNz(iNumNz), pRowMap(pRowMap)
{
     pMatrix = NM_create(eStorageType, iNumRows, iNumCols);

     if (!pMatrix) {
          throw ErrMemory(MBDYN_EXCEPT_ARGS);
     }

     if (iNumNz >= 0) {
          NM_triplet_alloc(pMatrix, iNumNz);
     }

#ifdef DEBUG
     IsValid();
#endif
}

SiconosMatrixHandler::~SiconosMatrixHandler()
{
     NM_free(pMatrix);
}

#ifdef DEBUG
void SiconosMatrixHandler::IsValid() const
{
     ASSERT(NM_check(pMatrix) == 0);

     pRowMap->IsValid();

     ASSERT(iGetNumRows() == pRowMap->iGetSize());
     ASSERT(iGetNumCols() == pRowMap->iGetSize());
}
#endif

void SiconosMatrixHandler::Resize(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SiconosMatrixHandler::ResizeReset(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SiconosMatrixHandler::Reset()
{
     switch (pMatrix->storageType) {
     case NM_SPARSE:
          NM_clearSparseStorage(pMatrix);

          if (iNumNz >= 0) {
               NM_triplet_alloc(pMatrix, iNumNz);
          }
          break;

     case NM_DENSE:
          NM_scal(0., pMatrix);
          break;

     default:
          ASSERT(0);
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

void
SiconosMatrixHandler::PutCoef(integer iRow, integer iCol, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iCol >= 1);
     ASSERT(iRow <= iGetNumRows());
     ASSERT(iCol <= iGetNumCols());

     if (pMatrix->storageType != NM_DENSE) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     ASSERT(pRowMap->iGetSize() == iGetNumRows());
     ASSERT(pRowMap->iGetIndex(iRow) >= 1);
     ASSERT(pRowMap->iGetIndex(iRow) <= iGetNumRows());

     pMatrix->matrix0[pRowMap->iGetIndex(iRow) - 1 + (iCol - 1) * pMatrix->size0] = dCoef;

#ifdef DEBUG
     IsValid();
#endif
}

void
SiconosMatrixHandler::IncCoef(integer iRow, integer iCol, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iCol >= 1);
     ASSERT(iRow <= iGetNumRows());
     ASSERT(iCol <= iGetNumCols());
     ASSERT(pRowMap->iGetSize() == iGetNumRows());

     DEBUGCERR("IncCoef(" << iRow << ", " << iCol << ", " << dCoef << ") -> " << pRowMap->iGetIndex(iRow) << "\n") ;

     switch (pMatrix->storageType) {
     case NM_SPARSE:
          NM_entry(pMatrix, pRowMap->iGetIndex(iRow) - 1, iCol - 1, dCoef);
          break;

     case NM_DENSE:
          // FIXME: NM_entry is equivalent to PutCoef in case dense matrices
          pMatrix->matrix0[pRowMap->iGetIndex(iRow) - 1 + (iCol - 1) * pMatrix->size0] += dCoef;
          break;

     default:
          ASSERT(0);
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

#ifdef DEBUG
     IsValid();
#endif
}

void
SiconosMatrixHandler::DecCoef(integer iRow, integer iCol, const doublereal& dCoef)
{
     IncCoef(iRow, iCol, -dCoef);
}

const doublereal&
SiconosMatrixHandler::operator()(integer iRow, integer iCol) const
{
     ASSERT(iRow >= 1);
     ASSERT(iCol >= 1);
     ASSERT(iRow <= iGetNumRows());
     ASSERT(iCol <= iGetNumCols());
     ASSERT(pRowMap->iGetSize() == iGetNumRows());

     if (pMatrix->storageType != NM_DENSE) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     return pMatrix->matrix0[pRowMap->iGetIndex(iRow) - 1 + (iCol - 1) * pMatrix->size0];
}

doublereal&
SiconosMatrixHandler::operator()(integer iRow, integer iCol)
{
     ASSERT(iRow >= 1);
     ASSERT(iCol >= 1);
     ASSERT(iRow <= iGetNumRows());
     ASSERT(iCol <= iGetNumCols());
     ASSERT(pRowMap->iGetSize() == iGetNumRows());

     if (pMatrix->storageType != NM_DENSE) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     return pMatrix->matrix0[pRowMap->iGetIndex(iRow) - 1 + (iCol - 1) * pMatrix->size0];
}

integer SiconosMatrixHandler::iGetNumRows() const
{
     ASSERT(pMatrix);

     return pMatrix->size0;
}

integer SiconosMatrixHandler::iGetNumCols() const
{
     ASSERT(pMatrix);

     return pMatrix->size1;
}

std::ostream& SiconosMatrixHandler::Print(std::ostream& os, MatPrintFormat eFormat) const
{
#ifdef __GNUC__
     char* pBuffer = nullptr;
     size_t nSize = 0u;

     FILE* pStream = open_memstream(&pBuffer, &nSize);

     if (!pStream) {
          return os;
     }

     NM_write_in_file(pMatrix, pStream);

     fclose(pStream);

     os.write(pBuffer, nSize);

     free(pBuffer);
#endif

     return os;
}

void SiconosMatrixHandler::EnumerateNz(const std::function<EnumerateNzCallback>& func) const
{
     switch (pMatrix->storageType) {
     case NM_DENSE:
          for (int j = 0; j < pMatrix->size1; ++j) {
               for (int i = 0; i < pMatrix->size0; ++i) {
                    const doublereal dCoef = pMatrix->matrix0[pRowMap->iGetIndex(i + 1) - 1 + j * pMatrix->size0];

                    if (dCoef) {
                         func(i + 1, j + 1, dCoef);
                    }
               }
          }
     default:
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

VectorHandler&
SiconosMatrixHandler::MatVecMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumRows());
     ASSERT(x.iGetSize() == iGetNumCols());

     y.Reset(); // Without this, valgrind complains about use of not initialized memory

     NM_gemv(1., pMatrix, x.pdGetVec(), 0., y.pdGetVec());

     return y;
}

VectorHandler&
SiconosMatrixHandler::MatTVecMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumCols());
     ASSERT(x.iGetSize() == iGetNumRows());

     y.Reset(); // Without this, valgrind complains about use of not initialized memory

     NM_tgemv(1., pMatrix, x.pdGetVec(), 0., y.pdGetVec());

     return y;
}

VectorHandler&
SiconosMatrixHandler::MatVecIncMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumRows());
     ASSERT(x.iGetSize() == iGetNumCols());

     NM_gemv(1., pMatrix, x.pdGetVec(), 1., y.pdGetVec());

     return y;
}

VectorHandler&
SiconosMatrixHandler::MatTVecIncMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumCols());
     ASSERT(x.iGetSize() == iGetNumRows());

     NM_tgemv(1., pMatrix, x.pdGetVec(), 1., y.pdGetVec());

     return y;
}

VectorHandler&
SiconosMatrixHandler::MatVecDecMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumRows());
     ASSERT(x.iGetSize() == iGetNumCols());

     NM_gemv(-1., pMatrix, x.pdGetVec(), 1., y.pdGetVec());

     return y;
}

VectorHandler&
SiconosMatrixHandler::MatTVecDecMul(VectorHandler& y, const VectorHandler& x) const
{
     ASSERT(y.iGetSize() == iGetNumCols());
     ASSERT(x.iGetSize() == iGetNumRows());

     NM_tgemv(-1., pMatrix, x.pdGetVec(), 1., y.pdGetVec());

     return y;
}

void SiconosMatrixHandler::SetRowMap(const SiconosRowMap* pRowMapNew)
{
     ASSERT(pRowMapNew->iGetSize() == iGetNumRows());

#ifdef DEBUG
     IsValid();
#endif

     pRowMap = pRowMapNew;

#ifdef DEBUG
     IsValid();
#endif
}

VectorHandler&
SiconosMatrixHandler::MatVecMul_base(void (VectorHandler::*op)(integer iRow,
                                                               const doublereal& dCoef),
                                     VectorHandler& out, const VectorHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

VectorHandler&
SiconosMatrixHandler::MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
                                               const doublereal& dCoef),
                     VectorHandler& out, const VectorHandler& in) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal SiconosMatrixHandler::ConditionNumber(enum Norm_t eNorm) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal SiconosMatrixHandler::Norm(enum Norm_t eNorm) const
{
     switch (eNorm) {
     case NORM_1:
          return NM_norm_1(pMatrix);

     case NORM_INF:
          return NM_norm_inf(pMatrix);

     default:
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

void SiconosMatrixHandler::Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

MatrixHandler* SiconosMatrixHandler::Copy() const
{
     SiconosMatrixHandler* pMH = nullptr;

     SAFENEWWITHCONSTRUCTOR(pMH,
                            SiconosMatrixHandler,
                            SiconosMatrixHandler(*this));

     return pMH;
}

bool SiconosMatrixHandler::AddItem(integer iRow, const sp_grad::SpGradient& oItem)
{
     for (const auto& oDer: oItem) {
          IncCoef(iRow, oDer.iDof, oDer.dDer);
     }

     return true;
}

#endif
