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

#ifndef __SICONOS_SPARSE_MATRIX_HANDLER_H__INCLUDED__
#define __SICONOS_SPARSE_MATRIX_HANDLER_H__INCLUDED__

#ifdef USE_SICONOS
#include <vector>

#include "mh.h"
#include <numerics/NumericsMatrix.h>

class SiconosIndexMap {
public:
     explicit SiconosIndexMap(integer iSize);

#ifdef DEBUG
     void IsValid() const;
#endif

     void SetIndex(std::vector<integer>&& rgIndex);

     integer iGetIndex(integer iIndex) const {
          ASSERT(iIndex >= 1);
          ASSERT(iIndex <= iGetSize());
          ASSERT(rgRowMap[iIndex - 1] >= 1);
          ASSERT(rgRowMap[iIndex - 1] <= iGetSize());
          
          return rgRowMap[iIndex - 1];
     }

     integer iGetSize() const {
          return rgRowMap.size();
     }
private:
     std::vector<integer> rgRowMap;
};

class SiconosVectorHandler: public MyVectorHandler {
public:
     SiconosVectorHandler(integer iSize = 0, doublereal* pdTmpVec = nullptr, const SiconosIndexMap* pIndexMap = nullptr);
     virtual ~SiconosVectorHandler();

     virtual SiconosVectorHandler& operator=(const VectorHandler& VH) override final;

     virtual SiconosVectorHandler& ScalarAddMul(const VectorHandler& VH, const doublereal& d) override final;
     
     virtual SiconosVectorHandler& ScalarAddMul(const VectorHandler& VH1, const VectorHandler& VH2, const doublereal& d) override final;
     
     virtual void Resize(integer iNewSize) override final;

     virtual void ResizeReset(integer) override final;

     virtual void PutCoef(integer iRow, const doublereal& dCoef) override final;

     virtual void IncCoef(integer iRow, const doublereal& dCoef) override final;

     virtual void DecCoef(integer iRow, const doublereal& dCoef) override final;

     virtual const doublereal& dGetCoef(integer iRow) const override final;

     virtual const doublereal& operator()(integer iRow) const override final;

     virtual doublereal& operator()(integer iRow) override final;

     void SetIndexMap(const SiconosIndexMap* pIndexMapNew);

     const SiconosIndexMap* pGetIndexMap() const { return pIndexMap; }
private:
     const SiconosIndexMap* pIndexMap;
};

class SiconosMatrixHandler: public MatrixHandler {
private:
     SiconosMatrixHandler(const SiconosMatrixHandler& oMH);

public:
     SiconosMatrixHandler(NM_types eStorageType, integer iNumRows, integer iNumCols, integer iNumNz, const SiconosIndexMap* pIndexMap);
     virtual ~SiconosMatrixHandler();

#ifdef DEBUG
     virtual void IsValid() const override final;
#endif

     virtual void Resize(integer, integer) override final;

     virtual void ResizeReset(integer, integer) override final;

     virtual void Reset() override final;

     virtual void
     PutCoef(integer iRow, integer iCol, const doublereal& dCoef) override final;

     virtual void
     IncCoef(integer iRow, integer iCol, const doublereal& dCoef) override final;

     virtual void
     DecCoef(integer iRow, integer iCol, const doublereal& dCoef) override final;

     virtual const doublereal&
     operator() (integer iRow, integer iCol) const override final;

     virtual doublereal&
     operator() (integer iRow, integer iCol) override final;

     virtual integer iGetNumRows() const override final;
     virtual integer iGetNumCols() const override final;

     virtual std::ostream& Print(std::ostream& os, MatPrintFormat eFormat) const override final;

     typedef void EnumerateNzCallback(integer, integer, doublereal);

     virtual void EnumerateNz(const std::function<EnumerateNzCallback>& func) const override final;

     virtual doublereal ConditionNumber(enum Norm_t eNorm = NORM_1) const override final;

     virtual doublereal Norm(enum Norm_t eNorm = NORM_1) const override final;

     virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale) override final;

     virtual MatrixHandler* Copy() const override final;

     virtual bool AddItem(integer iRow, const sp_grad::SpGradient& oItem) override final;

     virtual VectorHandler&
     MatVecMul(VectorHandler& y, const VectorHandler& x) const final;

     virtual VectorHandler&
     MatTVecMul(VectorHandler& y, const VectorHandler& x) const final;

     virtual VectorHandler&
     MatVecIncMul(VectorHandler& y, const VectorHandler& x) const final;

     virtual VectorHandler&
     MatTVecIncMul(VectorHandler& y, const VectorHandler& x) const final;

     virtual VectorHandler&
     MatVecDecMul(VectorHandler& y, const VectorHandler& x) const final;

     virtual VectorHandler&
     MatTVecDecMul(VectorHandler& y, const VectorHandler& x) const final;

     NumericsMatrix* pGetMatrix() { return pMatrix; }

     void SetIndexMap(const SiconosIndexMap* pIndexMapNew);

     const SiconosIndexMap* pGetIndexMap() const { return pIndexMap; }

private:
     virtual VectorHandler&
     MatVecMul_base(void (VectorHandler::*op)(integer iRow,
                                              const doublereal& dCoef),
                    VectorHandler& out, const VectorHandler& in) const override;
     virtual VectorHandler&
     MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
                                               const doublereal& dCoef),
                     VectorHandler& out, const VectorHandler& in) const override;
private:
     const integer iNumNz;
     NumericsMatrix* pMatrix;
     const SiconosIndexMap* pIndexMap;
};

#endif
#endif
