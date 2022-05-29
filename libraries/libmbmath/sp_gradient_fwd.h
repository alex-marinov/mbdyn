/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2020(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __SP_GRADIENT_FWD_H__INCLUDED__
#define __SP_GRADIENT_FWD_H__INCLUDED__

#include <initializer_list>
#include <vector>

#include "sp_gradient_base.h"
#include "sp_matrix_base_fwd.h"

namespace sp_grad {
     template <typename T>
     class SpGradientTraits;

     template <>
     struct SpGradientTraits<doublereal> {
          static inline void ResizeReset(doublereal& g, doublereal dVal, index_type iSize);
          inline static void InsertDof(doublereal, SpGradExpDofMap&);
          static inline void AddDeriv(doublereal f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) {}

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(doublereal g, const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          static constexpr inline doublereal
          dGetValue(doublereal a);
          static constexpr inline index_type
          iGetSize(doublereal a);
          static inline void InsertDeriv(const doublereal& f, doublereal& g, doublereal dCoef) noexcept {}
          static inline void Sort(doublereal);
          static bool bIsUnique(doublereal) { return true; }
          inline static void GetDofStat(doublereal, SpGradDofStat&);
          inline static constexpr doublereal dGetDeriv(doublereal, index_type);
          inline static void ZeroInit(doublereal& g) { g = 0.; }
     };

     template <>
     struct SpGradientTraits<SpGradient> {
          static inline void ResizeReset(SpGradient& g, doublereal dVal, index_type iSize);

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const SpGradient& g, const SpMatrixBase<SpGradient, NumRows, NumCols>& A);

          static inline void InsertDof(const SpGradient& g, SpGradExpDofMap& oDofMap);

          static inline void AddDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap);

          template <typename Expr>
          static constexpr inline doublereal
          dGetValue(const SpGradBase<Expr>& a);

          template <typename Expr>
          static inline index_type
          iGetSize(const SpGradBase<Expr>& a);

          static inline index_type
          iGetSize(const SpGradient& a);

          static inline void InsertDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef);
          static inline void InsertDeriv(const doublereal& f, SpGradient& g, doublereal dCoef) noexcept {}

          inline static void Sort(SpGradient& g);

          inline static void GetDofStat(const SpGradient& g, SpGradDofStat& s);

          inline static doublereal dGetDeriv(const SpGradient&g, index_type iDof);

          inline static bool bIsUnique(const SpGradient& g);

          inline static void ZeroInit(SpGradient&) {}
     };

     template <>
     struct SpGradientTraits<GpGradProd> {
          static inline doublereal
          dGetValue(const GpGradProd& a);

          static inline index_type
          iGetSize(const GpGradProd& a) { return 1; }

          static inline void InsertDeriv(const GpGradProd& f, GpGradProd& g, doublereal dCoef);
          static inline void InsertDeriv(const doublereal& f, GpGradProd& g, doublereal dCoef) noexcept {}
          static inline void ResizeReset(GpGradProd& g, doublereal dVal, index_type);

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const GpGradProd& g, const SpMatrixBase<GpGradProd, NumRows, NumCols>& A);

          inline static bool bIsUnique(const GpGradProd& g) { return true; }
          inline static void ZeroInit(GpGradProd&) {}
     };

     class SpGradient: public SpGradBase<SpGradient> {
     public:
          friend util::SpMatrixDataTraits<SpGradient>;

          static constexpr ExprEvalFlags eExprEvalFlags = ExprEvalDuplicate;

          inline SpGradient();

          inline explicit SpGradient(doublereal d);

          inline SpGradient(const SpGradient& g);

          inline SpGradient(SpGradient&& g);

          inline SpGradient(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer);

          inline SpGradient(doublereal dVal, const std::vector<SpDerivRec>& rgDer);

          template <typename Expr>
          inline SpGradient(const SpGradBase<Expr>& g);

          inline ~SpGradient();

          inline SpGradient& operator=(const SpGradient& g);

          inline SpGradient& operator=(SpGradient&& g);

          template <typename Expr>
          inline SpGradient& operator=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator+=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator-=(const SpGradBase<Expr>& g);

          inline SpGradient& operator+=(doublereal b);

          inline SpGradient& operator-=(doublereal b);

          template <typename Expr>
          inline SpGradient& operator*=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator/=(const SpGradBase<Expr>& g);

          inline SpGradient& operator*=(doublereal b);

          inline SpGradient& operator/=(doublereal b);

          inline void Reset(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer);

          inline void Reset(doublereal dVal, const std::vector<SpDerivRec>& rgDer);

          inline void Reset(doublereal dVal, index_type iDof, doublereal dDer);

          inline void ResetNumeric();

          inline void ResizeReset(doublereal dVal, index_type iSize);

          inline void Scale(doublereal dRowScale, const std::vector<doublereal>& oColScale);

          template <typename Expr>
          inline constexpr bool bHaveRefTo(const SpGradBase<Expr>&) const;

          inline bool bHaveRefTo(const SpGradBase<SpGradient>& g) const;

          template <index_type NumRows, index_type NumCols>
          inline bool bHaveRefTo(const SpMatrixBase<SpGradient, NumRows, NumCols>& A) const;

          inline doublereal dGetValue() const;

          inline doublereal dGetDeriv(index_type iDof) const;

          inline void InsertDeriv(SpGradient& g, doublereal dCoef) const;

          inline void InsertDof(SpGradExpDofMap& oExpDofMap) const;

          inline void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) const;

          inline const SpDerivRec* begin() const;

          inline const SpDerivRec* end() const;

          inline index_type iGetSize() const;

          inline void GetDofStat(SpGradDofStat& s) const;

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset);

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset,
                          const SpGradExpDofMap& oDofMap);

          template <typename AITER, typename BITER>
          inline void
          InnerProduct(AITER pAFirst,
                       AITER pALast,
                       index_type iAOffset,
                       BITER pBFirst,
                       BITER pBLast,
                       index_type iBOffset);

          template <typename Expr>
          inline void Assign(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline void MapAssign(const SpGradBase<Expr>& g);

          template <typename Func, typename Expr>
          inline void AssignOper(const SpGradBase<Expr>& g);

          template <typename Func, typename Expr>
          inline void MapAssignOper(const SpGradBase<Expr>& g);

          template <typename Func>
          inline void InitDerivAssign(doublereal f, doublereal df_du, const SpGradExpDofMap& oExpDofMap);

          template <typename Func>
          inline void InitDerivAssign(doublereal f, doublereal df_du, index_type iSizeRes);

          inline void InitDeriv(const SpGradExpDofMap& oExpDofMap);

          void Sort();

          inline bool bIsSorted() const;
          inline bool bIsUnique() const;

          void MakeUnique();

#ifdef SP_GRAD_DEBUG
          bool bValid() const;
          bool bCheckUnique() const;
          void PrintValue(std::ostream& os) const;
          void PrintDeriv(std::ostream& os, doublereal dCoef) const;
          static index_type iGetRefCntNullData() { return pGetNullData()->iRefCnt; }
#endif
     private:
          inline void UniqueOwner();

          inline explicit SpGradient(SpDerivData* pData);

          inline void SetValuePreserve(doublereal dVal);

          inline static size_t uGetAllocSize(index_type iSizeRes);

          inline static SpDerivData* pAllocMem(SpDerivData* ptr, index_type iSize);

          void Allocate(index_type iSizeRes, index_type iSizeInit, unsigned uFlags);

          inline void Free();
          inline void Cleanup();

          inline constexpr static bool bRecCompareWithDof(const SpDerivRec& a, index_type b);

          inline constexpr static doublereal AssignMulConst(doublereal a, doublereal b);

          inline constexpr static doublereal AssignDivConst(doublereal a, doublereal b);

          inline void MaybeSort() const;

          inline SpDerivRec* pFindRec(index_type iDof) const;

          inline SpDerivRec* pInsertRec(index_type iDof, doublereal dDer);

          template <typename CONT_TYPE>
          inline void CopyDeriv(doublereal dVal, const CONT_TYPE& rgDer);

          template <doublereal AssOpFn(doublereal, doublereal, doublereal&, doublereal&), typename Expr>
          inline void AssignOper(const SpGradBase<Expr>& g);

          template<doublereal AssOpFn(doublereal, doublereal)>
          inline void AssignOper(doublereal b);

          inline void InnerProductAddDer(const SpGradient& g, const doublereal dVal);

          inline void InnerProductAddDer(doublereal, doublereal) {}

          template <typename ITER>
          inline static index_type InnerProductSize(ITER pFirst,
                                                    ITER pLast,
                                                    index_type iOffset);

          template <typename ITER>
          inline static void InnerProductInsertDof(ITER pFirst,
                                                   ITER pLast,
                                                   index_type iOffset,
                                                   SpGradExpDofMap& oDofMap);

          inline void InnerProductAddDer(const SpGradient& g,
                                         doublereal dVal,
                                         const SpGradExpDofMap& oDofMap);

          inline static void InnerProductAddDer(doublereal,
                                                doublereal,
                                                const SpGradExpDofMap&) {}

          template <typename ITER>
          inline static void InnerProductDofStat(ITER pFirst,
                                                 ITER pLast,
                                                 index_type iOffset,
                                                 SpGradDofStat& s);

          inline static SpDerivData* pGetNullData();

          SpDerivData* pData;
          static SpDerivData oNullData;
     };

     class GpGradProd: public GpGradProdBase<GpGradProd> {
     public:
          constexpr explicit GpGradProd(doublereal dVal = 0., doublereal dDer = 0.)
               :dVal(dVal), dDer(dDer) {
          }

          void Reset(doublereal dNewVal = 0., doublereal dNewDer = 0.) {
               dVal = dNewVal;
               dDer = dNewDer;
          }

          template <typename Expr>
          inline constexpr GpGradProd(const GpGradProdBase<Expr>& oExpr);

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset);

          constexpr doublereal dGetValue() const {
               return dVal;
          }

          static inline constexpr doublereal dGetValue(doublereal dVal) {
               return dVal;
          }

          static inline doublereal dGetValue(const GpGradProd& g) {
               return g.dGetValue();
          }

          constexpr doublereal dGetDeriv() const {
               return dDer;
          }

          GpGradProd& operator+=(doublereal dCoef) {
               dVal += dCoef;

               return *this;
          }

          GpGradProd& operator-=(doublereal dCoef) {
               dVal -= dCoef;

               return *this;
          }

          GpGradProd& operator*=(doublereal dCoef) {
               dVal *= dCoef;
               dDer *= dCoef;

               return *this;
          }

          GpGradProd& operator/=(doublereal dCoef) {
               dVal /= dCoef;
               dDer /= dCoef;

               return *this;
          }

          template <typename Expr>
          inline GpGradProd& operator+=(const GpGradProdBase<Expr>& oExpr);

          template <typename Expr>
          inline GpGradProd& operator-=(const GpGradProdBase<Expr>& oExpr);

          template <typename Expr>
          inline GpGradProd& operator*=(const GpGradProdBase<Expr>& oExpr);

          template <typename Expr>
          inline GpGradProd& operator/=(const GpGradProdBase<Expr>& oExpr);

          static constexpr bool bIsScalarConst = false;

          inline void InsertDeriv(GpGradProd& g, doublereal dCoef) const;

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const SpMatrixBase<GpGradProd, NumRows, NumCols>& A) {
               return false;
          }

#ifdef SP_GRAD_DEBUG
          inline bool bValid() const;
#endif
     private:
          template <typename BinFunc, typename Expr>
          inline GpGradProd& AssignOper(const GpGradProdBase<Expr>& oExpr);

          inline void InnerProductAddDer(const GpGradProd& g, const doublereal dCoef);

          inline void InnerProductAddDer(doublereal, doublereal) {}

          doublereal dVal;
          doublereal dDer;
     };
}
#endif
