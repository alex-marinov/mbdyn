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
        Copyright (C) 2020(-2020) all rights reserved.

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
     class SpGradient: public SpGradBase<SpGradient> {
     public:
	  friend util::SpMatrixDataTraits<SpGradient>;
                
	  static constexpr ExprEvalFlags eExprEvalFlags = ExprEvalUncompressed;
                
	  inline SpGradient();

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

	  inline void ResizeReset(doublereal dVal, index_type iSize);

	  inline void SetValuePreserve(doublereal dVal);

	  static inline void SetValuePreserve(SpGradient& g, doublereal dVal);

	  static inline void SetValuePreserve(doublereal& g, doublereal dVal);

	  static inline void ResizeReset(SpGradient& g, doublereal dVal, index_type iSize);
                
	  static inline void ResizeReset(doublereal& g, doublereal dVal, index_type iSize);

	  template <typename Expr>
	  inline constexpr bool bHaveRefTo(const SpGradBase<Expr>&) const;

	  inline bool bHaveRefTo(const SpGradBase<SpGradient>& g) const;

	  template <typename Expr, index_type NumRows, index_type NumCols>
	  static constexpr inline bool bHaveRefTo(const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

	  template <index_type NumRows, index_type NumCols>
	  inline bool bHaveRefTo(const SpMatrixBase<SpGradient, NumRows, NumCols>& A) const;

	  template <index_type NumRows, index_type NumCols>
	  static inline bool bHaveRefTo(const SpGradient& g, const SpMatrixBase<SpGradient, NumRows, NumCols>& A) { return g.bHaveRefTo(A); }

	  template <typename Expr, index_type NumRows, index_type NumCols>
	  static constexpr inline bool bHaveRefTo(doublereal g, const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

	  inline doublereal dGetValue() const;

	  inline doublereal dGetDeriv(index_type iDof) const;

	  inline void InsertDeriv(SpGradient& g, doublereal dCoef) const;

	  inline void InsertDof(SpGradExpDofMap& oExpDofMap) const;

	  inline static void InsertDof(const SpGradient& g, SpGradExpDofMap& oDofMap);

	  inline static void InsertDof(doublereal, SpGradExpDofMap&);

	  inline void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) const;

	  static inline void AddDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap);

	  static inline void AddDeriv(doublereal f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) {}

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

	  inline void MakeUnique();

	  template <typename Expr>
	  static constexpr inline doublereal
	  dGetValue(const SpGradBase<Expr>& a);

	  static constexpr inline doublereal
	  dGetValue(doublereal a);

	  template <typename Expr>
	  static constexpr inline index_type
	  iGetSize(const SpGradBase<Expr>& a);

	  static constexpr inline index_type
	  iGetSize(doublereal a);

	  static inline void InsertDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef);
	  static inline void InsertDeriv(const doublereal& f, SpGradient& g, doublereal dCoef);
	  static inline void InsertDeriv(const doublereal& f, doublereal& g, doublereal dCoef);

	  static inline void Compress(doublereal);

	  static inline void Compress(SpGradient& g);

	  inline static void GetDofStat(const SpGradient& g, SpGradDofStat& s);

	  inline static void GetDofStat(doublereal, SpGradDofStat&);
                
	  inline static doublereal dGetDeriv(const SpGradient&g, index_type iDof);

	  inline static constexpr doublereal dGetDeriv(doublereal, index_type);

	  template <typename Expr>
	  inline void Assign(const SpGradBase<Expr>& g);

	  template <typename Expr>
	  inline void MapAssign(const SpGradBase<Expr>& g);

	  inline void InitDeriv(const SpGradExpDofMap& oExpDofMap);
                
	  void Compress();
#ifdef DEBUG
	  bool bValid() const;
	  void PrintValue(std::ostream& os) const;
	  void PrintDeriv(std::ostream& os, doublereal dCoef) const;
#endif
     private:
	  enum AllocFlags {
	       ALLOC_RESIZE,
	       ALLOC_UNIQUE
	  };
                
	  inline explicit SpGradient(SpDerivData* pData);

	  inline static constexpr size_t uGetAllocSize(index_type iSizeRes);

	  inline static SpDerivData* pAllocMem(SpDerivData* ptr, index_type iSize);

	  void Allocate(index_type iSizeRes, index_type iSizeInit, AllocFlags eAllocFlags);

	  inline void Free();

	  inline constexpr static bool bRecCompareWithDof(const SpDerivRec& a, index_type b);

	  inline static doublereal AssignAdd(doublereal a, doublereal b, doublereal& df_db);

	  inline static doublereal AssignSub(doublereal a, doublereal b, doublereal& df_db);

	  inline static doublereal AssignMul(doublereal a, doublereal b, doublereal& df_da, doublereal& df_db);

	  inline static doublereal AssignDiv(doublereal a, doublereal b, doublereal& df_da, doublereal& df_db);

	  inline constexpr static doublereal AssignMulConst(doublereal a, doublereal b);

	  inline constexpr static doublereal AssignDivConst(doublereal a, doublereal b);

	  inline void MaybeCompress() const;

	  inline SpDerivRec* pFindRec(index_type iDof) const;

	  inline SpDerivRec* pInsertRec(index_type iDof, doublereal dDer);

	  template <typename CONT_TYPE>
	  inline void CopyDeriv(doublereal dVal, const CONT_TYPE& rgDer);

	  template <doublereal AssOpFn(doublereal, doublereal, doublereal&), typename Expr>
	  inline void AssignOper(const SpGradBase<Expr>& g);

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
	  static SP_GRAD_THREAD_LOCAL SpDerivData oNullData;
     };
}
#endif
