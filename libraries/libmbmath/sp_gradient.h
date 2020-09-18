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

#ifndef __SP_GRADIENT_H__INCLUDED__
#define __SP_GRADIENT_H__INCLUDED__

#include <algorithm>

#include "sp_gradient_base.h"
#include "sp_exp_dof_map.h"
#include "sp_gradient_fwd.h"
#include "sp_gradient_util.h"
#include "sp_gradient_expr.h"
#include "sp_gradient_func.h"
#include "sp_gradient_op.h"
#include "sp_matrix_base.h"

namespace sp_grad {
     SpGradient::SpGradient(SpDerivData* p)
	  :pData(p)
     {
	  SP_GRAD_ASSERT(p != nullptr);

	  ++pData->iRefCnt;

	  if (pData->pOwner) {
	       pData->pOwner->Attach(this);
	  }

	  SP_GRAD_ASSERT(bValid());
     }

     SpGradient::SpGradient()
	  :SpGradient(pGetNullData()) {

	  SP_GRAD_ASSERT(bValid());
     }

     SpGradient::SpGradient(const SpGradient& g)
	  :SpGradient(g.pData) {

	  SP_GRAD_ASSERT(g.bValid());
     }

     SpGradient::SpGradient(SpGradient&& g)
	  :SpGradient() {

	  SP_GRAD_ASSERT(g.bValid());

	  *this = std::move(g);

	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());
     }

     SpGradient::SpGradient(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer)
	  :SpGradient() {

	  CopyDeriv(dVal, rgDer);

	  SP_GRAD_ASSERT(bValid());
     }

     SpGradient::SpGradient(doublereal dVal, const std::vector<SpDerivRec>& rgDer)
	  :SpGradient() {

	  CopyDeriv(dVal, rgDer);

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Expr>
     SpGradient::SpGradient(const SpGradBase<Expr>& g)
	  :SpGradient() {

	  // Determine at compile time if we are using a compressed or uncompressed evaluation
	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;

	  EvalHelp::Eval(*this, g);

	  SP_GRAD_ASSERT(bValid());
     }

     SpGradient::~SpGradient() {
	  SP_GRAD_ASSERT(bValid());

	  Free();
     }

     SpGradient& SpGradient::operator=(const SpGradient& g) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());

	  Assign(g);

	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator=(SpGradient&& g) {
	  SP_GRAD_ASSERT(g.bValid());
	  SP_GRAD_ASSERT(bValid());

	  std::swap(pData, g.pData);

	  if (pData->pOwner) {
	       pData->pOwner->Attach(this);
	       pData->pOwner->Detach(&g);
	  }

	  if (g.pData->pOwner) {
	       g.pData->pOwner->Attach(&g);
	       g.pData->pOwner->Detach(this);
	  }

	  SP_GRAD_ASSERT(g.bValid());
	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  // Determine at compile time if we are using a compressed or uncompressed evaluation
	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;

	  EvalHelp::Eval(*this, g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator+=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  AssignOper<AssignAdd>(g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator-=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  AssignOper<AssignSub>(g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator+=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  MakeUnique();

	  pData->dVal += b;

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator-=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  MakeUnique();

	  pData->dVal -= b;

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator*=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  AssignOper<AssignMul>(g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator/=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  AssignOper<AssignDiv>(g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator*=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  MakeUnique();

	  AssignOper<AssignMulConst>(b);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator/=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  MakeUnique();

	  AssignOper<AssignDivConst>(b);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     void SpGradient::Reset(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer) {
	  SP_GRAD_ASSERT(bValid());

	  CopyDeriv(dVal, rgDer);

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::Reset(doublereal dVal, const std::vector<SpDerivRec>& rgDer) {
	  SP_GRAD_ASSERT(bValid());

	  CopyDeriv(dVal, rgDer);

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::Reset(doublereal dVal, index_type iDof, doublereal dDer) {
	  SP_GRAD_ASSERT(bValid());

	  Allocate(1, 1, ALLOC_RESIZE);
	  pData->dVal = dVal;
	  pData->rgDer[0].iDof = iDof;
	  pData->rgDer[0].dDer = dDer;
	  pData->bCompressed = true;

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::ResizeReset(doublereal dVal, index_type iSize) {
	  SP_GRAD_ASSERT(bValid());

	  Allocate(iSize, 0, ALLOC_RESIZE);
	  pData->dVal = dVal;

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::ResizeReset(SpGradient& g, doublereal dVal, index_type iSize) {
	  g.ResizeReset(dVal, iSize);
     }

     void SpGradient::ResizeReset(doublereal& g, doublereal dVal, index_type) {
	  g = dVal;
     }

     void SpGradient::SetValuePreserve(doublereal dVal) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  pData->dVal = dVal;
	  pData->bCompressed = false;
     }

     void SpGradient::SetValuePreserve(SpGradient& g, doublereal dVal) {
	  g.SetValuePreserve(dVal);
     }

     void SpGradient::SetValuePreserve(doublereal& g, doublereal dVal) {
	  g = dVal;
     }

     template <typename Expr>
     constexpr bool SpGradient::bHaveRefTo(const SpGradBase<Expr>&) const {
	  SP_GRAD_ASSERT(bValid());

	  return false;
     }

     bool SpGradient::bHaveRefTo(const SpGradBase<SpGradient>& g) const {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.pGetRep()->bValid());

	  return g.pGetRep()->pData == pData && pData != pGetNullData();
     }

     template <index_type NumRows, index_type NumCols>
     bool SpGradient::bHaveRefTo(const SpMatrixBase<SpGradient, NumRows, NumCols>& A) const {
	  return A.bIsOwnerOf(pData->pOwner);
     }

     doublereal SpGradient::dGetValue() const {
	  SP_GRAD_ASSERT(bValid());

	  return pData->dVal;
     }

     doublereal SpGradient::dGetDeriv(index_type iDof) const {
	  SP_GRAD_ASSERT(bValid());

	  MaybeCompress();

	  auto pRec = pFindRec(iDof);

	  SP_GRAD_ASSERT(bValid());

	  return pRec ? pRec->dDer : 0.;
     }

     void SpGradient::InsertDeriv(SpGradient& g, doublereal dCoef) const {
	  for (const auto& r: *this) {
	       g.pInsertRec(r.iDof, r.dDer * dCoef);
	  }
     }

     void SpGradient::InsertDof(SpGradExpDofMap& oExpDofMap) const {
	  for (const auto& r: *this) {
	       oExpDofMap.InsertDof(r.iDof);
	  }
     }

     void SpGradient::InsertDof(const SpGradient& g, SpGradExpDofMap& oDofMap) {
	  g.InsertDof(oDofMap);
     }

     void SpGradient::InsertDof(doublereal, SpGradExpDofMap&) {
     }

     void SpGradient::AddDeriv(SpGradient& g, const doublereal dCoef, const SpGradExpDofMap& oDofMap) const {
	  SpDerivRec* const p = g.pData->rgDer;

	  SP_GRAD_ASSERT(g.pData != pGetNullData());
	  SP_GRAD_ASSERT(g.pData->iSizeCurr == oDofMap.iGetLocalSize());
	  
	  for (const auto& r: *this) {
	       const index_type i = oDofMap.iGetLocalIndex(r.iDof);

	       SP_GRAD_ASSERT(i >= 0);
	       SP_GRAD_ASSERT(i < g.pData->iSizeCurr);
	       SP_GRAD_ASSERT(p[i].iDof == r.iDof);

	       p[i].dDer += dCoef * r.dDer;
	  }
     }

     void SpGradient::AddDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) {
	  f.AddDeriv(g, dCoef, oDofMap);
     }

     const SpDerivRec* SpGradient::begin() const {
	  SP_GRAD_ASSERT(bValid());

	  return pData->rgDer;
     }

     const SpDerivRec* SpGradient::end() const {
	  SP_GRAD_ASSERT(bValid());

	  return pData->rgDer + pData->iSizeCurr;
     }

     index_type SpGradient::iGetSize() const {
	  SP_GRAD_ASSERT(bValid());

	  return pData->iSizeCurr;
     }

     void SpGradient::GetDofStat(SpGradDofStat& s) const {
	  for (const auto& r: *this) {
	       s.iMinDof = std::min(s.iMinDof, r.iDof);
	       s.iMaxDof = std::max(s.iMaxDof, r.iDof);
	       ++s.iNumNz;
	  }
     }

     template <typename AITER, typename BITER>
     void SpGradient::MapInnerProduct(AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset) {
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);
	  SP_GRAD_ASSERT((pALast - pAFirst) % iAOffset == 0);
	  SP_GRAD_ASSERT((pBLast - pBFirst) % iBOffset == 0);

	  SpGradDofStat s;

	  InnerProductDofStat(pAFirst, pALast, iAOffset, s);
	  InnerProductDofStat(pBFirst, pBLast, iBOffset, s);

	  SpGradExpDofMap oDofMap(s);

	  InnerProductInsertDof(pAFirst, pALast, iAOffset, oDofMap);
	  InnerProductInsertDof(pBFirst, pBLast, iBOffset, oDofMap);

	  oDofMap.InsertDone();

	  MapInnerProduct(pAFirst, pALast, iAOffset, pBFirst, pBLast, iBOffset, oDofMap);
     }

     template <typename AITER, typename BITER>
     void SpGradient::MapInnerProduct(AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset, const SpGradExpDofMap& oDofMap) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);
	  SP_GRAD_ASSERT((pALast - pAFirst) % iAOffset == 0);
	  SP_GRAD_ASSERT((pBLast - pBFirst) % iBOffset == 0);

	  Allocate(oDofMap.iGetLocalSize(), 0, ALLOC_RESIZE);

	  SetValuePreserve(0.);

	  SP_GRAD_ASSERT(bValid());

	  InitDeriv(oDofMap);

	  while (pAFirst < pALast) {
	       const auto& ai = *pAFirst;
	       const auto& bi = *pBFirst;
	       const doublereal aiv = dGetValue(ai);
	       const doublereal biv = dGetValue(bi);

	       pData->dVal += aiv * biv;

	       InnerProductAddDer(ai, biv, oDofMap);
	       InnerProductAddDer(bi, aiv, oDofMap);

	       pAFirst += iAOffset;
	       pBFirst += iBOffset;
	  }

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename AITER, typename BITER>
     void SpGradient::InnerProduct(AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);

	  index_type iNumNz;

	  iNumNz = InnerProductSize(pAFirst, pALast, iAOffset);
	  iNumNz += InnerProductSize(pBFirst, pBLast, iBOffset);

	  Allocate(iNumNz, 0, ALLOC_RESIZE);

	  SetValuePreserve(0.);

	  while (pAFirst != pALast) {
	       const auto& ai = *pAFirst;
	       const auto& bi = *pBFirst;
	       const doublereal aiv = dGetValue(ai);
	       const doublereal biv = dGetValue(bi);

	       pData->dVal += aiv * biv;

	       InnerProductAddDer(ai, biv);
	       InnerProductAddDer(bi, aiv);

	       pAFirst += iAOffset;
	       pBFirst += iBOffset;
	  }

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::MakeUnique() {
	  Allocate(pData->iSizeRes, pData->iSizeCurr, ALLOC_UNIQUE);
     }


     template <typename Expr>
     constexpr doublereal
     SpGradient::dGetValue(const SpGradBase<Expr>& a) {
	  return a.dGetValue();
     }

     constexpr doublereal
     SpGradient::dGetValue(doublereal a) {
	  return a;
     }

     template <typename Expr>
     constexpr index_type
     SpGradient::iGetSize(const SpGradBase<Expr>& a) {
	  return a.iGetSize();
     }

     constexpr index_type
     SpGradient::iGetSize(doublereal a) {
	  return 0;
     }

     void SpGradient::InsertDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef) {
	  f.InsertDeriv(g, dCoef);
     }

     void SpGradient::InsertDeriv(const doublereal&, SpGradient&, doublereal) {
     }

     void SpGradient::InsertDeriv(const doublereal&, doublereal&, doublereal) {
     }

     void SpGradient::Compress(doublereal) {
     }

     void SpGradient::Compress(SpGradient& g) {
	  g.Compress();
     }

     void SpGradient::GetDofStat(const SpGradient& g, SpGradDofStat& s) {
	  g.GetDofStat(s);
     }

     void SpGradient::GetDofStat(doublereal, SpGradDofStat&) {
     }

     doublereal SpGradient::dGetDeriv(const SpGradient&g, index_type iDof) {
	  return g.dGetDeriv(iDof);
     }

     constexpr doublereal SpGradient::dGetDeriv(doublereal, index_type) {
	  return 0.;
     }

     constexpr size_t SpGradient::uGetAllocSize(index_type iSizeRes) {
	  return offsetof(SpDerivData, rgDer[iSizeRes]);
     }

     void SpGradient::Free() {
	  --pData->iRefCnt;

	  SP_GRAD_ASSERT(pData->iRefCnt >= 0);

	  if (pData->pOwner) {
	       pData->pOwner->Detach(this);
	  } else if (!pData->iRefCnt) {
	       std::free(pData);
	  }

	  pData = pGetNullData();

	  SP_GRAD_ASSERT(pData->iRefCnt >= 0);

	  ++pData->iRefCnt;
     }

     constexpr  bool SpGradient::bRecCompareWithDof(const SpDerivRec& a, index_type b) {
	  return a.iDof < b;
     }

     doublereal SpGradient::AssignAdd(doublereal a, doublereal b, doublereal& df_db) {
	  df_db = 1.;
	  return a + b;
     }

     doublereal SpGradient::AssignSub(doublereal a, doublereal b, doublereal& df_db) {
	  df_db = -1.;
	  return a - b;
     }

     doublereal SpGradient::AssignMul(doublereal a, doublereal b, doublereal& df_da, doublereal& df_db) {
	  df_da = b;
	  df_db = a;
	  return a * b;
     }

     doublereal SpGradient::AssignDiv(doublereal a, doublereal b, doublereal& df_da, doublereal& df_db) {
	  df_da = 1. / b;
	  df_db = -a / (b * b);
	  return a / b;
     }

     constexpr  doublereal SpGradient::AssignMulConst(doublereal a, doublereal b) {
	  return a * b;
     }

     constexpr  doublereal SpGradient::AssignDivConst(doublereal a, doublereal b) {
	  return a / b;
     }

     void SpGradient::MaybeCompress() const {
	  if (!bIsCompressed()) {
	       const_cast<SpGradient*>(this)->Compress();
	  }

	  SP_GRAD_ASSERT(bIsCompressed());
     }

     bool SpGradient::bIsCompressed() const {
	  return pData->bCompressed;
     }

     SpDerivRec* SpGradient::pFindRec(index_type iDof) const {
	  SP_GRAD_ASSERT(bIsCompressed());

	  auto pBegin = pData->rgDer;
	  auto pEnd = pData->rgDer + pData->iSizeCurr;
	  auto pRec = std::lower_bound(pBegin, pEnd, iDof, bRecCompareWithDof);

	  if (pRec == pEnd || pRec->iDof != iDof) {
	       return nullptr;
	  }

	  return pRec;
     }

     SpDerivRec* SpGradient::pInsertRec(index_type iDof, doublereal dDer) {
	  SP_GRAD_ASSERT(pData->iSizeCurr < pData->iSizeRes);

	  auto p = pData->rgDer + pData->iSizeCurr++;

	  p->iDof = iDof;
	  p->dDer = dDer;

	  SP_GRAD_ASSERT(pData->iSizeCurr <= pData->iSizeRes);

	  return p;
     }

     template <typename CONT_TYPE>
     void SpGradient::CopyDeriv(doublereal dVal, const CONT_TYPE& rgDer) {
	  index_type iSize = rgDer.size();

	  if (iSize > 0) {
	       Allocate(iSize, iSize, ALLOC_RESIZE);

	       pData->dVal = dVal;

	       std::uninitialized_copy(std::begin(rgDer), std::end(rgDer), pData->rgDer);

	       SP_GRAD_ASSERT(!bIsCompressed());
	  }
     }

     template <typename Expr>
     void SpGradient::Assign(const SpGradBase<Expr>& g) {
	  SpGradient f;

	  if (!g.bHaveRefTo(*this)) {
	       f = std::move(*this);
	  }

	  f.Allocate(g.iGetSize(), 0, ALLOC_RESIZE);

	  f.pData->dVal = g.dGetValue();

	  g.InsertDeriv(f, 1.);

	  *this = std::move(f);

	  SP_GRAD_ASSERT(!bIsCompressed());
     }

     template <typename Expr>
     void SpGradient::MapAssign(const SpGradBase<Expr>& g) {
	  SpGradient f;

	  if (!g.bHaveRefTo(*this)) {
	       f = std::move(*this);
	  }

	  SpGradDofStat s;

	  g.GetDofStat(s);

	  SpGradExpDofMap oDofMap(s);

	  g.InsertDof(oDofMap);

	  oDofMap.InsertDone();

	  f.Allocate(oDofMap.iGetLocalSize(), 0, ALLOC_RESIZE);

	  f.pData->dVal = g.dGetValue();

	  f.InitDeriv(oDofMap);

	  g.AddDeriv(f, 1., oDofMap);

	  *this = std::move(f);

	  SP_GRAD_ASSERT(!bIsCompressed());
	  SP_GRAD_ASSERT(bIsUnique());
     }

     void SpGradient::InitDeriv(const SpGradExpDofMap& oExpDofMap) {
	  SP_GRAD_ASSERT(pData != pGetNullData());       
	  
	  for (index_type i = 0; i < oExpDofMap.iGetLocalSize(); ++i) {
	       pInsertRec(oExpDofMap.iGetGlobalIndex(i), 0.);
	  }
     }

     template <doublereal AssOpFn(doublereal, doublereal, doublereal&), typename Expr>
     void SpGradient::AssignOper(const SpGradBase<Expr>& g) {
	  const doublereal a = dGetValue();
	  const doublereal b = g.dGetValue();
	  // Assume df_da == 1
	  doublereal df_db;
	  const doublereal f = AssOpFn(a, b, df_db);
	  const index_type iSizePrev = iGetSize();

	  SpGradient r;

	  if (!g.bHaveRefTo(*this)) {
	       r = std::move(*this);
	  } else {
	       r = *this;
	  }

	  r.Allocate(iSizePrev + g.iGetSize(), iSizePrev, ALLOC_RESIZE);

	  r.pData->dVal = f;

	  g.InsertDeriv(r, df_db);

	  *this = std::move(r);
     }

     template <doublereal AssOpFn(doublereal, doublereal, doublereal&, doublereal&), typename Expr>
     void SpGradient::AssignOper(const SpGradBase<Expr>& g) {
	  const doublereal a = dGetValue();
	  const doublereal b = g.dGetValue();
	  doublereal df_da, df_db;
	  const doublereal f = AssOpFn(a, b, df_da, df_db);

	  SpGradient r;

	  r.Allocate(iGetSize() + g.iGetSize(), 0, ALLOC_RESIZE);

	  r.pData->dVal = f;

	  InsertDeriv(r, df_da);
	  g.InsertDeriv(r, df_db);

	  *this = std::move(r);
     }

     template<doublereal AssOpFn(doublereal, doublereal)>
     void SpGradient::AssignOper(doublereal b) {
	  if (pData == pGetNullData()) {
	       // f = 0 * b
	       // f' = 0' * b + 0 * b'
	       return;
	  }

	  const doublereal a = pData->dVal;
	  const doublereal f = AssOpFn(a, b);
	  const doublereal df_da = b;

	  pData->dVal = f;

	  for (SpDerivRec* p = pData->rgDer; p < pData->rgDer + pData->iSizeCurr; ++p) {
	       p->dDer = AssOpFn(p->dDer, df_da);
	  }
     }

     void SpGradient::InnerProductAddDer(const SpGradient& g, const doublereal dVal) {
	  auto pFirst = g.begin();
	  auto pLast = g.end();

	  while (pFirst < pLast) {
	       pInsertRec(pFirst->iDof, dVal * pFirst->dDer);
	       ++pFirst;
	  }
     }

     template <typename ITER>
     void SpGradient::InnerProductInsertDof(ITER pFirst, ITER pLast, index_type iOffset, SpGradExpDofMap& oDofMap) {
	  while (pFirst < pLast) {
	       InsertDof(*pFirst, oDofMap);
	       pFirst += iOffset;
	  }
     }

     template <typename ITER>
     index_type SpGradient::InnerProductSize(ITER pFirst, ITER pLast, index_type iOffset) {
	  index_type iNumNz = 0;

	  while (pFirst < pLast) {
	       iNumNz += iGetSize(*pFirst);
	       pFirst += iOffset;
	  }

	  return iNumNz;
     }

     void SpGradient::InnerProductAddDer(const SpGradient& g, doublereal dVal, const SpGradExpDofMap& oDofMap) {
	  g.AddDeriv(*this, dVal, oDofMap);
     }

     template <typename ITER>
     void SpGradient::InnerProductDofStat(ITER pFirst, ITER pLast, index_type iOffset, SpGradDofStat& s) {
	  while (pFirst < pLast) {
	       GetDofStat(*pFirst, s);
	       pFirst += iOffset;
	  }
     }

     SpDerivData* SpGradient::pGetNullData() {
	  SP_GRAD_ASSERT(oNullData.dVal == 0.);
	  SP_GRAD_ASSERT(oNullData.iSizeCurr == 0);
	  SP_GRAD_ASSERT(oNullData.iSizeRes == 0);
	  SP_GRAD_ASSERT(oNullData.iRefCnt >= 1);
	  SP_GRAD_ASSERT(oNullData.bCompressed);
	  SP_GRAD_ASSERT(oNullData.pOwner == nullptr);
	  
	  return &oNullData;
     }
}

#endif
