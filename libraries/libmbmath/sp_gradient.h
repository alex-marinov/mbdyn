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

     SpGradient::SpGradient(doublereal d)
	  :SpGradient(pGetNullData()) {

	  ResizeReset(d, 0);
	  
	  SP_GRAD_ASSERT(bValid());
     }

     SpGradient::SpGradient(const SpGradient& g)
	  :SpGradient(g.pData) {

	  SP_GRAD_ASSERT(bValid());
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

	  Cleanup(); // Avoid incrementing oNullData.iRefCnt here, because it could overflow
     }

     SpGradient& SpGradient::operator=(const SpGradient& g) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());

	  if (&g == this) {
	       return *this;
	  }
	  
	  Cleanup();

	  pData = g.pData;
	  ++pData->iRefCnt;

	  if (pData->pOwner) {
	       pData->pOwner->Attach(this);
	  }

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

	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;
	  
	  EvalHelp::template AssignOper<SpGradBinPlus, Expr>(*this, g);

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator-=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;
	  
	  EvalHelp::template AssignOper<SpGradBinMinus, Expr>(*this, g);	 

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator+=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  UniqueOwner();

	  pData->dVal += b;

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator-=(doublereal b) {
	  SP_GRAD_ASSERT(bValid());

	  UniqueOwner();

	  pData->dVal -= b;

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator*=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;
	  
	  EvalHelp::template AssignOper<SpGradBinMult, Expr>(*this, g);	 

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     template <typename Expr>
     SpGradient& SpGradient::operator/=(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());

	  typedef util::ExprEvalHelper<std::remove_reference<Expr>::type::eExprEvalFlags> EvalHelp;
	  
	  EvalHelp::template AssignOper<SpGradBinDiv, Expr>(*this, g);	 	  

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator*=(const doublereal v) {
	  SP_GRAD_ASSERT(bValid());

	  UniqueOwner();

	  pData->dVal *= v;

	  SpDerivRec* pCurr = pData->rgDer;
	  SpDerivRec* pEnd = pCurr + pData->iSizeCurr;

	  while (pCurr < pEnd) {
	       pCurr->dDer *= v;
	       ++pCurr;
	  }

	  SP_GRAD_ASSERT(bValid());

	  return *this;
     }

     SpGradient& SpGradient::operator/=(const doublereal v) {
	  SP_GRAD_ASSERT(bValid());

	  UniqueOwner();

	  pData->dVal /= v;

	  SpDerivRec* pCurr = pData->rgDer;
	  SpDerivRec* pEnd = pCurr + pData->iSizeCurr;

	  while (pCurr < pEnd) {
	       pCurr->dDer /= v;
	       ++pCurr;
	  }

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

	  Allocate(1, 1, SpDerivData::DER_SORTED | SpDerivData::DER_UNIQUE);
	  pData->dVal = dVal;
	  pData->rgDer[0].iDof = iDof;
	  pData->rgDer[0].dDer = dDer;

	  SP_GRAD_ASSERT(bValid());
     }

     // Preserves the nonzero pattern (needed for different sparse solvers)
     void SpGradient::ResetNumeric()
     {
	  SP_GRAD_ASSERT(bValid());

	  UniqueOwner();
	  
	  pData->dVal = 0.;

	  SpDerivRec* pCurr = pData->rgDer;
	  SpDerivRec* pEnd = pCurr + pData->iSizeCurr;
	  
	  while (pCurr != pEnd) {
	       pCurr->dDer = 0.;
	       ++pCurr;
	  }	       
	  
	  SP_GRAD_ASSERT(bValid());	  
     }

     void SpGradient::ResizeReset(doublereal dVal, index_type iSize) {
	  SP_GRAD_ASSERT(bValid());

	  Allocate(iSize, 0, SpDerivData::DER_GENERAL);
	  pData->dVal = dVal;

	  SP_GRAD_ASSERT(bValid());
     }

     void SpGradient::ResizeReset(SpGradient& g, doublereal dVal, index_type iSize) {
	  g.ResizeReset(dVal, iSize);
     }

     void SpGradient::ResizeReset(doublereal& g, doublereal dVal, index_type) {
     	  g = dVal;
     }

     void SpGradient::Scale(doublereal dRowScale, const std::vector<doublereal>& oColScale) {
	  UniqueOwner();

	  SpDerivRec* pCurr = pData->rgDer;
	  SpDerivRec* pEnd = pCurr + pData->iSizeCurr;
	  
	  while (pCurr < pEnd) {
	       SP_GRAD_ASSERT(pCurr->iDof >= 1);
	       SP_GRAD_ASSERT(oColScale.size() >= static_cast<size_t>(pCurr->iDof));
	       
	       pCurr->dDer *= dRowScale * oColScale[pCurr->iDof - 1];
	       ++pCurr;
	  }
     }

     void SpGradient::SetValuePreserve(doublereal dVal) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(pData != pGetNullData());

	  pData->dVal = dVal;

	  SP_GRAD_ASSERT(bValid());
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

	  MaybeSort();

	  auto pRec = pFindRec(iDof);

	  SP_GRAD_ASSERT(bValid());

	  return pRec ? pRec->dDer : 0.;
     }

     void SpGradient::InsertDeriv(SpGradient& g, doublereal dCoef) const {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());
	  
	  for (const auto& r: *this) {
	       g.pInsertRec(r.iDof, r.dDer * dCoef);
	  }

	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(g.bValid());	  
     }

     void SpGradient::InsertDof(SpGradExpDofMap& oExpDofMap) const {
	  SP_GRAD_ASSERT(bValid());

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
	  SP_GRAD_ASSERT(g.bValid());
	  SP_GRAD_ASSERT(oDofMap.bValid());
	  
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

	  SP_GRAD_ASSERT(g.bValid());
	  SP_GRAD_ASSERT(oDofMap.bValid());
     }

     void SpGradient::AddDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) {
	  f.AddDeriv(g, dCoef, oDofMap);
     }

     const SpDerivRec* SpGradient::begin() const {
	  return pData->rgDer;
     }

     const SpDerivRec* SpGradient::end() const {
	  return pData->rgDer + pData->iSizeCurr;
     }

     index_type SpGradient::iGetSize() const {
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
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);
	  SP_GRAD_ASSERT((pALast - pAFirst) % iAOffset == 0);
	  SP_GRAD_ASSERT((pBLast - pBFirst) % iBOffset == 0);

	  SpGradDofStat s;

	  InnerProductDofStat(pAFirst, pALast, iAOffset, s);
	  InnerProductDofStat(pBFirst, pBLast, iBOffset, s);

	  SpGradExpDofMap oDofMap(s);

	  SP_GRAD_ASSERT(oDofMap.bValid());

	  InnerProductInsertDof(pAFirst, pALast, iAOffset, oDofMap);
	  InnerProductInsertDof(pBFirst, pBLast, iBOffset, oDofMap);
	  
	  oDofMap.InsertDone();

	  SP_GRAD_ASSERT(oDofMap.bValid());

	  MapInnerProduct(pAFirst, pALast, iAOffset, pBFirst, pBLast, iBOffset, oDofMap);

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename AITER, typename BITER>
     void SpGradient::MapInnerProduct(AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset, const SpGradExpDofMap& oDofMap) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);
	  SP_GRAD_ASSERT((pALast - pAFirst) % iAOffset == 0);
	  SP_GRAD_ASSERT((pBLast - pBFirst) % iBOffset == 0);

	  Allocate(oDofMap.iGetLocalSize(), 0, SpDerivData::DER_UNIQUE);

	  SetValuePreserve(0.);

	  SP_GRAD_ASSERT(bValid());

	  SP_GRAD_ASSERT(oDofMap.bValid());
	  
	  InitDeriv(oDofMap);

	  SP_GRAD_ASSERT(oDofMap.bValid());

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

	  SP_GRAD_ASSERT(oDofMap.bValid());
	  SP_GRAD_ASSERT(bValid());
     }

     template <typename AITER, typename BITER>
     void SpGradient::InnerProduct(AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT((pBLast - pBFirst) / iBOffset == (pALast - pAFirst) / iAOffset);

	  index_type iNumNz;

	  iNumNz = InnerProductSize(pAFirst, pALast, iAOffset);
	  iNumNz += InnerProductSize(pBFirst, pBLast, iBOffset);

	  Allocate(iNumNz, 0, SpDerivData::DER_GENERAL);

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

     void SpGradient::UniqueOwner() {
	  Allocate(pData->iSizeRes, pData->iSizeCurr, pData->uFlags);
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

     void SpGradient::Sort(doublereal) {
     }

     void SpGradient::Sort(SpGradient& g) {
	  g.Sort();
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

     void SpGradient::Cleanup() {
	  SP_GRAD_ASSERT(bValid());
	  
	  --pData->iRefCnt;

	  SP_GRAD_ASSERT(pData->iRefCnt >= 0);

	  if (pData->pOwner) {
	       pData->pOwner->Detach(this);
	  } else if (!pData->iRefCnt && pData != &oNullData) {
	       std::free(pData);
	  }	  
     }
     
     void SpGradient::Free() {	  
	  Cleanup();

	  pData = pGetNullData();

	  SP_GRAD_ASSERT(pData->iRefCnt >= 0);

	  ++pData->iRefCnt;

	  SP_GRAD_ASSERT(bValid());
     }

     constexpr  bool SpGradient::bRecCompareWithDof(const SpDerivRec& a, index_type b) {
	  return a.iDof < b;
     }

     void SpGradient::MaybeSort() const {
	  if (!bIsSorted()) {
	       const_cast<SpGradient*>(this)->Sort();
	  }

	  SP_GRAD_ASSERT(bIsSorted());
     }

     bool SpGradient::bIsSorted() const {
	  return (pData->uFlags & SpDerivData::DER_SORTED) != 0u;
     }

     bool SpGradient::bIsUnique() const {
	  return (pData->uFlags & SpDerivData::DER_UNIQUE) != 0u;
     }

     SpDerivRec* SpGradient::pFindRec(index_type iDof) const {
	  SP_GRAD_ASSERT(bIsSorted());

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
	  index_type iSize = static_cast<index_type>(rgDer.size());

	  if (iSize > 0) {
	       Allocate(iSize, iSize, SpDerivData::DER_GENERAL);

	       pData->dVal = dVal;

	       std::uninitialized_copy(std::begin(rgDer), std::end(rgDer), pData->rgDer);

	       SP_GRAD_ASSERT(!bIsSorted());
	       SP_GRAD_ASSERT(!bIsUnique());
	  }
     }

     template <typename Expr>
     void SpGradient::Assign(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());
	  
	  SpGradient f;

	  if (!g.bHaveRefTo(*this)) {
	       f = std::move(*this);
	  }

	  f.Allocate(g.iGetSize(), 0, SpDerivData::DER_GENERAL);

	  f.pData->dVal = g.dGetValue();

	  g.InsertDeriv(f, 1.);

	  *this = std::move(f);

	  SP_GRAD_ASSERT(!bIsUnique());
	  SP_GRAD_ASSERT(!bIsSorted());
	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Expr>
     void SpGradient::MapAssign(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());
	  
	  SpGradient f;

	  if (!g.bHaveRefTo(*this)) {
	       f = std::move(*this);
	  }

	  SpGradDofStat s;

	  g.GetDofStat(s);

	  SpGradExpDofMap oDofMap(s);

	  g.InsertDof(oDofMap);

	  oDofMap.InsertDone();

	  f.Allocate(oDofMap.iGetLocalSize(), 0, SpDerivData::DER_UNIQUE);

	  f.pData->dVal = g.dGetValue();

	  f.InitDeriv(oDofMap);

	  g.AddDeriv(f, 1., oDofMap);
	  
	  *this = std::move(f);

	  SP_GRAD_ASSERT(!bIsSorted());
	  SP_GRAD_ASSERT(bIsUnique());
	  SP_GRAD_ASSERT(bValid());
     }
     
     void SpGradient::InitDeriv(const SpGradExpDofMap& oExpDofMap) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(pData != pGetNullData());       

	  pData->uFlags |= SpDerivData::DER_UNIQUE;
	  
	  for (index_type i = 0; i < oExpDofMap.iGetLocalSize(); ++i) {
	       pInsertRec(oExpDofMap.iGetGlobalIndex(i), 0.);
	  }

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Func, typename Expr>
     void SpGradient::AssignOper(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());
	  
	  const doublereal u = dGetValue();
	  const doublereal v = g.dGetValue();
	  const doublereal f = Func::f(u, v);
	  const doublereal df_du = Func::df_du(u, v);
	  const doublereal df_dv = Func::df_dv(u, v);
	  const index_type iSizePrev = iGetSize();

	  SpGradient r;

	  if (!g.bHaveRefTo(*this)) {
	       r = std::move(*this);
	  } else {
	       r = *this;
	  }

	  r.Allocate(iSizePrev + g.iGetSize(), iSizePrev, SpDerivData::DER_GENERAL);

	  r.pData->dVal = f;

	  Func::update_u(df_du, r.pData->rgDer, r.pData->rgDer + r.pData->iSizeCurr);
	  
	  g.InsertDeriv(r, df_dv);

	  r.pData->uFlags = 0u;
	  
	  *this = std::move(r);

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Func, typename Expr>
     void SpGradient::MapAssignOper(const SpGradBase<Expr>& g) {
	  SP_GRAD_ASSERT(bValid());
	  
	  const doublereal u = dGetValue();
	  const doublereal v = g.dGetValue();
	  const doublereal f = Func::f(u, v);
	  const doublereal df_du = Func::df_du(u, v);
	  const doublereal df_dv = Func::df_dv(u, v);
	  
	  SpGradient r;

	  if (!g.bHaveRefTo(*this)) {
	       r = std::move(*this);
	  } else {
	       r = *this;
	  }

	  if (!r.bIsUnique()) {
	       r.MakeUnique();
	  }
	  
	  SpGradDofStat s;

	  r.GetDofStat(s);
	  g.GetDofStat(s);

	  SpGradExpDofMap oDofMap(s);

	  r.InsertDof(oDofMap);
	  g.InsertDof(oDofMap);
	  
	  oDofMap.InsertDone();

	  r.template InitDerivAssign<Func>(f, df_du, oDofMap);	     

	  g.AddDeriv(r, df_dv, oDofMap);

	  *this = std::move(r);

	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Func>
     void SpGradient::InitDerivAssign(const doublereal f, const doublereal df_du, const SpGradExpDofMap& oDofMap) {
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(bIsUnique());

	  const index_type iNumNz = oDofMap.iGetLocalSize();
	  
	  Allocate(iNumNz, iGetSize(), SpDerivData::DER_UNIQUE);
	  
	  pData->dVal = f; // Must execute after Allocate

	  SP_GRAD_ASSERT(bValid());
	  
	  SP_GRAD_ASSERT(iNumNz <= pData->iSizeRes);
	  SP_GRAD_ASSERT(oDofMap.bValid());
	  
	  SpDerivRec* const pBegin = pData->rgDer;
	  SpDerivRec* pCurrEnd = pBegin + pData->iSizeCurr;
	  SpDerivRec* const pEndOfStorage = pBegin + iNumNz;

	  Func::update_u(df_du, pBegin, pCurrEnd);

#ifdef SP_GRAD_DEBUG
	  for (const SpDerivRec* pCurr = pBegin; pCurr < pCurrEnd; ++pCurr) {
	       SP_GRAD_ASSERT(pCurr->iDof == oDofMap.iGetGlobalIndex(pCurr - pBegin));
	  }
#endif

	  SP_GRAD_ASSERT(bValid());
	  
	  while (pCurrEnd < pEndOfStorage) {
	       pCurrEnd->iDof = oDofMap.iGetGlobalIndex(pCurrEnd - pBegin);
	       pCurrEnd->dDer = 0.;
	       ++pCurrEnd;
	  }

#ifdef SP_GRAD_DEBUG
	  for (const SpDerivRec* pCurr = pBegin; pCurr < pEndOfStorage; ++pCurr) {
	       SP_GRAD_ASSERT(pCurr->iDof == oDofMap.iGetGlobalIndex(pCurr - pBegin));
	  }
#endif

	  pData->iSizeCurr = iNumNz;

	  SP_GRAD_ASSERT(bIsUnique());
	  SP_GRAD_ASSERT(oDofMap.bValid());
	  SP_GRAD_ASSERT(bValid());
     }

     template <typename Func>
     void SpGradient::InitDerivAssign(const doublereal f, const doublereal df_du, index_type iSizeRes) {
	  SP_GRAD_ASSERT(bValid());
	  
	  Allocate(iSizeRes, pData->iSizeCurr, SpDerivData::DER_GENERAL);

	  pData->dVal = f; // Must execute after Allocate
	  
	  SpDerivRec* pBegin = pData->rgDer;
	  SpDerivRec* pCurrEnd = pBegin + pData->iSizeCurr;

	  Func::update_u(df_du, pBegin, pCurrEnd);

	  SP_GRAD_ASSERT(bValid());
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
	  SP_GRAD_ASSERT(oDofMap.bValid());
	  
	  while (pFirst < pLast) {
	       InsertDof(*pFirst, oDofMap);
	       pFirst += iOffset;
	  }

	  SP_GRAD_ASSERT(oDofMap.bValid());
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
	  SP_GRAD_ASSERT(bValid());
	  SP_GRAD_ASSERT(oDofMap.bValid());
	  
	  g.AddDeriv(*this, dVal, oDofMap);

	  SP_GRAD_ASSERT(bValid());
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
	  SP_GRAD_ASSERT(oNullData.uFlags & SpDerivData::DER_UNIQUE);
	  SP_GRAD_ASSERT(oNullData.uFlags & SpDerivData::DER_SORTED);
	  SP_GRAD_ASSERT(oNullData.pOwner == nullptr);
	  
	  return &oNullData;
     }
}

#endif
