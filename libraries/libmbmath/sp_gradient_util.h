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

#ifndef __SP_GRADIENT_UTIL_H__INCLUDED__
#define __SP_GRADIENT_UTIL_H__INCLUDED__

#include <iterator>
#include <type_traits>

#include "sp_gradient_base.h"

namespace sp_grad {
     namespace util {

          template <typename Expr>
          void ExprEvalHelper<SpGradCommon::ExprEvalDuplicate>::Eval(SpGradient& g, const SpGradBase<Expr>& f) {
               g.Assign(f);
          }

          template <typename Func, typename Expr>
          void ExprEvalHelper<SpGradCommon::ExprEvalDuplicate>::AssignOper(SpGradient& g1, const SpGradBase<Expr>& g2) {
               g1.AssignOper<Func, Expr>(g2);
          }
	  
          template <typename Expr>
          void ExprEvalHelper<SpGradCommon::ExprEvalUnique>::Eval(SpGradient& g, const SpGradBase<Expr>& f) {
               g.MapAssign(f);
          }

	  template <typename Func, typename Expr>
          void ExprEvalHelper<SpGradCommon::ExprEvalUnique>::AssignOper(SpGradient& g1, const SpGradBase<Expr>& g2) {
               g1.MapAssignOper<Func, Expr>(g2);
          }
	  
          template <typename AVAL, typename BVAL>
          struct InnerProductHelper {
               template <typename AITER, typename BITER>
               static void
               MapEval(SpGradient& g,
                       AITER pAFirst,
                       AITER pALast,
                       index_type iAOffset,
                       BITER pBFirst,
                       BITER pBLast,
                       index_type iBOffset) {
                    g.MapInnerProduct(pAFirst,
                                      pALast,
                                      iAOffset,
                                      pBFirst,
                                      pBLast,
                                      iBOffset);
               }

               template <typename AITER, typename BITER>
               static void
               MapEval(SpGradient& g,
                       AITER pAFirst,
                       AITER pALast,
                       index_type iAOffset,
                       BITER pBFirst,
                       BITER pBLast,
                       index_type iBOffset,
                       const SpGradExpDofMap& oDofMap) {
                    g.MapInnerProduct(pAFirst,
                                      pALast,
                                      iAOffset,
                                      pBFirst,
                                      pBLast,
                                      iBOffset,
                                      oDofMap);
               }

               template <typename AITER, typename BITER>
               static void
               Eval(SpGradient& g,
                    AITER pAFirst,
                    AITER pALast,
                    index_type iAOffset,
                    BITER pBFirst,
                    BITER pBLast,
                    index_type iBOffset) {
                    g.InnerProduct(pAFirst,
                                   pALast,
                                   iAOffset,
                                   pBFirst,
                                   pBLast,
                                   iBOffset);
               }
          };

          template <>
          struct InnerProductHelper<doublereal, doublereal> {
               template <typename AITER, typename BITER>
               static inline void
               Eval(doublereal& d,
                    AITER pAFirst,
                    AITER pALast,
                    index_type iAOffset,
                    BITER pBFirst,
                    BITER pBLast,
                    index_type iBOffset) noexcept {
                    SP_GRAD_ASSERT((pALast - pAFirst) / iAOffset == (pBLast - pBFirst) / iBOffset);
                    SP_GRAD_ASSERT((pALast - pAFirst) % iAOffset == 0);
                    SP_GRAD_ASSERT((pBLast - pBFirst) % iBOffset == 0);

                    d = 0.;

                    while (pAFirst < pALast) {
                         d += (*pAFirst) * (*pBFirst);
                         pAFirst += iAOffset;
                         pBFirst += iBOffset;
                    }
               }

               template <typename AITER, typename BITER>
               static void
               MapEval(doublereal& d, AITER pAFirst, AITER pALast, index_type iAOffset, BITER pBFirst, BITER pBLast, index_type iBOffset) {
                    Eval(d, pAFirst, pALast, iAOffset, pBFirst, pBLast, iBOffset);
               }

               template <typename AITER, typename BITER>
               static void
               MapEval(doublereal& d,
                       AITER pAFirst,
                       AITER pALast,
                       index_type iAOffset,
                       BITER pBFirst,
                       BITER pBLast,
                       index_type iBOffset,
                       const SpGradExpDofMap&) {
                    Eval(d, pAFirst, pALast, iAOffset, pBFirst, pBLast, iBOffset);
               }
          };

          template <typename AITER, typename BITER>
          struct IterResultType {
               typedef typename ResultType<typename std::iterator_traits<AITER>::value_type,
                                           typename std::iterator_traits<BITER>::value_type>::Type Type;
          };

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(typename IterResultType<AITER, BITER>::Type& g,
                          AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset) {
               typedef InnerProductHelper<typename std::iterator_traits<AITER>::value_type,
                                          typename std::iterator_traits<BITER>::value_type> IPH;
               IPH::MapEval(g,
                            pAFirst,
                            pALast,
                            iAOffset,
                            pBFirst,
                            pBLast,
                            iBOffset);
          }

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(typename IterResultType<AITER, BITER>::Type& g,
                          AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset,
                          const SpGradExpDofMap& oDofMap) {
               typedef InnerProductHelper<typename std::iterator_traits<AITER>::value_type,
                                          typename std::iterator_traits<BITER>::value_type> IPH;
               IPH::MapEval(g,
                            pAFirst,
                            pALast,
                            iAOffset,
                            pBFirst,
                            pBLast,
			    iBOffset,
			    oDofMap);
	  }

	  template <typename AITER, typename BITER>
	  inline void
	  InnerProduct(typename IterResultType<AITER, BITER>::Type& g,
		       AITER pAFirst,
		       AITER pALast,
		       index_type iAOffset,
		       BITER pBFirst,
		       BITER pBLast,
		       index_type iBOffset) {
	       typedef InnerProductHelper<typename std::iterator_traits<AITER>::value_type,
					  typename std::iterator_traits<BITER>::value_type> IPH;
	       IPH::Eval(g,
			 pAFirst,
			 pALast,
			 iAOffset,
			 pBFirst,
			 pBLast,
			 iBOffset);
	  }
     }
}
#endif
