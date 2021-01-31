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

#include <cstdlib>
#include <algorithm>

#ifdef SP_GRAD_DEBUG
#include <cmath>
#include <vector>
#endif

#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_gradient_op.h"

namespace sp_grad {
     SP_GRAD_THREAD_LOCAL SpDerivData SpGradient::oNullData{0., 0, 0, SpDerivData::DER_UNIQUE | SpDerivData::DER_SORTED, 1, nullptr};

     SpDerivData* SpGradient::pAllocMem(SpDerivData* ptr, index_type iSize) {
	  ptr = reinterpret_cast<SpDerivData*>(std::realloc(ptr, uGetAllocSize(iSize)));

	  if (!ptr) {
	       throw std::bad_alloc();
	  }

	  return ptr;
     }

     void SpGradient::Allocate(index_type iSizeRes, index_type iSizeInit, unsigned uFlags) {
	  SP_GRAD_ASSERT(iSizeRes >= 0);
	  SP_GRAD_ASSERT(iSizeInit >= 0);
	  SP_GRAD_ASSERT(iSizeRes >= iSizeInit);
	  SP_GRAD_ASSERT(pData->iRefCnt >= 1);
	  
	  const doublereal dVal = pData->dVal;

	  SpDerivData* pMem;

	  const bool bNeedToGrow = iSizeRes > pData->iSizeRes;
	  const bool bCannotReuse = pData->iRefCnt > 1 || (pData->pOwner && (bNeedToGrow || !pData->pOwner->bIsOwnerOf(this)));	  

	  if (bCannotReuse) {
	       pMem = pAllocMem(nullptr, iSizeRes);
	       index_type iSizeCopy = std::min(pData->iSizeCurr, iSizeInit);
	       std::uninitialized_copy(pData->rgDer,
				       pData->rgDer + iSizeCopy,
				       pMem->rgDer);
	       Cleanup();
	  } else if (bNeedToGrow) {
	       SP_GRAD_ASSERT(!pData->pOwner);
	       SP_GRAD_ASSERT(pData != pGetNullData());
	       
	       pMem = pAllocMem(pData, iSizeRes);
	  } else {
	       SP_GRAD_ASSERT(((pData->pOwner && pData->pOwner->bIsOwnerOf(this))) ? (pData->iSizeRes == 0 || pData->iSizeRes >= iSizeRes) : true);
	       SP_GRAD_ASSERT(pData != pGetNullData());
	       pData->iSizeCurr = iSizeInit;
	       pData->uFlags = uFlags;
	       return;
	  }
	  
	  pData = new(pMem) SpDerivData(dVal, iSizeRes, iSizeInit, uFlags, 1, nullptr);
     }

     void SpGradient::Sort() {
	  if (!(pData->uFlags & SpDerivData::DER_UNIQUE)) {
	       MakeUnique();
	  }

	  if (!(pData->uFlags & SpDerivData::DER_SORTED)) {
	       std::sort(pData->rgDer, pData->rgDer + pData->iSizeCurr);

	       pData->uFlags |= SpDerivData::DER_SORTED;
	  }
     }

     void SpGradient::MakeUnique() {
	  *this = EvalUnique(*this);
     }
     
#ifdef SP_GRAD_DEBUG
     bool SpGradient::bValid() const {
	  SP_GRAD_ASSERT(pData != nullptr);

	  if (pData == pGetNullData()) {
	       SP_GRAD_ASSERT(pData->dVal == 0.);
	       SP_GRAD_ASSERT(pData->iSizeCurr == 0);
	       SP_GRAD_ASSERT(pData->iSizeRes == 0);
	       SP_GRAD_ASSERT(pData->iRefCnt >= 1);
	       SP_GRAD_ASSERT(pData->uFlags & SpDerivData::DER_SORTED);
	       SP_GRAD_ASSERT(pData->uFlags & SpDerivData::DER_UNIQUE);
	  }

	  SP_GRAD_ASSERT(std::isfinite(pData->dVal));
	  SP_GRAD_ASSERT(pData->iSizeRes >= 0);
	  SP_GRAD_ASSERT(pData->iSizeCurr <= pData->iSizeRes);
	  SP_GRAD_ASSERT(pData->iRefCnt > 0);
	  SP_GRAD_ASSERT(!bIsUnique() || bCheckUnique());
	  
	  for (index_type i = 0; i < pData->iSizeCurr; ++i) {
	       SP_GRAD_ASSERT(std::isfinite(pData->rgDer[i].dDer));
	       SP_GRAD_ASSERT(pData->rgDer[i].iDof > 0);
	       
	       if (i > 0 && bIsSorted()) {
	       	    SP_GRAD_ASSERT(pData->rgDer[i].iDof > pData->rgDer[i - 1].iDof);
	       }
	  }

	  return true;
     }

     bool SpGradient::bCheckUnique() const {
	  SpGradDofStat s;
	  
	  GetDofStat(s);
	  
	  std::vector<bool> v(s.iMaxDof - s.iMinDof + 1, false);

	  for (const auto&r: *this) {
	       size_t i = r.iDof - s.iMinDof;

	       SP_GRAD_ASSERT(i < v.size());
	       
	       if (v[i]) {
		    return false;
	       }
	       
	       v[i] = true;
	  }

	  return true;
     }
     
     void SpGradient::PrintValue(std::ostream& os) const {
	  os << dGetValue() << ' ';
     }

     void SpGradient::PrintDeriv(std::ostream& os, doublereal dCoef) const {
	  for (const auto& d: *this) {
	       os << d.iDof << ':' << d.dDer << " * " << dCoef << ' ';
	  }
     }
#endif
}
