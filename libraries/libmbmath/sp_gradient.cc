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

#ifdef SP_GRAD_DEBUG
#include <cmath>
#include <vector>

#endif

#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_gradient_op.h"

namespace sp_grad {
     SpDerivData SP_GRAD_THREAD_LOCAL SpGradient::oNullData{0., 0, 0, true, 1, nullptr};

     SpDerivData* SpGradient::pAllocMem(SpDerivData* ptr, index_type iSize) {
	  ptr = reinterpret_cast<SpDerivData*>(std::realloc(ptr, uGetAllocSize(iSize)));

	  if (!ptr) {
	       throw std::bad_alloc();
	  }

	  return ptr;
     }

     void SpGradient::Allocate(index_type iSizeRes, index_type iSizeInit, AllocFlags eAllocFlags) {
	  SP_GRAD_ASSERT(iSizeRes >= 0);
	  SP_GRAD_ASSERT(iSizeInit >= 0);
	  SP_GRAD_ASSERT(iSizeRes >= iSizeInit);

	  const doublereal dVal = pData->dVal;

	  SpDerivData* pMem;

	  const bool bHaveOwner = pData->pOwner && pData->pOwner->bIsOwnerOf(this);

	  SP_GRAD_ASSERT((bHaveOwner && eAllocFlags == ALLOC_RESIZE) ? (pData->iSizeRes == 0 || pData->iSizeRes >= iSizeRes) : true);

	  if (bHaveOwner && pData->iSizeRes >= iSizeRes && eAllocFlags == ALLOC_RESIZE) {
	       pData->iSizeCurr = iSizeInit;
	       pData->bCompressed = false;
	       return;
	  } else if (pData->pOwner || pData->iRefCnt > 1 || pData->iSizeRes < iSizeRes) {
	       pMem = pAllocMem(nullptr, iSizeRes);
	       index_type iSizeCopy = std::min(pData->iSizeCurr, iSizeInit);
	       std::uninitialized_copy(pData->rgDer,
				       pData->rgDer + iSizeCopy,
				       pMem->rgDer);
	       Free();
	  } else {
	       SP_GRAD_ASSERT(!pData->pOwner);

	       pMem = pAllocMem(pData, iSizeRes);
	  }

	  pData = new(pMem) SpDerivData(dVal, iSizeRes, iSizeInit, false, 1, nullptr);
     }

     void SpGradient::Compress() {
	  auto ibeg = pData->rgDer;
	  auto iend = pData->rgDer + pData->iSizeCurr;

	  index_type iMaxDof = std::numeric_limits<index_type>::min();
	  index_type iMinDof = std::numeric_limits<index_type>::max();

	  for (auto i = ibeg; i < iend; ++i) {
	       iMaxDof = std::max(iMaxDof, i->iDof);
	       iMinDof = std::min(iMinDof, i->iDof);
	  }

	  index_type iSizeFlat = iMaxDof - iMinDof + 1;

	  std::vector<doublereal> vtmp(iSizeFlat, 0.);

	  index_type iSizeCompr = 0;

	  for (auto i = ibeg; i < iend; ++i) {
	       auto& v = vtmp[i->iDof - iMinDof];

	       if (!v) {
		    ++iSizeCompr;
	       }

	       v += i->dDer;
	  }

	  Allocate(iSizeCompr, 0, ALLOC_RESIZE);

	  for (index_type i = 0; i < iSizeFlat; ++i) {
	       if (auto d = vtmp[i]) {
		    pInsertRec(i + iMinDof, d);
	       }
	  }

	  pData->bCompressed = true;
     }

#ifdef SP_GRAD_DEBUG
     bool SpGradient::bValid() const {
	  SP_GRAD_ASSERT(pData != nullptr);

	  if (pData == pGetNullData()) {
	       SP_GRAD_ASSERT(pData->dVal == 0.);
	       SP_GRAD_ASSERT(pData->iSizeCurr == 0);
	       SP_GRAD_ASSERT(pData->iSizeRes == 0);
	       SP_GRAD_ASSERT(pData->iRefCnt >= 1);
	       SP_GRAD_ASSERT(pData->bCompressed);
	  }

	  SP_GRAD_ASSERT(std::isfinite(pData->dVal));
	  SP_GRAD_ASSERT(pData->iSizeRes >= 0);
	  SP_GRAD_ASSERT(pData->iSizeCurr <= pData->iSizeRes);
	  SP_GRAD_ASSERT(pData->iRefCnt > 0);

	  for (index_type i = 0; i < pData->iSizeCurr; ++i) {
	       SP_GRAD_ASSERT(std::isfinite(pData->rgDer[i].dDer));
	       SP_GRAD_ASSERT(pData->rgDer[i].iDof > 0);
	       
	       if (i > 0 && pData->bCompressed) {
	       	    SP_GRAD_ASSERT(pData->rgDer[i].iDof > pData->rgDer[i - 1].iDof);
	       }
	  }

	  return true;
     }

     bool SpGradient::bIsUnique() const {
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
