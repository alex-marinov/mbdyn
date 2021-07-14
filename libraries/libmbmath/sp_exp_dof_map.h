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

#ifndef __SP_EXP_DOF_MAP_H__INCLUDED__
#define __SP_EXP_DOF_MAP_H__INCLUDED__

#include <vector>

#include "sp_gradient_base.h"

namespace sp_grad {
     class SpGradExpDofMap {
     public:
	  inline SpGradExpDofMap();
	  inline SpGradExpDofMap(const SpGradExpDofMap&)=delete;
	  inline explicit SpGradExpDofMap(const SpGradDofStat& s);
	  inline ~SpGradExpDofMap();
	  
	  inline SpGradExpDofMap& operator=(const SpGradExpDofMap)=delete;
	  
	  inline bool InsertDof(index_type iDof);

	  inline void InsertDone();

	  inline index_type iGetLocalSize() const;

	  inline index_type iGetLocalIndex(index_type iDof) const;

	  inline index_type iGetGlobalIndex(index_type iSlot) const;

	  inline void Reset(const SpGradDofStat& s);
	  
#ifdef SP_GRAD_DEBUG
	  inline bool bValid() const;
#endif
     private:
	  index_type iMinDof, iMaxDof, iNumNzCurr, iNumNzMax;
	  index_type* pIdx;
	  index_type* pPerm;
	  std::vector<index_type> rgData;
     };

     SpGradExpDofMap::SpGradExpDofMap()
	  :iMinDof(std::numeric_limits<index_type>::max()),
	   iMaxDof(std::numeric_limits<index_type>::min()),
	   iNumNzCurr(0),
	   iNumNzMax(0),
	   pIdx(nullptr),
	   pPerm(nullptr)
     {
	  SP_GRAD_ASSERT(bValid());
     }

     SpGradExpDofMap::SpGradExpDofMap(const SpGradDofStat& s)
	  :SpGradExpDofMap()
     {
	  Reset(s);
     }

     SpGradExpDofMap::~SpGradExpDofMap()
     {
	  SP_GRAD_ASSERT(bValid());
     }

     bool SpGradExpDofMap::InsertDof(index_type iDof)
     {
	  SP_GRAD_ASSERT(iDof != SpGradCommon::iInvalidDof);
	  SP_GRAD_ASSERT(iDof != SpGradCommon::iDeletedDof);
	  SP_GRAD_ASSERT(iDof >= iMinDof);
	  SP_GRAD_ASSERT(iDof <= iMaxDof);

	  SP_GRAD_ASSERT(static_cast<size_t>(iNumNzMax + iMaxDof - iMinDof + 1) == rgData.size());
	  SP_GRAD_ASSERT(pIdx + iDof >= &rgData.front() + iNumNzMax);
	  SP_GRAD_ASSERT(pIdx + iDof <= &rgData.back());

	  index_type& p = pIdx[iDof];

	  if (p != SpGradCommon::iInvalidDof) {
	       SP_GRAD_ASSERT(p < iNumNzCurr);
	       SP_GRAD_ASSERT(pPerm[p] == iDof);
	       
	       return false;
	  } else {
	       p = iNumNzCurr;

	       SP_GRAD_ASSERT(iNumNzCurr < iNumNzMax);
	       SP_GRAD_ASSERT(pPerm[iNumNzCurr] == SpGradCommon::iInvalidDof);

	       pPerm[iNumNzCurr++] = iDof;

	       SP_GRAD_ASSERT(iNumNzCurr <= iNumNzMax);
	       
	       return true;
	  }
     }

     void SpGradExpDofMap::InsertDone()
     {
	  SP_GRAD_ASSERT(bValid());
     }

     index_type SpGradExpDofMap::iGetLocalSize() const {
	  SP_GRAD_ASSERT(iNumNzCurr <= iNumNzMax);
	  
	  return iNumNzCurr;
     }

     index_type SpGradExpDofMap::iGetGlobalIndex(index_type iSlot) const {
	  SP_GRAD_ASSERT(iSlot >= 0);
	  SP_GRAD_ASSERT(iSlot < iGetLocalSize());
	  SP_GRAD_ASSERT(rgData.size() == static_cast<size_t>(iNumNzMax + iMaxDof - iMinDof + 1));
	  SP_GRAD_ASSERT(pPerm[iSlot] >= iMinDof);
	  SP_GRAD_ASSERT(pPerm[iSlot] <= iMaxDof);
	  SP_GRAD_ASSERT(pIdx[pPerm[iSlot]] == iSlot);
	  
	  return pPerm[iSlot];
     }

     index_type SpGradExpDofMap::iGetLocalIndex(index_type iDof) const {
	  SP_GRAD_ASSERT(iDof >= iMinDof);
	  SP_GRAD_ASSERT(iDof <= iMaxDof);
	  SP_GRAD_ASSERT(&rgData.front() + iNumNzMax == pIdx + iMinDof);
	  SP_GRAD_ASSERT(pIdx + iDof >= &rgData.front() + iNumNzMax);
	  SP_GRAD_ASSERT(pIdx + iDof <= &rgData.back());
	  SP_GRAD_ASSERT(pIdx[iDof] >= 0);
	  SP_GRAD_ASSERT(pIdx[iDof] < iGetLocalSize());
	  SP_GRAD_ASSERT(pPerm[pIdx[iDof]] == iDof);
	  return pIdx[iDof];
     }

     void SpGradExpDofMap::Reset(const SpGradDofStat& s) {
	  SP_GRAD_ASSERT(bValid());
	  
	  iMinDof = s.iMinDof;
	  iMaxDof = s.iMaxDof;
	  iNumNzCurr = 0;
	  iNumNzMax = s.iNumNz;
	  rgData.clear();
	  rgData.resize(iNumNzMax + iMaxDof - iMinDof + 1, SpGradCommon::iInvalidDof);
	  pPerm = &rgData.front();
	  pIdx = &rgData.front() + iNumNzMax - iMinDof;

	  SP_GRAD_ASSERT(bValid());
     }

#ifdef SP_GRAD_DEBUG
     bool SpGradExpDofMap::bValid() const {
	  SP_GRAD_ASSERT(iNumNzCurr >= 0);
	  SP_GRAD_ASSERT(iNumNzMax >= 0);
	  SP_GRAD_ASSERT(iNumNzCurr <= iNumNzMax);
	  
	  if (iNumNzCurr > 0) {
	       SP_GRAD_ASSERT(iMinDof <= iMaxDof);
	       SP_GRAD_ASSERT(iMinDof >= 1);
	  }
	  
	  SP_GRAD_ASSERT(iNumNzMax == 0 || rgData.size() == static_cast<size_t>(iMaxDof - iMinDof + 1 + iNumNzMax));

	  std::vector<bool> v(iMaxDof - iMinDof + 1, false);
	  
	  for (index_type i = 0; i < iNumNzCurr; ++i) {
	       SP_GRAD_ASSERT(iGetLocalIndex(iGetGlobalIndex(i)) == i);

	       index_type idx = pPerm[i] - iMinDof;

	       SP_GRAD_ASSERT(!v[idx]);
	       
	       v[idx] = true;
	  }

	  for (index_type i = iNumNzCurr; i < iNumNzMax; ++i) {
	       SP_GRAD_ASSERT(pPerm[i] == SpGradCommon::iInvalidDof);
	  }

	  for (index_type i = iMinDof; i <= iMaxDof; ++i) {
	       if (pIdx[i] != SpGradCommon::iInvalidDof) {
		    index_type p = iGetLocalIndex(i);	       
		    SP_GRAD_ASSERT(iGetGlobalIndex(p) == i);
	       }
	  }
	  
	  return true;
     }
#endif
}
#endif
