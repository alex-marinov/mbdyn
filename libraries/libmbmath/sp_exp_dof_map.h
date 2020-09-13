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

#define SP_GRAD_USE_FLAT_VECTOR

#if defined(SP_GRAD_USE_DENSE_HASH_MAP) && defined(USE_DENSE_HASH_MAP)
#include <dense_hash_map>
#elif defined(SP_GRAD_USE_SPARSE_HASH_MAP) && defined(USE_DENSE_HASH_MAP)
#include <sparse_hash_map>
#elif defined(SP_GRAD_USE_UNORDERED_MAP)
#include <unordered_map>
#endif

#include <vector>
#include "sp_gradient_base.h"

namespace sp_grad {
     class SpGradExpDofMap {
     public:
	  inline SpGradExpDofMap();
	  inline explicit SpGradExpDofMap(const SpGradDofStat& s);

	  inline void InsertDof(index_type iDof);

	  inline void InsertDone();

	  inline index_type iGetLocalSize() const;

	  inline index_type iGetLocalIndex(index_type iDof) const;

	  inline index_type iGetGlobalIndex(index_type iSlot) const;

	  inline void Reset(const SpGradDofStat& s);
     private:
#ifdef SP_GRAD_USE_DENSE_HASH_MAP
	  typedef google::dense_hash_map<index_type, index_type> MapType;
#elif defined(SP_GRAD_USE_SPARSE_HASH_MAP)
	  typedef google::sparse_hash_map<index_type, index_type> MapType;
#elif defined(SP_GRAD_USE_UNORDERED_MAP)
	  typedef std::unordered_map<index_type, index_type> MapType;
#elif defined(SP_GRAD_USE_FLAT_VECTOR)
	  typedef std::vector<index_type> MapType;
	  index_type iMinDof, iMaxDof;
#endif
	  MapType rgIdx;
	  std::vector<index_type> rgPerm;
     };

     SpGradExpDofMap::SpGradExpDofMap()
#ifdef SP_GRAD_USE_FLAT_VECTOR
	  :iMinDof(std::numeric_limits<index_type>::max()),
	   iMaxDof(std::numeric_limits<index_type>::min())
#endif
     {
#if defined(SP_GRAD_USE_DENSE_HASH_MAP)
	  rgIdx.set_empty_key(SpGradCommon::iInvalidDof);
#endif
#if defined(SP_GRAD_USE_DENSE_HASH_MAP) || defined(SP_GRAD_USE_SPARSE_HASH_MAP)
	  rgIdx.set_deleted_key(SpGradCommon::iDeletedDof);
#endif
     }

     SpGradExpDofMap::SpGradExpDofMap(const SpGradDofStat& s)
	  :SpGradExpDofMap()
     {
	  Reset(s);
     }

     void SpGradExpDofMap::InsertDof(index_type iDof)
     {
	  SP_GRAD_ASSERT(iDof != SpGradCommon::iInvalidDof);
	  SP_GRAD_ASSERT(iDof != SpGradCommon::iDeletedDof);
	  SP_GRAD_ASSERT(iDof >= iMinDof);
	  SP_GRAD_ASSERT(iDof <= iMaxDof);

#ifdef SP_GRAD_USE_FLAT_VECTOR
	  SP_GRAD_ASSERT(static_cast<size_t>(iMaxDof - iMinDof + 1) == rgIdx.size());
	  SP_GRAD_ASSERT(static_cast<size_t>(iDof - iMinDof) < rgIdx.size());

	  index_type& p = rgIdx[iDof - iMinDof];

	  if (p != SpGradCommon::iInvalidDof) {
	       SP_GRAD_ASSERT(static_cast<size_t>(p) < rgPerm.size());
	       SP_GRAD_ASSERT(rgPerm[p] == iDof);
	  } else {
	       p = rgPerm.size();

	       SP_GRAD_ASSERT(rgPerm.capacity() > rgPerm.size());

	       rgPerm.push_back(iDof);
	  }
#else
	  auto p = rgIdx.find(iDof);

	  if (p == rgIdx.end()) {
	       rgIdx.insert(std::make_pair(iDof, rgIdx.size()));
	  }
#endif
     }

     void SpGradExpDofMap::InsertDone()
     {
#if !defined(SP_GRAD_USE_FLAT_VECTOR)
	  rgPerm.resize(iGetLocalSize());

	  for (auto idx: rgIdx) {
	       rgPerm[idx.second] = idx.first;
	  }
#endif
     }

     index_type SpGradExpDofMap::iGetLocalSize() const {
#ifdef SP_GRAD_USE_FLAT_VECTOR
	  return rgPerm.size();
#else
	  return rgIdx.size();
#endif
     }

     index_type SpGradExpDofMap::iGetGlobalIndex(index_type iSlot) const {
	  SP_GRAD_ASSERT(iSlot >= 0);
	  SP_GRAD_ASSERT(iSlot < iGetLocalSize());
#ifdef SP_GRAD_USE_FLAT_VECTOR
	  SP_GRAD_ASSERT(rgPerm[iSlot] >= iMinDof);
	  SP_GRAD_ASSERT(rgPerm[iSlot] <= iMaxDof);
	  SP_GRAD_ASSERT(rgIdx.size() == static_cast<size_t>(iMaxDof - iMinDof + 1));
	  SP_GRAD_ASSERT(rgIdx[rgPerm[iSlot] - iMinDof] == iSlot);
#else
	  SP_GRAD_ASSERT(rgPerm.size() == rgIdx.size());
#endif
	  return rgPerm[iSlot];
     }

     index_type SpGradExpDofMap::iGetLocalIndex(index_type iDof) const
     {
#ifdef SP_GRAD_USE_FLAT_VECTOR
	  SP_GRAD_ASSERT(iDof >= iMinDof);
	  SP_GRAD_ASSERT(iDof <= iMaxDof);
	  SP_GRAD_ASSERT(static_cast<size_t>(iDof - iMinDof) < rgIdx.size());
	  SP_GRAD_ASSERT(rgIdx[iDof - iMinDof] >= 0);
	  SP_GRAD_ASSERT(rgIdx[iDof - iMinDof] < iGetLocalSize());
	  SP_GRAD_ASSERT(rgPerm[rgIdx[iDof - iMinDof]] == iDof);

	  return rgIdx[iDof - iMinDof];
#else
	  auto p = rgIdx.find(iDof);

	  if (p != rgIdx.end()) {
	       return p->second;
	  }

	  SP_GRAD_ASSERT(false);

	  return SpGradCommon::iInvalidDof;
#endif
     }

     void SpGradExpDofMap::Reset(const SpGradDofStat& s) {
#ifdef SP_GRAD_USE_FLAT_VECTOR
	  iMinDof = s.iMinDof;
	  iMaxDof = s.iMaxDof;
	  rgIdx.clear();
	  rgIdx.resize(s.iMaxDof - s.iMinDof + 1, SpGradCommon::iInvalidDof);
	  rgPerm.clear();
	  rgPerm.reserve(s.iNumNz);
#else
	  rgIdx.clear();
#endif
     }
}
#endif
