/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2013(-2013) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___GRADIENT_H_INCLUDED___
#define ___GRADIENT_H_INCLUDED___

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "except.h"
#include "ac/f2c.h"

#ifndef GRADIENT_DEBUG
    #ifdef DEBUG
        #define GRADIENT_DEBUG 2
    #else
	    #define GRADIENT_DEBUG 0
    #endif
#endif

#if GRADIENT_DEBUG == 0 || defined(DEBUG)
	#define GRADIENT_ASSERT(expr) ASSERT(expr)
#elif GRADIENT_DEBUG > 0
	#define GRADIENT_ASSERT(expr) assert(expr)
#endif

#if GRADIENT_DEBUG >= 2
#define GRADIENT_TRACE(expr) \
		static_cast<void>(std::cerr << __FILE__ << ":" \
				  << __LINE__ << ":" \
				  << __FUNCTION__ << ": " << expr)
#else
#define GRADIENT_TRACE(expr) static_cast<void>(0)
#endif

#ifndef GRADIENT_MEMORY_STAT
    #ifdef DEBUG
        #define GRADIENT_MEMORY_STAT 1
    #else
        #define GRADIENT_MEMORY_STAT 0
    #endif
#endif

#if ! HAVE_COPYSIGN
inline doublereal copysign(doublereal x, doublereal y) {
	return y > 0. ? std::abs(x) : -std::abs(x);
}
#endif

namespace grad {

typedef integer index_type;

template <index_type N_SIZE>
class Gradient;

/*
 * This structure checks if the data type
 * at the left hand side of an assignment
 * can hold the number of derivatives
 * provided at the right hand side
*/
template <bool bSizeOK>
struct MaxSizeCheck {
private:
	/*
	 * For an expression like
	 * Gradient<3> a = Gradient<4>(5.);
	 * one will get an compilation error here!
	 */
	enum CheckType {};
};

template <>
struct MaxSizeCheck<true> {
	enum CheckType {};
};

template <typename T>
class GradientAllocator: public std::allocator<T> {
public:
    typedef typename std::allocator<T>::pointer pointer;
    typedef typename std::allocator<T>::size_type size_type;

    typedef typename std::allocator<T>::difference_type difference_type;
    typedef typename std::allocator<T>::const_pointer const_pointer;
    typedef typename std::allocator<T>::reference reference;
    typedef typename std::allocator<T>::const_reference const_reference;
    typedef typename std::allocator<T>::value_type value_type;

    template <typename U>
    struct rebind {
        typedef GradientAllocator<U> other; };

    GradientAllocator() throw() { }

    GradientAllocator(const GradientAllocator& a) throw()
    : std::allocator<T>(a) { }

    template <typename U>
    GradientAllocator(const GradientAllocator<U>& a) throw()
    : std::allocator<T>(a) { }

    pointer
    allocate(size_type n, const void* p=0) {
#if GRADIENT_MEMORY_STAT > 0
    	sMemUsage.Inc(n);
#endif
        return std::allocator<T>::allocate(n, p);
    }

    void
    deallocate(pointer p, size_type n) {
    	std::allocator<T>::deallocate(p, n);
#if GRADIENT_MEMORY_STAT > 0
    	sMemUsage.Dec(n);
#endif
    }

#if GRADIENT_MEMORY_STAT > 0
    static class MemStat {
    public:
    	MemStat()
    	:iCurrMem(0), iMaxMem(0), iNumAlloc(0), iNumDealloc(0) {

    	}

    	~MemStat() {
    		std::cerr << "GradientAllocator<" << typeid(T).name()
    				<< ">\n\tiMaxMem = "
    				<< std::setprecision(3)
    				<< doublereal(iMaxMem) / 1024 << "KB" << std::endl
    				<< "\tiCurrMem=" << doublereal(iCurrMem) / 1024 << "KB" << std::endl
    				<< "\tiNumAlloc=" << iNumAlloc << std::endl
    				<< "\tiNumDealloc=" << iNumDealloc << std::endl;
    	}

    	void Inc(size_t iSize) {
    		++iNumAlloc;
    		iCurrMem += iSize * sizeof(T);
    		if (iCurrMem > iMaxMem) {
    			iMaxMem = iCurrMem;
    		}
    	}

    	void Dec(size_t iSize) {
    		++iNumDealloc;
    		iCurrMem -= iSize * sizeof(T);
    	}

    private:
    	size_t iCurrMem;
    	size_t iMaxMem;
    	size_t iNumAlloc;
    	size_t iNumDealloc;
    } sMemUsage;
#endif
};

#if GRADIENT_MEMORY_STAT > 0
template <typename T>
typename GradientAllocator<T>::MemStat GradientAllocator<T>::sMemUsage;
#endif

template <typename T>
struct RangeVectorTraits;

template <>
struct RangeVectorTraits<doublereal> {
	static doublereal Zero(){
		return 0.;
	}

	static doublereal Invalid(){
		return NAN;
	}
};

template <>
struct RangeVectorTraits<bool> {
	static bool Zero(){
		return false;
	}

	static bool Invalid(){
		return false;
	}
};

template <typename T, index_type N_SIZE>
class RangeVector {    
public:
	static const index_type iMaxSize = N_SIZE;

    RangeVector()
        :iStart(0), iEnd(0) {

#if GRADIENT_DEBUG > 0
        Initialize(RangeVectorTraits<T>::Invalid());
#endif
    }
    
    template <typename T2, index_type N_SIZE2>
    RangeVector(const RangeVector<T2, N_SIZE2>& v)
    	:iStart(v.iGetStartIndex()), iEnd(v.iGetEndIndex())
    {
    	Copy(v);
    }

    template <typename T2, index_type N_SIZE2>
    RangeVector& operator=(const RangeVector<T2, N_SIZE2>& v)
    {
    	iStart = v.iGetStartIndex();
    	iEnd = v.iGetEndIndex();

    	Copy(v);

    	return *this;
    }

    explicit RangeVector(index_type iStart, index_type iEnd, const T& dVal)
        :iStart(iStart), iEnd(iEnd) {

        Initialize(dVal);
    }

    void ResizeReset(index_type iStartNew, index_type iEndNew, const T& dVal) {
        iStart = iStartNew;
        iEnd = iEndNew;

        Initialize(dVal);
    }

    void ResizePreserve(index_type iStartNew, index_type iEndNew) {
    	iStart = iStartNew;
    	iEnd = iEndNew;

        GRADIENT_ASSERT(iEnd <= N_SIZE);
        GRADIENT_ASSERT(iStart >= 0);
        GRADIENT_ASSERT(iStart <= iEnd);

    	for (index_type i = 0; i < iStart; ++i) {
    		rgVec[i] = RangeVectorTraits<T>::Zero();
    	}

    	for (index_type i = iEnd; i < N_SIZE; ++i) {
    		rgVec[i] = RangeVectorTraits<T>::Zero();
    	}
    }

    void Reset() {
    	// Reset the data but preserve iStart and iEnd

    	for (index_type i = iGetStartIndex(); i < iGetEndIndex(); ++i) {
    		SetValue(i, RangeVectorTraits<T>::Zero());
    	}

#if GRADIENT_DEBUG > 0
    	for (index_type i = 0; i < N_SIZE; ++i) {
    		// All elements must be zero
    		GRADIENT_ASSERT(rgVec[i] == RangeVectorTraits<T>::Zero());
    	}
#endif
    }

    void Reserve(index_type iMaxSize) {
    	GRADIENT_ASSERT(iMaxSize < N_SIZE);
    }

    index_type iGetStartIndex() const { return iStart; }
    index_type iGetEndIndex() const { return iEnd; }
    index_type iGetSize() const { return iEnd - iStart; }
    static index_type iGetMaxSize() { return iMaxSize; }
    
    T GetValue(index_type i) const {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < N_SIZE);
        GRADIENT_ASSERT((i >= iStart && i < iEnd) || rgVec[i] == 0.);
        return rgVec[i];
    }
    
    void SetValue(index_type i, const T& d) {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < N_SIZE);
        GRADIENT_ASSERT((i >= iStart && i < iEnd)/* || rgVec[i] == 0.*/);
        rgVec[i] = d;
    }

    static bool bUseDynamicMem() { return false; }
    
private:
    void Initialize(const T& dVal) {
        GRADIENT_ASSERT(iEnd <= N_SIZE);
        GRADIENT_ASSERT(iStart >= 0);
        GRADIENT_ASSERT(iStart <= iEnd);

        for (index_type i = 0; i < iStart; ++i) {
            rgVec[i] = RangeVectorTraits<T>::Zero(); // This allows fast access without frequent index checking
        }
        
        for (index_type i = iStart; i < iEnd; ++i) {
            rgVec[i] = dVal;
        }
        
        for (index_type i = iEnd; i < N_SIZE; ++i) {
            rgVec[i] = RangeVectorTraits<T>::Zero(); // This allows fast access without frequent index checking
        }
    }
    
    template <typename T2, index_type N_SIZE2>
    void Copy(const RangeVector<T2, N_SIZE2>& v) {
    	typedef typename MaxSizeCheck<iMaxSize >= RangeVector<T2, N_SIZE2>::iMaxSize>::CheckType check_iMaxSize;

		GRADIENT_ASSERT(iEnd <= N_SIZE);
		GRADIENT_ASSERT(iStart >= 0);
		GRADIENT_ASSERT(iStart <= iEnd);

		for (index_type i = 0; i < iStart; ++i) {
			rgVec[i] = RangeVectorTraits<T>::Zero(); // This allows fast access without frequent index checking
		}

		for (index_type i = iStart; i < iEnd; ++i) {
			rgVec[i] = v.GetValue(i);
		}

		for (index_type i = iEnd; i < N_SIZE; ++i) {
			rgVec[i] = RangeVectorTraits<T>::Zero(); // This allows fast access without frequent index checking
		}
    }
private:
    T rgVec[N_SIZE];
    index_type iStart, iEnd;
};

template <typename T>
class RangeVector<T, 0> {
	typedef std::vector<T, GradientAllocator<T> > vector_type;
	typedef typename vector_type::size_type size_type;

public:
	static const index_type iMaxSize = ((index_type(-1) < 0)
											? ~index_type(0) & ~(index_type(1) << (sizeof(index_type) * CHAR_BIT - 1))
											: ~index_type(0)
										 ) / sizeof(T);

    RangeVector()
        :iStart(0) {

    }
    
    RangeVector(const RangeVector& v)
    	:oVec(v.oVec), iStart(v.iStart) {

    }

    ~RangeVector() {

    }

    template <typename T2, index_type N_SIZE2>
    RangeVector(const RangeVector<T2, N_SIZE2>& v)
    	:oVec(v.iGetEndIndex() - v.iGetStartIndex()),
    	 iStart(v.iGetStartIndex())
    {
    	Copy(v);
    }

    template <typename T2, index_type N_SIZE2>
    RangeVector& operator=(const RangeVector<T2, N_SIZE2>& v)
    {
    	iStart = v.iGetStartIndex();

    	oVec.resize(v.iGetEndIndex() - iStart);

    	Copy(v);

    	return *this;
    }

    explicit RangeVector(index_type iStart, index_type iEnd, const T& dVal)
        :oVec(iEnd - iStart, dVal), iStart(iStart) {

    }
    
    void ResizeReset(index_type iStartNew, index_type iEnd, const T& dVal) {
        iStart = iStartNew;

        oVec.resize(iEnd - iStart);

        std::fill(oVec.begin(), oVec.end(), dVal);
    }

    void ResizePreserve(index_type iStartNew, index_type iEndNew) {

    	if (iStartNew == iStart) {
    		oVec.resize(iEndNew - iStartNew, RangeVectorTraits<T>::Zero());
    	} else {
    		vector_type oVecNew(iEndNew - iStartNew);

    		for (index_type i = 0; size_t(i) < oVecNew.size(); ++i) {
    			const index_type iOld = i + iStartNew - iStart;

    			if (iOld >= 0 && size_t(iOld) < oVec.size()) {
    				oVecNew[i] = oVec[iOld];
    			}
    		}

    		oVec = oVecNew;
    		iStart = iStartNew;
    	}
    }

    void Reset() {
    	// Reset the data but preserve memory and iStart
    	std::fill(oVec.begin(), oVec.end(), RangeVectorTraits<T>::Zero());
    }

    void Reserve(index_type iMaxSize) {
    	oVec.reserve(iMaxSize);
    }

    index_type iGetStartIndex() const { return iStart; }
    index_type iGetEndIndex() const { return iStart + oVec.size(); }
    index_type iGetSize() const { return oVec.size(); }
    index_type iGetMaxSize() const { return iMaxSize; }

    T GetValue(index_type i) const {
    	GRADIENT_ASSERT(i >= 0);
    	GRADIENT_ASSERT(i < iGetMaxSize());
        i -= iStart;
        return (i >= 0 && size_type(i) < oVec.size())
        		? oVec[i]
        		: RangeVectorTraits<T>::Zero();
    }
    
    void SetValue(index_type i, const T& d) {
    	GRADIENT_ASSERT(i >= 0);
    	GRADIENT_ASSERT(i < iGetMaxSize());

        i -= iStart;

        GRADIENT_ASSERT(i >= 0 && size_type(i) < oVec.size());
        oVec[i] = d;
    }
    
    static bool bUseDynamicMem() { return true; }
    
private:
    template <typename T2, index_type N_SIZE2>
    void Copy(const RangeVector<T2, N_SIZE2>& v)
    {
    	for (index_type i = iGetStartIndex(); i < iGetEndIndex(); ++i) {
    		SetValue(i,  v.GetValue(i));
    	}
    }

private:
    vector_type oVec;
    index_type iStart;
};

enum FunctionCall {
	// FIXME: There should be a flag for the initial derivatives phase
	// 		  However this information is not available for elements at the moment
	//		  The prototype for Element::AssRes and Element::AssJac should be changed like
	//		  AssRes(..., enum FunctionCall func);
	//		  AssJac(..., enum FunctionCall func);
	STATE_MASK 			= 0x0F,
	FUNCTION_MASK 		= 0xF0,
	INITIAL_ASS_FLAG 	= 0x01,
	INITIAL_DER_FLAG 	= 0x02,
	REGULAR_FLAG 		= 0x04,
	RESIDUAL_FLAG 		= 0x10,
	JACOBIAN_FLAG 		= 0x20,
	UNKNOWN_FUNC 	= 0x0,
	INITIAL_ASS_RES = INITIAL_ASS_FLAG | RESIDUAL_FLAG,
	INITIAL_ASS_JAC = INITIAL_ASS_FLAG | JACOBIAN_FLAG,
	INITIAL_DER_RES = INITIAL_DER_FLAG | RESIDUAL_FLAG,
	INITIAL_DER_JAC = INITIAL_DER_FLAG | JACOBIAN_FLAG,
	REGULAR_RES 	= REGULAR_FLAG 	   | RESIDUAL_FLAG,
	REGULAR_JAC 	= REGULAR_FLAG	   | JACOBIAN_FLAG
};

class LocalDofMap {
public:
    typedef std::vector<index_type, GradientAllocator<index_type> > VectorType;
    typedef VectorType::size_type size_type;
    typedef VectorType::const_iterator LocalIterator;
    typedef std::map<index_type,
    				 index_type,
    				 std::less<index_type>,
    				 GradientAllocator<std::pair<index_type, index_type> > >
    				 MapType;
    typedef MapType::const_iterator GlobalIterator;
    static const index_type INVALID_INDEX = -1;

    explicit LocalDofMap(index_type iMaxSize=0)
    	:eLastCall(UNKNOWN_FUNC) {

        if (iMaxSize > 0) {
            oLocalToGlobal.reserve(iMaxSize);
        }
    }

    index_type iGetGlobalDof(index_type iLocal) const {
    	GRADIENT_ASSERT(iLocal >= 0);
    	GRADIENT_ASSERT(size_type(iLocal) < oLocalToGlobal.size());

        index_type iGlobal = oLocalToGlobal[iLocal];

        GRADIENT_ASSERT(oGlobalToLocal.find(iGlobal)->second == iLocal);
        GRADIENT_ASSERT(oLocalToGlobal[iLocal] == iGlobal);
        GRADIENT_ASSERT(oLocalToGlobal.size() == oGlobalToLocal.size());

        return iGlobal;
    }

    index_type iGetLocalIndex(index_type iGlobal) const {
        MapType::const_iterator i = oGlobalToLocal.find(iGlobal);

        if (i == oGlobalToLocal.end()) {
            return INVALID_INDEX;
        }

        index_type iLocal = i->second;

        GRADIENT_ASSERT(iLocal >= 0);
        GRADIENT_ASSERT(size_type(iLocal) < oLocalToGlobal.size());
        GRADIENT_ASSERT(oLocalToGlobal[iLocal] == iGlobal);
        GRADIENT_ASSERT(oLocalToGlobal.size() == oGlobalToLocal.size());

        return iLocal;
    }

    index_type AllocateLocalDof(index_type iGlobal) {
        MapType::const_iterator i = oGlobalToLocal.find(iGlobal);

        if (i != oGlobalToLocal.end()) {
            GRADIENT_ASSERT(i->first == iGlobal);
            return i->second;
        }

        index_type iLocal = oLocalToGlobal.size();
        oLocalToGlobal.push_back(iGlobal);
        oGlobalToLocal[iGlobal] = iLocal;

        GRADIENT_ASSERT(oGlobalToLocal[iGlobal] == iLocal);
        GRADIENT_ASSERT(oLocalToGlobal[iLocal] == iGlobal);
        GRADIENT_ASSERT(oLocalToGlobal.size() == oGlobalToLocal.size());

        return iLocal;
    }

    void Reset(enum FunctionCall func = UNKNOWN_FUNC) {
    	bool bReset = (func & STATE_MASK) != (eLastCall & STATE_MASK)
    					 || func == UNKNOWN_FUNC;

    	eLastCall = func;

    	if (bReset) {
			oLocalToGlobal.clear();
			oGlobalToLocal.clear();
    	}
    }

    enum FunctionCall GetLastCall() const { return eLastCall; }

    LocalIterator BeginLocal() const { return oLocalToGlobal.begin(); }
    LocalIterator EndLocal() const { return oLocalToGlobal.end(); }
    GlobalIterator BeginGlobal() const { return oGlobalToLocal.begin(); }
    GlobalIterator EndGlobal() const { return oGlobalToLocal.end(); }
    index_type Size() const {
    	GRADIENT_ASSERT(oGlobalToLocal.size() == oLocalToGlobal.size());
    	return oLocalToGlobal.size();
    }

private:
    VectorType oLocalToGlobal;
    MapType oGlobalToLocal;
    enum FunctionCall eLastCall;
};

inline std::ostream& operator<<(std::ostream& os, const LocalDofMap& dof) {
    for (LocalDofMap::LocalIterator i = dof.BeginLocal(); i != dof.EndLocal(); ++i) {
        os << i - dof.BeginLocal() << "->" << *i << std::endl;
    }
    
    return os;
}

class MapVectorBase {
public:
    enum LocalScope { LOCAL };
    enum GlobalScope { GLOBAL };
};

template <index_type N_SIZE>
class MapVector: public MapVectorBase {
public:
	static const index_type iMaxSize = RangeVector<doublereal, N_SIZE>::iMaxSize;

	template <index_type N_SIZE2>
	MapVector(const MapVector<N_SIZE2>& v)
		:pDofMap(v.pGetDofMap()), oRange(v.GetLocalVector()) {

	}

    MapVector(LocalDofMap* pMap=0, index_type iStartLocal=0, index_type iEndLocal=0, LocalScope= LOCAL, doublereal dVal=0.)
        :pDofMap(pMap), oRange(iStartLocal, iEndLocal, dVal) {
    }

    MapVector(LocalDofMap* pMap, index_type iLocal, LocalScope, doublereal dVal)
        :pDofMap(pMap), oRange(iLocal, iLocal + 1, dVal) {
    }    

    MapVector(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, GlobalScope s, doublereal dVal) {
        ResizeReset(pMap, iStartGlobal, iEndGlobal, s, dVal);
    }

    MapVector(LocalDofMap* pMap, index_type iGlobal, GlobalScope s, doublereal dVal) {
        ResizeReset(pMap, iGlobal, iGlobal + 1, s, dVal);
    }  
    
    void ResizeReset(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, LocalScope, doublereal dVal) {
        pDofMap = pMap;
        oRange.ResizeReset(iStartLocal, iEndLocal, dVal);
    }

    void ResizePreserve(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, LocalScope) {
        pDofMap = pMap;
        oRange.ResizePreserve(iStartLocal, iEndLocal);
    }

    void ResizeReset(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, GlobalScope, doublereal dVal) {
        pDofMap = pMap;
        
        GRADIENT_ASSERT(pDofMap != 0);
        GRADIENT_ASSERT(iEndGlobal >= iStartGlobal);

        index_type iFirstIndex = std::numeric_limits<index_type>::max(), iLastIndex = 0;

        for (index_type iGlobal = iStartGlobal; iGlobal < iEndGlobal; ++iGlobal) {
            index_type iLocal = pDofMap->AllocateLocalDof(iGlobal);

            if (iLocal < iFirstIndex) {
                iFirstIndex = iLocal;
            }

            if (iLocal > iLastIndex) {
                iLastIndex = iLocal;
            }
        }

        GRADIENT_ASSERT(iLastIndex >= iFirstIndex);

        oRange.ResizeReset(iFirstIndex, iLastIndex + 1, dVal);
    }
    
    void Reset() {
    	oRange.Reset();
    }

    void Reserve(index_type iSize) {
    	oRange.Reserve(iSize);
    }

    doublereal dGetGlobalVector(index_type iGlobalDof) const {
        GRADIENT_ASSERT(pDofMap != 0);
        const index_type iLocalDof = pDofMap->iGetLocalIndex(iGlobalDof);

        if (iLocalDof == LocalDofMap::INVALID_INDEX) {
        	GRADIENT_ASSERT(0);
        	throw ErrOutOfRange(MBDYN_EXCEPT_ARGS);
        }

        GRADIENT_ASSERT(iLocalDof >= 0);
        GRADIENT_ASSERT(iLocalDof < iGetMaxSize());
        return oRange.GetValue(iLocalDof);
    }
    
    void SetGlobalVector(index_type iGlobalDof, doublereal dValue) {
        GRADIENT_ASSERT(pDofMap != 0);
        const index_type iLocalDof = pDofMap->iGetLocalIndex(iGlobalDof);

        if (iLocalDof == LocalDofMap::INVALID_INDEX) {
        	GRADIENT_ASSERT(0);
        	throw ErrOutOfRange(MBDYN_EXCEPT_ARGS);
        }

        GRADIENT_ASSERT(iLocalDof >= iGetStartIndexLocal());
        GRADIENT_ASSERT(iLocalDof < iGetEndIndexLocal());
        oRange.SetValue(iLocalDof, dValue);
    }    

    index_type iGetGlobalDof(index_type iLocalDof) const {
        GRADIENT_ASSERT(pDofMap != 0);
        GRADIENT_ASSERT(iLocalDof >= iGetStartIndexLocal());
        GRADIENT_ASSERT(iLocalDof < iGetEndIndexLocal());
        return pDofMap->iGetGlobalDof(iLocalDof);
    }
    
    index_type iGetSize() const { return oRange.iGetSize(); }
    index_type iGetMaxSize() const { return oRange.iGetMaxSize(); }
    doublereal dGetLocalVector(index_type i) const { return oRange.GetValue(i); }
    const RangeVector<doublereal, N_SIZE>& GetLocalVector() const { return oRange; }
    void SetLocalVector(index_type i, doublereal dValue){ oRange.SetValue(i, dValue); }
    index_type iGetStartIndexLocal() const { return oRange.iGetStartIndex(); }
    index_type iGetEndIndexLocal() const { return oRange.iGetEndIndex(); }
    LocalDofMap* pGetDofMap() const { return pDofMap; }
    bool bUseDynamicMem() const { return oRange.bUseDynamicMem(); }
    static const MapVector Zero;
    
private:
    LocalDofMap* pDofMap;
    RangeVector<doublereal, N_SIZE> oRange;
};

template <typename Expression>
class GradientExpression: public Expression {
public:
    GradientExpression(const Expression& u)
        :Expression(u) {

    }

    template <typename Expr>
    GradientExpression(const Expr& u)
        :Expression(u) {
            
    }
    
    template <typename LhsExpr, typename RhsExpr>
    GradientExpression(const LhsExpr& u, const RhsExpr& v)
        :Expression(u, v) {
            
    }
};


template <typename LhsExpr, typename RhsExpr>
struct MaxDerivatives {
	static const index_type iMaxDerivatives = LhsExpr::iMaxDerivatives > RhsExpr::iMaxDerivatives
			? LhsExpr::iMaxDerivatives : RhsExpr::iMaxDerivatives;
};

template <typename BinFunc, typename LhsExpr, typename RhsExpr>
class BinaryExpr {
public:
	static const index_type iMaxDerivatives = MaxDerivatives<LhsExpr, RhsExpr>::iMaxDerivatives;

	typedef LhsExpr LhsExprType;
	typedef RhsExpr RhsExprType;

    BinaryExpr(const LhsExpr& u, const RhsExpr& v)
        :oU(u), oV(v), bComputed(false) {
    }
    
    doublereal dGetValue() const {
        Compute();

        return f;
    }
    
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
        Compute();

        const doublereal du_dX = oU.dGetDerivativeLocal(iLocalDof);
        const doublereal dv_dX = oV.dGetDerivativeLocal(iLocalDof);

        return df_du * du_dX + df_dv * dv_dX;
    }
    
    index_type iGetStartIndexLocal() const {
        return std::min(oU.iGetStartIndexLocal(), oV.iGetStartIndexLocal());
    }
    
    index_type iGetEndIndexLocal() const {
        return std::max(oU.iGetEndIndexLocal(), oV.iGetEndIndexLocal());
    }
    
    LocalDofMap* pGetDofMap() const {
        LocalDofMap* pDofMap = oU.pGetDofMap();
        
        if (pDofMap == 0) {
            pDofMap = oV.pGetDofMap();
        } else {
            GRADIENT_ASSERT(oV.pGetDofMap() == 0 || pDofMap == oV.pGetDofMap());
        }
        
        return pDofMap;
    }
    
    bool bHaveReferenceTo(const void* p) const {
        return oU.bHaveReferenceTo(p) || oV.bHaveReferenceTo(p);
    }

    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

private:
    void Compute() const {
        if (bComputed) {
            return;
        }
        
        const doublereal u = oU.dGetValue();
        const doublereal v = oV.dGetValue();
        
        f = BinFunc::f(u, v);
        df_du = BinFunc::df_du(u, v);
        df_dv = BinFunc::df_dv(u, v);
        
        bComputed = true;
    }
    
private:
    const LhsExpr oU;
    const RhsExpr oV;
    mutable doublereal f, df_du, df_dv;
    mutable bool bComputed;
};

template <typename UnFunc, typename Expr>
class UnaryExpr {
public:
	static const index_type iMaxDerivatives = Expr::iMaxDerivatives;

    UnaryExpr(const Expr& u)
        :oU(u), bComputed(false) {
            
    }
        
    doublereal dGetValue() const {
        Compute();
        
        return f;
    }
        
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
        Compute();
        
        const doublereal du_dX = oU.dGetDerivativeLocal(iLocalDof);
        
        return df_du * du_dX;
    }
        
    index_type iGetStartIndexLocal() const {
        return oU.iGetStartIndexLocal();
    }
    
    index_type iGetEndIndexLocal() const {
        return oU.iGetEndIndexLocal();
    }
    
    LocalDofMap* pGetDofMap() const {
        return oU.pGetDofMap();
    }
    
    bool bHaveReferenceTo(const void* p) const {
        return oU.bHaveReferenceTo(p);
    }
    
    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

private:
    void Compute() const {
        if (bComputed) {
            return;
        }
        
        const doublereal u = oU.dGetValue();
        
        f = UnFunc::f(u);
        df_du = UnFunc::df_du(u);
        
        bComputed = true;
    }
    
private:
    const Expr oU;
    mutable doublereal f, df_du;
    mutable bool bComputed;
};

template <typename T>
class DirectExpr {
public:
	static const index_type iMaxDerivatives = T::iMaxDerivatives;

    DirectExpr(const T& g)
        :oG(g) {
            
    }
    
    doublereal dGetValue() const {
        return oG.dGetValue();
    }
        
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
        return oG.dGetDerivativeLocal(iLocalDof);
    }
        
    index_type iGetStartIndexLocal() const {
        return oG.iGetStartIndexLocal();
    }

    index_type iGetEndIndexLocal() const {
        return oG.iGetEndIndexLocal();
    }

    LocalDofMap* pGetDofMap() const {
        return oG.pGetDofMap();
    }
    
    bool bHaveReferenceTo(const void* p) const {
        return p == &oG;
    }
    
    index_type iGetMaxDerivatives() const {
    	return oG.iGetMaxDerivatives();
    }

private:
    const T& oG;
};

class ConstExpr {
public:
	static const index_type iMaxDerivatives = 0;

    ConstExpr(doublereal a)
        :dConst(a) {
            
    }
        
    doublereal dGetValue() const {
        return dConst;
    }
        
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
        return 0.;
    }
        
    index_type iGetStartIndexLocal() const {
        // Note: we must not combine two ConstExpr objects in a BinaryExpr!
        // If that happens, we would run out of memory!
        // However that would be meaningless at all.
        return std::numeric_limits<index_type>::max();
    }
    
    index_type iGetEndIndexLocal() const {
        return 0;
    }
    
    LocalDofMap* pGetDofMap() const {
        return 0;
    }
    
    bool bHaveReferenceTo(const void* p) const {
        return false;
    }
    
    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

private:
    const doublereal dConst;
};

template <typename BoolFunc, typename LhsExpr, typename RhsExpr>
class BoolExpr {
public:
	static const index_type iMaxDerivatives = MaxDerivatives<LhsExpr, RhsExpr>::iMaxDerivatives;

    BoolExpr(const LhsExpr& u, const RhsExpr& v)
        :oU(u), oV(v), bComputed(false) {
    }
    
    bool dGetValue() const {
        Compute();
        
        return f;
    }
    
    operator bool() const {
        return dGetValue();
    }
    
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
        return 0.;
    }
    
    index_type iGetStartIndexLocal() const {
        return std::numeric_limits<index_type>::max();
    }
    
    index_type iGetEndIndexLocal() const {
        return 0;
    }
    
    LocalDofMap* pGetDofMap() const {
        LocalDofMap* pDofMap = oU.pGetDofMap();
        
        if (pDofMap == 0) {
            pDofMap = oV.pGetDofMap();
        } else {
            GRADIENT_ASSERT(oV.pGetDofMap() == 0 || pDofMap == oV.pGetDofMap());
        }
        
        return pDofMap;
    }
    
    static index_type iGetMaxDerivatives() {
    	return iMaxDerivatives;
    }

private:
    void Compute() const {
        if (bComputed) {
            return;
        }
        
        const doublereal u = oU.dGetValue();
        const doublereal v = oV.dGetValue();
        
        f = BoolFunc::f(u, v);
        
        bComputed = true;
    }
    
    bool bHaveReferenceTo(const void* p) const {
        return oU.bHaveReferenceTo(p) || oV.bHaveReferenceTo(p);
    }
    
private:
    const LhsExpr oU;
    const RhsExpr oV;
    mutable bool f, bComputed;
};

class FuncPlus {
public:
    static doublereal f(doublereal u, doublereal v) {
        return u + v;
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return 1.;
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return 1.;
    }
};

class FuncMinus {
public:
    static doublereal f(doublereal u, doublereal v) {
        return u - v;
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return 1.;
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return -1.;
    }
};

class FuncMult {
public:
    static doublereal f(doublereal u, doublereal v) {
        return u * v;
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return v;
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return u;
    }
};

class FuncDiv {
public:
    static doublereal f(doublereal u, doublereal v) {
        return u / v;
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return 1. / v;
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return -u / (v * v);
    }
};

class FuncPow {
public:
    static doublereal f(doublereal u, doublereal v) {
        return pow(u, v);
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return v * pow(u, v - 1.);
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return pow(u, v) * log(u);
    }
};

class FuncAtan2 {
public:
    static doublereal f(doublereal u, doublereal v) {
        return atan2(u, v);
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return v / (v * v + u * u);
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return -u / (v * v + u * u);
    }
};

class FuncCopysign {
public:
    static doublereal f(doublereal u, doublereal v) {
        return copysign(u, v);
    }

    static doublereal df_du(doublereal u, doublereal v) {
        return copysign(1., u) * copysign(1., v);
    }

    static doublereal df_dv(doublereal u, doublereal v) {
        return 0.;
    }
};

class FuncFabs {
public:
    static doublereal f(doublereal u) {
        return fabs(u);
    }
    
    static doublereal df_du(doublereal u) {
        return copysign(1., u);
    }
};

class FuncSqrt {
public:
    static doublereal f(doublereal u) {
        return sqrt(u);
    }
    
    static doublereal df_du(doublereal u) {
        return 1. / (2. * sqrt(u));
    }
};

class FuncExp {
public:
    static doublereal f(doublereal u) {
        return exp(u);
    }
    
    static doublereal df_du(doublereal u) {
        return exp(u);
    }
};

class FuncLog {
public:
    static doublereal f(doublereal u) {
        return log(u);
    }
    
    static doublereal df_du(doublereal u) {
        return 1. / u;
    }
};

class FuncSin {
public:
    static doublereal f(doublereal u) {
        return sin(u);
    }

    static doublereal df_du(doublereal u) {
        return cos(u);
    }
};

class FuncCos {
public:
    static doublereal f(doublereal u) {
        return cos(u);
    }

    static doublereal df_du(doublereal u) {
        return -sin(u);
    }
};

class FuncTan {
public:
    static doublereal f(doublereal u) {
        return tan(u);
    }

    static doublereal df_du(doublereal u) {
        const doublereal tan_u = tan(u);
        return 1. + tan_u * tan_u;
    }
};

class FuncSinh {
public:
    static doublereal f(doublereal u) {
        return sinh(u);
    }

    static doublereal df_du(doublereal u) {
        return cosh(u);
    }
};

class FuncCosh {
public:
    static doublereal f(doublereal u) {
        return cosh(u);
    }

    static doublereal df_du(doublereal u) {
        return sinh(u);
    }
};

class FuncTanh {
public:
    static doublereal f(doublereal u) {
        return tanh(u);
    }

    static doublereal df_du(doublereal u) {
        const doublereal tanh_u = tanh(u);
        return 1. - tanh_u * tanh_u;
    }
};

class FuncAsin {
public:
    static doublereal f(doublereal u) {
        return asin(u);
    }

    static doublereal df_du(doublereal u) {
        return 1. / sqrt(1 - u * u);
    }
};

class FuncAcos {
public:
    static doublereal f(doublereal u) {
        return acos(u);
    }

    static doublereal df_du(doublereal u) {
        return -1. / sqrt(1 - u * u);
    }
};

class FuncAtan {
public:
    static doublereal f(doublereal u) {
        return atan(u);
    }

    static doublereal df_du(doublereal u) {
        return 1. / (1. + u * u);
    }
};

#if HAVE_ASINH
class FuncAsinh {
public:
    static doublereal f(doublereal u) {
        return asinh(u);
    }

    static doublereal df_du(doublereal u) {
        return 1. / sqrt(1. + u * u);
    }
};
#endif

#if HAVE_ACOSH
class FuncAcosh {
public:
    static doublereal f(doublereal u) {
        return acosh(u);
    }

    static doublereal df_du(doublereal u) {
        return 1. / sqrt(u * u - 1.);
    }
};
#endif

#if HAVE_ATANH
class FuncAtanh {
public:
    static doublereal f(doublereal u) {
        return atanh(u);
    }

    static doublereal df_du(doublereal u) {
        return 1. / (1. - u * u);
    }
};
#endif

class FuncUnaryMinus {
public:
    static doublereal f(doublereal u) {
        return -u;
    }

    static doublereal df_du(doublereal u) {
        return -1.;
    }
};

class FuncLessThan {
public:
    static bool f(doublereal u, doublereal v) {
        return u < v;
    }
};

class FuncLessEqual {
public:
    static bool f(doublereal u, doublereal v) {
        return u <= v;
    }
};

class FuncGreaterThan {
public:
    static bool f(doublereal u, doublereal v) {
        return u > v;
    }
};

class FuncGreaterEqual {
public:
    static bool f(doublereal u, doublereal v) {
        return u >= v;
    }
};

class FuncEqualTo {
public:
    static bool f(doublereal u, doublereal v) {
        return u == v;
    }
};

class FuncNotEqualTo {
public:
    static bool f(doublereal u, doublereal v) {
        return u != v;
    }
};

template <index_type N_SIZE>
class Gradient {
public:
	static const index_type iMaxDerivatives = MapVector<N_SIZE>::iMaxSize;

    explicit Gradient(doublereal a = 0., LocalDofMap* pDofMap = 0)
        :a(a), ad(pDofMap){
    }

    Gradient(doublereal a, const MapVector<N_SIZE>& da)
        :a(a), ad(da) {
            
    }

    template <index_type N_SIZE2>
    Gradient(const Gradient<N_SIZE2>& g)
    	:a(g.a), ad(g.ad) {

    }

    template <index_type N_SIZE2>
    Gradient(const Gradient<N_SIZE2>& g, LocalDofMap* pDofMap) {
    	Copy(g, pDofMap);
    }

    template <typename Expression>
    Gradient(const GradientExpression<Expression>& f)
        :a(f.dGetValue()),
         ad(f.pGetDofMap(), f.iGetStartIndexLocal(), f.iGetEndIndexLocal(), MapVector<N_SIZE>::LOCAL, 0.) {

        ApplyDerivative(f);
    }

    template <index_type N_SIZE2>
    void Copy(const Gradient<N_SIZE2>& g, LocalDofMap* pDofMap) {
    	LocalDofMap* pDofMap2 = g.pGetDofMap();

    	index_type iFirstLocal = std::numeric_limits<index_type>::max();
    	index_type iLastLocal = 0;

    	GRADIENT_TRACE("g=" << g << std::endl);
    	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
    	GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);

    	for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
    		if (g.dGetDerivativeLocal(i) == 0.) {
    			// FIXME: Checking for zero could cause problems with the KLU linear solver
    			// 		  because the matrix structure could change.
    			// But we have to omit zeros here, otherwise we could allocate
    			// too much local degrees of freedom and cause a buffer overflow
    			// in case of stack based Gradients.
    			continue;
    		}

    		const index_type iGlobal = pDofMap2->iGetGlobalDof(i);
    		const index_type iLocal = pDofMap->AllocateLocalDof(iGlobal);

    		if (iLocal < iFirstLocal) {
    			iFirstLocal = iLocal;
    		}

    		if (iLocal > iLastLocal) {
    			iLastLocal = iLocal;
    		}
    	}

    	++iLastLocal;

    	if (iLastLocal <= iFirstLocal) {
    		iFirstLocal = iLastLocal = 0;
    	}

    	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);

    	a = g.dGetValue();
		ad.ResizeReset(pDofMap, iFirstLocal, iLastLocal, MapVectorBase::LOCAL, 0.);

		GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);
		GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);

		for (index_type i = iFirstLocal; i < iLastLocal; ++i) {
			const index_type iGlobal = pDofMap->iGetGlobalDof(i);
			const index_type iLocal2 = pDofMap2->iGetLocalIndex(iGlobal);

			if (iLocal2 == LocalDofMap::INVALID_INDEX) {
				// Note: It can happen that there is no
				//		 corresponding entry in pDofMap2.
				//		 However it should be safe to ignore such entries.
				continue;
			}

			ad.SetLocalVector(i, g.dGetDerivativeLocal(iLocal2));
		}

#if GRADIENT_DEBUG > 0
		for (index_type i = ad.iGetStartIndexLocal(); i < ad.iGetEndIndexLocal(); ++i) {
			const index_type iGlobal = ad.iGetGlobalDof(i);
			const index_type iLocal2 = g.pGetDofMap()->iGetLocalIndex(iGlobal);
			if (iLocal2 != LocalDofMap::INVALID_INDEX) {
				GRADIENT_ASSERT(g.dGetDerivativeGlobal(iGlobal) == ad.dGetLocalVector(i));
			}
		}

		for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
			if (g.dGetDerivativeLocal(i) != 0.) {
				const index_type iGlobal = g.iGetGlobalDof(i);
				GRADIENT_ASSERT(ad.dGetGlobalVector(iGlobal) == g.dGetDerivativeLocal(i));
			}
		}
#endif
    }
    
    void Reset() {
    	a = 0.;
    	ad.Reset();
    }

    void Reserve(index_type iSize) {
    	ad.Reserve(iSize);
    }

    Gradient& operator=(doublereal d) {
    	// This operator is needed for matrix/vector expressions with different base types
    	SetValue(d);
    	return *this;
    }

    template <typename Expression>
    Gradient& operator=(const GradientExpression<Expression>& f) {
    	GRADIENT_ASSERT(f.iGetMaxDerivatives() <= iGetMaxDerivatives());

        if (!f.bHaveReferenceTo(this)) {
            a = f.dGetValue();

            ad.ResizeReset(f.pGetDofMap(), f.iGetStartIndexLocal(), f.iGetEndIndexLocal(), MapVector<N_SIZE>::LOCAL, 0.);

            ApplyDerivative(f);
        } else {
            // Attention: If we have an expression like a = (a + x)
            // we must not overwrite the contents of a 
            // until the expression (a + x) has been evaluated
            const Gradient tmp(f);
            *this = tmp;
        }
        
        return *this;
    }
    
    template <index_type N_SIZE2>
    Gradient& operator=(const Gradient<N_SIZE2>& g) {
    	a = g.dGetValue();
    	ad = g.GetDerivativeLocal();

    	return *this;
    }

    Gradient& operator=(const Gradient& g) {
    	if (&g != this) {
    		a = g.a;
    		ad = g.ad;
    	}

    	return *this;
    }

    inline Gradient& operator+=(const Gradient& g) {
    	ApplyBinaryFunction<FuncPlus>(GradientExpression<DirectExpr<Gradient> >(g));
    	return *this;
    }

    inline Gradient& operator-=(const Gradient& g) {
    	ApplyBinaryFunction<FuncMinus>(GradientExpression<DirectExpr<Gradient> >(g));
    	return *this;
    }

    inline Gradient& operator*=(const Gradient& g) {
    	ApplyBinaryFunction<FuncMult>(GradientExpression<DirectExpr<Gradient> >(g));
    	return *this;
    }

    inline Gradient& operator/=(const Gradient& g) {
    	ApplyBinaryFunction<FuncDiv>(GradientExpression<DirectExpr<Gradient> >(g));
    	return *this;
    }

    inline Gradient& operator++() {
        ++a;
        return *this;
    }

    inline Gradient& operator--(){
        --a;
        return *this;
    }

    inline Gradient operator++(int){
        const Gradient<N_SIZE> tmp(*this);
        a++;
        return tmp;
    }

    inline Gradient operator--(int) {
        const Gradient<N_SIZE> tmp(*this);
        a--;
        return tmp;
    }
    
    template <typename Expression>
    inline Gradient& operator+=(const GradientExpression<Expression>& f) {
    	ApplyBinaryFunction<FuncPlus>(f);
    	return *this;
    }
    
    template <typename Expression>
    inline Gradient& operator-=(const GradientExpression<Expression>& f) {
    	ApplyBinaryFunction<FuncMinus>(f);
    	return *this;
    }
    
    template <typename Expression>
    inline Gradient& operator*=(const GradientExpression<Expression>& f) {
    	ApplyBinaryFunction<FuncMult>(f);
    	return *this;
    }
    
    template <typename Expression>
    inline Gradient& operator/=(const GradientExpression<Expression>& f) {
    	ApplyBinaryFunction<FuncDiv>(f);
    	return *this;
    }
    
    inline Gradient& operator+=(doublereal d){
        a += d;
        return *this;
    }

    inline Gradient& operator-=(doublereal d){
        a -= d;
        return *this;
    }

    inline Gradient& operator*=(doublereal d) {
    	a *= d;

    	for (index_type i = iGetStartIndexLocal(); i < iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, ad.dGetLocalVector(i) * d);
    	}

        return *this;
    }

    inline Gradient& operator/=(doublereal d) {
    	a /= d;

    	for (index_type i = iGetStartIndexLocal(); i < iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, ad.dGetLocalVector(i) / d);
    	}

        return *this;
    }
    
    doublereal dGetValue() const {
        return a;
    }

    void SetValue(doublereal dVal) {
        a = dVal;
        ad.ResizeReset(0, 0, 0, MapVector<N_SIZE>::LOCAL, 0.);
    }
    
    void SetValuePreserve(doublereal dVal) {
    	// Keep the same derivatives - needed for the Node class for example
        a = dVal;
    }
    
    doublereal dGetDerivativeLocal(index_type iLocalDof) const {
    	GRADIENT_ASSERT(iLocalDof >= 0);
    	GRADIENT_ASSERT(iLocalDof < iGetMaxDerivatives());
        return ad.dGetLocalVector(iLocalDof);
    }
    
    const MapVector<N_SIZE>& GetDerivativeLocal() const {
    	return ad;
    }

    doublereal dGetDerivativeGlobal(index_type iGlobalDof) const {
    	GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) >= 0);
    	GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) < iGetMaxDerivatives());
    	return ad.dGetGlobalVector(iGlobalDof);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, MapVectorBase::GlobalScope s, doublereal dVal) {
    	ad.ResizeReset(pMap, iStartGlobal, iEndGlobal, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, MapVectorBase::LocalScope s, doublereal dVal) {
    	ad.ResizeReset(pMap, iStartLocal, iEndLocal, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iGlobal, MapVectorBase::GlobalScope s, doublereal dVal) {
    	DerivativeResizeReset(pMap, iGlobal, iGlobal + 1, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iLocal, MapVectorBase::LocalScope s, doublereal dVal) {
    	DerivativeResizeReset(pMap, iLocal, iLocal + 1, s, dVal);
    }

    void SetDerivativeLocal(index_type iLocalDof, doublereal dCoef) {
    	GRADIENT_ASSERT(iLocalDof >= iGetStartIndexLocal());
    	GRADIENT_ASSERT(iLocalDof < iGetEndIndexLocal());
        ad.SetLocalVector(iLocalDof, dCoef);
    }

    void SetDerivativeGlobal(index_type iGlobalDof, doublereal dCoef) {
        GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) >= iGetStartIndexLocal());
        GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) < iGetEndIndexLocal());
        ad.SetGlobalVector(iGlobalDof, dCoef);
    }

    index_type iGetGlobalDof(index_type iLocalDof) const {
    	return ad.iGetGlobalDof(iLocalDof);
    }

    index_type iGetStartIndexLocal() const {
        return ad.iGetStartIndexLocal();
    }
    
    index_type iGetEndIndexLocal() const {
        return ad.iGetEndIndexLocal();
    }
    
    index_type iGetLocalSize() const {
    	return iGetEndIndexLocal() - iGetStartIndexLocal();
    }

    index_type iGetMaxDerivatives() const {
    	return ad.iGetMaxSize();
    }

    bool bIsEqual(const Gradient& g) const {
        if (dGetValue() != g.dGetValue()) {
            return false;
        }
        
        const index_type iStart = std::min(iGetStartIndexLocal(), g.iGetStartIndexLocal());
        const index_type iEnd = std::max(iGetEndIndexLocal(), g.iGetEndIndexLocal());
        
        for (index_type i = iStart; i < iEnd; ++i) {
            if (dGetDerivativeLocal(i) != g.dGetDerivativeLocal(i)) {
                return false;
            }
        }
        
        return true;
    }
    
    template <typename Expression>
    bool bIsEqual(const GradientExpression<Expression>& g) const {
    	return bIsEqual(Gradient(g));
    }

    LocalDofMap* pGetDofMap() const {
        return ad.pGetDofMap();
    }

    bool bHaveReferenceTo(const void* p) const {
        return this == p;
    }
    
private:
    template <typename Expression>
    void ApplyDerivative(const GradientExpression<Expression>& f) {
    	// compile time check for the maximum number of derivatives
    	typedef typename MaxSizeCheck<iMaxDerivatives >= Expression::iMaxDerivatives>::CheckType check_iMaxDerivatives;

    	GRADIENT_ASSERT(f.iGetMaxDerivatives() <= iGetMaxDerivatives());

        for (index_type i = ad.iGetStartIndexLocal(); i < ad.iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, f.dGetDerivativeLocal(i));
        }
    }
    
    template <typename BinFunc, typename Expression>
    void ApplyBinaryFunction(const GradientExpression<Expression>& f) {
    	typedef typename MaxSizeCheck<iMaxDerivatives >= Expression::iMaxDerivatives>::CheckType check_iMaxDerivatives;

    	if (!f.bHaveReferenceTo(this)) {
        	LocalDofMap* pDofMap = pGetDofMap();
        	LocalDofMap* const pDofMap2 = f.pGetDofMap();

        	if (pDofMap == 0) {
        		pDofMap = pDofMap2;
        	}

    		const doublereal u = a;
    		const doublereal v = f.dGetValue();

    		a = BinFunc::f(u, v);
    		const doublereal df_du = BinFunc::df_du(u, v);
    		const doublereal df_dv = BinFunc::df_dv(u, v);

    		const index_type iStartFunc = f.iGetStartIndexLocal();
    		const index_type iEndFunc = f.iGetEndIndexLocal();

          	if (pDofMap == f.pGetDofMap()) {
				index_type iStartLocal = std::min(iGetStartIndexLocal(), iStartFunc);
				index_type iEndLocal = std::max(iGetEndIndexLocal(), iEndFunc);

				ad.ResizePreserve(pDofMap, iStartLocal, iEndLocal, MapVector<N_SIZE>::LOCAL);

				if (df_du == 1.) {
					// Optimize for operator+=() and operator-=()
					iStartLocal = iStartFunc;
					iEndLocal = iEndFunc;
				}

				for (index_type i = iStartLocal; i < iEndLocal; ++i) {
					const doublereal ud = ad.dGetLocalVector(i);
					const doublereal vd = f.dGetDerivativeLocal(i);
					ad.SetLocalVector(i, df_du * ud + df_dv * vd);
				}
          	} else {
          		index_type iFirstLocal = std::numeric_limits<index_type>::max();
          		index_type iLastLocal = 0;

            	for (index_type i = iStartFunc; i < iEndFunc; ++i) {
            		if (f.dGetDerivativeLocal(i) == 0.) {
            			// FIXME: In case of complex expressions, this extra function call might be expensive

            			// FIXME: Checking for zero could cause problems with the KLU linear solver
            			// 		  because the matrix structure could change.
            			// But we have to omit zeros here, otherwise we could allocate
            			// too much local degrees of freedom and cause a buffer overflow
            			// in case of stack based Gradients.
            			continue;
            		}

            		const index_type iGlobal = pDofMap2->iGetGlobalDof(i);
            		const index_type iLocal = pDofMap->AllocateLocalDof(iGlobal);

            		if (iLocal < iFirstLocal) {
            			iFirstLocal = iLocal;
            		}

            		if (iLocal > iLastLocal) {
            			iLastLocal = iLocal;
            		}
            	}

            	++iLastLocal;

            	if (iLastLocal <= iFirstLocal) {
            		// empty vector
              		iFirstLocal = std::numeric_limits<index_type>::max();
              		iLastLocal = 0;
            	}

            	GRADIENT_TRACE("iFirstLocal=" << iFirstLocal << std::endl);
            	GRADIENT_TRACE("iLastLocal=" << iLastLocal << std::endl);
            	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);

        		ad.ResizePreserve(pDofMap,
        						  std::min(ad.iGetStartIndexLocal(), iFirstLocal),
        						  std::max(ad.iGetEndIndexLocal(), iLastLocal),
        						  MapVectorBase::LOCAL);

        		GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);
        		GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
        		GRADIENT_TRACE("ad=" << ad.iGetStartIndexLocal() << ":" << ad.iGetEndIndexLocal() << std::endl);

        		if (df_du == 1.) {
        			// optimized loop for operator+= and operator-=
        			for (index_type i = iStartFunc; i < iEndFunc; ++i) {
        				const doublereal vd = f.dGetDerivativeLocal(i);

        				if (vd == 0.) {
        					// no effect for the current index
        					continue;
        				}

        				const index_type iGlobal = pDofMap2->iGetGlobalDof(i);
        				const index_type iLocal = pDofMap->iGetLocalIndex(iGlobal);

        				GRADIENT_ASSERT(iLocal != LocalDofMap::INVALID_INDEX); // because vd != 0

        				const doublereal ud = ad.dGetLocalVector(iLocal);

        				GRADIENT_TRACE("i=" << i << " iLocal=" << iLocal << " ud=" << ud << " vd=" << vd << std::endl);

        				ad.SetLocalVector(iLocal, ud + df_dv * vd);
        			}
        		} else {
        			// generic loop for operator*= and operator/=
					for (index_type i = ad.iGetStartIndexLocal(); i < ad.iGetEndIndexLocal(); ++i) {
						const doublereal ud = ad.dGetLocalVector(i);
						const index_type iGlobal = pDofMap->iGetGlobalDof(i);
						const index_type iLocal2 = pDofMap2->iGetLocalIndex(iGlobal);
						doublereal vd;

						if (iLocal2 == LocalDofMap::INVALID_INDEX) {
							// Note: It can happen that there is no
							//		 corresponding entry in pDofMap2.
							//		 However it should be safe to ignore such entries.
							vd = 0.;
						} else {
							vd = f.dGetDerivativeLocal(iLocal2);
						}

        				GRADIENT_TRACE("i=" << i << " iLocal2=" << iLocal2 << " ud=" << ud << " vd=" << vd << std::endl);

						ad.SetLocalVector(i, df_du * ud + df_dv * vd);
					}
        		}
          	}
    	} else {
    		const Gradient tmp(GradientExpression<BinaryExpr<BinFunc, DirectExpr<Gradient>, Expression> >(*this, f));
    		*this = tmp;
    	}
    }

private:
    doublereal a;
    MapVector<N_SIZE> ad;
};

// helper functions needed for templates (e.g. Matrix<T>, Vector<T>)
inline void Copy(doublereal& d1, const doublereal& d2, LocalDofMap*) {
	d1 = d2;
}

template <index_type N_SIZE1, index_type N_SIZE2>
inline void Copy(Gradient<N_SIZE1>& g1, const Gradient<N_SIZE2>& g2, LocalDofMap* pDofMap) {
	g1.Copy(g2, pDofMap);
}

inline void Reset(doublereal& d) {
	d = 0.;
}

template <index_type N_SIZE>
inline void Reset(Gradient<N_SIZE>& g) {
	g.Reset();
}

inline doublereal dGetValue(doublereal d) {
	return d;
}

template <index_type N_SIZE>
inline doublereal dGetValue(const Gradient<N_SIZE>& g) {
	return g.dGetValue();
}

template <index_type N_SIZE>
inline void Convert(doublereal& d, const Gradient<N_SIZE>& g) {
	// Attention: This operation must be explicit!
	d = g.dGetValue();
}

template <typename T1, typename T2>
inline void Convert(T1& g1, const T2& g2) {
	g1 = g2;
}

template <index_type N_SIZE>
inline std::ostream& operator<<(std::ostream& os, const Gradient<N_SIZE>& f) {
    os << std::setw(12) << f.dGetValue();
    const LocalDofMap* pDofMap = f.pGetDofMap();
    
    os << " [" << f.iGetStartIndexLocal() << "->" << f.iGetEndIndexLocal() << "]";

    if (pDofMap) {    
        os << " (";
        
        for (LocalDofMap::GlobalIterator i = pDofMap->BeginGlobal(); i != pDofMap->EndGlobal(); ++i) {
        	const index_type iGlobal = i->first;
        	const index_type iLocal = i->second;
        	os << std::setw(3) << iGlobal << ":" << std::setw(12) << f.dGetDerivativeLocal(iLocal) << " ";
        }

        os << ")";
    
    }
    
    return os;
}    

#define GRADIENT_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
    template <typename LhsExpr, typename RhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, LhsExpr, RhsExpr> > \
    FunctionName(const GradientExpression<LhsExpr>& u, const GradientExpression<RhsExpr>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, LhsExpr, RhsExpr> >(u, v); \
    } \
    \
    template <index_type N_SIZE, typename LhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, LhsExpr, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(const GradientExpression<LhsExpr>& u, const Gradient<N_SIZE>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, LhsExpr, DirectExpr<Gradient<N_SIZE> > > >(u, v); \
    } \
    \
    template <typename LhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, LhsExpr, ConstExpr> > \
    FunctionName(const GradientExpression<LhsExpr>& u, doublereal v) { \
        return GradientExpression<ExpressionName<FunctionClass, LhsExpr, ConstExpr> >(u, v); \
    } \
    template <index_type N_SIZE, typename RhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, RhsExpr> > \
    FunctionName(const Gradient<N_SIZE>& u, const GradientExpression<RhsExpr>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, RhsExpr> >(u, v); \
    } \
    template <typename RhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, ConstExpr, RhsExpr> > \
    FunctionName(doublereal u, const GradientExpression<RhsExpr>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, ConstExpr, RhsExpr> >(u, v); \
    } \
    template <index_type N_SIZE> \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(const Gradient<N_SIZE>& u, const Gradient<N_SIZE>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > (u, v); \
    } \
    template <index_type N_SIZE> \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, ConstExpr> > \
    FunctionName(const Gradient<N_SIZE>& u, doublereal v) { \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, ConstExpr> >(u, v); \
    } \
    \
    template <index_type N_SIZE> \
    inline GradientExpression<ExpressionName<FunctionClass, ConstExpr, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(doublereal u, const Gradient<N_SIZE>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, ConstExpr, DirectExpr<Gradient<N_SIZE> > > >(u, v); \
    }  

#define GRADIENT_DEFINE_UNARY_FUNCTION(FunctionName, FunctionClass) \
    template <typename Expr> \
    inline GradientExpression<UnaryExpr<FunctionClass, Expr> > \
    FunctionName(const GradientExpression<Expr>& u) { \
        return GradientExpression<UnaryExpr<FunctionClass, Expr> >(u); \
    } \
    \
    template <index_type N_SIZE> \
    inline GradientExpression<UnaryExpr<FunctionClass, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(const Gradient<N_SIZE>& u) { \
        return GradientExpression<UnaryExpr<FunctionClass, DirectExpr<Gradient<N_SIZE> > > >(u); \
    }

GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, operator +, FuncPlus)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, operator -, FuncMinus)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, operator *, FuncMult)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, operator /, FuncDiv)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, pow, FuncPow)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, atan2, FuncAtan2)
GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, copysign, FuncCopysign)

GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator <, FuncLessThan)
GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator <=, FuncLessEqual)
GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator >, FuncGreaterThan)
GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator >=, FuncGreaterEqual)
GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator ==, FuncEqualTo)
GRADIENT_DEFINE_BINARY_FUNCTION(BoolExpr, operator !=, FuncNotEqualTo)

GRADIENT_DEFINE_UNARY_FUNCTION(operator-, FuncUnaryMinus)
GRADIENT_DEFINE_UNARY_FUNCTION(fabs, FuncFabs)
GRADIENT_DEFINE_UNARY_FUNCTION(sqrt, FuncSqrt)
GRADIENT_DEFINE_UNARY_FUNCTION(exp, FuncExp)
GRADIENT_DEFINE_UNARY_FUNCTION(log, FuncLog)
GRADIENT_DEFINE_UNARY_FUNCTION(sin, FuncSin)
GRADIENT_DEFINE_UNARY_FUNCTION(cos, FuncCos)
GRADIENT_DEFINE_UNARY_FUNCTION(tan, FuncTan)
GRADIENT_DEFINE_UNARY_FUNCTION(sinh, FuncSinh)
GRADIENT_DEFINE_UNARY_FUNCTION(cosh, FuncCosh)
GRADIENT_DEFINE_UNARY_FUNCTION(tanh, FuncTanh)
GRADIENT_DEFINE_UNARY_FUNCTION(asin, FuncAsin)
GRADIENT_DEFINE_UNARY_FUNCTION(acos, FuncAcos)
GRADIENT_DEFINE_UNARY_FUNCTION(atan, FuncAtan)

#if HAVE_ASINH
GRADIENT_DEFINE_UNARY_FUNCTION(asinh, FuncAsinh)
#endif

#if HAVE_ACOSH
GRADIENT_DEFINE_UNARY_FUNCTION(acosh, FuncAcosh)
#endif

#if HAVE_ATANH
GRADIENT_DEFINE_UNARY_FUNCTION(atanh, FuncAtanh)
#endif

#undef GRADIENT_DEFINE_BINARY_FUNCTION
#undef GRADIENT_DEFINE_UNARY_FUNCTION

}
#endif
