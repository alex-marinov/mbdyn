/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 AUTHOR: Reinhard Resch <r.resch@a1.net>
        Copyright (C) 2013(-2017) all rights reserved.

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
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#if USE_DENSE_HASH_MAP == 1
#include <dense_hash_map>
#elif __cplusplus >= 201103L
#include <unordered_map>
#else
#include <map>
#endif
#include <new>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "myassert.h"
#include "except.h"
#include "ac/f2c.h"

#ifndef GRADIENT_DEBUG
    #ifdef DEBUG
#define GRADIENT_DEBUG 1
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
    struct IntegerTypeTraits
    {
        // FIXME: We should use std::numeric_limits<T>::max() here!
        // FIXME: Unfortunately that value is not accessible at compile time
        static const T iMaxValue = (T(-1) < 0)
            ? ~T(0) & ~(T(1) << (sizeof(T) * CHAR_BIT - 1))
                                : ~T(0);
    };

    static const index_type DYNAMIC_SIZE = IntegerTypeTraits<index_type>::iMaxValue - 1;
    
#ifndef GRADIENT_VECTOR_REGISTER_SIZE
	#define GRADIENT_VECTOR_REGISTER_SIZE 0
#endif

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
            typedef GradientAllocator<U> other;
        };

    GradientAllocator() noexcept { }

    GradientAllocator(const GradientAllocator& a) noexcept
            :std::allocator<T>(a) {
        }

    template <typename U>
    GradientAllocator(const GradientAllocator<U>& a) noexcept
            :std::allocator<T>(a) {
        }

    pointer
        allocate(size_type n, const void* p=0) {
#if GRADIENT_MEMORY_STAT > 0
    	sMemUsage.Inc(n);
#endif
            const pointer pMem = allocate_aligned(GRADIENT_VECTOR_REGISTER_SIZE, n);

#if GRADIENT_DEBUG > 0
            std::memset(pMem, 0xFF, n);
#endif      
            return pMem;
    }

    void
        deallocate(pointer p, size_type n) {
            deallocate_aligned(p, n);

#if GRADIENT_MEMORY_STAT > 0
    	sMemUsage.Dec(n);
#endif
    }

        static T* allocate_aligned(size_type alignment, size_t n, size_type extra_bytes = 0u) {
#if GRADIENT_MEMORY_STAT > 0
    	sMemUsage.Inc(n);
#endif
    	void* p;

    	const size_type byte_size = sizeof(T) * n + extra_bytes;

    	GRADIENT_ASSERT(alignment % sizeof(void*) == 0);

#if defined(HAVE_POSIX_MEMALIGN)
        if (0 != posix_memalign(&p, alignment, byte_size)) {
    		p = 0;
    	}
#elif defined(HAVE_MEMALIGN)
    	p = memalign(alignment, byte_size);
#elif defined(HAVE_ALIGNED_MALLOC)
        p = _aligned_malloc(byte_size, alignment);
#else
        p = malloc(byte_size);
#endif
            if (p == 0) {
    		throw std::bad_alloc();
    	}

#if GRADIENT_DEBUG > 0
            const ptrdiff_t alignment_curr = (reinterpret_cast<const char*>(p) - reinterpret_cast<const char*>(0)) % sizeof(alignment);
        
            if (alignment_curr != 0) {
                silent_cout("address " << p << " has invalid alignment " << alignment_curr << std::endl);
            }

            GRADIENT_ASSERT(alignment_curr == 0);

            std::memset(p, 0xFF, byte_size);
#endif        
    	return reinterpret_cast<T*>(p);
    }


    static void
        deallocate_aligned(pointer p, size_type n) {
#if defined(HAVE_ALIGNED_MALLOC)
        _aligned_free(p);
#else
    	free(p);
#endif

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
                silent_cerr("GradientAllocator<" << typeid(T).name()
    				<< ">\n\tiMaxMem = "
    				<< std::setprecision(3)
    				<< iMaxMem << std::endl
    				<< "\tiCurrMem=" << iCurrMem << std::endl
    				<< "\tiNumAlloc=" << iNumAlloc << std::endl
                            << "\tiNumDealloc=" << iNumDealloc << std::endl);
    	}

    	void Inc(size_t iSize) {
    		++iNumAlloc;
    		iCurrMem += iSize;
    		if (iCurrMem > iMaxMem) {
    			iMaxMem = iCurrMem;
    		}
    	}

    	void Dec(size_t iSize) {
    		++iNumDealloc;
    		iCurrMem -= iSize;
    	}

    private:
    	size_t iCurrMem;
    	size_t iMaxMem;
    	size_t iNumAlloc;
    	size_t iNumDealloc;
    } sMemUsage;
#endif
};

    struct AlignedAlloc {
#if USE_AUTODIFF > 0 && GRADIENT_VECTOR_REGISTER_SIZE > 0
        void* operator new(size_t size) {
		const size_t nBound = GRADIENT_VECTOR_REGISTER_SIZE > sizeof(void*) ? GRADIENT_VECTOR_REGISTER_SIZE : sizeof(void*);
		return GradientAllocator<char>::allocate_aligned(nBound, 0, size);
	}
        void operator delete(void *p) {
		GradientAllocator<char>::deallocate_aligned(reinterpret_cast<char*>(p), 0);
	}
#endif
};

#if GRADIENT_MEMORY_STAT > 0
template <typename T>
typename GradientAllocator<T>::MemStat GradientAllocator<T>::sMemUsage;
#endif

template <typename T>
    inline void array_fill(T* first, T* const last, const T& val) {
        while (first < last) {
            *first++ = val;
        }
    }

    template <typename T>
    inline T* array_copy(const T* first, const T* const last, T* result) {
        while (first < last) {
            *result++ = *first++;
        }
        
        return result;
    }

    template <typename T>
struct RangeVectorTraits {
	static T Zero(){
		return T();
	}

	static T Invalid(){
		return T(NAN);
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

typedef doublereal scalar_func_type; // data type for the value of a function
typedef doublereal scalar_deriv_type; // data type for the value of a functions derivative

template <typename T, index_type N_SIZE=0>
    class RangeVectorBase {
public:

#if GRADIENT_VECTOR_REGISTER_SIZE > 0
	static const int iVectorSize = GRADIENT_VECTOR_REGISTER_SIZE / sizeof(T);
	typedef T __attribute__((vector_size(iVectorSize * sizeof(T)))) vector_type;
#else
	#warning "vectorization not supported for this compiler"
	static const int iVectorSize = 1;
	typedef T vector_type;
#endif

	typedef T scalar_type;

        static index_type iRoundStartIndexVector(index_type iStart) {
    	return iStart / iVectorSize;
    }

        static index_type iRoundEndIndexVector(index_type iEnd) {
            return iEnd / iVectorSize + (iEnd % iVectorSize ? 1 : 0);
    }
};

typedef doublereal scalar_func_type; // data type for the value of a function
typedef doublereal scalar_deriv_type; // data type for the value of a functions derivative
typedef typename RangeVectorBase<scalar_func_type>::vector_type vector_deriv_type;

template <typename T, index_type N_SIZE>
class RangeVector: public RangeVectorBase<T, N_SIZE> {
public:
	static const index_type iVectorSize = RangeVectorBase<T, N_SIZE>::iVectorSize;
	typedef typename RangeVectorBase<T, N_SIZE>::scalar_type scalar_type;
	typedef typename RangeVectorBase<T, N_SIZE>::vector_type vector_type;

	static const index_type iMaxSizeVector = N_SIZE / iVectorSize + (N_SIZE % iVectorSize ? 1 : 0);
	static const index_type iMaxSize = N_SIZE;

        RangeVector() {
#if GRADIENT_DEBUG > 0
        Initialize(0, iMaxSize, RangeVectorTraits<scalar_type>::Invalid());
#endif
    	GRADIENT_ASSERT(bIsAligned());

        iStart = 0;
        iEnd = 0;
        iStartVec = 0;
        iEndVec = 0;

        GRADIENT_ASSERT(bInvariant());
    }
    
        RangeVector(const RangeVector& v) {
        GRADIENT_ASSERT(bIsAligned());

    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());
    }
    
    template <typename T2, index_type N_SIZE2>
        RangeVector(const RangeVector<T2, N_SIZE2>& v) {
        GRADIENT_ASSERT(bIsAligned());

    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());
    }

        explicit RangeVector(index_type iStartNew, index_type iEndNew, const scalar_type& dVal) {
        GRADIENT_ASSERT(bIsAligned());

        Initialize(iStartNew, iEndNew, dVal);

        GRADIENT_ASSERT(bInvariant());
    }

        RangeVector& operator=(const RangeVector& v) {
        GRADIENT_ASSERT(bIsAligned());
        
    	GRADIENT_ASSERT(bInvariant());

    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());

    	return *this;
    }

    template <typename T2, index_type N_SIZE2>
        RangeVector& operator=(const RangeVector<T2, N_SIZE2>& v) {
        GRADIENT_ASSERT(bIsAligned());
        
    	GRADIENT_ASSERT(bInvariant());

    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());

    	return *this;
    }

        void ResizeReset(index_type iStartNew, index_type iEndNew, const T& dVal) {
        GRADIENT_ASSERT(bIsAligned());

    	GRADIENT_ASSERT(bInvariant());

        Initialize(iStartNew, iEndNew, dVal);

        GRADIENT_ASSERT(bInvariant());
    }

        void ResizePreserve(index_type iStartNew, index_type iEndNew) {
        GRADIENT_ASSERT(bIsAligned());
        
    	GRADIENT_ASSERT(bInvariant());

    	iStart = iStartNew;
    	iEnd = iEndNew;
    	iStartVec = iRoundStartIndexVector(iStartNew);
    	iEndVec = iRoundEndIndexVector(iEndNew);

            array_fill(rgArray, rgArray + iStart, RangeVectorTraits<scalar_type>::Zero());
            array_fill(rgArray + iEnd, rgArray + iMaxSize, RangeVectorTraits<scalar_type>::Zero());

    	GRADIENT_ASSERT(bInvariant());
    	}

        void Reset() {
        GRADIENT_ASSERT(bIsAligned());

    	GRADIENT_ASSERT(bInvariant());
    	// Reset the data but preserve iStart and iEnd
            array_fill(rgArrayVec + iStartVec, rgArrayVec + iEndVec, RangeVectorTraits<vector_type>::Zero());

    	GRADIENT_ASSERT(bInvariant());
    	}

        void Reserve(index_type iMaxSizeNew) {
        GRADIENT_ASSERT(bIsAligned());
    	GRADIENT_ASSERT(bInvariant());
    	GRADIENT_ASSERT(iMaxSizeNew <= iMaxSize);
    }

    index_type iGetStartIndex() const { return iStart; }
    index_type iGetEndIndex() const { return iEnd; }
    index_type iGetSize() const { return iEnd - iStart; }
    static index_type iGetMaxSize() { return iMaxSize; }
    index_type iGetStartIndexVector() const { return iStartVec; }
    index_type iGetEndIndexVector() const { return iEndVec; }
    index_type iGetSizeVector() const { return iEndVec - iStartVec; }
    static index_type iGetMaxSizeVector() { return iMaxSizeVector; }
    
        scalar_type GetValue(index_type i) const {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSize);
        GRADIENT_ASSERT((i >= iGetStartIndex() && i < iGetEndIndex()) || rgArray[i] == RangeVectorTraits<scalar_type>::Zero());
        GRADIENT_ASSERT(bIsAligned());

        return rgArray[i];
    }
    
        void SetValue(index_type i, const scalar_type& d) {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSize);
        GRADIENT_ASSERT(i >= iGetStartIndex() && i < iGetEndIndex());
        GRADIENT_ASSERT(bIsAligned());
        
        rgArray[i] = d;
    }

    vector_type GetVectorValue(index_type i) const {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSizeVector);
        GRADIENT_ASSERT(bIsAligned());
        //GRADIENT_ASSERT(i >= iStartVec && i < iEndVec || rgArrayVec[i] == RangeVectorTraits<vector_type>::Zero());

        return rgArrayVec[i];
    }
    
    void SetVectorValue(index_type i, const vector_type& d) {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSizeVector);
        GRADIENT_ASSERT(i >= iGetStartIndexVector() && i < iGetEndIndexVector());
        GRADIENT_ASSERT(bIsAligned());
        
        rgArrayVec[i] = d;
    }

    static bool bUseDynamicMem() { return false; }
    
private:
    using RangeVectorBase<T, N_SIZE>::iRoundStartIndexVector;
    using RangeVectorBase<T, N_SIZE>::iRoundEndIndexVector;

        void Initialize(index_type iStartNew, index_type iEndNew, const scalar_type& dVal) {
    	iStart = iStartNew;
    	iEnd = iEndNew;
    	iStartVec = iRoundStartIndexVector(iStartNew);
    	iEndVec = iRoundEndIndexVector(iEndNew);
        
            array_fill(rgArrayVec, rgArrayVec + iMaxSizeVector, RangeVectorTraits<vector_type>::Zero());
        
        if (dVal != RangeVectorTraits<scalar_type>::Zero())
        {
                array_fill(rgArray + iStart, rgArray + iEnd, dVal);
        }
    }
    
    template <typename T2, index_type N_SIZE2>
        void Copy(const RangeVector<T2, N_SIZE2>& v) {
    	typedef typename MaxSizeCheck<iMaxSize >= RangeVector<T2, N_SIZE2>::iMaxSize>::CheckType check_iMaxSize;

    	iStart = v.iGetStartIndex();
    	iEnd = v.iGetEndIndex();
    	iStartVec = iRoundStartIndexVector(iStart);
    	iEndVec = iRoundEndIndexVector(iEnd);

		// This allows fast access without frequent index checking
            array_fill(rgArray, rgArray + iStart, RangeVectorTraits<scalar_type>::Zero());
            array_fill(rgArray + iEnd, rgArray + iMaxSize, RangeVectorTraits<scalar_type>::Zero());

		CopyData(v);
		}

    template <typename T2, index_type N_SIZE2>
        void CopyData(const RangeVector<T2, N_SIZE2>& v) {
		for (index_type i = iStart; i < iEnd; ++i) {
			SetValue(i, v.GetValue(i));
		}
    }

        void CopyData(const RangeVector& v) {
            array_copy(v.rgArrayVec + v.iStartVec, v.rgArrayVec + v.iEndVec, rgArrayVec + iStartVec);
    }

#if GRADIENT_DEBUG > 0
        bool bInvariant() const {
		GRADIENT_ASSERT(iEndVec <= iMaxSizeVector);
		GRADIENT_ASSERT(iStartVec >= 0);
		GRADIENT_ASSERT(iStartVec <= iEndVec);

		const index_type iStartOffset = iStart - iStartVec * iVectorSize;

		GRADIENT_ASSERT(iStartOffset >= 0);
		GRADIENT_ASSERT(iStartOffset < iVectorSize);

		const index_type iEndOffset = iEndVec * iVectorSize - iEnd;

		GRADIENT_ASSERT(iEndOffset >= 0);
		GRADIENT_ASSERT(iEndOffset < iVectorSize);

            if (iGetSize() > 0) {
                for (index_type i = 0; i < iStart; ++i) {
				GRADIENT_ASSERT(GetValue(i) == RangeVectorTraits<scalar_type>::Zero());
		}

                for (index_type i = iEnd; i < iMaxSize; ++i) {
				GRADIENT_ASSERT(GetValue(i) == RangeVectorTraits<scalar_type>::Zero());
		}
    }

		return true;
    }
    
        bool bIsAligned() const {
#ifdef __GNUC__
    	GRADIENT_ASSERT(__alignof__(rgArrayVec) >= sizeof(vector_type));
#endif        
        const ptrdiff_t alignment = (reinterpret_cast<const char*>(rgArrayVec) - reinterpret_cast<const char*>(0)) % sizeof(vector_type);
        
            if (alignment != 0) {
                silent_cout("address " << rgArrayVec << " has invalid alignment " << alignment << std::endl);
        }
        
#if defined(WIN32) || defined(__CYGWIN__)
            // FIXME: Stack is aligned at 16 byte boundaries on 64bit Windows systems
            // FIXME: Therefore this test will always fail on Windows
            return true;
#else
        GRADIENT_ASSERT(alignment == 0);
        return alignment == 0;
#endif
    }
#endif

    union
    {
    	scalar_type rgArray[iMaxSize];
    	vector_type rgArrayVec[iMaxSizeVector];
};

    index_type iStart, iEnd, iStartVec, iEndVec;
};

template <typename T>
    class RangeVector<T, 0>: public RangeVectorBase<T, 0> {
public:
	static const index_type iVectorSize = RangeVectorBase<T, 0>::iVectorSize;
	typedef typename RangeVectorBase<T, 0>::scalar_type scalar_type;
	typedef typename RangeVectorBase<T, 0>::vector_type vector_type;

        static const index_type iMaxSize = IntegerTypeTraits<index_type>::iMaxValue / sizeof(scalar_type);

	static const index_type iMaxSizeVector = iMaxSize / iVectorSize;

    RangeVector()
            :pData(pNullData()) {
    	GRADIENT_ASSERT(bInvariant());
    }
    
    RangeVector(const RangeVector& v)
            :pData(pNullData()) {
    	ReserveMem(v.iGetStartIndex(), v.iGetEndIndex(), RESIZE);
    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());
    }

    template <typename T2, index_type N_SIZE2>
    RangeVector(const RangeVector<T2, N_SIZE2>& v)
            :pData(pNullData()) {
    	ReserveMem(v.iGetStartIndex(), v.iGetEndIndex(), RESIZE);
    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());
    }

    explicit RangeVector(index_type iStart, index_type iEnd, const scalar_type& dVal)
            :pData(pNullData()) {
    	ReserveMem(iStart, iEnd, RESIZE);
    	Initialize(iStart, iEnd, dVal);

    	GRADIENT_ASSERT(bInvariant());
    }

        ~RangeVector() {
    	GRADIENT_ASSERT(bInvariant());

    	FreeMem();
    }

        RangeVector& operator=(const RangeVector& v) {
    	GRADIENT_ASSERT(bInvariant());

    	ReserveMem(v.iGetStartIndex(), v.iGetEndIndex(), RESIZE);
    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());

    	return *this;
    }

    template <typename T2, index_type N_SIZE2>
        RangeVector& operator=(const RangeVector<T2, N_SIZE2>& v) {
    	GRADIENT_ASSERT(bInvariant());

    	ReserveMem(v.iGetStartIndex(), v.iGetEndIndex(), RESIZE);
    	Copy(v);

    	GRADIENT_ASSERT(bInvariant());

    	return *this;
    }

        void ResizeReset(index_type iStartNew, index_type iEndNew, const scalar_type& dVal) {
    	GRADIENT_ASSERT(bInvariant());

    	ReserveMem(iStartNew, iEndNew, RESIZE);
    	Initialize(iStartNew, iEndNew, dVal);

    	GRADIENT_ASSERT(bInvariant());
    }
    
        void ResizePreserve(index_type iStartNew, index_type iEndNew) {
    	GRADIENT_ASSERT(bInvariant());

            const index_type iStartNewVec = iRoundStartIndexVector(iStartNew);
            const index_type iEndNewVec = iRoundEndIndexVector(iEndNew);

            if (iStartNewVec == iGetStartIndexVector() && iEndNewVec - iStartNewVec <= iGetCapacityVector()) {
                const index_type iEndCurr = iGetEndIndex();

    		ReserveMem(iStartNew, iEndNew, RESIZE);

                const index_type iEndNewRound = iGetEndIndexVector() * iVectorSize;

                for (index_type i = iEndCurr; i < iEndNewRound; ++i) {
                    SetValueUnchecked(i, RangeVectorTraits<T>::Zero());
    }
            } else {
			RangeVector<T, 0> oTmpVec;
			oTmpVec.ReserveMem(iStartNew, iEndNew, RESIZE);
			oTmpVec.Reset();
                const index_type iComStart = std::max(iGetStartIndex(), oTmpVec.iGetStartIndex());
                const index_type iComEnd = std::min(iGetEndIndex(), oTmpVec.iGetEndIndex());

                GRADIENT_ASSERT(iComStart - iGetStartIndex() >= 0);
                GRADIENT_ASSERT(iComEnd - iGetStartIndex() <= iGetSize());
                GRADIENT_ASSERT(iComEnd - oTmpVec.iGetStartIndex() <= oTmpVec.iGetSize());
                GRADIENT_ASSERT(iComStart - oTmpVec.iGetStartIndex() < oTmpVec.iGetSize());
                GRADIENT_ASSERT(iComEnd - iComStart <= oTmpVec.iGetSize());
                GRADIENT_ASSERT(iComEnd - iComStart <= iGetSize());

                for (index_type i = iComStart; i < iComEnd; ++i) {
                    oTmpVec.SetValue(i, GetValue(i));
                }

#if GRADIENT_DEBUG > 0
                for (index_type i = 0; i < iComStart; ++i) {
                    GRADIENT_ASSERT(oTmpVec.GetValue(i) == RangeVectorTraits<T>::Zero());
                }
                
                for (index_type i = iComStart; i < iComEnd; ++i) {
                    GRADIENT_ASSERT(GetValue(i) == oTmpVec.GetValue(i));
                }

                for (index_type i = iComEnd; i < oTmpVec.iGetEndIndex(); ++i) {
                    GRADIENT_ASSERT(oTmpVec.GetValue(i) == RangeVectorTraits<T>::Zero());
                }
#endif

			std::swap(pData, oTmpVec.pData);
    	}

    	GRADIENT_ASSERT(bInvariant());
    }

        void Reset() {
    	GRADIENT_ASSERT(bInvariant());
    	// Reset the data but preserve memory and iStart
            array_fill(beginVec(), endVec(), RangeVectorTraits<vector_type>::Zero());

    	GRADIENT_ASSERT(bInvariant());
    			}

        void Reserve(index_type iMaxSize) {
    	GRADIENT_ASSERT(bInvariant());

    	const index_type iCurrCap = iGetCapacity();

            if (iCurrCap >= iMaxSize) {
    		return;
    		}

            if (iGetSizeVector() > 0) {
    		RangeVector oTmpVec;

    		GRADIENT_ASSERT(oTmpVec.iGetSizeVector() == 0);

    		oTmpVec.Reserve(iMaxSize);
    		oTmpVec = *this;

    		std::swap(pData, oTmpVec.pData);
            } else {
    		ReserveMem(0, iMaxSize, RESERVE);
    }

    	GRADIENT_ASSERT(bInvariant());
    }

        index_type iGetStartIndex() const {
    	return pData->iStart;
    }

        index_type iGetEndIndex() const {
    	return pData->iEnd;
    }

        index_type iGetStartIndexVector() const {
    	return pData->iStartVec;
    }

        index_type iGetEndIndexVector() const {
    	return pData->iEndVec;
    }

        index_type iGetSizeVector() const {
    	return pData->iEndVec - pData->iStartVec;
    }

        index_type iGetSize() const {
    	return pData->iEnd - pData->iStart;
    }

    static index_type iGetMaxSize() { return iMaxSize; }
    static index_type iGetMaxSizeVector() { return iMaxSizeVector; }

        scalar_type GetValue(index_type i) const {
    	GRADIENT_ASSERT(i >= 0);
    	GRADIENT_ASSERT(i < iGetMaxSize());

        return i >= iGetStartIndex() && i < iGetEndIndex()
        		? begin()[i - iGetStartIndexVector() * iVectorSize]
        		: RangeVectorTraits<scalar_type>::Zero();
    }
    
        void SetValue(index_type i, const scalar_type& d) {
        GRADIENT_ASSERT(i >= iGetStartIndex() && i < iGetEndIndex());

            SetValueUnchecked(i, d);
    }
    
        vector_type GetVectorValue(index_type i) const {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSizeVector);

        return i >= iGetStartIndexVector() && i < iGetEndIndexVector()
        		? beginVec()[i - iGetStartIndexVector()]
        		: RangeVectorTraits<vector_type>::Zero();
    }

    void SetVectorValue(index_type i, const vector_type& d) {
        GRADIENT_ASSERT(i >= 0);
        GRADIENT_ASSERT(i < iMaxSizeVector);
        GRADIENT_ASSERT(i >= iGetStartIndexVector() && i < iGetEndIndexVector());

        beginVec()[i - iGetStartIndexVector()] = d;
    }

    static bool bUseDynamicMem() { return true; }
    
private:
        void SetValueUnchecked(index_type i, const scalar_type& d) {
            GRADIENT_ASSERT(i >= 0);
            GRADIENT_ASSERT(i < iGetMaxSize());
            GRADIENT_ASSERT(i >= iGetStartIndexVector() * iVectorSize && i < iGetEndIndexVector() * iVectorSize);
            begin()[i - iGetStartIndexVector() * iVectorSize] = d;
        }

    using RangeVectorBase<T, 0>::iRoundStartIndexVector;
    using RangeVectorBase<T, 0>::iRoundEndIndexVector;

    struct Data
	{
		index_type iStart;
		index_type iEnd;
		index_type iStartVec;
		index_type iEndVec;
		index_type iCapacityVec;

		union
		{
                scalar_type rgArray[1];
                vector_type rgArrayVec[1];
		};
	};

        index_type iGetCapacity() const {
    	return iGetCapacityVector() * iVectorSize;
    }

        index_type iGetCapacityVector() const {
    	return pData->iCapacityVec;
    }

    template <typename T2, index_type N_SIZE2>
        void Copy(const RangeVector<T2, N_SIZE2>& v) {
    	for (index_type i = iGetStartIndex(); i < iGetEndIndex(); ++i) {
    		SetValue(i,  v.GetValue(i));
    	}
    }

        void Copy(const RangeVector& v) {
            array_copy(v.beginVec(), v.endVec(), beginVec());
    }

        void Initialize(index_type iStart, index_type iEnd, const scalar_type& dVal) {
    	Reset();

    	if (dVal != RangeVectorTraits<scalar_type>::Zero())
    	{
                array_fill(begin() + iStart - iGetStartIndexVector() * iVectorSize, begin() + iEnd - iGetStartIndexVector() * iVectorSize, dVal);
    	}
    }

    enum MemFlags
    {
    	RESIZE,
    	RESERVE
    };

        scalar_type* begin() {
    	return pData->rgArray;
    }

        scalar_type* end() {
    	return begin() + iGetSize();
    }

        const scalar_type* begin() const {
    	return pData->rgArray;
    }

        const scalar_type* end() const {
    	return begin() + iGetSize();
    }

        vector_type* beginVec() {
    	return pData->rgArrayVec;
    }

        vector_type* endVec() {
    	return beginVec() + iGetSizeVector();
    }

        const vector_type* beginVec() const {
    	return pData->rgArrayVec;
    }

        const vector_type* endVec() const {
    	return beginVec() + iGetSizeVector();
    }

        void ReserveMem(index_type iStartNew, index_type iEndNew, MemFlags eFlags = RESIZE) {
    	const index_type iCapCurr = iGetCapacityVector();
    	index_type iStartVec = iRoundStartIndexVector(iStartNew);
    	index_type iEndVec = iRoundEndIndexVector(iEndNew);
    	const index_type iSizeVec = iEndVec - iStartVec;

            if (iSizeVec > iCapCurr) {
    		FreeMem();

                GRADIENT_ASSERT(iSizeVec > 0);
                
                pData = GradientAllocator<Data>::allocate_aligned(sizeof(vector_type), 1u, sizeof(vector_type) * (iSizeVec - 1));

    		pData->iCapacityVec = iSizeVec;
    	}

    	GRADIENT_ASSERT(((char*)(pData->rgArrayVec) - (char*)0) % sizeof(vector_type) == 0);

            if (pData == pNullData()) {
    		GRADIENT_ASSERT(iSizeVec == 0);

    		return;
    	}

            if (eFlags == RESERVE) {
    		iStartNew = iEndNew = iStartVec = iEndVec = 0;
    	}

		pData->iStart = iStartNew;
		pData->iEnd = iEndNew;
		pData->iStartVec = iStartVec;
		pData->iEndVec = iEndVec;
    }

        void FreeMem() {
            if (pData != pNullData()) {
    		GradientAllocator<Data>::deallocate_aligned(pData, 1u);
    		pData = pNullData();
    	}
    }

        static Data* pNullData() {
    	GRADIENT_ASSERT(sNullData.iStart == 0);
    	GRADIENT_ASSERT(sNullData.iEnd == 0);
    	GRADIENT_ASSERT(sNullData.iStartVec == 0);
    	GRADIENT_ASSERT(sNullData.iEndVec == 0);
    	GRADIENT_ASSERT(sNullData.iCapacityVec == 0);

    	return const_cast<Data*>(&sNullData);
    }
    
#if GRADIENT_DEBUG > 0
        bool bInvariant() const {
    	GRADIENT_ASSERT(pData != 0);
    	GRADIENT_ASSERT(pData->iStart >= 0);
    	GRADIENT_ASSERT(pData->iStart <= pData->iEnd);
    	GRADIENT_ASSERT(pData->iEnd - pData->iStart <= iVectorSize * pData->iCapacityVec);
    	GRADIENT_ASSERT(pData->iStartVec >= 0);
    	GRADIENT_ASSERT(pData->iStartVec <= pData->iEndVec);
    	GRADIENT_ASSERT(pData->iEndVec - pData->iStartVec <= pData->iCapacityVec);

    	const index_type iStartOffset = pData->iStart - pData->iStartVec * iVectorSize;

    	GRADIENT_ASSERT(iStartOffset >= 0);
    	GRADIENT_ASSERT(iStartOffset < iVectorSize);

    	const index_type iEndOffset = pData->iEndVec * iVectorSize - pData->iEnd;
    
    	GRADIENT_ASSERT(iEndOffset >= 0);
    	GRADIENT_ASSERT(iEndOffset < iVectorSize);
    
    	return true;
    }
#endif

    Data* pData;
    static const Data sNullData;
};

template <typename T>
const typename RangeVector<T, 0>::Data RangeVector<T, 0>::sNullData = {0};

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
#if USE_DENSE_HASH_MAP == 1
        typedef google::dense_hash_map<index_type,
                                       index_type,
                                       std::hash<index_type>,
                                       std::equal_to<index_type>,
                                       GradientAllocator<std::pair<const index_type, index_type> > > MapType;
#elif __cplusplus >= 201103L
        typedef std::unordered_map<index_type,
                                   index_type,
                                   std::hash<index_type>,
                                   std::equal_to<index_type>,
                                   GradientAllocator<std::pair<const index_type, index_type> > > MapType;
                
#else            
    typedef std::map<index_type,
    				 index_type,
    				 std::less<index_type>,
                         GradientAllocator<std::pair<index_type, index_type> > > MapType;
#endif
    typedef MapType::const_iterator GlobalIterator;
    static const index_type INVALID_INDEX = -1;

    explicit LocalDofMap(index_type iMaxSize=0)
    	:eLastCall(UNKNOWN_FUNC) {
#if USE_DENSE_HASH_MAP == 1
            oGlobalToLocal.set_empty_key(-1);
            oGlobalToLocal.set_deleted_key(-2);
#endif
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
	typedef RangeVector<scalar_deriv_type, N_SIZE> RangeVectorType;
	static const index_type iMaxSize = RangeVectorType::iMaxSize;
	typedef typename RangeVectorType::scalar_type scalar_type;
	typedef typename RangeVectorType::vector_type vector_type;

	template <index_type N_SIZE2>
	MapVector(const MapVector<N_SIZE2>& v)
		:pDofMap(v.pGetDofMap()), oRange(v.GetLocalVector()) {

	}

    MapVector(LocalDofMap* pMap=0, index_type iStartLocal=0, index_type iEndLocal=0, LocalScope= LOCAL, scalar_func_type dVal=0.)
        :pDofMap(pMap), oRange(iStartLocal, iEndLocal, dVal) {
    }

    MapVector(LocalDofMap* pMap, index_type iLocal, LocalScope, scalar_func_type dVal)
        :pDofMap(pMap), oRange(iLocal, iLocal + 1, dVal) {
    }    

    MapVector(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, GlobalScope s, scalar_func_type dVal) {
        ResizeReset(pMap, iStartGlobal, iEndGlobal, s, dVal);
    }

    MapVector(LocalDofMap* pMap, index_type iGlobal, GlobalScope s, scalar_func_type dVal) {
        ResizeReset(pMap, iGlobal, iGlobal + 1, s, dVal);
    }  
    
    void ResizeReset(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, LocalScope, scalar_func_type dVal) {
        pDofMap = pMap;
        oRange.ResizeReset(iStartLocal, iEndLocal, dVal);
    }

    void ResizePreserve(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, LocalScope) {
        pDofMap = pMap;
        oRange.ResizePreserve(iStartLocal, iEndLocal);
    }

    void ResizeReset(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, GlobalScope, scalar_func_type dVal) {
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

    scalar_deriv_type dGetGlobalVector(index_type iGlobalDof) const {
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
    
    void SetGlobalVector(index_type iGlobalDof, scalar_deriv_type dValue) {
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
    index_type iGetSizeVector() const { return oRange.iGetSizeVector(); }
    index_type iGetMaxSize() const { return oRange.iGetMaxSize(); }
    index_type iGetMaxSizeVector() const { return oRange.iGetMaxSizeVector(); }
    scalar_type dGetLocalVector(index_type i) const { return oRange.GetValue(i); }
    vector_type dGetLocalVectorVector(index_type i) const { return oRange.GetVectorValue(i); }
    const RangeVectorType& GetLocalVector() const { return oRange; }
    void SetLocalVector(index_type i, scalar_type dValue){ oRange.SetValue(i, dValue); }
    void SetLocalVectorVector(index_type i, vector_type dVector) { oRange.SetVectorValue(i, dVector); }
    index_type iGetStartIndexLocal() const { return oRange.iGetStartIndex(); }
    index_type iGetEndIndexLocal() const { return oRange.iGetEndIndex(); }
    index_type iGetStartIndexLocalVector() const { return oRange.iGetStartIndexVector(); }
    index_type iGetEndIndexLocalVector() const { return oRange.iGetEndIndexVector(); }
    LocalDofMap* pGetDofMap() const { return pDofMap; }
    bool bUseDynamicMem() const { return oRange.bUseDynamicMem(); }
    static const MapVector Zero;
    
private:
    LocalDofMap* pDofMap;
    RangeVectorType oRange;
};

template <typename Expression>
class GradientExpression: public Expression {
public:
	typedef typename Expression::GradientType GradientType;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename Expression::scalar_func_type scalar_func_type;
	typedef typename Expression::scalar_deriv_type scalar_deriv_type;
	typedef typename Expression::vector_deriv_type vector_deriv_type;

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

template <index_type N_SIZE1, index_type N_SIZE2>
    struct GradientSizeHelper {
	static const index_type iDimension = (N_SIZE1 == 0 || N_SIZE2 == 0) ? 0 : (N_SIZE1 > N_SIZE2 ? N_SIZE1 : N_SIZE2);
	typedef Gradient<iDimension> GradientType;
};

template <typename BinFunc, typename LhsExpr, typename RhsExpr>
class BinaryExpr {
public:
	static const bool bAlias = LhsExpr::bAlias || RhsExpr::bAlias;
	static const index_type iMaxDerivatives = MaxDerivatives<LhsExpr, RhsExpr>::iMaxDerivatives;
	static const bool bVectorize = sizeof(typename LhsExpr::vector_deriv_type) == sizeof(typename RhsExpr::vector_deriv_type)
									&& LhsExpr::bVectorize && RhsExpr::bVectorize && BinFunc::bVectorize;
	typedef typename GradientSizeHelper<LhsExpr::iDimension, RhsExpr::iDimension>::GradientType GradientType;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename GradientType::scalar_func_type scalar_func_type;
	typedef typename GradientType::scalar_deriv_type scalar_deriv_type;
	typedef typename GradientType::vector_deriv_type vector_deriv_type;

	typedef LhsExpr LhsExprType;
	typedef RhsExpr RhsExprType;

    BinaryExpr(const LhsExpr& u, const RhsExpr& v)
        :oU(u), oV(v) {
#if GRADIENT_DEBUG > 0
    	f = df_du = df_dv = NAN;
#endif
    }
    
    scalar_func_type dGetValue() const {
    	GRADIENT_ASSERT(!std::isnan(f));
        return f;
    }
    
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
    	GRADIENT_ASSERT(!std::isnan(f));
    	GRADIENT_ASSERT(!std::isnan(df_du));
    	GRADIENT_ASSERT(!std::isnan(df_dv));

        const scalar_deriv_type du_dX = oU.dGetDerivativeLocal(iLocalDof);
        const scalar_deriv_type dv_dX = oV.dGetDerivativeLocal(iLocalDof);

        return EvalDeriv(du_dX, dv_dX);
    }
    
    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
    	GRADIENT_ASSERT(!std::isnan(f));
    	GRADIENT_ASSERT(!std::isnan(df_du));
    	GRADIENT_ASSERT(!std::isnan(df_dv));

        const vector_deriv_type du_dX = oU.dGetDerivativeLocalVector(iLocalVecDof);
        const vector_deriv_type dv_dX = oV.dGetDerivativeLocalVector(iLocalVecDof);

        return EvalDeriv(du_dX, dv_dX);
    }
    
    index_type iGetStartIndexLocal() const {
        return std::min(oU.iGetStartIndexLocal(), oV.iGetStartIndexLocal());
    }
    
    index_type iGetEndIndexLocal() const {
        return std::max(oU.iGetEndIndexLocal(), oV.iGetEndIndexLocal());
    }
    
    index_type iGetStartIndexLocalVector() const {
        return std::min(oU.iGetStartIndexLocalVector(), oV.iGetStartIndexLocalVector());
    }

    index_type iGetEndIndexLocalVector() const {
        return std::max(oU.iGetEndIndexLocalVector(), oV.iGetEndIndexLocalVector());
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

    void Compute() const {
    	GRADIENT_ASSERT(std::isnan(f));
    	GRADIENT_ASSERT(std::isnan(df_du));
    	GRADIENT_ASSERT(std::isnan(df_dv));

    	oU.Compute();
    	oV.Compute();
        
            const auto u = oU.dGetValue();
            const auto v = oV.dGetValue(); // Use auto here in order to support pow(doublereal, integer)
        
        f = BinFunc::f(u, v);
        df_du = BinFunc::df_du(u, v);
        df_dv = BinFunc::df_dv(u, v);
    }
        
private:
    template <typename T>
    T EvalDeriv(T du_dX, T dv_dX) const {
    	return df_du * du_dX + df_dv * dv_dX;
    }
    
private:
    const LhsExpr oU;
    const RhsExpr oV;
    mutable scalar_func_type f;
    mutable scalar_deriv_type df_du, df_dv;
};

template <typename UnFunc, typename Expr>
class UnaryExpr {
public:
	static const bool bAlias = Expr::bAlias;
	static const index_type iMaxDerivatives = Expr::iMaxDerivatives;
	static const bool bVectorize = Expr::bVectorize && UnFunc::bVectorize;
	typedef typename Expr::GradientType GradientType;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename GradientType::scalar_func_type scalar_func_type;
	typedef typename GradientType::scalar_deriv_type scalar_deriv_type;
	typedef typename GradientType::vector_deriv_type vector_deriv_type;

    UnaryExpr(const Expr& u)
        :oU(u) {
#if GRADIENT_DEBUG > 0
    	f = df_du = NAN;
#endif
    }
        
    scalar_func_type dGetValue() const {
    	GRADIENT_ASSERT(!std::isnan(f));
    	GRADIENT_ASSERT(!std::isnan(df_du));
        
        return f;
    }
        
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
    	GRADIENT_ASSERT(!std::isnan(f));
    	GRADIENT_ASSERT(!std::isnan(df_du));
        
        const scalar_deriv_type du_dX = oU.dGetDerivativeLocal(iLocalDof);
        
        return EvalDeriv(du_dX);
    }

    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
    	GRADIENT_ASSERT(!std::isnan(f));
    	GRADIENT_ASSERT(!std::isnan(df_du));

    	const vector_deriv_type du_dX = oU.dGetDerivativeLocalVector(iLocalVecDof);

    	return EvalDeriv(du_dX);
    }
        
    index_type iGetStartIndexLocal() const {
        return oU.iGetStartIndexLocal();
    }
    
    index_type iGetEndIndexLocal() const {
        return oU.iGetEndIndexLocal();
    }
    
    index_type iGetStartIndexLocalVector() const {
        return oU.iGetStartIndexLocalVector();
    }

    index_type iGetEndIndexLocalVector() const {
        return oU.iGetEndIndexLocalVector();
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

    void Compute() const {
    	GRADIENT_ASSERT(std::isnan(f));
    	GRADIENT_ASSERT(std::isnan(df_du));

    	oU.Compute();
        
        const scalar_func_type u = oU.dGetValue();
        
        f = UnFunc::f(u);
        df_du = UnFunc::df_du(u);
    }

        
private:
    template <typename T>
    T EvalDeriv(T du_dX) const
    {
        return df_du * du_dX;
    }
    
private:
    const Expr oU;
    mutable scalar_func_type f;
    mutable scalar_deriv_type df_du;
};

template <typename T, bool ALIAS=false>
class DirectExpr {
public:
	static const bool bAlias = ALIAS;
	typedef typename T::GradientType GradientType;
	static const index_type iMaxDerivatives = GradientType::iMaxDerivatives;
	static const bool bVectorize = true;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename GradientType::scalar_func_type scalar_func_type;
	typedef typename GradientType::scalar_deriv_type scalar_deriv_type;
	typedef typename GradientType::vector_deriv_type vector_deriv_type;

    DirectExpr(const T& g)
        :oG(g) {
            
    }
    
    scalar_func_type dGetValue() const {
        return oG.dGetValue();
    }
        
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
        return oG.dGetDerivativeLocal(iLocalDof);
    }
        
    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
        return oG.dGetDerivativeLocalVector(iLocalVecDof);
    }

    index_type iGetStartIndexLocal() const {
        return oG.iGetStartIndexLocal();
    }

    index_type iGetEndIndexLocal() const {
        return oG.iGetEndIndexLocal();
    }

    index_type iGetStartIndexLocalVector() const {
        return oG.iGetStartIndexLocalVector();
    }

    index_type iGetEndIndexLocalVector() const {
        return oG.iGetEndIndexLocalVector();
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

    void Compute() const {}

private:
    const GradientType& oG;
};

    template <typename T, typename C = scalar_func_type>
class ConstExpr {
public:
	static const bool bAlias = false;
	static const index_type iMaxDerivatives = 0;
	static const bool bVectorize = true;
	typedef T GradientType;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename GradientType::scalar_func_type scalar_func_type;
	typedef typename GradientType::scalar_deriv_type scalar_deriv_type;
	typedef typename GradientType::vector_deriv_type vector_deriv_type;

        ConstExpr(C a)
        :dConst(a) {
            
    }
        
        C dGetValue() const {
        return dConst;
    }
        
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
        return scalar_deriv_type();
    }

    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
    	return vector_deriv_type();
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
    
    index_type iGetStartIndexLocalVector() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocalVector() const {
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

    void Compute() const {}

private:
        const C dConst;
};

template <typename BoolFunc, typename LhsExpr, typename RhsExpr>
class BoolExpr {
public:
	static const bool bAlias = LhsExpr::bAlias || RhsExpr::bAlias;
	static const index_type iMaxDerivatives = MaxDerivatives<LhsExpr, RhsExpr>::iMaxDerivatives;
	static const bool bVectorize = true;
	typedef typename GradientSizeHelper<LhsExpr::iDimension, RhsExpr::iDimension>::GradientType GradientType;
	static const index_type iDimension = GradientType::iDimension;
	typedef typename GradientType::scalar_func_type scalar_func_type;
	typedef typename GradientType::scalar_deriv_type scalar_deriv_type;
	typedef typename GradientType::vector_deriv_type vector_deriv_type;

    BoolExpr(const LhsExpr& u, const RhsExpr& v)
        :oU(u), oV(v) {
    }
    
    bool dGetValue() const {
    	oU.Compute();
    	oV.Compute();
        
        const scalar_func_type u = oU.dGetValue();
        const scalar_func_type v = oV.dGetValue();

        return BoolFunc::f(u, v);
    }
    
    operator bool() const {
        return dGetValue();
    }
    
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
        return scalar_deriv_type();
    }
    
    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalDof) const {
        return vector_deriv_type();
    }
    
    index_type iGetStartIndexLocal() const {
        return std::numeric_limits<index_type>::max();
    }
    
    index_type iGetEndIndexLocal() const {
        return 0;
    }
    
    index_type iGetStartIndexLocalVector() const {
        return std::numeric_limits<index_type>::max();
    }

    index_type iGetEndIndexLocalVector() const {
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

    void Compute() const {
    }
    
private:
    bool bHaveReferenceTo(const void* p) const {
        return oU.bHaveReferenceTo(p) || oV.bHaveReferenceTo(p);
    }
    
private:
    const LhsExpr oU;
    const RhsExpr oV;
};

class FuncPlus {
public:
	static const bool bVectorize = true;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return u + v;
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return 1.;
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return 1.;
    }
};

class FuncMinus {
public:
	static const bool bVectorize = true;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return u - v;
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return 1.;
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return -1.;
    }
};

class FuncMult {
public:
	static const bool bVectorize = true;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return u * v;
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return v;
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return u;
    }
};

class FuncDiv {
public:
	static const bool bVectorize = true;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return u / v;
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return 1. / v;
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return -u / (v * v);
    }
};

class FuncPow {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return pow(u, v);
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return v * pow(u, v - 1.);
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return pow(u, v) * log(u);
    }
};

    class FuncPowInt {
    public:
        static const bool bVectorize = false;

        static scalar_func_type f(scalar_func_type u, integer v) {
            return pow(u, v);
        }

        static scalar_deriv_type df_du(scalar_func_type u, integer v) {
            return v * pow(u, v - 1);
        }

        static scalar_deriv_type df_dv(scalar_func_type u, integer v) {
            return 0.;
        }
    };
    
class FuncAtan2 {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return atan2(u, v);
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return v / (v * v + u * u);
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return -u / (v * v + u * u);
    }
};

class FuncCopysign {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
        return copysign(u, v);
    }

    static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
        return copysign(1., u) * copysign(1., v);
    }

    static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
        return 0.;
    }
};

class FuncFabs {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return fabs(u);
    }
    
    static scalar_deriv_type df_du(scalar_func_type u) {
        return copysign(1., u);
    }
};

class FuncSqrt {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return sqrt(u);
    }
    
    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / (2. * sqrt(u));
    }
};

class FuncExp {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return exp(u);
    }
    
    static scalar_deriv_type df_du(scalar_func_type u) {
        return exp(u);
    }
};

class FuncLog {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return log(u);
    }
    
    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / u;
    }
};

class FuncSin {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return sin(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return cos(u);
    }
};

class FuncCos {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return cos(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return -sin(u);
    }
};

class FuncTan {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return tan(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        const scalar_deriv_type tan_u = tan(u);
        return 1. + tan_u * tan_u;
    }
};

class FuncSinh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return sinh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return cosh(u);
    }
};

class FuncCosh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return cosh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return sinh(u);
    }
};

class FuncTanh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return tanh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        const scalar_deriv_type tanh_u = tanh(u);
        return 1. - tanh_u * tanh_u;
    }
};

class FuncAsin {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return asin(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / sqrt(1 - u * u);
    }
};

class FuncAcos {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return acos(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return -1. / sqrt(1 - u * u);
    }
};

class FuncAtan {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return atan(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / (1. + u * u);
    }
};

#if HAVE_ASINH
class FuncAsinh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return asinh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / sqrt(1. + u * u);
    }
};
#endif

#if HAVE_ACOSH
class FuncAcosh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return acosh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / sqrt(u * u - 1.);
    }
};
#endif

#if HAVE_ATANH
class FuncAtanh {
public:
	static const bool bVectorize = false;

    static scalar_func_type f(scalar_func_type u) {
        return atanh(u);
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return 1. / (1. - u * u);
    }
};
#endif

    class FuncFmod {
    public:
        static const bool bVectorize = false;

        static scalar_func_type f(scalar_func_type u, scalar_func_type v) {
            return fmod(u, v);
        }

        static scalar_deriv_type df_du(scalar_func_type u, scalar_func_type v) {
            return 1.;
        }

        static scalar_deriv_type df_dv(scalar_func_type u, scalar_func_type v) {
            return -int(u / v);
        }
    };
        
class FuncUnaryMinus {
public:
	static const bool bVectorize = true;

    static scalar_func_type f(scalar_func_type u) {
        return -u;
    }

    static scalar_deriv_type df_du(scalar_func_type u) {
        return -1.;
    }
};

class FuncLessThan {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u < v;
    }
};

class FuncLessEqual {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u <= v;
    }
};

class FuncGreaterThan {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u > v;
    }
};

class FuncGreaterEqual {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u >= v;
    }
};

class FuncEqualTo {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u == v;
    }
};

class FuncNotEqualTo {
public:
	static const bool bVectorize = true;

    static bool f(scalar_func_type u, scalar_func_type v) {
        return u != v;
    }
};

template <bool bVectorize>
struct ApplyDerivativeHelper	{	};

template <>
    struct ApplyDerivativeHelper<false> {
    template <typename MapVectorType, typename Expression>
    static void ApplyDerivative(MapVectorType& ad, const GradientExpression<Expression>& f) {
        for (index_type i = ad.iGetStartIndexLocal(); i < ad.iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, f.dGetDerivativeLocal(i));
        }
    }

    template <typename MapVectorType, typename Expression>
    static void ApplyBinaryFunction1(MapVectorType& ad,
    								 const GradientExpression<Expression>& f,
    								 const index_type iStartLocal,
    								 const index_type iEndLocal,
    								 const scalar_deriv_type df_du,
    								 const scalar_deriv_type df_dv) {
		for (index_type i = iStartLocal; i < iEndLocal; ++i) {
			const scalar_deriv_type ud = ad.dGetLocalVector(i);
			const scalar_deriv_type vd = f.dGetDerivativeLocal(i);
			ad.SetLocalVector(i, df_du * ud + df_dv * vd);
		}
    }
};

template <>
    struct ApplyDerivativeHelper<true> {
    template <typename MapVectorType, typename Expression>
    static void ApplyDerivative(MapVectorType& ad, const GradientExpression<Expression>& f) {
        for (index_type i = ad.iGetStartIndexLocalVector(); i < ad.iGetEndIndexLocalVector(); ++i) {
    		ad.SetLocalVectorVector(i, f.dGetDerivativeLocalVector(i));
        }
    }

    template <typename MapVectorType, typename Expression>
    static void ApplyBinaryFunction1(MapVectorType& ad,
    								 const GradientExpression<Expression>& f,
    								 const index_type iStartLocal,
    								 const index_type iEndLocal,
    								 const typename GradientExpression<Expression>::scalar_deriv_type df_du,
    								 const typename GradientExpression<Expression>::scalar_deriv_type df_dv) {
    	typedef RangeVectorBase<typename MapVectorType::scalar_type> RangeVectorType;
    	typedef typename GradientExpression<Expression>::vector_deriv_type vector_deriv_type;
    	const index_type iStartLocalVec = RangeVectorType::iRoundStartIndexVector(iStartLocal);
    	const index_type iEndLocalVec = RangeVectorType::iRoundEndIndexVector(iEndLocal);

		for (index_type i = iStartLocalVec; i < iEndLocalVec; ++i) {
			const vector_deriv_type ud = ad.dGetLocalVectorVector(i);
			const vector_deriv_type vd = f.dGetDerivativeLocalVector(i);
			ad.SetLocalVectorVector(i, df_du * ud + df_dv * vd);
		}
    }
};

template <bool ALIAS>
struct ApplyAliasHelper {};

template <>
    struct ApplyAliasHelper<false> {
	template <typename GradientType, typename Expression>
        static inline void ApplyExpression(GradientType& g, const GradientExpression<Expression>& f) {
		g.ApplyNoAlias(f);
	}

    template <typename BinFunc, typename GradientType, typename Expression>
        static void ApplyBinaryFunction(GradientType& g, const GradientExpression<Expression>& f, const BinFunc& bfunc) {
    	g.ApplyBinaryFunctionNoAlias(f, bfunc);
    }
};

template <>
    struct ApplyAliasHelper<true> {
	template <typename GradientType, typename Expression>
        static inline void ApplyExpression(GradientType& g, const GradientExpression<Expression>& f) {
		g.ApplyWithAlias(f);
	}

    template <typename BinFunc, typename GradientType, typename Expression>
        static void ApplyBinaryFunction(GradientType& g, const GradientExpression<Expression>& f, const BinFunc& bfunc) {
    	g.ApplyBinaryFunctionWithAlias(f, bfunc);
    }
};

template <index_type N_SIZE>
class Gradient {
	typedef MapVector<N_SIZE> MapVectorType;

public:
	typedef Gradient GradientType;
	static const index_type iDimension = N_SIZE;
	static const index_type iMaxDerivatives = MapVectorType::iMaxSize;
	typedef grad::scalar_func_type scalar_func_type;
	typedef typename MapVectorType::scalar_type scalar_deriv_type;
	typedef typename MapVectorType::vector_type vector_deriv_type;

    explicit Gradient(scalar_func_type a = 0., LocalDofMap* pDofMap = 0)
        :a(a), ad(pDofMap){
    }

    Gradient(scalar_func_type a, const MapVectorType& da)
        :a(a), ad(da) {
            
    }

    Gradient(const Gradient& g)
            :a(g.a), ad(g.ad) {
        
    }

    template <index_type N_SIZE2>
    Gradient(const Gradient<N_SIZE2>& g)
    	:a(g.a), ad(g.ad) {

    }

    template <index_type N_SIZE2>
    Gradient(const Gradient<N_SIZE2>& g, LocalDofMap* pDofMap) {
#if GRADIENT_DEBUG > 0
    	a = NAN;
#endif
    	Copy(g, pDofMap);
    }

    template <typename Expression>
    Gradient(const GradientExpression<Expression>& f)
        :ad(f.pGetDofMap(), f.iGetStartIndexLocal(), f.iGetEndIndexLocal(), MapVector<N_SIZE>::LOCAL, 0.) {

    	// No aliases are possible because the object did not exist before
    	GRADIENT_ASSERT(!f.bHaveReferenceTo(this));

    	f.Compute();
    	a = f.dGetValue();

        ApplyDerivative(f);
    }

    template <index_type N_SIZE2>
    void Copy(const Gradient<N_SIZE2>& g, LocalDofMap* pDofMap) {
    	LocalDofMap* pDofMap2 = g.pGetDofMap();

    	index_type iFirstLocal = std::numeric_limits<index_type>::max();
    	index_type iLastLocal = 0;

#if GRADIENT_DEBUG > 0
    	GRADIENT_TRACE("g=" << g << std::endl);
            if (pDofMap) {
    	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
            }
            if (pDofMap2) {
    	GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);
            }
#endif

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

#if GRADIENT_DEBUG > 0
            if (pDofMap) {
    	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
            }
#endif
    	a = g.dGetValue();
		ad.ResizeReset(pDofMap, iFirstLocal, iLastLocal, MapVectorBase::LOCAL, 0.);

#if GRADIENT_DEBUG > 0
            if (pDofMap2) {
		GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);
            }
            if (pDofMap) {
		GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
            }
#endif
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
		for (index_type i = iFirstLocal; i < iLastLocal; ++i) {
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

    Gradient& operator=(scalar_func_type d) {
    	// This operator is needed for matrix/vector expressions with different base types
    	SetValue(d);
    	return *this;
    }

    template <typename Expression>
    Gradient& operator=(const GradientExpression<Expression>& f) {
    	ApplyAliasHelper<GradientExpression<Expression>::bAlias>::ApplyExpression(*this, f);
        
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
    
    inline Gradient& operator+=(scalar_func_type d){
        a += d;
        return *this;
    }

    inline Gradient& operator-=(scalar_func_type d){
        a -= d;
        return *this;
    }

    inline Gradient& operator*=(scalar_func_type d) {
    	a *= d;

    	for (index_type i = iGetStartIndexLocal(); i < iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, ad.dGetLocalVector(i) * d);
    	}

        return *this;
    }

    inline Gradient& operator/=(scalar_func_type d) {
    	a /= d;

    	for (index_type i = iGetStartIndexLocal(); i < iGetEndIndexLocal(); ++i) {
    		ad.SetLocalVector(i, ad.dGetLocalVector(i) / d);
    	}

        return *this;
    }
    
    scalar_func_type dGetValue() const {
        return a;
    }

    void SetValue(scalar_func_type dVal) {
        a = dVal;
        ad.ResizeReset(0, 0, 0, MapVector<N_SIZE>::LOCAL, 0.);
    }
    
    void SetValuePreserve(scalar_func_type dVal) {
    	// Keep the same derivatives - needed for the Node class for example
        a = dVal;
    }
    
    scalar_deriv_type dGetDerivativeLocal(index_type iLocalDof) const {
    	GRADIENT_ASSERT(iLocalDof >= 0);
    	GRADIENT_ASSERT(iLocalDof < iGetMaxDerivatives());
        return ad.dGetLocalVector(iLocalDof);
    }
    
    vector_deriv_type dGetDerivativeLocalVector(index_type iLocalVecDof) const {
    	GRADIENT_ASSERT(iLocalVecDof >= 0);
    	GRADIENT_ASSERT(iLocalVecDof < iGetMaxDerivativesVector());
    	return ad.dGetLocalVectorVector(iLocalVecDof);
    }

    const MapVector<N_SIZE>& GetDerivativeLocal() const {
    	return ad;
    }

    scalar_deriv_type dGetDerivativeGlobal(index_type iGlobalDof) const {
    	GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) >= 0);
    	GRADIENT_ASSERT(ad.pGetDofMap()->iGetLocalIndex(iGlobalDof) < iGetMaxDerivatives());
    	return ad.dGetGlobalVector(iGlobalDof);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iStartGlobal, index_type iEndGlobal, MapVectorBase::GlobalScope s, scalar_deriv_type dVal) {
    	ad.ResizeReset(pMap, iStartGlobal, iEndGlobal, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iStartLocal, index_type iEndLocal, MapVectorBase::LocalScope s, scalar_deriv_type dVal) {
    	ad.ResizeReset(pMap, iStartLocal, iEndLocal, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iGlobal, MapVectorBase::GlobalScope s, scalar_deriv_type dVal) {
    	DerivativeResizeReset(pMap, iGlobal, iGlobal + 1, s, dVal);
    }

    void DerivativeResizeReset(LocalDofMap* pMap, index_type iLocal, MapVectorBase::LocalScope s, scalar_deriv_type dVal) {
    	DerivativeResizeReset(pMap, iLocal, iLocal + 1, s, dVal);
    }

    void SetDerivativeLocal(index_type iLocalDof, scalar_deriv_type dCoef) {
    	GRADIENT_ASSERT(iLocalDof >= iGetStartIndexLocal());
    	GRADIENT_ASSERT(iLocalDof < iGetEndIndexLocal());
        ad.SetLocalVector(iLocalDof, dCoef);
    }

    void SetDerivativeLocalVector(index_type iLocalVecDof, vector_deriv_type dVec) {
    	GRADIENT_ASSERT(iLocalVecDof >= iGetStartIndexLocalVector());
    	GRADIENT_ASSERT(iLocalVecDof < iGetEndIndexLocalVector());
    	ad.SetLocalVectorVector(iLocalVecDof, dVec);
    }

    void SetDerivativeGlobal(index_type iGlobalDof, scalar_deriv_type dCoef) {
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
    
    index_type iGetStartIndexLocalVector() const {
    	return ad.iGetStartIndexLocalVector();
    }

    index_type iGetEndIndexLocalVector() const {
    	return ad.iGetEndIndexLocalVector();
    }

    index_type iGetLocalSize() const {
    	return iGetEndIndexLocal() - iGetStartIndexLocal();
    }

    index_type iGetMaxDerivatives() const {
    	return ad.iGetMaxSize();
    }

    index_type iGetMaxDerivativesVector() const {
    	return ad.iGetMaxSizeVector();
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
    friend struct ApplyAliasHelper<true>;
    friend struct ApplyAliasHelper<false>;

    template <typename Expression>
    void ApplyNoAlias(const GradientExpression<Expression>& f) {
    	GRADIENT_ASSERT(f.iGetMaxDerivatives() <= iGetMaxDerivatives());
		GRADIENT_ASSERT(!f.bHaveReferenceTo(this));

		f.Compute();

		a = f.dGetValue();

		ad.ResizeReset(f.pGetDofMap(), f.iGetStartIndexLocal(), f.iGetEndIndexLocal(), MapVector<N_SIZE>::LOCAL, 0.);

		ApplyDerivative(f);
    }

    template <typename Expression>
    void ApplyWithAlias(const GradientExpression<Expression>& f) {
        // Attention: If we have an expression like a = (a + x)
        // we must not overwrite the contents of a
        // until the expression (a + x) has been evaluated
        *this = Gradient(f);
    }

    template <typename Expression>
    void ApplyDerivative(const GradientExpression<Expression>& f) {
    	// compile time check for the maximum number of derivatives
    	typedef typename MaxSizeCheck<iMaxDerivatives >= Expression::iMaxDerivatives>::CheckType check_iMaxDerivatives;
    	const bool bVectorize = sizeof(vector_deriv_type) == sizeof(typename GradientExpression<Expression>::vector_deriv_type)
    							&& GradientExpression<Expression>::bVectorize;
    	GRADIENT_ASSERT(f.iGetMaxDerivatives() <= iGetMaxDerivatives());

    	ApplyDerivativeHelper<bVectorize>::ApplyDerivative(ad, f);
    }
    
    template <typename BinFunc, typename Expression>
    void ApplyBinaryFunction(const GradientExpression<Expression>& f) {
    	ApplyAliasHelper<GradientExpression<Expression>::bAlias>::ApplyBinaryFunction(*this, f, BinFunc());
    }

    template <typename BinFunc, typename Expression>
    void ApplyBinaryFunctionNoAlias(const GradientExpression<Expression>& f, const BinFunc&) {
    	typedef typename MaxSizeCheck<iMaxDerivatives >= Expression::iMaxDerivatives>::CheckType check_iMaxDerivatives;
    	const bool bVectorize = sizeof(vector_deriv_type) == sizeof(typename GradientExpression<Expression>::vector_deriv_type)
    							&& GradientExpression<Expression>::bVectorize;

    	GRADIENT_ASSERT(!f.bHaveReferenceTo(this));
        	LocalDofMap* pDofMap = pGetDofMap();
        LocalDofMap* pDofMap2 = f.pGetDofMap();

        if (pDofMap2 == 0) {
        	pDofMap2 = pDofMap;
        }

        	if (pDofMap == 0) {
        		pDofMap = pDofMap2;
        	}

        f.Compute();

        const scalar_func_type u = a;
        const scalar_func_type v = f.dGetValue();

    		a = BinFunc::f(u, v);
        const scalar_deriv_type df_du = BinFunc::df_du(u, v);
        const scalar_deriv_type df_dv = BinFunc::df_dv(u, v);

    		const index_type iStartFunc = f.iGetStartIndexLocal();
    		const index_type iEndFunc = f.iGetEndIndexLocal();

        if (pDofMap == pDofMap2) {
                index_type iStartLocal = iGetStartIndexLocal();
                index_type iEndLocal = iGetEndIndexLocal();

                if (iEndFunc > iStartFunc) {
                    iStartLocal = std::min(iStartLocal, iStartFunc);
                    iEndLocal = std::max(iEndLocal, iEndFunc);
                }

				ad.ResizePreserve(pDofMap, iStartLocal, iEndLocal, MapVector<N_SIZE>::LOCAL);

				if (df_du == 1.) {
					// Optimize for operator+=() and operator-=()
					iStartLocal = iStartFunc;
					iEndLocal = iEndFunc;
				}

            ApplyDerivativeHelper<bVectorize>::ApplyBinaryFunction1(ad, f, iStartLocal, iEndLocal, df_du, df_dv);
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

#if GRADIENT_DEBUG > 0
            	GRADIENT_TRACE("iFirstLocal=" << iFirstLocal << std::endl);
            	GRADIENT_TRACE("iLastLocal=" << iLastLocal << std::endl);
                if (pDofMap) {
            	GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
                }
#endif

        		ad.ResizePreserve(pDofMap,
        						  std::min(ad.iGetStartIndexLocal(), iFirstLocal),
        						  std::max(ad.iGetEndIndexLocal(), iLastLocal),
        						  MapVectorBase::LOCAL);

#if GRADIENT_DEBUG > 0
                if (pDofMap2) {
        		GRADIENT_TRACE("*pDofMap2=" << *pDofMap2 << std::endl);
                }
                if (pDofMap) {
        		GRADIENT_TRACE("*pDofMap=" << *pDofMap << std::endl);
                }
        		GRADIENT_TRACE("ad=" << ad.iGetStartIndexLocal() << ":" << ad.iGetEndIndexLocal() << std::endl);
#endif

        		if (df_du == 1.) {
        			// optimized loop for operator+= and operator-=
        			for (index_type i = iStartFunc; i < iEndFunc; ++i) {
                    const scalar_deriv_type vd = f.dGetDerivativeLocal(i);

        				if (vd == 0.) {
        					// no effect for the current index
        					continue;
        				}

        				const index_type iGlobal = pDofMap2->iGetGlobalDof(i);
        				const index_type iLocal = pDofMap->iGetLocalIndex(iGlobal);

        				GRADIENT_ASSERT(iLocal != LocalDofMap::INVALID_INDEX); // because vd != 0

                    const scalar_deriv_type ud = ad.dGetLocalVector(iLocal);

        				GRADIENT_TRACE("i=" << i << " iLocal=" << iLocal << " ud=" << ud << " vd=" << vd << std::endl);

        				ad.SetLocalVector(iLocal, ud + df_dv * vd);
        			}
        		} else {
        			// generic loop for operator*= and operator/=
					for (index_type i = ad.iGetStartIndexLocal(); i < ad.iGetEndIndexLocal(); ++i) {
                    const scalar_deriv_type ud = ad.dGetLocalVector(i);
						const index_type iGlobal = pDofMap->iGetGlobalDof(i);
						const index_type iLocal2 = pDofMap2->iGetLocalIndex(iGlobal);
                    scalar_deriv_type vd;

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
    	}

    template <typename BinFunc, typename Expression>
    void ApplyBinaryFunctionWithAlias(const GradientExpression<Expression>& f, const BinFunc&) {
    	typedef typename MaxSizeCheck<iMaxDerivatives >= Expression::iMaxDerivatives>::CheckType check_iMaxDerivatives;

        *this = Gradient(GradientExpression<BinaryExpr<BinFunc, DirectExpr<Gradient>, Expression> >(*this, f));
    }

private:
    scalar_func_type a;
    MapVectorType ad;
};

// helper functions needed for templates (e.g. Matrix<T>, Vector<T>)
inline void Copy(scalar_func_type& d1, const scalar_func_type& d2, LocalDofMap*) {
	d1 = d2;
}

template <index_type N_SIZE1, index_type N_SIZE2>
inline void Copy(Gradient<N_SIZE1>& g1, const Gradient<N_SIZE2>& g2, LocalDofMap* pDofMap) {
	g1.Copy(g2, pDofMap);
}

    template <typename T>
    inline void Reset(T& d) {
        d = T();
}

template <index_type N_SIZE>
inline void Reset(Gradient<N_SIZE>& g) {
	g.Reset();
}

inline scalar_func_type dGetValue(scalar_func_type d) {
	return d;
}

template <index_type N_SIZE>
inline scalar_func_type dGetValue(const Gradient<N_SIZE>& g) {
	return g.dGetValue();
}

template <index_type N_SIZE>
inline void Convert(scalar_func_type& d, const Gradient<N_SIZE>& g) {
	// Attention: This operation must be explicit!
	d = g.dGetValue();
}

template <typename T1, typename T2>
inline void Convert(T1& g1, const T2& g2) {
	g1 = g2;
}

template <index_type N_SIZE>
inline GradientExpression<DirectExpr<Gradient<N_SIZE>, true> > Alias(const Gradient<N_SIZE>& g)
{
	return GradientExpression<DirectExpr<Gradient<N_SIZE>, true> >(g);
}

inline scalar_func_type Alias(scalar_func_type d)
{
	return d;
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
                const scalar_deriv_type d = f.dGetDerivativeLocal(iLocal);
                
                if (d) {
                    os << iGlobal << ":" << d << " ";
                }
        }

        os << ")";
    
    }
    
    return os;
}    

#define GRADIENT_DEFINE_BINARY_FUNCTION_CONST_ARG_LHS(ExpressionName, FunctionName, FunctionClass, ScalarFuncTypeLhs) \
    template <typename RhsExpr>                                         \
    inline GradientExpression<ExpressionName<FunctionClass, ConstExpr<typename RhsExpr::GradientType, ScalarFuncTypeLhs>, RhsExpr> > \
    FunctionName(ScalarFuncTypeLhs u, const GradientExpression<RhsExpr>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, ConstExpr<typename RhsExpr::GradientType, ScalarFuncTypeLhs>, RhsExpr> >(u, v); \
    }                                                                   \
                                                                        \
    template <index_type N_SIZE>                                        \
    inline GradientExpression<ExpressionName<FunctionClass, ConstExpr<Gradient<N_SIZE>, ScalarFuncTypeLhs>, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(ScalarFuncTypeLhs u, const Gradient<N_SIZE>& v) {       \
        return GradientExpression<ExpressionName<FunctionClass, ConstExpr<Gradient<N_SIZE>, ScalarFuncTypeLhs>, DirectExpr<Gradient<N_SIZE> > > >(u, v); \
    }
    
#define GRADIENT_DEFINE_BINARY_FUNCTION_CONST_ARG_RHS(ExpressionName, FunctionName, FunctionClass, ScalarFuncTypeRhs) \
    template <typename LhsExpr>                                         \
    inline GradientExpression<ExpressionName<FunctionClass, LhsExpr, ConstExpr<typename LhsExpr::GradientType, ScalarFuncTypeRhs> > > \
    FunctionName(const GradientExpression<LhsExpr>& u, ScalarFuncTypeRhs v) { \
        return GradientExpression<ExpressionName<FunctionClass, LhsExpr, ConstExpr<typename LhsExpr::GradientType, ScalarFuncTypeRhs> > >(u, v); \
    }                                                                   \
                                                                        \
    template <index_type N_SIZE>                                        \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, ConstExpr<Gradient<N_SIZE>, ScalarFuncTypeRhs> > > \
    FunctionName(const Gradient<N_SIZE>& u, ScalarFuncTypeRhs v) {       \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, ConstExpr<Gradient<N_SIZE>, ScalarFuncTypeRhs> > >(u, v); \
    }                                                                   \

#define GRADIENT_DEFINE_BINARY_FUNCTION(ExpressionName, FunctionName, FunctionClass) \
    GRADIENT_DEFINE_BINARY_FUNCTION_CONST_ARG_LHS(ExpressionName, FunctionName, FunctionClass, scalar_func_type) \
    GRADIENT_DEFINE_BINARY_FUNCTION_CONST_ARG_RHS(ExpressionName, FunctionName, FunctionClass, scalar_func_type) \
                                                                        \
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
    template <index_type N_SIZE, typename RhsExpr> \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, RhsExpr> > \
    FunctionName(const Gradient<N_SIZE>& u, const GradientExpression<RhsExpr>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, RhsExpr> >(u, v); \
    } \
                                                                        \
    template <index_type N_SIZE> \
    inline GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > \
    FunctionName(const Gradient<N_SIZE>& u, const Gradient<N_SIZE>& v) { \
        return GradientExpression<ExpressionName<FunctionClass, DirectExpr<Gradient<N_SIZE> >, DirectExpr<Gradient<N_SIZE> > > > (u, v); \
    } \

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
    GRADIENT_DEFINE_BINARY_FUNCTION(BinaryExpr, fmod, FuncFmod)
    
    GRADIENT_DEFINE_BINARY_FUNCTION_CONST_ARG_RHS(BinaryExpr, pow, FuncPowInt, integer)

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
