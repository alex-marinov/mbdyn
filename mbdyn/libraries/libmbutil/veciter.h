/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

/* Iteratore per vettori */


#ifndef VECITER_H
#define VECITER_H


#include "myassert.h"

/* GetFirst ritorna 1 ed assegna a TReturn il primo item del vettore. 
 * Si assume che il vettore contenga almeno un termine;
 * GetNext ritorna 1 ed assegna a TReturn il termine successivo se esiste,
 * altrimenti ritorna 0 e l'assegnamento a TReturn e' unpredictable. */


template <class T> 
class Iter {
 public:
   virtual ~Iter(void) { NO_OP; };
   virtual inline flag fGetFirst(T& TReturn) const = 0;
   virtual inline flag fGetNext(T& TReturn) const = 0;
};

template<class T>
class VecIter : public Iter<T> {
 private:
   const T* pStart;
   T* pCount;
   long int iSize;
   
 public:
   VecIter(void) : pStart(NULL), pCount(NULL), iSize(0) { NO_OP; };
   VecIter(const T* p, long int i) : pStart(p), pCount((T*)p), iSize(i) 
     { 
	ASSERT(pStart != NULL);
	ASSERT(iSize > 0);
     };
   
   virtual ~VecIter(void) { NO_OP; };
   
   void Init(const T* p, long int i)
     {
	ASSERT(p != NULL);
	ASSERT(i > 0);
	
	pStart = pCount = (T*)p;
	iSize = i;
     };
   
   inline flag fGetFirst(T& TReturn) const
     {
	ASSERT(pStart != NULL);
	ASSERT(iSize > 0);
	
	TReturn = *((T*&)pCount = (T*)pStart);
	
	return flag(1);
     };
   
   inline flag fGetNext(T& TReturn) const
     {
	ASSERT(pStart != NULL);     
	ASSERT(iSize > 0);
	ASSERT(pCount >= (T*)pStart);
	
	if (++((T*&)pCount) == (T*)pStart+iSize)
	  return flag(0);
	
	TReturn = *pCount;
	
	return flag(1);
     };        
};

#endif
