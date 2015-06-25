/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/*****************************************************************************
 * 
 * Reference Counter
 * 
 * Used for shared objects. A class derived from this automatically
 * keeps a count of the number of references to a common resource
 * and switches from one to another when the object is assigned
 * to a new shared resource.
 * 
 * Suggested use:
 * 
 * class MyClass : protected RefCount {
 *    // ...
 *  protected: 
 *    int* pi;
 *  public:
 *    MyClass(unsigned int iSize) : RefCount(), pi(NULL) 
 *      {
 * 	   if (iSize > 0) {
 *	      SAFENEWARR(pi, int, iSize);
 *	      if(pi == NULL) { NULL; } // ....
 *	   }
 *      };
 *    MyClass(const MyClass& m) : RefCount(m), pi(m.pi) { NULL; }
 *    ~MyClass(void)
 *      {
 *	   if(GetCount() == 1) {
 *	      SAFEDELETEARR(pi);
 *	   }
 *      };
 *    MyClass& operator = (const MyClass& m)
 *      {  
 *         // Create temporaries
 *         RefCount r = this->CreateTmp(*this);
 *         int* piTmp = pi;  // your data
 *
 *         // Assign 
 *         Assign(m);
 *         pi = m.pi; // your data
 *
 *         // Check of temporary
 *         if (r.GetCount() == 1) { // r is going to be deleted by destructor
 *            SAFEDELETEARR(piTmp);
 *         }
 *
 *         return *this;
 *      };
 *    //...
 * };
 * 
 ****************************************************************************/

#ifndef REFCOUNT_H
#define REFCOUNT_H

#include "mynewmem.h"
#include "except.h"

class RefCount {
 protected:
   unsigned int* piCount;
   
 protected:
   RefCount(void) : piCount(NULL) {
      SAFENEW(piCount, unsigned int);
      if (piCount == NULL) {
	 cerr << "Out of memory" << endl;
	 throw ErrMemory(MBDYN_EXCEPT_ARGS);
      }
      *piCount = 1;
   };
   
   RefCount(const RefCount& r) : piCount(r.piCount)
      {
	 ASSERT(r.piCount != NULL);
	 ASSERT(*r.piCount > 0);
	 (*piCount)++;
      };
   
   RefCount& operator = (const RefCount& r) {
      Assign(r);
      return *this;
   };
   
 public:
   virtual ~RefCount(void) {
      Decrease();
   };
   
 protected:
   int Decrease(void) {
      ASSERT(piCount != NULL);
      ASSERT(*piCount > 0);
      if (--(*piCount) == 0) {
	 SAFEDELETE(piCount);
	 return 1;
      }
      return 0;
   };   
   
   int Assign(const RefCount& r) {
      ASSERT(r.piCount != NULL);
      ASSERT(*r.piCount > 0);
      int i = Decrease();
      piCount = r.piCount;
      (*piCount)++;
      return i;
   };
   
 public:
   unsigned int GetCount(void) const { 
      ASSERT(piCount != NULL);
      ASSERT(*piCount > 0);
      return *piCount; 
   };
   
   RefCount CreateTmp(const RefCount& r) const {
      return RefCount(r);
   };
   
   RefCount& New(void) {
      if (GetCount() == 1) {
	 return *this;
      }
      
      Decrease();
      piCount = NULL;
      SAFENEW(piCount, unsigned int);
      if (piCount == NULL) {
	 cerr << "Out of memory" << endl;
	 throw ErrMemory(MBDYN_EXCEPT_ARGS);
      }
      *piCount = 1;
      return *this;
   }
};

#endif // REFCOUNT_H
