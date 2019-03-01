/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef HANDLE_H
#define HANDLE_H

#include "bool.h"

#include "myassert.h"
#include "refcount.h"
#include "destr.h"

template <class T>
class Handle : protected RefCount {
 protected:
   T* p;
   const Destructor<T>& d;

 private:
   Handle(void);
   
 public:
   Handle(const T* pnt, const Destructor<T>& dst)
     : RefCount(), p((T*)pnt), d(dst) { 
	NO_OP; 
     };
   
   Handle(const T& t, const Destructor<T>& dst)
     : RefCount(), p(&(T&)t), d(dst) { 
	NO_OP;
     };
   
   Handle(const Handle<T>& h) : RefCount(h), p(h.p), d(h.d) { 
      NO_OP; 
   };
   
   // Binding functions
   Handle<T>& Bind(const T* pnt, const Destructor<T>& dst) {
      // warning: this function fails if the handle is rebound to itself
      if (RefCount::GetCount() == 1) {
	 if (p != NULL) {	  
	    d.Destroy(p);
	 }
      }
      
      RefCount::New();     
      (Destructor<T>&)d = (Destructor<T>&)dst;
      p = (T*)pnt;
      
      return *this;           
   };
   
   Handle<T>& Bind(const T& t, const Destructor<T>& dst) {
      // warning: this function fails if the handle is rebound to itself
      if (RefCount::GetCount() == 1) {
	 if (p != NULL) {	   
	    d.Destroy(p);
	 }
      }
      
      RefCount::New();     
      (Destructor<T>&)d = dst;
      p = &(T&)t;
      
      return *this;           
   };
   
   Handle<T>& Bind(const Handle<T>& h) {
      // Create temporaries
      RefCount r = this->CreateTmp(*this);
      const Destructor<T>& dTmp = d;
      T* pTmp = p;  // your data
 
      // Assign 
      Assign(h);
      (Destructor<T>&)d = h.d;
      p = h.p; // your data
      
      // Check of temporary
      if (r.GetCount() == 1) { // r is going to be deleted by its destructor
	 if (pTmp != NULL) {
	    dTmp.Destroy(pTmp);
	 }
      }
      
      return *this;     
   };

   // Binding operators
   Handle<T>& operator = (const Handle<T>& h) {
      return Bind(h);
   };   
   
   // Utility
   inline bool IsValid(void) const { 
      return (p == NULL) ? FALSE : TRUE; 
   };
   
   inline T* GetPtr(void) { 
      return p;
   };
   
   inline T& GetRef(void) { 
      ASSERT(p != NULL); 
      return *p;
   };
   
   inline const T* GetPtr(void) const { 
      return p; 
   };
   
   inline const T& GetRef(void) const { 
      ASSERT(p != NULL); return *p; 
   };

   void Destroy(void) {
      if (GetCount() == 1) {
	 d.Destroy(p);
      }
   };   
   
   
   // Destructor.
   // Be sure the Destructor<T>& has been properly initialised, namely:
   // - LinkDestructor<T> if p refers to system allocated memory,
   // - HardDestructor<T> if p refers to user allocated memory,
   ~Handle(void) {
      Destroy();    
   };   
   
   // Referencing operators
   inline T* operator -> (void) { 
      ASSERT(p != NULL); 
      return p;
   };
   
   inline T& operator * (void) { 
      ASSERT(p != NULL); 
      return *p;
   };
   
   inline const T* operator -> (void) const { 
      ASSERT(p != NULL); 
      return p;
   };
   
   inline const T& operator * (void) const { 
      ASSERT(p != NULL); 
      return *p;
   };
};

#endif // HANDLE_H
