/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifndef DESTR_H
#define DESTR_H

#include "myassert.h"
#include "mynewmem.h"

template <class T>
class Destructor {
protected:   
#if defined(DEBUG_MEMMANAGER)
   	clMemMan& mm;
#endif   
   
public:
#if defined(DEBUG_MEMMANAGER)
   	Destructor(clMemMan& m) : mm(m) { 
      		NO_OP;
   	};
#else /* !DEBUG_MEMMANAGER */
   	Destructor(void) { 
      		NO_OP; 
   	};
#endif /* !DEBUG_MEMMANAGER */
   
   	virtual ~Destructor(void) { 
      		NO_OP;
   	};
   
   	virtual void Destroy(T*&) const = 0;
   
#if defined(DEBUG_MEMMANAGER)
   	clMemMan& GetMemMan(void) const {
      		return mm;
   	};
#endif /* DEBUG_MEMMANAGER */
};


template <class T>
class LinkDestructor : public Destructor<T> {
public:
#if defined(DEBUG_MEMMANAGER)
   	LinkDestructor(clMemMan& m) : Destructor<T>(m) { 
      		NO_OP;
   	};
#else /* !DEBUG_MEMMANAGER */
   	LinkDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};
#endif /* !DEBUG_MEMMANAGER */
   
   	~LinkDestructor(void) { 
      		NO_OP;
   	};
   
   	void Destroy(T*& pnt) const { 
      		NO_OP;
   	};
};


template <class T>
class HardDestructor : public Destructor<T> {
public:
#if defined(DEBUG_MEMMANAGER)
   	HardDestructor(clMemMan& m) : Destructor<T>(m) { 
      		NO_OP;
   	};
#else /* !DEBUG_MEMMANAGER */
   	HardDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};
#endif /* !DEBUG_MEMMANAGER */
   
   	~HardDestructor(void) { 
      		NO_OP;
   	};
   
   	void Destroy(T*& pnt) const {
      		if (pnt != NULL) {
	 		SAFEDELETE(pnt, mm);
	 		pnt = NULL;
      		}
   	};   
};

template <class T>
class ArrayHardDestructor : public Destructor<T> {
public:
#if defined(DEBUG_MEMMANAGER)
   	ArrayHardDestructor(clMemMan& m) : Destructor<T>(m) { 
      		NO_OP;
   	};
#else /* !DEBUG_MEMMANAGER */
   	ArrayHardDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};
#endif /* !DEBUG_MEMMANAGER */

   	~ArrayHardDestructor(void) { 
      		NO_OP;
   	};
   
   	void Destroy(T*& pnt) const {
      		if (pnt != NULL) {
	 		SAFEDELETEARR(pnt, mm);
	 		pnt = NULL;
      		}
   	};
};


#if defined(DEBUG_MEMMANAGER)

#define LINKDESTRUCTOR(d, type, m) LinkDestructor<type> d((m))
#define HARDDESTRUCTOR(d, type, m) HardDestructor<type> d((m))
#define ARRAYHARDDESTRUCTOR(d, type, m) ArrayHardDestructor<type> d((m))

#else /* !DEBUG_MEMMANAGER */

#define LINKDESTRUCTOR(d, type, m) LinkDestructor<type> d
#define HARDDESTRUCTOR(d, type, m) HardDestructor<type> d
#define ARRAYHARDDESTRUCTOR(d, type, m) ArrayHardDestructor<type> d

#endif /* !DEBUG_MEMMANAGER */

#endif /* DESTR_H */

