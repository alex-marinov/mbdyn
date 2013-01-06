/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef DESTR_H
#define DESTR_H

#include "myassert.h"
#include "mynewmem.h"

template <class T>
class Destructor {
public:
   	Destructor(void) { 
      		NO_OP; 
   	};
   
   	virtual ~Destructor(void) { 
      		NO_OP;
   	};
   
   	virtual void Destroy(T*&) const = 0;
};


template <class T>
class LinkDestructor : public Destructor<T> {
public:
   	LinkDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};
   
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
   	HardDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};
   
   	~HardDestructor(void) { 
      		NO_OP;
   	};
   
   	void Destroy(T*& pnt) const {
      		if (pnt != NULL) {
	 		SAFEDELETE(pnt);
	 		pnt = NULL;
      		}
   	};   
};

template <class T>
class ArrayHardDestructor : public Destructor<T> {
public:
   	ArrayHardDestructor(void) : Destructor<T>() { 
      		NO_OP;
   	};

   	~ArrayHardDestructor(void) { 
      		NO_OP;
   	};
   
   	void Destroy(T*& pnt) const {
      		if (pnt != NULL) {
	 		SAFEDELETEARR(pnt);
	 		pnt = NULL;
      		}
   	};
};

#define LINKDESTRUCTOR(d, type, m) LinkDestructor<type> d
#define HARDDESTRUCTOR(d, type, m) HardDestructor<type> d
#define ARRAYHARDDESTRUCTOR(d, type, m) ArrayHardDestructor<type> d

#endif /* DESTR_H */

