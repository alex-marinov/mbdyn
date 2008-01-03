/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifndef STACK_H
#define STACK_H

template <class T> class Stack {
 private: 
   struct Tlist {
      T t;
      Tlist* next;
      
      Tlist(const T& t) : t(t), next(NULL) {
	 NO_OP;
      };	 
   };   
   
   Tlist* p;
   
 public:
   Stack(void) : p(NULL) {
      NO_OP;
   }
   
   ~Stack(void) {
      while (p != NULL) {
	 Tlist* pp = p;
	 p = p->next;
	 SAFEDELETE(pp);
      }
   };
   
   void Push(const T& t) {
      Tlist* pp = NULL;
      SAFENEWWITHCONSTRUCTOR(pp, Tlist, Tlist(t));
      pp->next = p;
      p = pp;
   }
   
   int Pop(T& t) {
      if (p == NULL) {
	 return 0;
      }
      t = p->t;
      Tlist* pp = p;
      p = p->next;
      SAFEDELETE(pp);
      return 1;
   };
};

#endif
