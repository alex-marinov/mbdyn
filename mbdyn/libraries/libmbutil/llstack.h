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

/* Stack realizzata mediante linked list */

/* L'oggetto LLStack contiene una linked list di strutture
 * LLItem, che a loro volta contengono un puntatore ad un oggetto di tipo T
 * ed un puntatore all'item successivo.
 * Quando la struttura viene creata, viene allocato il primo elemento 
 * della lista, che cosi' come successivo ha NULL e contiene un puntatore
 * a NULL.
 * Gli elementi che vengono aggiunti puntano al primo della fila e diventano 
 * cosi' il top della stack.
 * Durante lo svuotmento, l'ultimo elemento si riconosce dal fatto che come 
 * successivo ha NULL.
 */

#ifndef LLSTACK_H
#define LLSTACK_H

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

template <class T>
class LLStack {
 private:
   struct LLItem {
      const T* pTItem;
      const LLItem* pnext;
      LLItem(const T* pT = NULL, const LLItem* pn = NULL)
	: pTItem(pT), pnext(pn) { NO_OP; };
   };
   
   LLItem* pTop;
   
 public:   
   LLStack(void) : pTop(NULL) {
      SAFENEW(pTop, LLItem);
   };
   
   ~LLStack(void) {
      ASSERT(pTop != NULL);
      LLItem* p = pTop;
      while(p != NULL) {
	 pTop = (LLItem*)(pTop->pnext);
	 SAFEDELETE(p);
	 p = pTop;
      }	
   };
   
   void Push(const T* pTIn) { 
      ASSERT(pTop != NULL);
      ASSERT(pTIn != NULL); /* Non necessario, solo per debug */
      LLItem* p = NULL;
      SAFENEWWITHCONSTRUCTOR(p, LLItem, LLItem(pTIn, pTop));
      pTop = p;
   };
   
   int iPop(T*& pTOut) {
      ASSERT(pTop != NULL);
      
      /* Empty stack ? */
      if(iIsEmpty()) {
#ifdef DEBUG
	 pTOut = NULL;
#endif	   
	 return 0;
      }
      
      /* No, still something to pop */
      pTOut = (T*)(pTop->pTItem);
      LLItem* p = pTop;
      pTop = (LLItem*)(pTop->pnext);
      SAFEDELETE(p);
      return 1;
   };
   
   inline int iIsEmpty(void) const { 
      return int(pTop->pnext == NULL);
   };
};

#endif
