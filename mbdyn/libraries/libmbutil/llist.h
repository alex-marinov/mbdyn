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

/* Linked list generica
 * Il tipo <T> deve essere derivato da <WithLabel> 
 * e deve avere il costruttore nullo T(0) e l'operatore di assegnazione "="
 * 
 * esempio:
 * class MyRand : public WithLabel {
 *  public:
 *    MyRand(unsigned int uL) : WithLabel(uL) { NO_OP; };
 *    integer iRand;
 * };
 * HARDDESTRUCTOR(D, MyRand, mmDMY);
 * MyLList<MyRand> LL(D);
 * MyRand* pmrIn = NULL;
 * SAFENEWWITHCONSTRUCTOR(pmrIn, MyRand, MyRand(0));
 * if(pmrIn == NULL) {
 *    // error ...
 * }
 * if(!LL.iAdd(pmrIn)) {
 *    // error ...
 * }
 * MyRand* pmrOut = NULL;
 * if(!LL.iGetFirst(pmrOut)) {
 *    // error ...
 * }
 * if(pmrOut->GetLabel() != pmrIn->GetLabel()) {
 *    // error ...
 * }
 * 
 * // and so on ...
 */

#ifndef LLIST_H
#define LLIST_H

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "destr.h"

template <class T>
class MyLList {
 private:
   Destructor<T>& D;
   
   struct LLItem {
      const LLItem* pNext;
      const T* pTItem;
      LLItem(const LLItem* p, const T* pTIn) 
	: pNext(p), pTItem(pTIn) { 
	   NO_OP; 
	};
   };      
   
   LLItem* pBase;
   LLItem* pCurr;
   
   long int iSize;
   
   LLItem* pFindPrev(unsigned int uLabel) const
     {       
	// ASSERT(uLabel > 0);
	ASSERT(pBase != NULL);
	LLItem* p = pBase;
	do {	
	   if ((p->pNext->pTItem == NULL)
	       || (p->pNext->pTItem->GetLabel() >= uLabel)) {
	      return p;
	   }
	   p = (LLItem*)(p->pNext);
	} while (p != pBase);
	return pBase;
     };
   
   
   LLItem* pFind(unsigned int uLabel) const
     {
	// ASSERT(uLabel > 0);
	LLItem* p = pFindPrev(uLabel);
	if ((p != NULL) && (p->pNext->pTItem) != NULL 
	    && (p->pNext->pTItem->GetLabel() == uLabel)) {
	   return (LLItem*)(p->pNext);
	}
	return NULL;
     };
   
 public:
   /* Costruzione/distruzione */
   MyLList(Destructor<T>& d);
   ~MyLList(void);
   
   /* Gestione lista */
   int iAdd(const T* pTIn);
   int iRemove(unsigned int uLabel);
      
   /* Ricerca per chiave */
   const T* Get(unsigned int uLabel) const;

   /* Iteratore built-in */
   int iGetFirst(T*& pTOut) const;
   int iGetNext(T*& pTOut) const;
   
   /* Numero di elementi nella lista */
   long int iGetSize(void) const { return iSize; };
};

template <class T>
MyLList<T>::MyLList(Destructor<T>& d)
: D(d), pBase(NULL), pCurr(NULL), iSize(0)
{
   SAFENEWWITHCONSTRUCTOR(pBase, LLItem, LLItem(NULL, NULL));
   (LLItem*&)(pBase->pNext) = pBase;
   
   pCurr = pBase;
}


template <class T>
MyLList<T>::~MyLList(void)
{
   ASSERT(pBase != NULL);
   LLItem* p = pBase;
   LLItem* pd = NULL;
   do {           
      pd = p;
      p = (LLItem*)(p->pNext);
      
      D.Destroy((T*&)(pd->pTItem));
      SAFEDELETE(pd);
   } while (p != pBase);
}


/*
template <class T>
MyLList<T>::LLItem* MyLList<T>::pFindPrev(unsigned int uLabel)
{
   // ASSERT(uLabel > 0);
   ASSERT(pBase != NULL);
   LLItem* p = pBase;
   while (p->pNext != pBase) {
      if (p->pNext->TItem.GetLabel() >= uLabel) {
	 return p;
      }
      p = p->pNext;
   }
   return NULL;
}


template <class T>
MyLList<T>::LLItem* MyLList<T>::pFind(unsigned int uLabel)
{
   // ASSERT(uLabel > 0);
   LLItem* p = pFindPrev(uLabel);
   if (p != NULL && p->pNext.GetLabel() == uLabel) {
      return p->pNext;
   }
   return NULL;
}
 */


template <class T>
int MyLList<T>::iAdd(const T* pTIn)
{
   ASSERT(pTIn != NULL);
      
   LLItem* p = pFindPrev(pTIn->GetLabel());
   if ((p->pNext->pTItem != NULL) 
      && (p->pNext->pTItem->GetLabel() == pTIn->GetLabel())) {
      return 1;
   }
   
   LLItem* pn = NULL;
   SAFENEWWITHCONSTRUCTOR(pn, LLItem, LLItem(p->pNext, pTIn));
   (LLItem*&)(p->pNext) = pn;
   iSize++;
   return 0;
}


template <class T>
int MyLList<T>::iRemove(unsigned int uLabel)
{
   // ASSERT(uLabel > 0);
   LLItem* p = pFindPrev(uLabel);
   if ((p != NULL) && (p->pNext->pTItem != NULL) 
      && (p->pNext->pTItem->GetLabel() == uLabel)) {
      LLItem* pd = (LLItem*)(p->pNext);
      (LLItem*&)(p->pNext) = (LLItem*)(pd->pNext);
      
      D.Destroy((T*&)(pd->pTItem));
      SAFEDELETE(pd);
      
      ASSERT(iSize > 0);
      iSize--;
      
      return 1;
   }
   return 0;
}

  
template <class T>
const T* MyLList<T>::Get(unsigned int uLabel) const
{
   LLItem* p = pFind(uLabel);
   if (p != NULL) {
      return p->pTItem;
   }
   return pBase->pTItem;
}


template <class T>
int MyLList<T>::iGetFirst(T*& pTOut) const
{
   if (pBase->pNext == pBase) {
      return 0;
   }
   (LLItem*&)pCurr = (LLItem*)pBase->pNext;
   pTOut = (T*)(pCurr->pTItem);
   return 1;
}


template <class T>
int MyLList<T>::iGetNext(T*& pTOut) const
{
   (LLItem*&)pCurr = (LLItem*)(pCurr->pNext);
   if (pCurr == pBase) {
      return 0;
   }
   pTOut = (T*)(pCurr->pTItem);
   return 1;
}

#endif
