/*
 * MBDyn (C) is a multibody analysis code.
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

#include <loadable.h>
#include <dataman.h>
#include <strnode.h>

struct spring {
   StructNode* pNode1;
   StructNode* pNode2;
   integer c;
   doublereal k;
};

/* funzioni di default */
void* read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
   DEBUGCOUT("Spring Elem: " << __PRETTY_FUNCTION__ << endl);

   spring* s = NULL;
   SAFENEW(s, spring);
   
   unsigned int ul = HP.GetInt();
   DEBUGCOUT("Linked to Node " << ul << endl);
       
   /* verifica di esistenza del nodo */  
   if ((s->pNode1 = pDM->pFindStructNode(ul)) == NULL) {
      cerr << endl << __PRETTY_FUNCTION__
	<< " at line " << HP.GetLineData() 
	<< ": structural node " << ul
	<< " not defined" << endl;
      
      THROW(DataManager::ErrGeneric());
   }                  
   
   ul = HP.GetInt();
   DEBUGCOUT("Linked to Node " << ul << endl);
       
   /* verifica di esistenza del nodo */  
   if ((s->pNode2 = pDM->pFindStructNode(ul)) == NULL) {
      cerr << endl << __PRETTY_FUNCTION__
	<< " at line " << HP.GetLineData() 
	<< ": structural node " << ul
	<< " not defined" << endl;
      
      THROW(DataManager::ErrGeneric());
   }  

   s->c = HP.GetInt();
   if (s->c < 1 || s->c > 3) {
      cerr << __PRETTY_FUNCTION__ << "illegal component" << endl;
      THROW(DataManager::ErrGeneric());
   }

   s->k = HP.GetReal();
   
   DEBUGCOUT("Component: " << s->c << ", stiffness: " << s->k << endl);
   
   return (void *)s;
}

void work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUT("Spring Elem: " << __PRETTY_FUNCTION__ << endl);
   *piNumRows = 2;
   *piNumCols = 2;
}

VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Spring Elem: " << __PRETTY_FUNCTION__ << endl);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   spring* s = (spring*)pEl->pGetData();
   
   WM.fPutRowIndex(1, s->pNode1->iGetFirstRowIndex()+s->c);
   WM.fPutRowIndex(2, s->pNode2->iGetFirstRowIndex()+s->c);
   WM.fPutColIndex(1, s->pNode1->iGetFirstColIndex()+s->c);
   WM.fPutColIndex(2, s->pNode2->iGetFirstColIndex()+s->c);
   
   doublereal k = s->k*dCoef;
   WM.fPutCoef(1, 1, k);
   WM.fPutCoef(1, 2, -k);
   WM.fPutCoef(2, 1, -k);
   WM.fPutCoef(2, 2, k);   
   
   return WorkMat;
}

SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Spring Elem: " << __PRETTY_FUNCTION__ << endl);
   WorkVec.Resize(2);
   
   spring* s = (spring*)pEl->pGetData();
   
   WorkVec.fPutRowIndex(1, s->pNode1->iGetFirstRowIndex()+s->c);
   WorkVec.fPutRowIndex(2, s->pNode2->iGetFirstRowIndex()+s->c);
   
   doublereal f = s->k*(s->pNode2->GetXCurr().dGet(s->c)
			-s->pNode1->GetXCurr().dGet(s->c));
   
   WorkVec.fPutCoef(1, f);
   WorkVec.fPutCoef(2, -f);
   
   return WorkVec;
}

void destroy(LoadableElem* pEl)
{
   DEBUGCOUT("Spring Elem: " << __PRETTY_FUNCTION__ << endl);
   
   spring* s = (spring*)pEl->pGetData();
   SAFEDELETE(s);
}
