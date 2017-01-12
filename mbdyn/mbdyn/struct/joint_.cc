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

#include <joint_.h>

/* Joint_1Node - begin */

SubVectorHandler& 
Joint_1Node::AssRes(SubVectorHandler& WorkVec,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   integer iNodeFirstRowIndex = pNode->iGetFirstRowIndex();
   for (int i = 1; i <= 6; i++) {
      WorkVec.PutRowIndex(i, iNodeFirstRowIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WorkVec.PutRowIndex(6+i, iFirstIndex+i);
   }

   AssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_1Node::AssJac(VariableSubVectorHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);
   
   integer iNodeFirstRowIndex = pNode->iGetFirstRowIndex();
   integer iNodeFirstColIndex = pNode->iGetFirstColIndex();
   for (int i = 1; i <= 6; i++) {
      WM.PutRowIndex(i, iNodeFirstRowIndex+i);
      WM.PutColIndex(i, iNodeFirstColIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WM.PutRowIndex(6+i, iFirstIndex+i);
      WM.PutColIndex(6+i, iFirstIndex+i);
   }

   AssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}

SubVectorHandler& 
Joint_1Node::InitialAssRes(SubVectorHandler& WorkVec,
			   const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   integer iNodeFirstIndex = pNode->iGetFirstIndex();   
   for (int i = 1; i <= 12; i++) {
      WorkVec.PutRowIndex(i, iNodeFirstIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WorkVec.PutRowIndex(12+i, iFirstIndex+i);
   }

   InitialAssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_1Node::InitialAssJac(VariableSubVectorHandler& WorkMat,
			   const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);
   
   integer iNodeFirstIndex = pNode->iGetFirstIndex();  
   for (int i = 1; i <= 12; i++) {
      WM.PutRowIndex(i, iNodeFirstIndex+i);
      WM.PutColIndex(i, iNodeFirstIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WM.PutRowIndex(12+i, iFirstIndex+i);
      WM.PutColIndex(12+i, iFirstIndex+i);
   }

   InitialAssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}

/* Joint_1Node - end */


/* Joint_2Nodes - begin */

SubVectorHandler& 
Joint_2Nodes::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   integer iNode1FirstRowIndex = pNode1->iGetFirstRowIndex();
   integer iNode2FirstRowIndex = pNode2->iGetFirstRowIndex();
   for (int i = 1; i <= 6; i++) {
      WorkVec.PutRowIndex(i, iNode1FirstRowIndex+i);
      WorkVec.PutRowIndex(6+i, iNode2FirstRowIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WorkVec.PutRowIndex(12+i, iFirstIndex+i);
   }

   AssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_2Nodes::AssJac(VariableSubVectorHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);
   
   integer iNode1FirstRowIndex = pNode1->iGetFirstRowIndex();
   integer iNode1FirstColIndex = pNode1->iGetFirstColIndex();
   integer iNode2FirstRowIndex = pNode2->iGetFirstRowIndex();
   integer iNode2FirstColIndex = pNode2->iGetFirstColIndex();
   for (int i = 1; i <= 6; i++) {
      WM.PutRowIndex(i, iNode1FirstRowIndex+i);
      WM.PutColIndex(i, iNode1FirstColIndex+i);
      WM.PutRowIndex(6+i, iNode2FirstRowIndex+i);
      WM.PutColIndex(6+i, iNode2FirstColIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WM.PutRowIndex(12+i, iFirstIndex+i);
      WM.PutColIndex(12+i, iFirstIndex+i);
   }

   AssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}


SubVectorHandler& 
Joint_2Nodes::InitialAssRes(SubVectorHandler& WorkVec,
			    const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   integer iNode1FirstIndex = pNode1->iGetFirstIndex();
   integer iNode2FirstIndex = pNode2->iGetFirstIndex();
   for (int i = 1; i <= 6; i++) {
      WorkVec.PutRowIndex(i, iNode1FirstIndex+i);
      WorkVec.PutRowIndex(12+i, iNode2FirstIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WorkVec.PutRowIndex(24+i, iFirstIndex+i);
   }

   InitialAssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_2Nodes::AssJac(VariableSubVectorHandler& WorkMat,
		    const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);
   
   integer iNode1FirstIndex = pNode1->iGetFirstIndex();
   integer iNode2FirstIndex = pNode2->iGetFirstIndex();
   for (int i = 1; i <= 6; i++) {
      WM.PutRowIndex(i, iNode1FirstIndex+i);
      WM.PutColIndex(i, iNode1FirstIndex+i);
      WM.PutRowIndex(12+i, iNode2FirstIndex+i);
      WM.PutColIndex(12+i, iNode2FirstIndex+i);
   }
   
   integer iFirstIndex = iGetFirstIndex();
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WM.PutRowIndex(24+i, iFirstIndex+i);
      WM.PutColIndex(24+i, iFirstIndex+i);
   }

   InitialAssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}

/* Joint_2Nodes - end */


/* Joint_NNodes - begin */

SubVectorHandler& 
Joint_NNodes::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   int iNNod = iGetNNodes();
   for (int j = 0; j < iNNod; j++) {
      ASSERT(pNodes[j] != NULL);
      ASSERT(pNodes[j]->GetNodeType() == Node::STRUCTURAL);
      integer iNodeFirstRowIndex = pNodes[j]->iGetFirstRowIndex();    
   
      for (int i = 1; i <= 6; i++) {
	 WorkVec.PutRowIndex(6*j+i, iNodeFirstRowIndex+i);
      }
   }
   
   integer iFirstIndex = iGetFirstIndex();
   integer iOff = 6*iNNod;
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WorkVec.PutRowIndex(iOff+i, iFirstIndex+i);
   }

   AssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_NNodes::AssJac(VariableSubVectorHandler& WorkMat,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);

   int iNNod = iGetNNodes();
   for (int j = 0; j < iNNod; j++) {
      integer iNodeFirstRowIndex = pNodes[j]->iGetFirstRowIndex();
      integer iNodeFirstColIndex = pNodes[j]->iGetFirstColIndex();
   
      for (int i = 1; i <= 6; i++) {
	 WM.PutRowIndex(6*j+i, iNodeFirstRowIndex+i);
	 WM.PutColIndex(6*j+i, iNodeFirstColIndex+i);
      }
   }
   
   integer iFirstIndex = iGetFirstIndex();
   integer iOff = 6*iNNod;
   for (int i = 1; i <= iGetNumDofs(); i++) {
      WM.PutRowIndex(iOff+i, iFirstIndex+i);
      WM.PutColIndex(iOff+i, iFirstIndex+i);
   }

   AssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}


SubVectorHandler& 
Joint_NNodes::InitialAssRes(SubVectorHandler& WorkVec,
			    const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   
   int iNNod = iGetNNodes();
   for (int j = 0; j < iNNod; j++) {
      ASSERT(pNodes[j] != NULL);
      ASSERT(pNodes[j]->GetNodeType() == Node::STRUCTURAL);
      integer iNodeFirstIndex = pNodes[j]->iGetFirstIndex();
   
      for (int i = 1; i <= 12; i++) {
	 WorkVec.PutRowIndex(12*j+i, iNodeFirstIndex+i);
      }
   }
   
   integer iFirstIndex = iGetFirstIndex();
   integer iOff = 12*iNNod;
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WorkVec.PutRowIndex(iOff+i, iFirstIndex+i);
   }

   InitialAssRes_(WorkVec, dCoef, XCurr, XPrimeCurr);
   return WorkVec;
}


VariableSubVectorHandler& 
Joint_NNodes::InitialAssJac(VariableSubVectorHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WorkMat.Resize(iNumRows, iNumCols);

   int iNNod = iGetNNodes();
   for (int j = 0; j < iNNod; j++) {
      integer iNodeFirstIndex = pNodes[j]->iGetFirstIndex();  
   
      for (int i = 1; i <= 12; i++) {
	 WM.PutRowIndex(12*j+i, iNodeFirstIndex+i);
	 WM.PutColIndex(12*j+i, iNodeFirstIndex+i);
      }
   }
   
   integer iFirstIndex = iGetFirstIndex();
   integer iOff = 12*iNNod;
   for (int i = 1; i <= iGetInitialNumDofs(); i++) {
      WM.PutRowIndex(iOff+i, iFirstIndex+i);
      WM.PutColIndex(iOff+i, iFirstIndex+i);
   }

   InitialAssJac_(WM, dCoef, XCurr, XPrimeCurr);
   return WorkMat;
}

/* Joint_NNodes - end */
