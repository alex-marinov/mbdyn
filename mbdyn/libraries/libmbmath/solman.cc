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

/* solution manager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <solman.h>
#include <submat.h>
#include <matvec3.h>

/* VectorHandler - begin */

/* Somma un Vec3 nella posizione desiderata */
void VectorHandler::Add(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif
   
   fIncCoef(iRow, v.dGet(1));
   fIncCoef(++iRow, v.dGet(2));
   fIncCoef(++iRow, v.dGet(3));
}
   
/* Sottrae un Vec3 nella posizione desiderata */
void VectorHandler::Sub(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif
   
   fDecCoef(iRow, v.dGet(1));
   fDecCoef(++iRow, v.dGet(2));
   fDecCoef(++iRow, v.dGet(3));
}

/* Scrive un Vec3 nella posizione desiderata */
void VectorHandler::Put(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif
   
   fPutCoef(iRow, v.dGet(1));
   fPutCoef(++iRow, v.dGet(2));
   fPutCoef(++iRow, v.dGet(3));
}

/* Somma e moltiplica per uno scalare */
VectorHandler& VectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d) {      
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif
   
   for (integer i = iGetSize(); i > 0; ) {
      fIncCoef(i, d*VH.dGetCoef(i));
      i--;
   }
   
   return *this;
}

/* Moltiplica per uno scalare */
VectorHandler& VectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif
   
   for (integer i = iGetSize(); i > 0; ) {
      fPutCoef(i, d*VH.dGetCoef(i));
      i--;
   }
   
   return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler& VectorHandler::operator += (const VectorHandler& VH) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif
   
   for (integer i = iGetSize(); i > 0; ) {
      fIncCoef(i, VH.dGetCoef(i));
      i--;
   }
   
   return *this;
}
   
/* Overload di += usato per l'assemblaggio del residuo */
VectorHandler& VectorHandler::operator += (const SubVectorHandler& SubVH)
{
#ifdef DEBUG
   SubVH.IsValid();
#endif

   if (SubVH.iGetSize() == 0) {
      return *this;
   }

   return SubVH.AddTo(*this);
}

/* Overload di -= */
VectorHandler& VectorHandler::operator -= (const VectorHandler& VH) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif
   
   for (integer i = iGetSize(); i > 0; ) {
      fDecCoef(i, VH.dGetCoef(i));
      i--;
   }
   
   return *this;
}

/* Assegnazione che copia il contenuto della memoria di due handlers */
VectorHandler& VectorHandler::operator = (const VectorHandler& VH) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif
   
   for (integer i = iGetSize(); i > 0; ) {
      fPutCoef(i, VH.dGetCoef(i));
      i--;
   }
   
   return *this;
}
   
/* Norma 2 del vettore */
doublereal VectorHandler::Dot(void) const {
#ifdef DEBUG
   IsValid();
#endif
   
   doublereal d2 = 0.;
   
   for (integer i = iGetSize(); i > 0;) {
      doublereal d = dGetCoef(i);
      d2 += d*d;
      i--;
   }
   
   return d2;
}

/* Norma del vettore */
doublereal VectorHandler::Norm(void) const {
   return sqrt(Dot());
}

/* VectorHandler - end */


/* MyVectorHandler - begin */

MyVectorHandler::MyVectorHandler(integer iSize) 
: fOwnsMemory(1), 
iMaxSize(iSize), iCurSize(iSize), pdVec(NULL), pdVecm1(NULL) 
{
   if (iSize > 0) {
      Resize(iSize);
   }
}

MyVectorHandler::MyVectorHandler(integer iSize, doublereal* pdTmpVec) 
: fOwnsMemory(0), 
iMaxSize(iSize), iCurSize(iSize), pdVec(pdTmpVec), pdVecm1(pdVec-1) 
{
#ifdef DEBUG
   IsValid();
#else
   NO_OP;
#endif
}
   
void MyVectorHandler::Resize(integer iSize) 
{
   if (iSize < 0) {
      cerr << "Negative size!" << endl;
      THROW(ErrGeneric());
   }

   
   if (fOwnsMemory) {
      if (pdVec != NULL) {
	 if (iSize > iMaxSize) {
	    doublereal* pd = NULL;
	    SAFENEWARR(pd, doublereal, iSize, SMmm);
	    for (integer i = iCurSize; i-- > 0; ) {
	       pd[i] = pdVec[i];
	    }
	    SAFEDELETEARR(pdVec, SMmm);
	    pdVec = pd;
	    pdVecm1 = pdVec-1;
	    iMaxSize = iCurSize = iSize;	    
	 } else {
	    iCurSize = iSize;
	 }
      } else {
	 SAFENEWARR(pdVec, doublereal, iSize, SMmm);
	 pdVecm1 = pdVec-1;
	 iMaxSize = iCurSize = iSize;
      }
   } else {
      if (pdVec != NULL) {
	 if (iSize > iMaxSize) {
	    cerr << "Can't resize to " << iSize 
	      << ": larger than max size " << iMaxSize << endl;
	    THROW(ErrGeneric());
	 } else {
	    iCurSize = iSize;
	 }
      } else {
	 cerr << "internal error!" << endl;
	 THROW(ErrGeneric());
      }
   }   
}

void MyVectorHandler::Detach(void) 
{
   if (fOwnsMemory) {
      if (pdVec != NULL) {
	 SAFEDELETEARR(pdVec, SMmm);
      }
      fOwnsMemory = 0;
   }
   iMaxSize = iCurSize = 0;
   pdVec = pdVecm1 = NULL;
}

void MyVectorHandler::Attach(integer iSize, doublereal* pd, integer iMSize) 
{
   if (fOwnsMemory || pdVec != NULL) {
      Detach();
      fOwnsMemory = 0;
   }
   iMaxSize = iCurSize = iSize;
   if (iMSize > 0) {
      if (iMSize >= iSize) {
	 iMaxSize = iMSize;
      } else if (iMSize < iSize) {
	 THROW(ErrGeneric());
      }
   }
   pdVec = pd;
   pdVecm1 = pdVec-1;
}

void MyVectorHandler::IsValid(void) const
{
   ASSERT(iMaxSize > 0);
   ASSERT(iCurSize >= 0 && iCurSize <= iMaxSize);
   ASSERT(pdVec != NULL);
   ASSERT(pdVecm1 < pdVec);
   ASSERT(pdVecm1+1 == pdVec);
   
#ifdef DEBUG_MEMMANAGER
   if (fOwnsMemory) {      
      ASSERT(SMmm.fIsBlock((void*)pdVec, iMaxSize*sizeof(doublereal)));
   } else {
      ASSERT(SMmm.fIsValid((void*)pdVec, iMaxSize*sizeof(doublereal)));
   }
#endif
}

void MyVectorHandler::Reset(doublereal dResetVal)
{
#ifdef DEBUG
   IsValid();
#endif
   
   ASSERT(iCurSize > 0);
   ASSERT(pdVec+iCurSize > pdVec);
   
   for (integer i = iGetSize(); i-- > 0; ) {
      pdVec[i] = dResetVal;
   }	
}

/* Somma e moltiplica per uno scalare */
VectorHandler& 
MyVectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iCurSize > 0);
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());
#endif

   if (d != 0.) {
      for (integer i = iGetSize(); i > 0; ) {
	 pdVecm1[i] += d*VH.dGetCoef(i);
	 i--;
      }
   }   

   return *this;
}

/* Moltiplica per uno scalare */
VectorHandler& 
MyVectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iCurSize > 0);
   ASSERT(VH.iGetSize() > 0);    
   ASSERT(iCurSize == VH.iGetSize());
#endif   

   if (d != 0.) {
      for (integer i = iGetSize(); i > 0; ) {	 
	 pdVecm1[i] = d*VH.dGetCoef(i);
	 i--;
      }
   } else {
      Reset(0.);
   }   
   
   return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler& MyVectorHandler::operator += (const VectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   for (integer i = iGetSize(); i > 0; i--) {
      pdVecm1[i] += VH.dGetCoef(i);
   }

   return *this;
}

/* Overload di += usato per la correzione della soluzione */
MyVectorHandler& MyVectorHandler::operator += (const MyVectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   doublereal* pdFrom = VH.pdGetVec();
   for (integer i = iGetSize(); i-- > 0; ) {
      pdVec[i] += pdFrom[i];
   }

   return *this;
}

/* Overload di -= usato per la correzione della soluzione */
VectorHandler& MyVectorHandler::operator -= (const VectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   for (integer i = iGetSize(); i > 0; i--) {
      pdVecm1[i] -= VH.dGetCoef(i);
   }

   return *this;
}

/* Overload di -= usato per la correzione della soluzione */
MyVectorHandler& MyVectorHandler::operator -= (const MyVectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   doublereal* pdFrom = VH.pdGetVec();
   for (integer i = iGetSize(); i-- > 0; ) {
      pdVec[i] -= pdFrom[i];
   }

   return *this;
}

/* Overload di = con copia termine per termine */
VectorHandler& MyVectorHandler::operator = (const VectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   for (integer i = iGetSize(); i > 0; i--) {
      pdVecm1[i] = VH.dGetCoef(i);
   }

   return *this;
}

/* Overload di = con copia termine per termine */
MyVectorHandler& MyVectorHandler::operator = (const MyVectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif
   
   ASSERT(VH.iGetSize() > 0);
   ASSERT(iCurSize == VH.iGetSize());

   doublereal* pdFrom = VH.pdGetVec();
   for (integer i = iGetSize(); i-- > 0; ) {
      pdVec[i] = pdFrom[i];
   }

   return *this;
}

void MyVectorHandler::Add(integer iRow, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif   

   pdVecm1[iRow] += v.dGet(1);
   pdVecm1[++iRow] += v.dGet(2);
   pdVecm1[++iRow] += v.dGet(3);
}

void MyVectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif   

   pdVecm1[iRow] -= v.dGet(1);
   pdVecm1[++iRow] -= v.dGet(2);
   pdVecm1[++iRow] -= v.dGet(3);
}

void MyVectorHandler::Put(integer iRow, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif   

   pdVecm1[iRow] = v.dGet(1);
   pdVecm1[++iRow] = v.dGet(2);
   pdVecm1[++iRow] = v.dGet(3);
}

/* MyVectorHandler - end */


/* MatrixHandler - begin */

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler& MatrixHandler::operator +=(const SubMatrixHandler& SubMH) {
   return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler& MatrixHandler::operator -=(const SubMatrixHandler& SubMH) {
   return SubMH.SubFrom(*this);
}

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler& MatrixHandler::operator +=(const VariableSubMatrixHandler& SubMH) {
   return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler& MatrixHandler::operator -=(const VariableSubMatrixHandler& SubMH) {
   return SubMH.SubFrom(*this);
}

MatrixHandler& MatrixHandler::ScalarMul(const doublereal& d)
{
   for (integer i = 1; i <= iGetNumRows(); i++) {
      for (integer j = 0; j <= iGetNumCols(); i++) {
	 fPutCoef(i, j, dGetCoef(i, j)*d);
      }
   }
   return *this;
}



ostream& operator << (ostream& out, const MatrixHandler& MH)
{
   for (integer i = 1; i <= MH.iGetNumRows(); i++) {
      for (integer j = 1; j <= MH.iGetNumCols(); j++) {
	 out << setw(16) << MH.dGetCoef(i, j);
      }
      out << endl;
   }
   return out;
}

/* MatrixHandler - end */
