/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* solution manager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <string.h>	/* for memset() */

#include <ac/iostream>
#include <ac/iomanip>

#include <solman.h>
#include <submat.h>
#include <matvec3.h>

/* Zero for sparse vector and matrix handlers */
const doublereal dZero = 0.;

/* VectorHandler - begin */

VectorHandler::~VectorHandler(void) {
	NO_OP;
}

/* Somma un Vec3 nella posizione desiderata */
void VectorHandler::Add(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif /* DEBUG */

   IncCoef(iRow, v.dGet(1));
   IncCoef(++iRow, v.dGet(2));
   IncCoef(++iRow, v.dGet(3));
}

/* Sottrae un Vec3 nella posizione desiderata */
void VectorHandler::Sub(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif /* DEBUG */

   DecCoef(iRow, v.dGet(1));
   DecCoef(++iRow, v.dGet(2));
   DecCoef(++iRow, v.dGet(3));
}

/* Scrive un Vec3 nella posizione desiderata */
void VectorHandler::Put(integer iRow, const Vec3& v) {
#ifdef DEBUG
   IsValid();
   ASSERT(iRow > 0);
   ASSERT(iGetSize() >= iRow+2);
#endif /* DEBUG */

   PutCoef(iRow, v.dGet(1));
   PutCoef(++iRow, v.dGet(2));
   PutCoef(++iRow, v.dGet(3));
}

/* Somma e moltiplica per uno scalare */
VectorHandler& VectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      IncCoef(i, d*VH.dGetCoef(i));
   }

   return *this;
}

/* Somma e moltiplica per uno scalare this = VH + d * VH1 */
VectorHandler& VectorHandler::ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
				const doublereal& d) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
   VH1.IsValid();
   ASSERT(iGetSize() == VH1.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      PutCoef(i, VH.dGetCoef(i) + d*VH1.dGetCoef(i));
   }

   return *this;
}

/* Moltiplica per uno scalare */
VectorHandler& VectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      PutCoef(i, d*VH.dGetCoef(i));
   }

   return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler& VectorHandler::operator += (const VectorHandler& VH) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      IncCoef(i, VH.dGetCoef(i));
   }

   return *this;
}

/* Overload di += usato per l'assemblaggio del residuo */
VectorHandler& VectorHandler::operator += (const SubVectorHandler& SubVH)
{
#ifdef DEBUG
   SubVH.IsValid();
#endif /* DEBUG */

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
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      DecCoef(i, VH.dGetCoef(i));
   }

   return *this;
}

/* Assegnazione che copia il contenuto della memoria di due handlers */
VectorHandler& VectorHandler::operator = (const VectorHandler& VH) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
      PutCoef(i, VH.dGetCoef(i));
   }

   return *this;
}

/* Norma 2 del vettore */
doublereal VectorHandler::Dot(void) const {
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   doublereal d2 = 0.;

   for (integer i = iGetSize(); i > 0; i--) {
      doublereal d = dGetCoef(i);
      d2 += d*d;
   }

   return d2;
}

/* Norma del vettore */
doublereal VectorHandler::Norm(void) const {
   return sqrt(Dot());
}

/* Prodotto Scalare fra due Vettori */
doublereal VectorHandler::InnerProd(const VectorHandler& VH) const {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

   doublereal d2 = 0.;

   for (integer i = iGetSize(); i > 0; i--) {
      d2 += dGetCoef(i)*VH.dGetCoef(i);
   }

   return d2;
}


/* VectorHandler - end */


/* MyVectorHandler - begin */

MyVectorHandler::MyVectorHandler(integer iSize)
: bOwnsMemory(true),
iMaxSize(iSize), iCurSize(iSize), pdVec(NULL), pdVecm1(NULL)
{
   if (iSize > 0) {
      Resize(iSize);
   }
}

MyVectorHandler::MyVectorHandler(integer iSize, doublereal* pdTmpVec)
: bOwnsMemory(false),
iMaxSize(iSize), iCurSize(iSize), pdVec(pdTmpVec), pdVecm1(pdVec-1)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */
}

MyVectorHandler::~MyVectorHandler(void) {
	Detach();
}

void MyVectorHandler::Resize(integer iSize)
{
   if (iSize < 0) {
      std::cerr << "Negative size!" << std::endl;
      THROW(ErrGeneric());
   }


   if (bOwnsMemory) {
      if (pdVec != NULL) {
	 if (iSize > iMaxSize) {
	    doublereal* pd = NULL;
	    SAFENEWARR(pd, doublereal, iSize);
	    for (integer i = iCurSize; i-- > 0; ) {
	       pd[i] = pdVec[i];
	    }
	    SAFEDELETEARR(pdVec);
	    pdVec = pd;
	    pdVecm1 = pdVec-1;
	    iMaxSize = iCurSize = iSize;
	 } else {
	    iCurSize = iSize;
	 }
      } else {
	 SAFENEWARR(pdVec, doublereal, iSize);
	 pdVecm1 = pdVec-1;
	 iMaxSize = iCurSize = iSize;
      }
   } else {
      if (pdVec != NULL) {
	 if (iSize > iMaxSize) {
	    std::cerr << "Can't resize to " << iSize
	      << ": larger than max size " << iMaxSize << std::endl;
	    THROW(ErrGeneric());
	 } else {
	    iCurSize = iSize;
	 }
      } else {
         std::cerr << "internal error!" << std::endl;
	 THROW(ErrGeneric());
      }
   }
}

void MyVectorHandler::Detach(void)
{
   if (bOwnsMemory) {
      if (pdVec != NULL) {
	 SAFEDELETEARR(pdVec);
      }
      bOwnsMemory = false;
   }
   iMaxSize = iCurSize = 0;
   pdVec = pdVecm1 = NULL;
}

void MyVectorHandler::Attach(integer iSize, doublereal* pd, integer iMSize)
{
   if (bOwnsMemory || pdVec != NULL) {
      Detach();
      bOwnsMemory = false;
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

#ifdef DEBUG
void MyVectorHandler::IsValid(void) const
{
   ASSERT(iMaxSize > 0);
   ASSERT(iCurSize >= 0 && iCurSize <= iMaxSize);
   ASSERT(pdVec != NULL);
   ASSERT(pdVecm1 < pdVec);
   ASSERT(pdVecm1+1 == pdVec);

#ifdef DEBUG_MEMMANAGER
   if (bOwnsMemory) {
      ASSERT(defaultMemoryManager.fIsBlock(pdVec, iMaxSize*sizeof(doublereal)));
   } else {
      ASSERT(defaultMemoryManager.fIsValid(pdVec, iMaxSize*sizeof(doublereal)));
   }
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */

void MyVectorHandler::Reset(doublereal dResetVal)
{
#ifdef DEBUG
   IsValid();
#endif /* DEBUG */

   ASSERT(iCurSize > 0);
   ASSERT(pdVec+iCurSize > pdVec);

#ifdef HAVE_MEMSET
   if (dResetVal == 0.) {
	  memset(pdVec, 0, iGetSize()*sizeof(doublereal));
   } else {
#endif /* HAVE_MEMSET */
   for (integer i = iGetSize(); i-- > 0; ) {
      pdVec[i] = dResetVal;
   }
#ifdef HAVE_MEMSET
   }
#endif /* HAVE_MEMSET */
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
#endif /* DEBUG */

   if (d != 0.) {
      for (integer i = iGetSize(); i > 0; i--) {
	 pdVecm1[i] += d*VH.dGetCoef(i);
      }
   }

   return *this;
}

/* Somma e moltiplica per uno scalare this = VH + d * VH1 */
VectorHandler& MyVectorHandler::ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
				const doublereal& d) {
#ifdef DEBUG
   IsValid();
   VH.IsValid();
   ASSERT(iGetSize() == VH.iGetSize());
   VH1.IsValid();
   ASSERT(iGetSize() == VH1.iGetSize());
#endif /* DEBUG */

   for (integer i = iGetSize(); i > 0; i--) {
	 pdVecm1[i] = VH.dGetCoef(i) + d*VH1.dGetCoef(i);
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
#endif /* DEBUG */

   if (d == 0.) {
      Reset(0.);
   } else {
      for (integer i = iGetSize(); i > 0; i--) {
	 pdVecm1[i] = d*VH.dGetCoef(i);
      }
   }

   return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler& MyVectorHandler::operator += (const VectorHandler& VH)
{
#ifdef DEBUG
   IsValid();
   VH.IsValid();
#endif /* DEBUG */

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
#endif /* DEBUG */

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
#endif /* DEBUG */

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
#endif /* DEBUG */

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
#endif /* DEBUG */

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
#endif /* DEBUG */

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
#endif /* DEBUG */

   pdVecm1[iRow] += v.dGet(1);
   pdVecm1[++iRow] += v.dGet(2);
   pdVecm1[++iRow] += v.dGet(3);
}

void MyVectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif /* DEBUG */

   pdVecm1[iRow] -= v.dGet(1);
   pdVecm1[++iRow] -= v.dGet(2);
   pdVecm1[++iRow] -= v.dGet(3);
}

void MyVectorHandler::Put(integer iRow, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif /* DEBUG */

   pdVecm1[iRow] = v.dGet(1);
   pdVecm1[++iRow] = v.dGet(2);
   pdVecm1[++iRow] = v.dGet(3);
}

/* Norma 2 del vettore */
doublereal
MyVectorHandler::Dot(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	doublereal d2 = 0.;

	for (integer i = iCurSize; --i > 0;) {
		d2 += pdVec[i]*pdVec[i];
	}

	return d2;
}

/* Norma del vettore */
doublereal
MyVectorHandler::Norm(void) const
{
	return sqrt(Dot());
}

/* MyVectorHandler - end */


/* MatrixHandler - begin */

MatrixHandler::~MatrixHandler(void)
{
	NO_OP;
}

void
MatrixHandler::Reset(const doublereal& dResetVal)
{
	Init(dResetVal);
}

/* Impacchetta la matrice; restituisce il numero di elementi 
 * diversi da zero */
integer
MatrixHandler::PacMat(void)
{
	return 0L;
}

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
      for (integer j = 0; j <= iGetNumCols(); j++) {
	 PutCoef(i, j, dGetCoef(i, j)*d);
      }
   }
   return *this;
}

VectorHandler& MatrixHandler::MatVecMul(
		VectorHandler& out,
		const VectorHandler& in
) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumRows(); r++) {
		doublereal d = 0.;
		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += dGetCoef(r, c)*in.dGetCoef(c);
		}
		out.PutCoef(r, d);
	}
	return out;
}

VectorHandler& MatrixHandler::MatTVecMul(
		VectorHandler& out,
		const VectorHandler& in
) const
{
	if (out.iGetSize() != iGetNumCols()
			|| in.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumCols(); r++) {
		doublereal d = 0.;
		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += dGetCoef(c, r)*in.dGetCoef(c);
		}
		out.PutCoef(r, d);
	}
	return out;
}

VectorHandler& MatrixHandler::MatVecIncMul(
		VectorHandler& out,
		const VectorHandler& in
) const
{
	if (out.iGetSize() != iGetNumRows()
			|| in.iGetSize() != iGetNumCols()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumRows(); r++) {
		doublereal d = 0.;
		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += dGetCoef(r, c)*in.dGetCoef(c);
		}
		out.IncCoef(r, d);
	}
	return out;
}

VectorHandler& MatrixHandler::MatTVecIncMul(
		VectorHandler& out,
		const VectorHandler& in
) const
{
	if (out.iGetSize() != iGetNumCols()
			|| in.iGetSize() != iGetNumRows()) {
		THROW(ErrGeneric());
	}

	for (integer r = 1; r <= iGetNumCols(); r++) {
		doublereal d = 0.;
		for (integer c = 1; c <= in.iGetSize(); c++) {
			d += dGetCoef(c, r)*in.dGetCoef(c);
		}
		out.IncCoef(r, d);
	}
	return out;
}

std::ostream&
operator << (std::ostream& out, const MatrixHandler& MH)
{
   for (integer i = 1; i <= MH.iGetNumRows(); i++) {
      for (integer j = 1; j <= MH.iGetNumCols(); j++) {
	 out << std::setw(16) << MH.dGetCoef(i, j);
      }
      out << std::endl;
   }
   return out;
}

/* MatrixHandler - end */


/* SolutionDataManager - begin */

SolutionDataManager::~SolutionDataManager(void)
{
	NO_OP;
}

/* SolutionDataManager - end */


/* SolutionManager - begin */

SolutionManager::~SolutionManager(void)
{
	NO_OP;
}

/* Inizializzatore "speciale" */
void
SolutionManager::MatrInitialize(const doublereal& d)
{
	MatrInit(d);
}

void
SolutionManager::LinkToSolution(const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) {
	NO_OP;
}

/* SolutionManager - end */

