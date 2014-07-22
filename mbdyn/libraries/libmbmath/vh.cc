/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>	/* for memset() */
#include <cmath>

#include <iostream>
#include <iomanip>

#include "matvec3.h"
#include "submat.h"

class SubVectorHandler;

/* VectorHandler - begin */

VectorHandler::~VectorHandler(void)
{
	NO_OP;
}

void
VectorHandler::ResizeReset(integer iNewRow)
{
	Resize(iNewRow);
	Reset();
}

/* Somma un Vec3 nella posizione desiderata */
void
VectorHandler::Add(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
	ASSERT(iRow > 0);
	ASSERT(iGetSize() >= iRow + 2);
#endif /* DEBUG */

	operator()(iRow) += v(1);
	operator()(++iRow) += v(2);
	operator()(++iRow) += v(3);
}

/* Sottrae un Vec3 nella posizione desiderata */
void
VectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
	ASSERT(iRow > 0);
	ASSERT(iGetSize() >= iRow + 2);
#endif /* DEBUG */

	operator()(iRow) -= v(1);
	operator()(++iRow) -= v(2);
	operator()(++iRow) -= v(3);
}

/* Scrive un Vec3 nella posizione desiderata */
void
VectorHandler::Put(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
	ASSERT(iRow > 0);
	ASSERT(iGetSize() >= iRow + 2);
#endif /* DEBUG */

	operator()(iRow) = v(1);
	operator()(++iRow) = v(2);
	operator()(++iRow) = v(3);
}

/* Somma e moltiplica per uno scalare */
VectorHandler&
VectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) += d*VH(i);
	}

	return *this;
}

/* Somma e moltiplica per uno scalare this = VH + d * VH1 */
VectorHandler&
VectorHandler::ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
		const doublereal& d)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
	VH1.IsValid();
	ASSERT(iGetSize() == VH1.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) = VH(i) + d*VH1(i);
	}

	return *this;
}

/* Moltiplica per uno scalare */
VectorHandler&
VectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) = d*VH(i);
	}

	return *this;
}


/* Overload di += usato per la correzione della soluzione */
VectorHandler&
VectorHandler::operator += (const VectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) += VH(i);
	}

	return *this;
}

/* Overload di += usato per l'assemblaggio del residuo */
VectorHandler&
VectorHandler::operator += (const SubVectorHandler& SubVH)
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
VectorHandler&
VectorHandler::operator -= (const VectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) -= VH(i);
	}

	return *this;
}

/* Overload di *= */
VectorHandler&
VectorHandler::operator *= (const doublereal& d)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		operator()(i) *= d;
	}

	return *this;
}

/* Assegnazione che copia il contenuto della memoria di due handlers */
VectorHandler&
VectorHandler::operator = (const VectorHandler& VH)
{
#ifdef DEBUG
	VH.IsValid();
#endif /* DEBUG */

	integer nr = VH.iGetSize();
	Resize(nr);
	for (integer i = nr; i > 0; i--) {
		operator()(i) = VH(i);
	}

	return *this;
}

/* Norma 2 del vettore */
doublereal
VectorHandler::Dot(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	doublereal d2 = 0.;

	for (integer i = iGetSize(); i > 0; i--) {
 		doublereal d = operator()(i);
		d2 += d*d;
	}

	return d2;
}

/* Norma del vettore */
doublereal
VectorHandler::Norm(void) const
{
	return std::sqrt(Dot());
}

/* Prodotto Scalare fra due Vettori */
doublereal
VectorHandler::InnerProd(const VectorHandler& VH) const
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
#endif /* DEBUG */

	doublereal d2 = 0.;

	for (integer i = iGetSize(); i > 0; i--) {
		d2 += operator()(i)*VH(i);
	}

	return d2;
}

std::ostream&
operator << (std::ostream& out, const VectorHandler& VH)
{
	for (integer i = 1; i <= VH.iGetSize(); i++) {
		out << std::setw(16) << VH(i) << std::endl;
	}
	return out;
}

/* VectorHandler - end */


/* MyVectorHandler - begin */

MyVectorHandler::MyVectorHandler(integer iSize, doublereal* pdTmpVec)
: bOwnsMemory(false),
iMaxSize(iSize), iCurSize(iSize), pdVecm1(0)
{
	if (iSize == 0) {
		ASSERT(pdVecm1 == NULL);

	} else {
		if (pdTmpVec == NULL) {
			bOwnsMemory = true;
			Resize(iSize);
			Reset();
		} else {
			pdVecm1 = pdTmpVec - 1;
		}
#ifdef DEBUG
		IsValid();
#endif /* DEBUG */
	}
}

MyVectorHandler::MyVectorHandler(const MyVectorHandler& VH)
: bOwnsMemory(false),
iMaxSize(VH.iCurSize), iCurSize(VH.iCurSize), pdVecm1(0)
{
	if (iCurSize == 0) {
		ASSERT(VH.pdVecm1 == 0);

	} else {
		bOwnsMemory = true;
		Resize(iCurSize);
#ifdef DEBUG
		IsValid();
#endif /* DEBUG */

		for (integer i = 1; i <= iCurSize; i++) {
			pdVecm1[i] = VH.pdVecm1[i];
		}
	}
}

MyVectorHandler::~MyVectorHandler(void)
{
	Detach();
}

void
MyVectorHandler::Resize(integer iSize)
{
	if (iSize < 0) {
		silent_cerr("Negative size!" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (pdVecm1 == NULL || bOwnsMemory) {
		if (pdVecm1 != NULL) {
			if (iSize > iMaxSize) {
				doublereal* pd = NULL;

				SAFENEWARR(pd, doublereal, iSize);
				pd--;
#ifdef HAVE_MEMMOVE
				memmove(pd + 1, pdVecm1 + 1, iCurSize*sizeof(doublereal));
#else /* ! HAVE_MEMMOVE */
				for (integer i = iCurSize; i > 0; i--) {
					pd[i] = pdVecm1[i];
				}
#endif /* ! HAVE_MEMMOVE */
				doublereal *pdv = pdVecm1 + 1;
				SAFEDELETEARR(pdv);
				pdVecm1 = pd;
				iMaxSize = iCurSize = iSize;

			} else {
				iCurSize = iSize;
			}

		} else {
			SAFENEWARR(pdVecm1, doublereal, iSize);
			pdVecm1--;
			iMaxSize = iCurSize = iSize;
			bOwnsMemory = true;
		}

	} else {
		if (pdVecm1 != NULL) {
			if (iSize > iMaxSize) {
				silent_cerr("Can't resize to " << iSize
					<< ": larger than "
					"max size " << iMaxSize << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			
			} else {
				iCurSize = iSize;
			}
		} else {
			silent_cerr("internal error!" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

void
MyVectorHandler::Detach(void)
{
	if (bOwnsMemory) {
		if (pdVecm1 != NULL) {
			doublereal *pd = pdVecm1 + 1;
			SAFEDELETEARR(pd);
		}
		bOwnsMemory = false;
	}

	iMaxSize = iCurSize = 0;
	pdVecm1 = NULL;
}

void
MyVectorHandler::Attach(integer iSize, doublereal* pd, integer iMSize)
{
	if (bOwnsMemory || pdVecm1 != NULL) {
		Detach();
		bOwnsMemory = false;
	}

	iMaxSize = iCurSize = iSize;
	if (iMSize > 0) {
		if (iMSize >= iSize) {
			iMaxSize = iMSize;

		} else if (iMSize < iSize) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	pdVecm1 = pd - 1;
}

#ifdef DEBUG
void
MyVectorHandler::IsValid(void) const
{
	ASSERT(iMaxSize > 0);
	ASSERT(iCurSize >= 0 && iCurSize <= iMaxSize);
	ASSERT(pdVecm1 != NULL);

#ifdef DEBUG_MEMMANAGER
	if (bOwnsMemory) {
		ASSERT(defaultMemoryManager.fIsBlock(pdVecm1 + 1,
					iMaxSize*sizeof(doublereal)));
	} else {
		ASSERT(defaultMemoryManager.fIsValid(pdVecm1 + 1,
					iMaxSize*sizeof(doublereal)));
	}
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */

void
MyVectorHandler::Reset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	ASSERT(iCurSize > 0);

#if defined HAVE_MEMSET
	memset(pdVecm1 + 1, 0, iGetSize()*sizeof(doublereal));
#else /* !HAVE_MEMSET */
	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] = 0.;
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
			pdVecm1[i] += d*VH(i);
		}
	}

	return *this;
}

/* Somma e moltiplica per uno scalare this = VH + d * VH1 */
VectorHandler&
MyVectorHandler::ScalarAddMul(const VectorHandler& VH,
		const VectorHandler& VH1, const doublereal& d)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
	ASSERT(iGetSize() == VH.iGetSize());
	VH1.IsValid();
	ASSERT(iGetSize() == VH1.iGetSize());
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] = VH(i) + d*VH1(i);
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
		Reset();
	} else {
		for (integer i = iGetSize(); i > 0; i--) {
			pdVecm1[i] = d*VH(i);
		}
	}

	return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler&
MyVectorHandler::operator += (const VectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] += VH(i);
	}

	return *this;
}

/* Overload di += usato per la correzione della soluzione */
MyVectorHandler&
MyVectorHandler::operator += (const MyVectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	doublereal* pdFrom = VH.pdGetVec() - 1;
	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] += pdFrom[i];
	}

	return *this;
}

/* Overload di -= usato per la correzione della soluzione */
VectorHandler&
MyVectorHandler::operator -= (const VectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] -= VH(i);
	}

	return *this;
}

/* Overload di *= */
VectorHandler&
MyVectorHandler::operator *= (const doublereal& d)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] *= d;
	}

	return *this;
}

/* Overload di -= usato per la correzione della soluzione */
MyVectorHandler&
MyVectorHandler::operator -= (const MyVectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	doublereal* pdFrom = VH.pdGetVec() - 1;
	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] -= pdFrom[i];
	}

	return *this;
}

/* Overload di = con copia termine per termine */
VectorHandler&
MyVectorHandler::operator = (const VectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] = VH(i);
	}

	return *this;
}

/* Overload di = con copia termine per termine */
MyVectorHandler&
MyVectorHandler::operator = (const MyVectorHandler& VH)
{
#ifdef DEBUG
	IsValid();
	VH.IsValid();
#endif /* DEBUG */

	ASSERT(VH.iGetSize() > 0);
	ASSERT(iCurSize == VH.iGetSize());

	doublereal* pdFrom = VH.pdGetVec() - 1;
	for (integer i = iGetSize(); i > 0; i--) {
		pdVecm1[i] = pdFrom[i];
	}

	return *this;
}

void
MyVectorHandler::Add(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif /* DEBUG */

	pdVecm1[iRow] += v.dGet(1);
	pdVecm1[++iRow] += v.dGet(2);
	pdVecm1[++iRow] += v.dGet(3);
}

void
MyVectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize-2));
#endif /* DEBUG */

	pdVecm1[iRow] -= v.dGet(1);
	pdVecm1[++iRow] -= v.dGet(2);
	pdVecm1[++iRow] -= v.dGet(3);
}

void
MyVectorHandler::Put(integer iRow, const Vec3& v)
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

	for (integer i = iCurSize; i > 0; i--) {
		d2 += pdVecm1[i]*pdVecm1[i];
	}

	return d2;
}

/* MyVectorHandler - end */
