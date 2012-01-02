/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <algorithm>

#include "stlvh.h"
#include "matvec3.h"

/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

/* STLVectorHandler - begin */

STLVectorHandler::STLVectorHandler(integer iSize)
: std::vector<doublereal>(iSize)
{
	ASSERT(iSize >= 0);
}

STLVectorHandler::~STLVectorHandler(void)
{
	NO_OP;
}

#ifdef DEBUG
/* Usata per il debug */
void
STLVectorHandler::IsValid(void) const
{
	ASSERT(!empty());
}
#endif // DEBUG

doublereal *
STLVectorHandler::pdGetVec(void) const
{
#ifdef DEBUG
	IsValid();
#endif // DEBUG

	return &((*const_cast<STLVectorHandler *>(this))[0]);
}

integer
STLVectorHandler::iGetSize(void) const
{
	return size();
}

void
STLVectorHandler::Reset(void)
{
#ifdef DEBUG
	IsValid();
#endif // DEBUG

	std::fill(begin(), end(), 0.);
}

void
STLVectorHandler::Resize(integer iNewSize)
{
	ASSERT(iNewSize >= 0);

	resize(iNewSize);
}

void
STLVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	(*this)[--iRow] = dCoef;
}

void
STLVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	(*this)[--iRow] += dCoef;
}

void
STLVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	(*this)[--iRow] -= dCoef;
}

const doublereal&
STLVectorHandler::dGetCoef(integer iRow) const
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	return (*this)[--iRow];
}

const doublereal&
STLVectorHandler::operator () (integer iRow) const
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	return (*this)[--iRow];
}

doublereal&
STLVectorHandler::operator () (integer iRow)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) <= size());

	return (*this)[--iRow];
}

/* Somma un Vec3 nella posizione desiderata */
void
STLVectorHandler::Add(integer iRow, const Vec3& v)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) + 2 <= size());

	(*this)[--iRow] += v(1);
	(*this)[++iRow] += v(2);
	(*this)[++iRow] += v(3);
}

/* Sottrae un Vec3 nella posizione desiderata */
void
STLVectorHandler::Sub(integer iRow, const Vec3& v)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) + 2 <= size());

	(*this)[--iRow] -= v(1);
	(*this)[++iRow] -= v(2);
	(*this)[++iRow] -= v(3);
}

/* Scrive un Vec3 nella posizione desiderata */
void
STLVectorHandler::Put(integer iRow, const Vec3& v)
{
	ASSERT(iRow > 0);
	ASSERT(unsigned(iRow) + 2 <= size());

	(*this)[--iRow] = v(1);
	(*this)[++iRow] = v(2);
	(*this)[++iRow] = v(3);
}

/* Somma e moltiplica per uno scalare */
VectorHandler&
STLVectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = d*VH(i);
		i--;
		(*this)[i] = dd;
	}

	return *this;
}

/* Somma e moltiplica per uno scalare v = VH + d * VH1 */
VectorHandler&
STLVectorHandler::ScalarAddMul(const VectorHandler& VH,
	const VectorHandler& VH1,
	const doublereal& d)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());
	ASSERT(VH1.iGetSize() > 0);
	ASSERT(unsigned(VH1.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = VH(i) + d*VH1(i);
		i--;
		(*this)[i] = dd;
	}

	return *this;
}

/* Moltiplica per uno scalare */
VectorHandler&
STLVectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0; ) {
		doublereal dd = d*VH(i);
		i--;
		(*this)[i] = dd;
	}

	return *this;
}

/* Overload di += usato per la correzione della soluzione */
VectorHandler&
STLVectorHandler::operator += (const VectorHandler& VH)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = VH(i);
		i--;
		(*this)[i] += dd;
	}

	return *this;
}

/* Overload di -= */
VectorHandler&
STLVectorHandler::operator -= (const VectorHandler& VH)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = VH(i);
		i--;
		(*this)[i] -= dd;
	}

	return *this;
}

/* Overload di *= */
VectorHandler&
STLVectorHandler::operator *= (const doublereal &d)
{
	for (std::vector<doublereal>::iterator i = begin(); i != end(); ++i) {
		*i *= d;
	}

	return *this;
}

/* Assegnazione che copia il contenuto della memoria di due handlers */
VectorHandler&
STLVectorHandler::operator = (const VectorHandler& VH)
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = VH(i);
		i--;
		(*this)[i] = dd;
	}

	return *this;
}

/* Norma 2 del vettore */
doublereal
STLVectorHandler::Dot(void) const
{
#ifdef DEBUG
	IsValid();
#endif // DEBUG

	doublereal d2 = 0.;

	for (std::vector<doublereal>::const_iterator i = begin();
		i != end(); ++i)
	{
		doublereal d = *i;
		d2 += d*d;
	}

	return d2;
}


/* Prodotto Scalare di due vettori */
doublereal
STLVectorHandler::InnerProd(const VectorHandler& VH) const
{
	ASSERT(VH.iGetSize() > 0);
	ASSERT(unsigned(VH.iGetSize()) == size());

	doublereal d = 0.;

	for (integer i = iGetSize(); i > 0;) {
		doublereal dd = VH(i);
		i--;
		d += (*this)[i]*dd;
	}

	return d;
}

/* STLVectorHandler - end */

