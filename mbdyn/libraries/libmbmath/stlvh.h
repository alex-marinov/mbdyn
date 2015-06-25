/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#ifndef STLVH_H
#define STLVH_H

#include <vector>
#include "vh.h"

class Vec3;
class SubVectorHandler;

/* STLVectorHandler - begin */

class STLVectorHandler : public VectorHandler, public std::vector<doublereal> {
public:
	STLVectorHandler(integer iSize = 0);

	virtual ~STLVectorHandler(void);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const;
#endif /* DEBUG */

	virtual doublereal* pdGetVec(void) const;

	virtual integer iGetSize(void) const;

	virtual void Reset(void);

	virtual void Resize(integer iNewSize);

	virtual void PutCoef(integer iRow, const doublereal& dCoef);

	virtual void IncCoef(integer iRow, const doublereal& dCoef);

	virtual void DecCoef(integer iRow, const doublereal& dCoef);

	virtual const doublereal& dGetCoef(integer iRow) const;

	virtual const doublereal& operator () (integer iRow) const;

	virtual doublereal& operator () (integer iRow);

	/* Somma un Vec3 nella posizione desiderata */
	virtual void Add(integer iRow, const Vec3& v);

	/* Sottrae un Vec3 nella posizione desiderata */
	virtual void Sub(integer iRow, const Vec3& v);

	/* Scrive un Vec3 nella posizione desiderata */
	virtual void Put(integer iRow, const Vec3& v);

	/* Somma e moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const doublereal& d);

	/* Somma e moltiplica per uno scalare v = VH + d * VH1 */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
		const doublereal& d);

	/* Moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarMul(const VectorHandler& VH, const doublereal& d);

	/* Overload di += usato per la correzione della soluzione */
	virtual VectorHandler& operator += (const VectorHandler& VH);

	/* Overload di -= */
	virtual VectorHandler& operator -= (const VectorHandler& VH);

	/* Overload di *= */
	virtual VectorHandler& operator *= (const doublereal &d);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual VectorHandler& operator = (const VectorHandler& VH);

	/* Norma 2 del vettore */
	virtual doublereal Dot(void) const;

	/* Prodotto Scalare di due vettori */
	virtual doublereal InnerProd(const VectorHandler& VH) const;
};

/* STLVectorHandler - end */

#endif // STLVH_H

