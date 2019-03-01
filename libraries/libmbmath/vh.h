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

/*****************************************************************************
 *                                                                           *
 *                            SOLUTION MANAGER                               *
 *                                                                           *
 *****************************************************************************/

/* Pierangelo Masarati */


#ifndef VH_H
#define VH_H

#include <iostream>
#include "ac/f2c.h"

/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

class Vec3;
class SubVectorHandler;

/* VectorHandler - begin */

/* Gestore di vettori. Usa spazio messo a disposizione da altri. */

/* questa classe non e' virtuale pura, ma viene definita nel caso banale di
 * vettore pieno in quanto e' di interesse pratico. Il pacchetto <harwrap>
 * si basa su questa per la gestione dei vettori residuo e soluzione */

class VectorHandler {
public:
	virtual ~VectorHandler(void);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const = 0;
#endif /* DEBUG */

	virtual doublereal* pdGetVec(void) const = 0;

	virtual integer iGetSize(void) const = 0;

	virtual void Reset(void) = 0;

	virtual void Resize(integer iNewSize) = 0;

	virtual void ResizeReset(integer);

	virtual void PutCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual void IncCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual void DecCoef(integer iRow, const doublereal& dCoef) = 0;

	virtual const doublereal& dGetCoef(integer iRow) const = 0;

	virtual const doublereal& operator () (integer iRow) const = 0;

	virtual doublereal& operator () (integer iRow) = 0;

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

	/* Overload di += usato per l'assemblaggio del residuo */
	virtual VectorHandler& operator += (const SubVectorHandler& SubVH);

	/* Overload di -= */
	virtual VectorHandler& operator -= (const VectorHandler& VH);

	/* Overload di *= */
	virtual VectorHandler& operator *= (const doublereal &d);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual VectorHandler& operator = (const VectorHandler& VH);

	/* Norma 2 del vettore */
	virtual doublereal Dot(void) const;

	/* Norma del vettore */
	virtual doublereal Norm(void) const;

	/* Prodotto Scalare di due vettori */
	virtual doublereal InnerProd(const VectorHandler& VH) const;
};

extern std::ostream&
operator << (std::ostream& out, const VectorHandler& VH);

/* VectorHandler - end */


/* MyVectorHandler - begin */

class MyVectorHandler : public VectorHandler {
	friend class FullMatrixHandler;
protected:
	bool bOwnsMemory;

protected:
	integer iMaxSize;
	integer iCurSize;

	doublereal* pdVecm1;

public:		// needed by Shell4 :(
	MyVectorHandler(const MyVectorHandler&);

public:
	MyVectorHandler(integer iSize = 0, doublereal* pdTmpVec = NULL);

	virtual ~MyVectorHandler(void);

	virtual void Resize(integer iSize);

	void Detach(void);

	void Attach(integer iSize, doublereal* pd, integer iMSize = 0);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const;
#endif /* DEBUG */

	virtual inline doublereal* pdGetVec(void) const;

	virtual inline integer iGetSize(void) const;

	virtual void Reset(void);

	virtual inline void PutCoef(integer iRow, const doublereal& dCoef);

	virtual inline void IncCoef(integer iRow, const doublereal& dCoef);

	virtual inline void DecCoef(integer iRow, const doublereal& dCoef);

	virtual inline const doublereal& dGetCoef(integer iRow) const;

	virtual inline const doublereal& operator () (integer iRow) const;

	virtual inline doublereal& operator () (integer iRow);

	/* Somma un Vec3 nella posizione desiderata */
	virtual void Add(integer iRow, const Vec3& v);

	/* Sottrae un Vec3 nella posizione desiderata */
	virtual void Sub(integer iRow, const Vec3& v);

	/* Scrive un Vec3 nella posizione desiderata */
	virtual void Put(integer iRow, const Vec3& v);

	/* Somma e moltiplica per uno scalare v = VH + d * VH1 */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
			const doublereal& d);

	/* Somma e moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarAddMul(const VectorHandler& VH, const doublereal& d);

	/* Moltiplica per uno scalare */
	virtual VectorHandler&
	ScalarMul(const VectorHandler& VH, const doublereal& d);

	/* Overload di += usato per la correzione della soluzione */
	virtual VectorHandler& operator += (const VectorHandler& VH);

	/* Overload di += usato per la correzione della soluzione */
	virtual MyVectorHandler& operator += (const MyVectorHandler& VH);

	/* Overload di -= */
	virtual VectorHandler& operator -= (const VectorHandler& VH);

	/* Overload di *= */
	virtual VectorHandler& operator *= (const doublereal& d);

	/* Overload di -= */
	virtual MyVectorHandler& operator -= (const MyVectorHandler& VH);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual VectorHandler& operator = (const VectorHandler& VH);

	/* Assegnazione che copia il contenuto della memoria di due handlers */
	virtual MyVectorHandler& operator = (const MyVectorHandler& VH);

	/* Norma 2 del vettore */
	doublereal Dot(void) const;
};

inline doublereal*
MyVectorHandler::pdGetVec(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	return &pdVecm1[1];
}

inline integer
MyVectorHandler::iGetSize(void) const
{
	return iCurSize;
}

inline void
MyVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
	/* Nota: il flag di ritorno e' pleonastico. Lo si e' messo per
	 * analogia con le matrici sparse, in cui l'aggiunta
	 * di un coefficiente puo' risultare in un errore.
	 * Qui, per motivi di efficienza, il controllo sulla validita'
	 * dell'indice e del vettore viene svolto solo in debug */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] = dCoef;
}

inline void
MyVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] += dCoef;
}

inline void
MyVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	pdVecm1[iRow] -= dCoef;
}

inline const doublereal&
MyVectorHandler::dGetCoef(integer iRow) const
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	return pdVecm1[iRow];
}

inline const doublereal&
MyVectorHandler::operator () (integer iRow) const
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	return pdVecm1[iRow];
}

inline doublereal&
MyVectorHandler::operator () (integer iRow)
{
	/* Vedi nota di PutCoef() */

#ifdef DEBUG
	IsValid();
	ASSERT((iRow > 0) && (iRow <= iCurSize));
#endif /* DEBUG */

	return pdVecm1[iRow];
}

/* MyVectorHandler - end */

#endif /* VH_H */

