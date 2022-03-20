/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifndef ___EPETRA_VECTOR_HANDLER__INCLUDED___
#define ___EPETRA_VECTOR_HANDLER__INCLUDED___

#ifdef USE_TRILINOS
#include "vh.h"
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

class EpetraVectorHandler: public VectorHandler {
public:
     EpetraVectorHandler(integer iSize, const Epetra_Comm& oComm);
     virtual ~EpetraVectorHandler(void);

#ifdef DEBUG
     virtual void IsValid(void) const;
#endif
     virtual doublereal* pdGetVec(void) const;

     virtual integer iGetSize(void) const;

     virtual void Reset(void);

     virtual void Resize(integer iNewSize);

     virtual void ResizeReset(integer);

     virtual void PutCoef(integer iRow, const doublereal& dCoef);

     virtual void IncCoef(integer iRow, const doublereal& dCoef);

     virtual void DecCoef(integer iRow, const doublereal& dCoef);

     virtual const doublereal& dGetCoef(integer iRow) const;

     virtual const doublereal& operator () (integer iRow) const;

     virtual doublereal& operator () (integer iRow);

     virtual void Add(integer iRow, const Vec3& v);

     virtual void Sub(integer iRow, const Vec3& v);

     virtual void Put(integer iRow, const Vec3& v);

     virtual VectorHandler&
     ScalarAddMul(const VectorHandler& VH, const doublereal& d);

     virtual VectorHandler&
     ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
                  const doublereal& d);

     virtual VectorHandler&
     ScalarMul(const VectorHandler& VH, const doublereal& d);

     virtual VectorHandler& operator += (const VectorHandler& VH);

     virtual VectorHandler& operator += (const SubVectorHandler& SubVH);

     virtual VectorHandler& operator -= (const VectorHandler& VH);

     virtual VectorHandler& operator *= (const doublereal &d);

     virtual VectorHandler& operator = (const VectorHandler& VH);

     virtual doublereal Dot(void) const;

     virtual doublereal Norm(void) const;

     virtual doublereal InnerProd(const VectorHandler& VH) const;

     const Epetra_Vector* pGetEpetraVector() const { return &oEPV; }
     Epetra_Vector* pGetEpetraVector() { return &oEPV; }
private:
     Epetra_Vector oEPV;
};

#endif
#endif
