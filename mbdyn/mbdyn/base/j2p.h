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

#ifndef J2P_H
#define J2P_H

#include "node.h"
#include "elem.h"


/* Elem2Param - begin */

class Elem2Param : public ParameterNode {
 protected:   
   Elem* pElem;  
   unsigned int iNum;
   
 public:
   Elem2Param(unsigned int uL, const DofOwner* pDO, flag fOut);
   virtual ~Elem2Param(void);
   
   virtual void Bind(const Elem* pEl, unsigned int i);
   
   /* Contributo del nodo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   /* Restituisce il valore del dof iDof;
    * se differenziale, iOrder puo' essere = 1 per la derivata */
   virtual inline const doublereal& 
     dGetDofValue(int iDof, int iOrder = 0) const;

   /* virtual void SetX(const doublereal& d); */
   virtual inline const doublereal& dGetX(void) const;
   
   /* Setta il valore del dof iDof a dValue;
    * se differenziale, iOrder puo' essere = 1 per la derivata */
   virtual void SetDofValue(const doublereal& dValue,
			    unsigned int iDof, 
			    unsigned int iOrder = 0);
};


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
#ifdef DEBUG
inline const doublereal& Elem2Param::dGetDofValue(int iDof, int iOrder) const
#else
inline const doublereal& Elem2Param::dGetDofValue(int /* iDof */ , int /* iOrder */ ) const
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0);
   return dGetX();
}

/* Restituisce il valore del dof */
inline const doublereal& Elem2Param::dGetX(void) const
{
   return ((doublereal&)dX = pElem->dGetPrivData(iNum));
}

/* Elem2Param - end */


/* StrainGageParam - begin */

class StrainGageParam : public Elem2Param {
 protected:  
   doublereal dY;
   doublereal dZ;
   
 public:
   StrainGageParam(unsigned int uL, const DofOwner* pDO,
		   doublereal dy, doublereal dz, flag fOut);
   virtual ~StrainGageParam(void);
   
   virtual void Bind(const Elem* pEl, unsigned int i);
   
   /* Contributo del nodo al file di restart */
   virtual ostream& Restart(ostream& out) const;

   /* Restituisce il valore del dof iDof;
    * se differenziale, iOrder puo' essere = 1 per la derivata */
   virtual inline const doublereal& 
     dGetDofValue(int iDof, int iOrder = 0) const;

   /* virtual void SetX(const doublereal& d); */
   virtual inline const doublereal& dGetX(void) const;   
   
   /* Setta il valore del dof iDof a dValue;
    * se differenziale, iOrder puo' essere = 1 per la derivata
    * (usa quella di Elem2Param)
   virtual void SetDofValue(const doublereal& dValue,
			    unsigned int iDof, 
			    unsigned int iOrder = 0);
    */
};


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
inline const doublereal&
#ifdef DEBUG
StrainGageParam::dGetDofValue(int iDof, int iOrder) const
#else
StrainGageParam::dGetDofValue(int /* iDof */ , int /* iOrder */ ) const
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0);
   return dGetX();
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
inline const doublereal& StrainGageParam::dGetX(void) const
{  
   unsigned int i = iNum-1;
   return ((doublereal&)dX = pElem->dGetPrivData(6*i+1)
           +dZ*pElem->dGetPrivData(6*i+5)
	   -dY*pElem->dGetPrivData(6*i+6));
}

/* StrainGageParam - end */

#endif // J2P_H
