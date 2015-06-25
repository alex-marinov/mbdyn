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

#ifndef FORGFACT_H
#define FORGFACT_H

#include <cfloat>
#include <cmath>
#include <ac/f2c.h>

#include <myassert.h>
#include <mynewmem.h>


/* Takes care of the forgetting factor */

class ForgettingFactor {
 protected:
   integer iNumErr;
 public:
   ForgettingFactor(integer i);
   virtual ~ForgettingFactor(void);
     
   virtual void Update(const doublereal* pErr = NULL) = 0;
   virtual doublereal dGet(void) const = 0;
};


class ConstForgettingFactor : public ForgettingFactor {
 protected:
   doublereal dk;
   
 public:
   ConstForgettingFactor(doublereal d);
   ~ConstForgettingFactor(void);
   
   void Update(const doublereal* /* pErr */ = NULL);
   inline doublereal dGet(void) const;
};


inline doublereal ConstForgettingFactor::dGet(void) const 
{
   return dk;
}


class DynamicForgettingFactor : public ForgettingFactor {
 protected:
   integer N1;
   integer N2;     
   
   doublereal dRho;
   doublereal dFact;
   doublereal dkRef;
   doublereal dkLim;
   
   doublereal dk;
   
   doublereal* pdErrM;
   doublereal* pdErrS;
   integer iRef1;
   integer iRef2;
   
   doublereal dErr1M;
   doublereal dErr1S;
   doublereal dErr2M;
   doublereal dErr2S;   
   
 public:
   DynamicForgettingFactor(integer n1, integer n2, integer i,
			    doublereal r, doublereal f,
			    doublereal kr, doublereal kl);
   ~DynamicForgettingFactor(void);
   
   void Update(const doublereal* pErr = NULL);
   inline doublereal dGet(void) const;
};


inline doublereal DynamicForgettingFactor::dGet(void) const 
{
   return dk;
}


class DynamicForgettingFactor2 : public ForgettingFactor {
 protected:
   integer N1;
   integer N2;     
   
   doublereal dRho;
   doublereal dFact;
   doublereal dkRef;
   doublereal dkLim;
   
   doublereal dk;
   
   doublereal* pdErr;
   doublereal** ppdErr;
   integer iRef1;
   integer iRef2;
          
   doublereal* pdErr1M;
   doublereal* pdErr1S;
   doublereal* pdErr2M;
   doublereal* pdErr2S; 
   
 public:
   DynamicForgettingFactor2(integer n1, integer n2, integer i,
			    doublereal r, doublereal f,
			    doublereal kr, doublereal kl);
   ~DynamicForgettingFactor2(void);
   
   void Update(const doublereal* pErr = NULL);
   inline doublereal dGet(void) const;
};


inline doublereal DynamicForgettingFactor2::dGet(void) const
{
   return dk;
}


#endif
