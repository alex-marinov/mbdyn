/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#ifndef ID_H
#define ID_H

#include "ac/f2c.h"
#include "myassert.h"
#include "ldl.h"
#include "forgfact.h"

const doublereal LDL_INIT = 1.e12;

/* Ident - begin */

class Ident {   
 protected:
   integer size;
   integer nout;
   
   doublereal* pdBase;
   doublereal** ppdBase;
   
   doublereal* ldl;
   doublereal** vldl;   
   
   doublereal* z;
   doublereal** vz;
   
   doublereal* theta;
   doublereal** vtheta;
   
   doublereal* phi;
   doublereal* y;
   doublereal* err;

   ForgettingFactor* pF;
   doublereal k;
   doublereal w;
   
   /* helper */
   void SetForgettingFactor(const doublereal& kk);
   
 public:
   Ident(integer size, integer nout, ForgettingFactor* pf,
	 const doublereal& ldl_init = LDL_INIT);
   virtual ~Ident(void);
   doublereal dGetForgettingFactor(void) const;   
   doublereal* pdGetTheta(void);
   void UpdateTheta(const doublereal* pd); // non usata
   doublereal* pdGetErr(void);
   void Update(const doublereal* pphi, const doublereal* yy);
};

/* Ident - end */


/* IdentProcess - begin */

class IdentProcess {
 protected:
   unsigned int iNumOutput;
   unsigned int iNumInput;
   unsigned int iOrdA;
   unsigned int iOrdB;

   Ident* pIdent; 
   
 public:
   IdentProcess(unsigned int iOut, unsigned int iIn, 
		unsigned int iA, unsigned int iB);
   
   virtual ~IdentProcess(void);
   
   void CreateIdent(integer size, integer nout, ForgettingFactor* pf);   
   virtual void Update(const doublereal* pdY, const doublereal* pdU) = 0;
   virtual void GetTheta(doublereal* pdTheta) = 0;
   virtual void GetErr(doublereal* pdE);
   
   virtual inline doublereal dGetForgettingFactor(void) {
      ASSERT(pIdent != NULL);
      return pIdent->dGetForgettingFactor();
   };   
   
   virtual integer iGetSize(void) const = 0;
   
   integer iGetNumOutput(void) const {
      return iNumOutput;
   };   
};

/* IdentProcess - end */


/* IdentARXProcess - begin */

class IdentARXProcess : public IdentProcess {
 protected:
   integer size;
   
   // e' il puntatore allo spazio di lavoro; usato solo per gestirlo
   doublereal* pdBase;
   
   doublereal* pdPhi;
   doublereal* pdY;
   
 public:
   IdentARXProcess(unsigned int iOut, unsigned int iIn,
		   unsigned int iA, unsigned int iB,
		   ForgettingFactor* pf);
   
   virtual ~IdentARXProcess(void);
   virtual void Update(const doublereal* pdYTmp, const doublereal* pdUTmp);
   virtual void GetTheta(doublereal* pdTheta);  

   virtual integer iGetSize(void) const {
      return size;
   };
};

/* IdentARXProcess - end */


/* IdentARMAXProcess - begin */

class IdentARMAXProcess : public IdentProcess {
 protected:
   integer size;
   
   // e' il puntatore allo spazio di lavoro; usato solo per gestirlo
   doublereal* pdBase;
   
   doublereal* pdPhi;
   doublereal* pdY;
   doublereal* pdErr;
   
   flag fCheckMA(doublereal* pdTheta);
   
 public:
   IdentARMAXProcess(unsigned int iOut, unsigned int iIn,
		     unsigned int iA, unsigned int iB,
		     ForgettingFactor* pf);
   
   virtual ~IdentARMAXProcess(void);
   virtual void Update(const doublereal* pdYTmp, const doublereal* pdUTmp);
   virtual void GetTheta(doublereal* pdTheta);
   void GetErr(doublereal* pdE);

   virtual integer iGetSize(void) const {
      return size;
   };
};

/* IdentARMAXProcess - end */

#endif
