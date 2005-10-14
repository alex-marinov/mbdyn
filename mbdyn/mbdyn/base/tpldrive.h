/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

#ifndef TPLDRIVE_H
#define TPLDRIVE_H

#include "drive.h"
//#include "mbpar.h"


/* TplDriveCaller - begin */

template <class T>
class TplDriveCaller {
 public:
   virtual ~TplDriveCaller(void) { 
      NO_OP;
   };
   
   /* copia */
   virtual TplDriveCaller<T>* pCopy(void) const = 0;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const = 0;
   
   /* Restituisce il valore del driver */
   virtual T Get(void) const = 0;
};

/* TplDriveCaller - end */


/* TplDriveOwner - begin */

template <class T>
class TplDriveOwner {
 protected:
   TplDriveCaller<T>* pTplDriveCaller;
   
 public:
   TplDriveOwner(const TplDriveCaller<T>* pDC = NULL)
     : pTplDriveCaller((TplDriveCaller<T>*)pDC) { 
	NO_OP;
     };
   
   virtual ~TplDriveOwner(void) { 
      NO_OP;
   };
   
   void Set(const TplDriveCaller<T>* pDC) {
      ASSERT(pDC != NULL);
      ASSERT(pTplDriveCaller == NULL);
      pTplDriveCaller = (TplDriveCaller<T>*)pDC;
   };
   
   TplDriveCaller<T>* pGetDriveCaller(void) const { 
      return pTplDriveCaller;
   };
   
   T Get(void) const { 
      return pTplDriveCaller->Get(); 
   };
};

/* TplDriveOwner - end */

#endif /* TPLDRIVE_H */

