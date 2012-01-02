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

#include "px.h"


/* PersistentExcitation - begin */

PersistentExcitation::PersistentExcitation(int i) 
: iNumDrives(i)
{
   NO_OP;
}


PersistentExcitation::~PersistentExcitation(void) 
{
   NO_OP;
}
   

int PersistentExcitation::iGetNumDrives(void) const 
{
   return iNumDrives;
}

/* PersistentExcitation - end */



NullPX::NullPX(void)
: PersistentExcitation(0)
{
   NO_OP;
}


NullPX::~NullPX(void)
{
   NO_OP;
}


void NullPX::AddInput(doublereal* pd) const
{
   NO_OP;
}


ScalarPX::ScalarPX(DriveCaller* p)
: PersistentExcitation(1), pDrive(NULL) 
{
   ASSERT(p != NULL);
   SAFENEWWITHCONSTRUCTOR(pDrive, DriveOwner, DriveOwner(p));
}


ScalarPX::~ScalarPX(void)
{
   SAFEDELETE(pDrive);
}


void ScalarPX::AddInput(doublereal* pd) const
{
   pd[0] += pDrive->dGet();
}


VectorPX::VectorPX(int i, DriveCaller** pp)
: PersistentExcitation(i), pvDrives(NULL) 
{
#ifdef DEBUG
   ASSERT(iNumDrives > 0);
   ASSERT(pp != NULL);
   for (int i = iNumDrives; i-- > 0; ) {
      ASSERT(pp[i] != NULL);
   }
#endif // DEBUG
   
   SAFENEWARR(pvDrives, DriveOwner*, iNumDrives);
   
   for (int i = iNumDrives; i-- > 0; ) {
      pvDrives[i] = NULL;
      SAFENEWWITHCONSTRUCTOR(pvDrives[i], DriveOwner, DriveOwner(pp[i]));
   }
   
   SAFEDELETEARR(pp);
}


VectorPX::~VectorPX(void)
{
   for (int i = iNumDrives; i-- > 0; ) {
      SAFEDELETE(pvDrives[i]);
   }
   SAFEDELETEARR(pvDrives);
}


void VectorPX::AddInput(doublereal* pd) const
{
   for (int i = iNumDrives; i-- > 0; ) {
      pd[i] += pvDrives[i]->dGet();
   }
}

