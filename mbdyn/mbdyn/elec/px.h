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

#ifndef PX_H
#define PX_H

#include "drive.h"

class PersistentExcitation {
 protected:
   int iNumDrives;
   
 public:
   PersistentExcitation(int i);
   virtual ~PersistentExcitation(void);
   
   virtual int iGetNumDrives(void) const;
   virtual void AddInput(doublereal* pd) const = 0;
};

class NullPX : public PersistentExcitation {
 public:
   NullPX(void);
   virtual ~NullPX(void);
   
   virtual void AddInput(doublereal* pd) const;
};

class ScalarPX : public PersistentExcitation {
 protected:
   DriveOwner* pDrive;
   
 public:
   ScalarPX(DriveCaller* p);
   virtual ~ScalarPX(void);
   
   virtual void AddInput(doublereal* pd) const;
};

class VectorPX : public PersistentExcitation {
 protected:
   DriveOwner** pvDrives;
   
 public:
   VectorPX(int i, DriveCaller** p);
   virtual ~VectorPX(void);
   
   virtual void AddInput(doublereal* pd) const;
};

#endif // PX_H
