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

/* Drive con Vec3 */

#ifndef VEC3DRV_H
#define VEC3DRV_H

#include "drive.h"
#include "matvec3.h"

class Vec3DriveOwner : public DriveOwner {
 protected:
   const Vec3 Dir;
   
 public:
   Vec3DriveOwner(const DriveCaller* pDC = NULL, const Vec3& V = Vec3())
     : DriveOwner(pDC), Dir(V) {
      NO_OP;
   };
   
   virtual ~Vec3DriveOwner(void) {
      NO_OP;
   };
   
   virtual ostream& Restart(ostream& out) const {
      out << "reference, local, eye, ", Dir.Write(out, ", ") << ", ";
      return pGetDriveCaller()->Restart(out);
   };
   
   virtual Vec3 GetVec(void) const {
      return Dir*dGet();
   };
   
   virtual Vec3 GetDir(void) const {
      return Dir;
   };
};

#endif
