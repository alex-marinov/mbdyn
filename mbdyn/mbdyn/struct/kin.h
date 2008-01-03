/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifndef KIN_H
# define KIN_H

# include <matvec3.h>
# include <drive.h>

/* Kinematics - begin */

class Kinematics {
 public:
   virtual ~Kinematics(void) {
      NO_OP;
   };
   
   virtual const Vec3& GetXCurr(void) const = 0;
   virtual const Mat3x3& GetRCurr(void) const = 0;
   virtual const Vec3& GetVCurr(void) const = 0;
   virtual const Vec3& GetWCurr(void) const = 0;
};

/* Kinematics - end */


/* KinematicsTest - begin */

class KinematicsTest : public Kinematics, public DriveOwner {
 protected:
   Vec3 X;
   Mat3x3 R;
   Vec3 V;
   Vec3 W;
   
 public:   
   KinematicsTest(const DriveCaller* pDC);   
   ~KinematicsTest(void);
   
   const Vec3& GetXCurr(void) const;
   const Mat3x3& GetRCurr(void) const;
   const Vec3& GetVCurr(void) const;
   const Vec3& GetWCurr(void) const;
};

/* KinematicsTest - end */

#endif /* KIN_H */
