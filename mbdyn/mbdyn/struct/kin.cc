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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <kin.h>

KinematicsTest::KinematicsTest(const DriveCaller* pDC)
: DriveOwner(pDC)
{
   NO_OP;
}
 

KinematicsTest::~KinematicsTest(void) 
{
   NO_OP;
}
   
const Vec3& 
KinematicsTest::GetXCurr(void) const 
{
   return (Vec3&)X = Vec3(dGet(), 0., 0.);
}


const Mat3x3& 
KinematicsTest::GetRCurr(void) const 
{
   doublereal c = cos(.01*dGet());
   doublereal s = sin(.01*dGet());
   return (Mat3x3&)R = Mat3x3(c, s, 0., -s, c, 0., 0., 0., 1.);
}
 

const Vec3& 
KinematicsTest::GetVCurr(void) const 
{
   return (Vec3&)V = Vec3(1., 0., 0.);
}


const Vec3& 
KinematicsTest::GetWCurr(void) const 
{
   return (Vec3&)W = Vec3(0., 0., .01);
}
