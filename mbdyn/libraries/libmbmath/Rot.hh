/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 * Teodoro Merlini  <merlini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * This code is a partial merge of HmFe and MBDyn.
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


#ifndef Rot_hh
#define Rot_hh

namespace RotManip {

Mat3x3 Rot(const Vec3 & phi);
 
Mat3x3 DRot(const Vec3 & phi);

void RotAndDRot(const Vec3 & phi, Mat3x3 & Phi, Mat3x3 & Ga);

Mat3x3 DRot_IT(const Vec3 & phi);

void RotAndDRot_IT(const Vec3 & phi, Mat3x3 & PhiIT, Mat3x3 & GaIT);

Vec3 VecRot(const Mat3x3 & Phi);

} //end of namespace RotManip



namespace RoTrManip {

MatExp RoTr(const VecExp & eta);

MatExp DRoTr(const VecExp & eta);

void RoTrAndDRoTr(const VecExp & eta, MatExp & H, MatExp & Th);

MatExp DRoTr_It (const VecExp & eta);

void RoTrAndDRoTr_It(const VecExp & eta, MatExp & HIt, MatExp & ThIt);

VecExp Helix (const MatExp & H);

} //end of namespace RoTrManip
#endif // Rot_hh

