/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2014
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
 * Copyright (C) 1996-2014
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


#ifndef Rot_hh
#define Rot_hh

#include "matvecexp.h"

namespace RotManip {

/**
 * Compute the rotation matrix Phi given Euler Rogriguez's parameters phi
 */
Mat3x3 Rot(const Vec3 & phi);

/**
 * Compute G matrix given Euler Rogriguez's parameters Phi
 * G defined in such a way that dPhi * PhiT = G * dphi
 */
 Mat3x3 DRot(const Vec3 & phi);

/**
 * Compute rotation matrix Phi and Ga matrix 
 * given Euler Rogriguez's parameters Phi
 */
void RotAndDRot(const Vec3 & phi, Mat3x3 & Phi, Mat3x3 & Ga);

/**
 * Compute the inverse transpose of G matrix given Euler Rogriguez's parameters Phi
 */
Mat3x3 DRot_IT(const Vec3 & phi);

/**
 * Compute the inverse of G matrix given Euler Rogriguez's parameters Phi
 */
Mat3x3 DRot_I(const Vec3 & phi);

/**
 * Compute inverse transpose ot rotation matrix Phi and Ga matrix 
 * given Euler Rogriguez's parameters Phi
 */
void RotAndDRot_IT(const Vec3 & phi, Mat3x3 & PhiIT, Mat3x3 & GaIT);

/**
 * Compute Euler Rogriguez's parameters phi given rotation matrix Phi
 */
Vec3 VecRot(const Mat3x3 & Phi);

/**
 * Compute, given Euler Rogriguez's parameters phi, L matrix such that
 * dG * a = L(phi, a) * dphi
 */
Mat3x3 Elle(const Vec3 & phi, const Vec3 & a);

} //end of namespace RotManip


/*
 * Euler-Rodrigues rotation manipulation namespace
 */
namespace ER_Rot {

class Param_Manip : public Vec3_Manip {
public:
	inline Vec3 operator << (const Mat3x3& m) const {
		return RotManip::VecRot(m);
  	};

	inline void Manipulate(Vec3& v, const Mat3x3& m) const {
		v = RotManip::VecRot(m);
	};
};

class MatR_Manip : public Mat3x3_Manip {   
public:
	inline Mat3x3 operator << (const Vec3& v) const {
		return RotManip::Rot(v);
	};
   
	inline void Manipulate(Mat3x3& m, const Vec3& v) const {
		m = RotManip::Rot(v);
	};
};

class MatG_Manip : public Mat3x3_Manip {
public:
	inline Mat3x3 operator << (const Vec3& v) const {
		return RotManip::DRot(v);
	};

	inline void Manipulate(Mat3x3& m, const Vec3& v) const {
		m = RotManip::DRot(v);
	};
};

class MatGm1_Manip : public Mat3x3_Manip {
public:
	inline Mat3x3 operator << (const Vec3& v) const {
		return RotManip::DRot_I(v);
	};
   
	inline void Manipulate(Mat3x3& m, const Vec3& v) const {
		m = RotManip::DRot_I(v);
	};
};

extern Param_Manip Param;
extern MatR_Manip MatR;
extern MatG_Manip MatG;
extern MatGm1_Manip MatGm1;

} // end of namespace ER_Rot



namespace RoTrManip {

MatExp Elle(const VecExp & phi, const VecExp & a);

MatExp RoTr(const VecExp & eta);

MatExp DRoTr(const VecExp & eta);

void RoTrAndDRoTr(const VecExp & eta, MatExp & H, MatExp & Th);

MatExp DRoTr_It (const VecExp & eta);

MatExp DRoTr_I (const VecExp & eta);

void RoTrAndDRoTr_It(const VecExp & eta, MatExp & HIt, MatExp & ThIt);

VecExp Helix (const MatExp & H);

} //end of namespace RoTrManip
#endif // Rot_hh

