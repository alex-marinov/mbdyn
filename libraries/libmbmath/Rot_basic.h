#ifndef Rot_basic_h
#define Rot_basic_h


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

#endif // Rot_basic_h
