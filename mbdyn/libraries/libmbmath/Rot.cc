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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "matvecexp.h"
#include "RotCoeff.hh"
#include "Rot.hh"

using namespace RotCoeff;

Mat3x3 RotManip::Rot(const Vec3 & phi) {
	doublereal coeff[COEFF_B];

	CoeffB(phi.Dot(),coeff);
#if 0
	Mat3x3 Phi(1.);
	Mat3x3 phi_cross(phi);

	Phi+=phi_cross*coeff[0];
	Phi+=phi_cross*(phi_cross)*coeff[1];
#else /* 0 */
	Mat3x3 Phi(1., phi*coeff[0]);		/* I + c[0] * phi x */
	Phi += Mat3x3(phi, phi*coeff[1]);	/* += c[1] * phi x phi x */
#endif /* 0 */

	return Phi;
};


Mat3x3 RotManip::DRot(const Vec3 & phi) {
	doublereal coeff[COEFF_C];

	CoeffC(phi.Dot(),coeff);
#if 0
	Mat3x3 Ga(1.), phi_cross(phi);

	Ga+=phi_cross*coeff[1];
	Ga+=phi_cross*(phi_cross)*coeff[2];
#else /* 0 */
	Mat3x3 Ga(1., phi*coeff[1]);		/* I + c[0] * phi x */
	Ga += Mat3x3(phi, phi*coeff[2]);	/* += c[1] * phi x phi x */
#endif /* 0 */

	return Ga;
};


void RotManip::RotAndDRot(const Vec3 & phi, Mat3x3 & Phi, Mat3x3 & Ga) {
	doublereal coeff[COEFF_C];

	CoeffC(phi.Dot(),coeff);
#if 0
	Mat3x3 phi_cross(phi);
	Phi=Mat3x3(1.);
	Ga=Mat3x3(1.);

	Phi+=phi_cross*coeff[0];
	Phi+=phi_cross*(phi_cross)*coeff[1];
	Ga+=phi_cross*coeff[1];
	Ga+=phi_cross*(phi_cross)*coeff[2];
#else /* 0 */
	Phi = Mat3x3(1., phi*coeff[0]);
	Phi += Mat3x3(phi, phi*coeff[1]);

	Ga = Mat3x3(1., phi*coeff[1]);
	Ga += Mat3x3(phi, phi*coeff[2]);
#endif /* 0 */

	return;
};


Mat3x3 RotManip::DRot_IT(const Vec3 & phi) {
	doublereal coeff[COEFF_D], coeffs[COEFF_C_STAR];
	
	CoeffCStar(phi.Dot(),coeff,coeffs);
#if 0
	Mat3x3 GaIT(1.);
	Mat3x3 phi_cross(phi);

	GaIT+=phi_cross*0.5;
	GaIT+=phi_cross*(phi_cross)*coeffs[0];
#else /* 0 */
	Mat3x3 GaIT(1., phi*.5);
	GaIT += Mat3x3(phi, phi*coeffs[0]);
#endif /* 0 */

	return GaIT;
};


void RotManip::RotAndDRot_IT(const Vec3 & phi, Mat3x3 & PhiIT, Mat3x3 & GaIT) {
	doublereal coeff[COEFF_D], coeffs[COEFF_C_STAR];

	CoeffCStar(phi.Dot(),coeff,coeffs);
#if 0
	PhiIT=Mat3x3(1.);
	GaIT=Mat3x3(1.);
	Mat3x3 phi_cross(phi);

	PhiIT+=phi_cross*coeff[0];
	PhiIT+=phi_cross*(phi_cross)*coeff[1];
	GaIT+=phi_cross*0.5;
	GaIT+=phi_cross*(phi_cross)*coeffs[0];
#else /* 0 */
	PhiIT = Mat3x3(1., phi*coeff[0]);
	PhiIT += Mat3x3(phi, phi*coeff[1]);

	GaIT = Mat3x3(1., phi*.5);
	GaIT += Mat3x3(phi, phi*coeffs[0]);
#endif /* 0 */

	return;
};

Vec3 RotManip::VecRot(const Mat3x3 & Phi) {
	doublereal phi, a, cosphi, sinphi;
	Vec3 unit;

	cosphi = .5*(Phi.Trace()-1.);
	if (cosphi > 0.) {
		unit = Phi.Ax();
		sinphi = unit.Norm();
		phi = atan2(sinphi, cosphi);
		CoeffA(phi*phi, &a);
		unit = unit/a;
	} else {
		Mat3x3 eet(Phi.Symm());
		eet -= Mat3x3(cosphi);
		eet /= (1.-cosphi);
		Int maxcol = 1;
		if ((eet.dGet(2, 2) > eet.dGet(1, 1))
				||(eet.dGet(3, 3)>eet.dGet(1, 1))) {
			maxcol = 2;
			if (eet.dGet(3, 3) > eet.dGet(2, 2)) {
				maxcol = 3;
			}
		}
		//M> spero che GetVec mi dia la colonna e non la riga!!!!
		//P> SI.
		unit = (eet.GetVec(maxcol)/sqrt(eet.dGet(maxcol, maxcol)));
		sinphi = (Mat3x3(unit)*Phi).Trace()*(-.5);
		unit *= atan2(sinphi, cosphi);
	}
	return unit;
};

MatExp RoTrManip::RoTr(const VecExp & eta) {
	doublereal phi2, coeff[COEFF_D];
	MatExp H(1., 0.);

	MatExp etaCross(eta.Cross());
	MatExp etaCross2(etaCross*etaCross);
	MatExp etaCross3(etaCross2*etaCross);
	MatExp etaCross4(etaCross3*etaCross);

	//P> Dot senza argomenti viene eseguito sempre su se stesso ...
	phi2 = eta.Vec().Dot();
	CoeffD(phi2, coeff);

	H += etaCross*(coeff[0]-.5*phi2*(coeff[2]-coeff[1]));
	H += etaCross2*(coeff[1]-.5*phi2*coeff[3]);
	H += etaCross3*(-.5*(coeff[2]-coeff[1]));
	H += etaCross4*(-.5*coeff[3]);

	return H;
};


MatExp RoTrManip::DRoTr(const VecExp & eta) {
	doublereal phi2, coeff[COEFF_E];
	MatExp Th(1., 0.);

	MatExp etaCross(eta.Cross());
	MatExp etaCross2(etaCross*etaCross);
	MatExp etaCross3(etaCross2*etaCross);
	MatExp etaCross4(etaCross3*etaCross);

	phi2 = eta.Vec().Dot();
	CoeffE(phi2, coeff);

	Th += etaCross*(coeff[1]-.5*phi2*coeff[3]);
	Th += etaCross2*(coeff[2]-.5*phi2*coeff[4]);
	Th += etaCross3*(-.5*coeff[3]);
	Th += etaCross4*(-.5*coeff[4]);
	return Th;
};


void RoTrManip::RoTrAndDRoTr(const VecExp & eta,
				MatExp & H,
				MatExp & Th) {
	doublereal phi2, coeff[COEFF_E];
	H = MatExp(1., 0.);
	Th = MatExp(1., 0.);

	MatExp etaCross(eta.Cross());
	MatExp etaCross2(etaCross*etaCross);
	MatExp etaCross3(etaCross2*etaCross);
	MatExp etaCross4(etaCross3*etaCross);

	phi2 = eta.Vec().Dot();
	CoeffE(phi2, coeff);

	H += etaCross*(coeff[0]-.5*phi2*(coeff[2]-coeff[1]));
	H += etaCross2*(coeff[1]-.5*phi2*coeff[3]);
	H += etaCross3*(-.5*(coeff[2]-coeff[1]));
	H += etaCross4*(-.5*coeff[3]);
	Th += etaCross*(coeff[1]-.5*phi2*coeff[3]);
	Th += etaCross2*(coeff[2]-.5*phi2*coeff[4]);
	Th += etaCross3*(-.5*coeff[3]);
	Th += etaCross4*(-.5*coeff[4]);
	return;
};


MatExp RoTrManip::DRoTr_It(const VecExp & eta) {
	doublereal phi2, coeff[COEFF_F], coeffs[COEFF_E_STAR];
	MatExp ThIt(1., 0.);

	MatExp etaCross(eta.Cross());
	MatExp etaCross2(etaCross*etaCross);
	MatExp etaCross4(etaCross2*etaCross2);

	phi2 = eta.Vec().Dot();
	CoeffEStar(phi2, coeff, coeffs);

	ThIt += etaCross*.5;
	ThIt += etaCross2*(coeffs[0]-.5*phi2*coeffs[1]);
	ThIt += etaCross4*(-.5*coeffs[1]);
	return ThIt;
};


void RoTrManip::RoTrAndDRoTr_It(const VecExp & eta,
					MatExp & HIt,
					MatExp & ThIt) {
	doublereal phi2, coeff[COEFF_F], coeffs[COEFF_E_STAR];
	HIt = MatExp(1., 0.);
	ThIt = MatExp(1., 0.);

	MatExp etaCross(eta.Cross());
	MatExp etaCross2(etaCross*etaCross);
	MatExp etaCross3(etaCross2*etaCross);
	MatExp etaCross4(etaCross3*etaCross);

	phi2 = eta.Vec().Dot();
	CoeffEStar(phi2, coeff, coeffs);

	HIt += etaCross*(coeff[0]-.5*phi2*(coeff[2]-coeff[1]));
	HIt += etaCross2*(coeff[1]-.5*phi2*coeff[3]);
	HIt += etaCross3*(-.5*(coeff[2]-coeff[1]));
	HIt += etaCross4*(-.5*coeff[3]);
	ThIt += etaCross*.5;
	ThIt += etaCross2*(coeffs[0]-.5*phi2*coeffs[1]);
	ThIt += etaCross4*(-.5*coeffs[1]);
	return;
};

VecExp RoTrManip::Helix(const MatExp & H)  {
	Vec3 phi(RotManip::VecRot(H.Vec()));
	return VecExp(phi,
			RotManip::DRot_IT(phi).Transpose()*(
				(H.Mom()*(H.Vec().Transpose())).Ax())
			);
};

