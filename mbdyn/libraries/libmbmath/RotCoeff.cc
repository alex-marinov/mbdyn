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

#include <myassert.h>
#include "RotCoeff.hh"

namespace RotCoeff {

    const doublereal SerCoeff[6][9]={
	{1.,	    // a
	 -6.,
	 120.,
	 -5040.,
	 362880.,
	 -39916800.,
	 6227020800.,
	 -1307674368000.,
	 355687428096000.},
	{2.,	    // b
	 -24.,
	 720.,
	 -40320.,
	 3628800.,
	 -479001600.,
	 87178291200.,
	 -20922789888000.,
	 6402373705728000.},
	{6.,	    // c
	 -120.,
	 5040.,
	 -362880.,
	 39916800.,
	 -6227020800.,
	 1307674368000.,
	 -355687428096000.,
	 121645100408832000.},
	{-12.,      // d
	 180.,
	 -6720.,
	 453600.,
	 -47900160.,
	 7264857600.,
	 -1494484992000.,
	 400148356608000.,
	 -135161222676480000.},
	{-60.,      // e
	 1260.,
	 -60480.,
	 4989600.,
	 -622702080.,
	 108972864000.,
	 -25406244864000.,
	 7602818775552000.,
	 -2838385676206080000.},
	{90.,	    // f
	 -1680.,
	 75600.,
	 -5987520.,
	 726485760.,
	 -124540416000.,
	 28582025472000.,
	 -8447576417280000.,
	 3122224243826688000.}};


    const doublereal SerTrunc[6]={
				9,    // for a: phi^16
				9,    // for b: phi^16
				9,    // for c: phi^16
				9,    // for d: phi^16
				9,    // for e: phi^16
				9};   // for f: phi^16



    const doublereal SerThrsh[6]={
				1.1,		  // for a = coeff[0]
				1.3,		  // for b = coeff[1]
				1.5,		  // for c = coeff[2]
				1.6,		  // for d = coeff[3]
				1.7,		  // for e = coeff[4]
				1.8};		  // for f = coeff[5]


void
RotCoeff(const Int cid, const doublereal &phi2, doublereal *const cf)
{
	ASSERT(cid >= 1 && cid <= 6);
	ASSERT(phi2 >= 0.);
	
	doublereal phi(sqrt(phi2)), phip[10];
	Int k, j;

	if (phi < SerThrsh[cid-1]) {
		phip[0] = 1.;
		for (j = 1; j <= 9; j++) {
			phip[j] = phip[j-1]*phi2;
		}
		for (k = 0; k < cid; k++) {
			cf[k] = 0;
			for (j = 0; j < SerTrunc[k]; j++) {
				cf[k] += phip[j]/SerCoeff[k][j];
			}
		}

		return;
	} 
	
	cf[0]=sin(phi)/phi;                 // a = sin(phi)/phi
	if (cid == 1) return;
	cf[1]=(1.-cos(phi))/phi2;           // b = (1.-cos(phi))/phi2
	if (cid == 2) return;
	cf[2]=(1.-cf[0])/phi2;              // c = (1.-a)/phi2
	if (cid == 3) return;
	cf[3]=(cf[0]-2.*cf[1])/phi2;        // d = (a-2*b)/phi2
	if (cid == 4) return;
	cf[4]=(cf[1]-3.*cf[2])/phi2;        // e = (b-3*c)/phi2
	if (cid == 5) return;
	cf[5]=(cf[2]-cf[1]-4.*cf[3])/phi2;  // f = (c-b-4*d)/phi2
	//if (cid == 6) return; inutile
	return;
};

// Coefficients:            up to a     (COEFF_A)
void CoeffA(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_A,phi2,coeff);
};

// Coefficients:            up to b     (COEFF_B)
void CoeffB(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_B,phi2,coeff);
};

// Coefficients:            up to c     (COEFF_C)
void CoeffC(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_C,phi2,coeff);
};

// Coefficients:            up to d     (COEFF_D)
void CoeffD(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_D,phi2,coeff);
};

// Coefficients:            up to e     (COEFF_E)
void CoeffE(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_E,phi2,coeff);
};

// Coefficients:            up to f     (COEFF_F)
void CoeffF(const doublereal &phi2, doublereal *const coeff) {
    RotCoeff(COEFF_F,phi2,coeff);
};

// Starred coefficients:    up to c*    (COEFF_C_STAR)
// Coefficients:            up to d     (COEFF_D)
void CoeffCStar(const doublereal &phi2, doublereal *const coeff, doublereal *const coeffs) {
    RotCoeff(COEFF_D,phi2,coeff);
    coeffs[0]=-coeff[3]/(2.*coeff[1]);
};

// Starred coefficients:    up to e*    (COEFF_E_STAR)
// Coefficients:            up to f     (COEFF_F)
void CoeffEStar(const doublereal &phi2, doublereal *const coeff, doublereal *const coeffs) {
    RotCoeff(COEFF_F,phi2,coeff);
    coeffs[0]=-coeff[3]/(2.*coeff[1]);
    coeffs[1]=-(coeff[4]+coeff[5])/(4*coeff[1]);
};



}
