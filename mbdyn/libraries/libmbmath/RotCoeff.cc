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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <myassert.h>
#include "RotCoeff.hh"

namespace RotCoeff {

    const doublereal SerCoeff[6][12]={
	{1.,	    // a
	 -6.,
	 120.,
	 -5040.,
	 362880.,
	 -39916800.,
	 6227020800.,
	 -1307674368000.,
	 355687428096000.,
	 0.,
	 0.,
	 0.},
	{2.,	    // b
	 -24.,
	 720.,
	 -40320.,
	 3628800.,
	 -479001600.,
	 87178291200.,
	 -20922789888000.,
	 6402373705728000.,
	 0.,
	 0.,
	 0.},
	{6.,	    // c
	 -120.,
	 5040.,
	 -362880.,
	 39916800.,
	 -6227020800.,
	 1307674368000.,
	 -355687428096000.,
	 121645100408832000.,
	 0.,
	 0.,
	 0.},
	{-12.,      // d
	 180.,
	 -6720.,
	 453600.,
	 -47900160.,
	 7264857600.,
	 -1494484992000.,
	 400148356608000.,
	 -135161222676480000.,
	 56200036388880384000.,
	 0.,
	 0.},
	{-60.,      // e
	 1260.,
	 -60480.,
	 4989600.,
	 -622702080.,
	 108972864000.,
	 -25406244864000.,
	 7602818775552000.,
	 -2838385676206080000.,
	 1292600836944248832000.,
	 0.,
	 0.},
	{90.,	    // f
	 -1680.,
	 75600.,
	 -5987520.,
	 726485760.,
	 -124540416000.,
	 28582025472000.,
	 -8447576417280000.,
	 3122224243826688000.,
	 0.,
	 0.,
	 0.}};


    const doublereal SerTrunc[10]={
				9,    // for a: phi^16
				9,    // for b: phi^16
				9,    // for c: phi^16
				9,    // for d: phi^16
				9,    // for e: phi^16
				9};   // for f: phi^16



    const doublereal SerThrsh[10]={
				1.1,		  // for a = coeff[0]
				1.3,		  // for b = coeff[1]
				1.5,		  // for c = coeff[2]
				1.6,		  // for d = coeff[3]
				1.7,		  // for e = coeff[4]
				1.8};		  // for f = coeff[5]

void
CoeffFun(
		const Int cid,
		const doublereal &phi,
		const doublereal &phi2,
		doublereal *const cf
		)
{
    ASSERT(phi != 0.);
    ASSERT(phi2 != 0.);

    switch (cid) {
	case 0:
            cf[0]=sin(phi)/phi;                 // a = sin(phi)/phi
            break;
        case 1:
            cf[1]=(1.-cos(phi))/phi2;           // b = (1.-cos(phi))/phi2
            break;
        case 2:
            cf[2]=(1.-cf[0])/phi2;              // c = (1.-a)/phi2
            break;
        case 3:
            cf[3]=(cf[0]-2.*cf[1])/phi2;        // d = (a-2*b)/phi2
            break;
        case 4:
            cf[4]=(cf[1]-3.*cf[2])/phi2;        // e = (b-3*c)/phi2
            break;
        case 5:
            cf[5]=(cf[2]-cf[1]-4.*cf[3])/phi2;  // f = (c-b-4*d)/phi2
            break;
        case 6:
            cf[6]=(cf[3]-5.*cf[4])/phi2;        // g = (d-5*e)/phi2
            break;
        case 7:
            cf[7]=(cf[4]-cf[3]-6.*cf[5])/phi2;  // h = (e-d-6*f)/phi2
            break;
        case 8:
            cf[8]=(cf[5]-7.*cf[6])/phi2;        // i = (f-7*g)/phi2
            break;
        case 9:
            cf[9]=(cf[6]-cf[5]-8.*cf[7])/phi2;  // j = (g-f-8*h)/phi2
            break;
	}
};

#if 1 /* molto piu' efficiente: evito n chiamate a unzione con switch ... */
void
RotCoeff(const Int cid, const doublereal &phi2, doublereal *const cf)
{
	ASSERT(cid >= 0 && cid <= 9);
	ASSERT(phi2 >= 0.);
	
	doublereal phi(sqrt(phi2)), phip[10];
	Int k, j;

	//P> Non si mette un 'if' in un ciclo ...
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
	// if (cid == 0) return; /* impossibile */
	cf[1]=(1.-cos(phi))/phi2;           // b = (1.-cos(phi))/phi2
	if (cid == 1) return;
	cf[2]=(1.-cf[0])/phi2;              // c = (1.-a)/phi2
	if (cid == 2) return;
	cf[3]=(cf[0]-2.*cf[1])/phi2;        // d = (a-2*b)/phi2
	if (cid == 3) return;
	cf[4]=(cf[1]-3.*cf[2])/phi2;        // e = (b-3*c)/phi2
	if (cid == 4) return;
	cf[5]=(cf[2]-cf[1]-4.*cf[3])/phi2;  // f = (c-b-4*d)/phi2
	if (cid == 5) return;
	cf[6]=(cf[3]-5.*cf[4])/phi2;        // g = (d-5*e)/phi2
	if (cid == 6) return;
	cf[7]=(cf[4]-cf[3]-6.*cf[5])/phi2;  // h = (e-d-6*f)/phi2
	if (cid == 7) return;
	cf[8]=(cf[5]-7.*cf[6])/phi2;        // i = (f-7*g)/phi2
	if (cid == 8) return;
	cf[9]=(cf[6]-cf[5]-8.*cf[7])/phi2;  // j = (g-f-8*h)/phi2
	// if (cid == 9) return; /* inutile */
};
#else
void
RotCoeff(const Int cid, const doublereal &phi2, doublereal *const coeff)
{
	ASSERT(cid >= 0 && cid <= 9);
	ASSERT(phi2 >= 0.);
	
	doublereal phi(sqrt(phi2)), phip[10];
	Int k, j;

	//P> Non si mette un 'if' in un ciclo ...
	if (phi < SerThrsh[cid-1]) {
		phip[0] = 1.;
		for (j = 1; j <= 9; j++) {
			phip[j] = phip[j-1]*phi2;
		}
		for (k = 0; k < cid; k++) {
			coeff[k] = 0;
			for (j = 0; j < SerTrunc[k]; j++) {
				coeff[k] += phip[j]/SerCoeff[k][j];
			}
		}
	} else {
		for (k = 0; k < cid; k++) {
			CoeffFun(k, phi, phi2, coeff);
		}
	}
};
#endif

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
