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

#ifndef RotCoeff_hh
#define RotCoeff_hh

#include "myf2c.h"
#include "mystddef.h"

namespace RotCoeff {

const Int COEFF_A = 1;
const Int COEFF_B = 2;
const Int COEFF_C = 3;
const Int COEFF_D = 4;
const Int COEFF_E = 5;
const Int COEFF_F = 6;

const Int COEFF_C_STAR = 1;
const Int COEFF_E_STAR = 2;


void CoeffA(const doublereal &phi2, doublereal *const coeff);
void CoeffB(const doublereal &phi2, doublereal *const coeff);
void CoeffC(const doublereal &phi2, doublereal *const coeff);
void CoeffD(const doublereal &phi2, doublereal *const coeff);
void CoeffE(const doublereal &phi2, doublereal *const coeff);
void CoeffF(const doublereal &phi2, doublereal *const coeff);

void CoeffCStar(const doublereal &phi2,doublereal *const coeff, doublereal *const coeffs);
void CoeffEStar(const doublereal &phi2,doublereal *const coeff, doublereal *const coeffs);


}

#endif // RotCoeff_hh

