/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2012
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
 * Copyright (C) 1996-2012
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

#ifndef RotCoeff_hh
#define RotCoeff_hh

#include <ac/f2c.h>

#include <matvec3.h>
#include <matvecexp.h>

namespace RotCoeff {

const Int COEFF_A = 1;
const Int COEFF_B = 2;
const Int COEFF_C = 3;
const Int COEFF_D = 4;
const Int COEFF_E = 5;
const Int COEFF_F = 6;

const Int COEFF_C_STAR = 1;
const Int COEFF_E_STAR = 2;


template<class T1, class T2> void CoeffA(const T1 &phi, const Vec3 &p, T2 *const coeff);
template<class T1, class T2> void CoeffB(const T1 &phi, const Vec3 &p, T2 *const coeff);
template<class T1, class T2> void CoeffC(const T1 &phi, const Vec3 &p, T2 *const coeff);
template<class T1, class T2> void CoeffD(const T1 &phi, const Vec3 &p, T2 *const coeff);
template<class T1, class T2> void CoeffE(const T1 &phi, const Vec3 &p, T2 *const coeff);
template<class T1, class T2> void CoeffF(const T1 &phi, const Vec3 &p, T2 *const coeff);

template<class T1, class T2> void CoeffCStar(const T1 &phi, const Vec3 &p,
				T2 *const coeff, T2 *const coeffs);
template<class T1, class T2> void CoeffEStar(const T1 &phi, const Vec3 &p,
				T2 *const coeff, T2 *const coeffs);


}

#include "RotCoeff.hc"

#endif // RotCoeff_hh

