/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2010-2017
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
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

/*
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "constltp_impl.h"
#include "tpldrive_impl.h"
#include "shell.h"
#include "mynewmem.h"


Shell::Shell(unsigned uLabel, const DofOwner* pDO, flag fOut)
: Elem(uLabel, fOut),
ElemGravityOwner(uLabel, fOut),
ElemWithDofs(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut)
{
	NO_OP;
}

Shell::~Shell(void)
{
	NO_OP;
}

int
ReadShellConstLaw(MBDynParser& HP, Shell::fmh& pD, Shell::vh& PreStress)
{
	ASSERT(pD.iGetNumRows() == 12);
	ASSERT(pD.iGetNumCols() == 12);

	if (HP.IsKeyWord("diag")) {
		pD.Reset();
		for (unsigned ir = 1; ir <= 12; ir++) {
			pD(ir, ir) = HP.GetReal();
		}

	} else if (HP.IsKeyWord("sym")) {
		for (unsigned ir = 1; ir <= 12; ir++) {
			pD(ir, ir) = HP.GetReal();
			for (unsigned ic = ir + 1; ic <= 12; ic++) {
				doublereal d = HP.GetReal();
				pD(ir, ic) = d;
				pD(ic, ir) = d;
			}
		}

	} else if (HP.IsKeyWord("isotropic")) {
		doublereal dE;
		doublereal dnu;
		doublereal dG;
		doublereal dh;
		doublereal das = 1.;
		doublereal dat = .01;
		bool bGot_E(false);
		bool bGot_nu(false);
		bool bGot_G(false);
		bool bGot_h(false);
		bool bGot_as(false);
		bool bGot_at(false);

		while (HP.IsArg()) {
			if (HP.IsKeyWord("E") || HP.IsKeyWord("Young" "modulus")) {
				if (bGot_E) {
					silent_cerr("Shell isotropic constitutive law: Young's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_E = true;
				dE = HP.GetReal();
				if (dE <= 0.) {
					silent_cerr("Shell isotropic constitutive law: invalid Young's modulus " << dE << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("nu") || HP.IsKeyWord("Poisson" "modulus")) {
				if (bGot_nu) {
					silent_cerr("Shell isotropic constitutive law: Poisson's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_nu = true;
				dnu = HP.GetReal();
				if (dnu < 0. || dnu >= .5) {
					silent_cerr("Shell isotropic constitutive law: invalid Poisson's modulus " << dnu << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("G") || HP.IsKeyWord("shear" "modulus")) {
				if (bGot_G) {
					silent_cerr("Shell isotropic constitutive law: shear modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_G = true;
				dG = HP.GetReal();
				if (dG <= 0.) {
					silent_cerr("Shell isotropic constitutive law: invalid shear modulus " << dG << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("thickness")) {
				if (bGot_h) {
					silent_cerr("Shell isotropic constitutive law: thickness already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_h = true;
				dh = HP.GetReal();
				if (dh <= 0.) {
					silent_cerr("Shell isotropic constitutive law: invalid thickness " << dh << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("as") /* better name? */ ) {
				if (bGot_as) {
					silent_cerr("Shell isotropic constitutive law: as (?) already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_as = true;
				das = HP.GetReal();
				if (das <= 0.) {
					silent_cerr("Shell isotropic constitutive law: invalid as " << das << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("at") /* better name? */ ) {
				if (bGot_at) {
					silent_cerr("Shell isotropic constitutive law: at (?) already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_at = true;
				dat = HP.GetReal();
				if (dat <= 0.) {
					silent_cerr("Shell isotropic constitutive law: invalid at " << dat << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else {
				break;
			}
		}

		int got = 0;
		if (bGot_E) got++;
		if (bGot_nu) got++;
		if (bGot_G) got++;

		if (got < 2) {
			silent_cerr("Shell isotropic constitutive law: incomplete material data (need at least two among Young's modulus, Poisson's modulus, shear modulus)" << std::endl);
			return -1;
		}

		if (got > 2) {
			if (std::abs(dE/(2.*(1 + dnu)) - dG) >
				std::numeric_limits<doublereal>::epsilon())
			{
				silent_cerr("Shell isotropic constitutive law: inconsistent material data" << std::endl);
				return -1;
			}
		} else {
			if (!bGot_G) {
				dG = dE/(2.*(1 + dnu));

			} else if (!bGot_E) {
				dE = 2.*(1 + dnu)*dG;

			} else if (!bGot_nu) {
				dnu = dE/(2.*dG) - 1.;
				if (dnu <= 0. || dnu >= .5) {
					silent_cerr("Shell isotropic constitutive law: inconsistent material data (computed Poisson's modulus: " << dnu << ")" << std::endl);
					return -1;
				}
			}
		}

		if (!bGot_h) {
			silent_cerr("Shell isotropic constitutive law: shell thickness missing" << std::endl);
			return -1;
		}

		doublereal C = dE/(1. - dnu*dnu)*dh;
		doublereal D = C*dh*dh/12;
		doublereal G = dG*dh;
		doublereal F = G*dh*dh/12;

		pD.Reset();

		pD(1, 1) = C;
		pD(1, 5) = C*dnu;
		pD(2, 2) = 2*G;
		pD(3, 3) = G*das;
		pD(4, 4) = 2*G;
		pD(5, 1) = C*dnu;
		pD(5, 5) = C;
		pD(6, 6) = G*das;

		pD(7, 7) = 2*F;
		pD(8, 8) = D;
		pD(8, 10) = -D*dnu;
		pD(9, 9) = F*dat;
		pD(10, 10) = D;
		pD(10, 8) = -D*dnu;
		pD(11, 11) = 2*F;
		pD(12, 12) = F*dat;

	} else if (HP.IsKeyWord("plane" "stress" "orthotropic")) {
/*
Eshell =
[              h/(1-nu_lt^2/E_l*E_t)*E_l,                                      0,                                      0,                                      0,        h/(1-nu_lt^2/E_l*E_t)*nu_lt*E_t,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[                                      0,                                  2*h*G,                                      0,                                  2*h*G,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[                                      0,                                      0,                                 h*G*as,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[                                      0,                                  2*h*G,                                      0,                                  2*h*G,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[        h/(1-nu_lt^2/E_l*E_t)*nu_lt*E_t,                                      0,                                      0,                                      0,              h/(1-nu_lt^2/E_l*E_t)*E_t,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                 h*G*as,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,       1/12*h^3/(1-nu_lt^2/E_l*E_t)*E_l,                                      0,                                      0,                                      0, 1/12*h^3/(1-nu_lt^2/E_l*E_t)*nu_lt*E_t,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                              1/6*h^3*G,                                      0,                              1/6*h^3*G,                                      0,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                          1/12*h^3*G*dat,                                      0,                                      0,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                              1/6*h^3*G,                                      0,                              1/6*h^3*G,                                      0,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0, 1/12*h^3/(1-nu_lt^2/E_l*E_t)*nu_lt*E_t,                                      0,                                      0,                                      0,       1/12*h^3/(1-nu_lt^2/E_l*E_t)*E_t,                                      0]
[                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                                      0,                          1/12*h^3*G*dat]
 
*/

		/* TODO:
		     - check singularity in shear coefficients
		     - add (optional) re-orientation?
		 */

		doublereal dE_l;
		doublereal dE_t;
		doublereal dnu_lt;
		doublereal dnu_tl;
		doublereal dG;
		doublereal dBeta;
		doublereal dh;
		doublereal das = 1.;
		doublereal dat = .01;
		bool bGot_E_l(false);
		bool bGot_E_t(false);
		bool bGot_nu_lt(false);
		bool bGot_nu_tl(false);
		bool bGot_G(false);
		bool bGot_Beta(false);
		bool bGot_h(false);
		bool bGot_as(false);
		bool bGot_at(false);

		while (HP.IsArg()) {
			if (HP.IsKeyWord("E_l") || HP.IsKeyWord("longitudinal" "Young" "modulus")) {
				if (bGot_E_l) {
					silent_cerr("Shell plane stress orthotropic constitutive law: longitudinal Young's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_E_l = true;
				dE_l = HP.GetReal();
				if (dE_l <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid longitudinal Young's modulus " << dE_l << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("E_t") || HP.IsKeyWord("transverse" "Young" "modulus")) {
				if (bGot_E_t) {
					silent_cerr("Shell plane stress orthotropic constitutive law: transverse Young's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_E_t = true;
				dE_t = HP.GetReal();
				if (dE_t <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid transverse Young's modulus " << dE_t << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("nu_lt") || HP.IsKeyWord("longitudinal" "transverse" "Poisson" "modulus")) {
				if (bGot_nu_lt) {
					silent_cerr("Shell plane stress orthotropic constitutive law: longitudinal transverse Poisson's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_nu_lt = true;
				dnu_lt = HP.GetReal();
#if 0
				if (dnu_lt <= 0. || dnu_lt >= .5) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid longitudinal transverse Poisson's modulus " << dnu_lt << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}
#endif

			} else if (HP.IsKeyWord("nu_tl") || HP.IsKeyWord("transverse" "longitudinal" "Poisson" "modulus")) {
				if (bGot_nu_lt) {
					silent_cerr("Shell plane stress orthotropic constitutive law: transverse longitudinal Poisson's modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_nu_tl = true;
				dnu_tl = HP.GetReal();
#if 0
				if (dnu_tl <= 0. || dnu_tl >= .5) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid transverse longitudinal Poisson's modulus " << dnu_lt << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}
#endif

			} else if (HP.IsKeyWord("G") || HP.IsKeyWord("shear" "modulus")) {
				if (bGot_G) {
					silent_cerr("Shell plane stress orthotropic constitutive law: shear modulus already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_G = true;
				dG = HP.GetReal();
				if (dG <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid shear modulus " << dG << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("beta")) {
				if (bGot_Beta) {
					silent_cerr("Shell plane stress orthotropic constitutive law: fiber angle \"beta\" already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_Beta = true;
				dBeta = HP.GetReal();

				// temporary
				silent_cerr("Shell plane stress orthotropic constitutive law: fiber angle \"beta\" not supported yet" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);

			} else if (HP.IsKeyWord("thickness")) {
				if (bGot_h) {
					silent_cerr("Shell plane stress orthotropic constitutive law: thickness already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_h = true;
				dh = HP.GetReal();
				if (dh <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid thickness " << dh << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("as") /* better name? */ ) {
				if (bGot_as) {
					silent_cerr("Shell plane stress orthotropic constitutive law: as (?) already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_as = true;
				das = HP.GetReal();
				if (das <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid as " << das << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else if (HP.IsKeyWord("at") /* better name? */ ) {
				if (bGot_at) {
					silent_cerr("Shell plane stress orthotropic constitutive law: at (?) already provided at line " << HP.GetLineData() << std::endl);
					return -1;
				}
				bGot_at = true;
				dat = HP.GetReal();
				if (dat <= 0.) {
					silent_cerr("Shell plane stress orthotropic constitutive law: invalid at " << dat << " at line " << HP.GetLineData() << std::endl);
					return -1;
				}

			} else {
				break;
			}
		}

		if (!bGot_E_l) {
				silent_cerr("Shell plane stress orthotropic constitutive law: longitudinal Young's modulus missing at line " << HP.GetLineData() << std::endl);
				return -1;
		}

		if (!bGot_E_t) {
				silent_cerr("Shell plane stress orthotropic constitutive law: transverse Young's modulus missing at line " << HP.GetLineData() << std::endl);
				return -1;
		}

		if (!bGot_nu_lt && !bGot_nu_tl) {
				silent_cerr("Shell plane stress orthotropic constitutive law: Poisson's modulus missing at line " << HP.GetLineData() << std::endl);
				return -1;

		} else if (bGot_nu_lt && bGot_nu_tl) {
			if (std::abs(dnu_tl*dE_l - dnu_lt*dE_t) >
				std::numeric_limits<doublereal>::epsilon())
			{
				silent_cerr("Shell plane stress orthotropic constitutive law: inconsistent material data" << std::endl);
				return -1;
			}

		} else if (bGot_nu_tl) {
			dnu_lt = dnu_tl*dE_l/dE_t;
		}

		if (!bGot_G) {
				silent_cerr("Shell plane stress orthotropic constitutive law: shear modulus missing at line " << HP.GetLineData() << std::endl);
				return -1;
		}

		doublereal C = dh/(1. - dnu_lt*dnu_lt/dE_l*dE_t);
		doublereal G = dh*dG;
		doublereal D = 1./12.*dh*dh*dh/(1. - dnu_lt*dnu_lt/dE_l*dE_t);
		doublereal F = 1./12.*dh*dh*dh*dG;

		pD(1, 1) = C*dE_l;
		pD(1, 5) = C*dnu_lt*dE_t;
		pD(2, 2) = 2.*G;
		// pD(2, 4) = 2.*G;
		pD(3, 3) = G*das;
		// pD(4, 2) = 2.*G;
		pD(4, 4) = 2.*G;
		pD(5, 1) = C*dnu_lt*dE_t;
		pD(5, 5) = C*dE_t;
		pD(6, 6) = G*das;

#if 0
		pD(7, 7) = D*dE_l;
		pD(7, 11) = D*dnu_lt*dE_t;
		pD(8, 8) = 2.*F;
		// pD(8, 10) = 2*F;
		pD(9, 9) = F*dat;
		// pD(10, 8) = 2*F;
		pD(10, 10) = 2*F;
		pD(11, 7) = D*dnu_lt*dE_t;
		pD(11, 11) = D*dE_t;
		pD(12, 12) = F*dat;
#endif

		pD(7, 7) = 2.*F;
		pD(8, 8) = D*dE_l;
		pD(8, 10) = -D*dnu_lt*dE_t;
		// pD(8, 10) = 2*F;
		pD(9, 9) = F*dat;
		// pD(10, 8) = 2*F;
		pD(10, 8) = -D*dnu_lt*dE_t;
		pD(10, 10) = D*dE_t;
		pD(11, 11) = 2*F;
		pD(12, 12) = F*dat;

	} else {
		if (HP.IsKeyWord("matr")) {
			// tolerate it
		}

		for (unsigned ir = 1; ir <= 12; ir++) {
			for (unsigned ic = 1; ic <= 12; ic++) {
				pD(ir, ic) = HP.GetReal();
			}
		}
	}

	if (HP.IsKeyWord("prestress")) {
		ASSERT(PreStress.iGetSize() == 12);
		for (unsigned ir = 1; ir <= 12; ir++) {
			PreStress(ir) = HP.GetReal();
		}
	}

	return 0;
}

