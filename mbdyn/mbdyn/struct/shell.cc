/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2010
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

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
ReadShellConstLaw(MBDynParser& HP, Shell::fmh& pD)
{
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
				if (dnu <= 0. || dnu >= .5) {
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
			silent_cerr("Shell isotropic constitutive law: incomplete material data (need two of Young's modulus, Poisson's modulus, shear modulus)" << std::endl);
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

		doublereal C = dE*dh/(1. - dnu*dnu);
		doublereal D = dE*dh*dh*dh/12/(1. - dnu*dnu);

		pD.Reset();

		pD(1, 1) = C;
		pD(1, 5) = C*dnu;
		pD(2, 2) = C*(1. - dnu);
		pD(3, 3) = C*(1. - dnu)/2*das;
		pD(4, 4) = C*(1. - dnu);
		pD(5, 1) = pD(1, 5);
		pD(5, 5) = C;
		pD(6, 6) = C*(1. - dnu)/2*das;

		pD(7, 7) = D;
		pD(7, 11) = D*dnu;
		pD(8, 8) = D*(1. - dnu);
		pD(9, 9) = D*(1. - dnu)/2*dat;
		pD(10, 10) = D*(1. - dnu);
		pD(11, 7) = pD(7, 11);
		pD(11, 11) = D;
		pD(12, 12) = D*(1. - dnu)/2*dat;

	} else {
		for (unsigned ir = 1; ir <= 12; ir++) {
			for (unsigned ic = 1; ic <= 12; ic++) {
				pD(ir, ic) = HP.GetReal();
			}
		}
	}

	return 0;
}

