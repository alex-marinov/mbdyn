/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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
/*
 * Author: Andrea Zanoni <>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"

class MusclePennestriCL
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	// define parameters
	doublereal Li;
	doublereal L0;
	doublereal F0;
	DriveOwner Activation;

public:
	MusclePennestriCL(doublereal Li, doublereal L0, doublereal F0, const DriveCaller *pAct)
	: Li(Li), L0(L0), F0(F0), Activation(pAct) {
		// pass parameters via constructor and initialize
	};

	virtual ~MusclePennestriCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef MusclePennestriCL cl;
		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(Li, L0, F0, Activation.pGetDriveCaller()));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "muscle, initial length, " << Li << ", " << L0 << ", " << F0 << ", ",
		       Activation.pGetDriveCaller()->Restart(out);
		return out;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal a = Activation.dGet();
		if (a < 0. || a > 1.) {
			silent_cerr("MusclePennestriCL: activation overflow (a=" << a << ")" << std::endl);
			// error out?
		}

		doublereal dxdEps = Li/L0;
		doublereal dvdEpsPrime = dxdEps/2.5; // ?!?
		doublereal x = (1. + ConstitutiveLaw<doublereal, doublereal>::Epsilon)*dxdEps;
		doublereal v = ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime*dvdEpsPrime;

		doublereal f1 = std::exp(std::pow(x - 0.95, 2) - 40*std::pow(x - 0.95, 4));
		doublereal f2 = 1.6 - 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1, 4));
		doublereal f3 = 1.3*std::atan(0.1*std::pow(x - 0.22, 10));

		doublereal df1dx = f1*(2*(x - 0.95) - 4*40.*std::pow(x - 0.95, 3));
		doublereal df2dv = 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1, 4))*(2*0.1/std::pow(v - 1., 3) - 4*1.1/std::pow(v - 1., 5));
		doublereal df3dx = 1.3*std::pow(x - 0.22, 9)/(0.01*std::pow(x - 0.22, 20) + 1);

		ConstitutiveLaw<doublereal, doublereal>::F = F0*(f1*f2*a + f3);
		ConstitutiveLaw<doublereal, doublereal>::FDE = F0*(df1dx*f2*a + df3dx)*dxdEps;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = F0*f1*df2dv*a*dvdEpsPrime;
	};
};

/* specific functional object(s) */
struct MusclePennestriCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		doublereal Li = -1.;
		if (HP.IsKeyWord("initial" "length")) {
			Li = HP.GetReal();
			if (Li <= 0.) {
				silent_cerr("null or negative initial length "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal L0 = HP.GetReal();
		if (L0 <= 0.) {
			silent_cerr("null or negative reference length "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal F0 = HP.GetReal();
		if (F0 <= 0.) {
			silent_cerr("null or negative reference force "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const DriveCaller *pAct = HP.GetDriveCaller();

		if (Li == -1.) {
			Li = L0;
		}

		typedef MusclePennestriCL L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(Li, L0, F0, pAct));

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new MusclePennestriCLR;
	if (!SetCL1D("muscle" "pennestri", rf1D)) {
		delete rf1D;

		silent_cerr("MusclePennestriCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

