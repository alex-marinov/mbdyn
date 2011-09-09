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

class MuscleConstitutiveLaw
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	// define parameters
	doublereal dStiffness;
	DriveOwner Activation;

public:
	MuscleConstitutiveLaw(doublereal dStiff, const DriveCaller *pAct)
	: dStiffness(dStiff), Activation(pAct) {
		// pass parameters via constructor and initialize
		// ConstitutiveLaw<doublereal, doublereal>::FDE = dStiffness;
	};

	virtual ~MuscleConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = NULL;

		typedef MuscleConstitutiveLaw cl;
		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(dStiffness, Activation.pGetDriveCaller()));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "muscle, " << dStiffness << ", ",
		       Activation.pGetDriveCaller()->Restart(out);
		return out;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps;

		doublereal a = Activation.dGet();
		if (a < 0. || a > 1.) {
			silent_cerr("MuscleConstitutiveLaw: activation overflow (a=" << a << ")" << std::endl);
			// error out?
		}
		ConstitutiveLaw<doublereal, doublereal>::FDE = dStiffness*(1. + a);
		ConstitutiveLaw<doublereal, doublereal>::F
			= ConstitutiveLaw<doublereal, doublereal>::FDE*ConstitutiveLaw<doublereal, doublereal>::Epsilon;
	};
};

/* specific functional object(s) */
struct MuscleCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		doublereal dS = HP.GetReal();
		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness "
				"at line " << HP.GetLineData() << std::endl);
		}

		const DriveCaller *pAct = HP.GetDriveCaller();

		typedef MuscleConstitutiveLaw L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(dS, pAct));

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

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new MuscleCLR;
	if (!SetCL1D("muscle", rf1D)) {
		delete rf1D;

		silent_cerr("MuscleConstitutiveLaw1D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

