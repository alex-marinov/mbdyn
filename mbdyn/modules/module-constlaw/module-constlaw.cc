/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"

template <class T, class Tder>
class DummyConstitutiveLaw
: public ConstitutiveLaw<T, Tder> {
private:
	doublereal dStiffness;

public:
	DummyConstitutiveLaw(doublereal dStiff)
	: dStiffness(dStiff) {
		ConstitutiveLaw<T, Tder>::FDE = mb_deye<Tder>(dStiffness);
	};

	virtual ~DummyConstitutiveLaw(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const {
		ConstitutiveLaw<T, Tder>* pCL = NULL;

		typedef DummyConstitutiveLaw<T, Tder> cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(dStiffness));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "dummy, " << dStiffness;
	};

	virtual void Update(const T& Eps, const T& /* EpsPrime */  = 0.) {
		ConstitutiveLaw<T, Tder>::Epsilon = Eps;
		ConstitutiveLaw<T, Tder>::F = ConstitutiveLaw<T, Tder>::Epsilon*dStiffness;
	};
};

/* specific functional object(s) */
template <class T, class Tder>
struct DummyCLR : public ConstitutiveLawRead<T, Tder> {
	virtual ConstitutiveLaw<T, Tder> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<T, Tder>* pCL = 0;

		CLType = ConstLawType::ELASTIC;

		doublereal dS = HP.GetReal();
		if (dS <= 0.) {
			silent_cerr("warning, null or negative stiffness "
				"at line " << HP.GetLineData() << std::endl);
		}

		typedef DummyConstitutiveLaw<T, Tder> L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(dS));

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

	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new DummyCLR<doublereal, doublereal>;
	if (!SetCL1D("dummy", rf1D)) {
		delete rf1D;

		silent_cerr("DummyConstitutiveLaw1D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new DummyCLR<Vec3, Mat3x3>;
	if (!SetCL3D("dummy", rf3D)) {
		delete rf3D;

		silent_cerr("DummyConstitutiveLaw3D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	ConstitutiveLawRead<Vec6, Mat6x6> *rf6D = new DummyCLR<Vec6, Mat6x6>;
	if (!SetCL6D("dummy", rf6D)) {
		delete rf6D;

		silent_cerr("DummyConstitutiveLaw6D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

