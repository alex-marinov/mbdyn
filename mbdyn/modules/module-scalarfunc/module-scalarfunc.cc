/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
#include "ScalarFunctions.h"

class DummySF : public DifferentiableScalarFunction {
private:
	doublereal dConst;

public:
	DummySF(const doublereal d) : dConst(d) { NO_OP; };
	virtual ~DummySF(void) { NO_OP; };
	virtual doublereal operator()(const doublereal x) const { return dConst; };
	virtual doublereal ComputeDiff(const doublereal t, const integer order = 1) const { return 0; };
};

/* prototype of the functional object: reads a drive caller */
struct DummySFR : public ScalarFunctionRead {
	const BasicScalarFunction *
	Read(DataManager* pDM, MBDynParser& HP) const {
		doublereal d = HP.GetReal();
		return new DummySF(d);
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	ScalarFunctionRead	*rf = new DummySFR;

	if (!SetSF("dummy", rf)) {
		delete rf;

		silent_cerr("DummyScalarFunction: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

