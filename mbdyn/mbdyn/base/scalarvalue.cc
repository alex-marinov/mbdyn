/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 * Michele Attolico <attolico@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "scalarvalue.h"

/* ScalarValue - begin */

ScalarValue::~ScalarValue(void)
{
	NO_OP;
}

ScalarDofValue::ScalarDofValue(const ScalarDof& sd)
: ScalarDof(sd)
{
	NO_OP;
}

doublereal
ScalarDofValue::dGetValue(void) const
{
	return ScalarDof::dGetValue();
}

ScalarDriveValue::ScalarDriveValue(const DriveCaller *pdc)
: pDC(pdc)
{
	NO_OP;
}

ScalarDriveValue::~ScalarDriveValue(void)
{
	if (pDC) {
		delete pDC;
	}
}

doublereal
ScalarDriveValue::dGetValue(void) const
{
	return pDC->dGet();
}

ScalarValue *
ReadScalarValue(DataManager *pDM, MBDynParser& HP)
{
	ScalarValue *svp = 0;

	if (HP.IsKeyWord("drive")) {
		svp = new ScalarDriveValue(HP.GetDriveCaller(false));

	} else {
		if (!HP.IsKeyWord("node" "dof")) {
			silent_cerr("Warning, missing keyword \"node dof\" "
				"at line " << HP.GetLineData() << std::endl);
		}

		svp = new ScalarDofValue(ReadScalarDof(pDM, HP, false, true));
	}

	return svp;
}

/* NOTE: the array must have the desired size */
void
ReadScalarValues(DataManager *pDM, MBDynParser& HP,
	std::vector<ScalarValue *>& Values)
{
	unsigned nch = Values.size();
	if (nch == 0) {
		silent_cerr("Request to read an empty ScalarValue vector "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (unsigned int i = 0; i < nch; i++) {
		Values[i] = ReadScalarValue(pDM, HP);
	}
}

/* ScalarValue - end */

