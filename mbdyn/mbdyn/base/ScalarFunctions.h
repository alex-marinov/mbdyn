/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2005
 *
 * Marco Morandini  <morandini@aero.polimi.it>
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
 * Copyright (C) 1996-2005
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
#ifndef ScalarFunctions_hh
#define ScalarFunctions_hh

#include "ac/f2c.h"

class BasicScalarFunction {
public:
	virtual ~BasicScalarFunction();
	virtual doublereal operator()(const doublereal x) const = 0;
};

class DifferentiableScalarFunction : public BasicScalarFunction {
public:
	virtual ~DifferentiableScalarFunction();
	virtual doublereal operator()(const doublereal x) const = 0;
	virtual doublereal ComputeDiff(const doublereal t, const integer order = 1) const = 0;
};

#include "mbpar.h"
#include "dataman.h"
//implemented in ScalarFunctionsImpl.cc
const BasicScalarFunction *const ParseScalarFunction(MBDynParser& HP,
	DataManager* const pDM);


#endif //ScalarFunctions_hh
