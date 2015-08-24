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

/*
 * With the contribution of Ankit Aggarwal <ankit.ankit.aggarwal@gmail.com>
 * during Google Summer of Code 2015
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <cfloat>
#include <cstdlib>
#include <climits>
#include <limits>
#include <sstream>

#include "mathp.h"
#include "parser.h"

#ifdef USE_EE

#include "evaluator_impl.h"

unsigned ExpressionElement::m_uEEFlags = ExpressionElement::EE_OPTIMIZE;

std::string
EEStrOut(const ExpressionElement *e)
{
	std::ostringstream out;
	e->Output(out);
	return out.str();
}

bool
EE_Eval(const ExpressionElement *ee, TypedValue& dst)
{
	if (ee == 0) {
		return false;
	}

	dst = ee->Eval();

	return true;
}

bool
EE_Eval(const ExpressionElement *ee, bool& dst)
{
	if (ee == 0) {
		return false;
	}

	dst = ee->Eval().GetBool();

	return true;
}

bool
EE_Eval(const ExpressionElement *ee, Int& dst)
{
	if (ee == 0) {
		return false;
	}

	dst = ee->Eval().GetInt();

	return true;
}

bool
EE_Eval(const ExpressionElement *ee, Real& dst)
{
	if (ee == 0) {
		return false;
	}

	dst = ee->Eval().GetReal();

	return true;
}

bool
EE_Eval(const ExpressionElement *ee, std::string& dst)
{
	if (ee == 0) {
		return false;
	}

	dst = ee->Eval().GetString();

	return true;
}

#endif // USE_EE
