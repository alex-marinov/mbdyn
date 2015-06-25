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

#include "dummypgin.h"

DummyPlugIn::DummyPlugIn(MathParser& mp, void *arg)
: MathParser::PlugIn(mp)
{
	ASSERT(arg != NULL);
	throw MathParser::ErrGeneric(&mp, MBDYN_EXCEPT_ARGS, "error: '",
		static_cast<const char *>(arg), "'");
}

DummyPlugIn::~DummyPlugIn(void)
{
}

const char *
DummyPlugIn::sName(void) const
{
	return NULL;
}

int
DummyPlugIn::Read(int /* argc */ , char * /* argv */ [])
{
	return 0;
}

TypedValue::Type
DummyPlugIn::GetType(void) const
{
	return TypedValue::VAR_UNKNOWN;
}

TypedValue 
DummyPlugIn::GetVal(void) const
{
	return TypedValue(0);
}

MathParser::PlugIn *
dummy_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = NULL;
	
	SAFENEWWITHCONSTRUCTOR(p, DummyPlugIn, DummyPlugIn(mp, arg));

	return p;
}

