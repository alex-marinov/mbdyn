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

#include "dataman.h"
#include "userelem.h"

// include here known loadable element parsers that go into InitUDE()
#include "loadable.h"

#ifdef STATIC_MODULES
#include "module-wheel2/module-wheel2.h"
#include "module-asynchronous_machine/module-asynchronous_machine.h"
#include "module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.h"
#include "module-inline_friction/module-inline_friction.h"
#include "module-cyclocopter/module-cyclocopter.h"
#ifdef HAVE_CHARM
#include "module-charm/mbcharm.h"
#endif // HAVE_CHARM
#endif // STATIC_MODULES

typedef std::map<std::string, UserDefinedElemRead *, ltstrcase> UDEMapType;
static UDEMapType UDEMap;

struct UDEWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::UDEMap.find(std::string(s)) != ::UDEMap.end();
	};
};
static UDEWordSetType UDEWordSet;

UserDefinedElem *
ParseUserDefinedElem(unsigned uLabel, DofOwner* pDO,
	DataManager* const pDM, MBDynParser& HP)
{
	DEBUGCOUTFNAME("ParseUserDefinedElem(" << uLabel << ")");

	const char *s = HP.IsWord(::UDEWordSet);
	if (s == 0) {
		silent_cerr("ParseUserDefinedElem(" << uLabel << "): "
			"unknown user-defined element type "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	UDEMapType::iterator i = ::UDEMap.find(std::string(s));
	if (i == ::UDEMap.end()) {
		silent_cerr("ParseUserDefinedElem(" << uLabel << "): "
			"unknown user-defined element type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return i->second->Read(uLabel, pDO, pDM, HP);
}

// legacy
Elem * 
ReadLoadable(DataManager* pDM, MBDynParser& HP, 
     const DofOwner* pDO, unsigned int uLabel)
{
	UDEMapType::iterator i = UDEMap.find("loadable");
	ASSERT(i != UDEMap.end());

	return i->second->Read(uLabel, pDO, pDM, HP);
}

bool
SetUDE(const std::string &s, UserDefinedElemRead *rude)
{
	pedantic_cout("registering user-defined element \"" << s << "\""
		<< std::endl );
	return UDEMap.insert(UDEMapType::value_type(s, rude)).second;
}

static unsigned done = 0;

void
InitUDE(void)
{
	if (::done++ > 0) {
		return;
	}

	bool b;

	b = SetUDE("loadable", new LoadableElemRead);
	ASSERT(b != false);
#ifdef STATIC_MODULES
	b = wheel2_set();
	ASSERT(b != false);
	b = asynchronous_machine_set();
	ASSERT(b != false);
	b = hydrodynamic_plain_bearing_set();
	ASSERT(b != false);
	b = inline_friction_set();
	ASSERT(b != false);
	b = mbdyn_cyclocopter_set();
	ASSERT(b != false);
#ifdef HAVE_CHARM
	b = mbcharm_set();
	ASSERT(b != false);
#endif // HAVE_CHARM
#endif // STATIC_MODULES
}

void
DestroyUDE(void)
{
	if (::done == 0) {
		silent_cerr("DestroyUDE() called once too many" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--done > 0) {
		return;
	}

	for (UDEMapType::iterator i = UDEMap.begin(); i != UDEMap.end(); ++i) {
		delete i->second;
	}
	UDEMap.clear();
}

// base class for user-defined elements
UserDefinedElem::UserDefinedElem(unsigned uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
InitialAssemblyElem(uLabel, flag(0)),
AerodynamicElem(uLabel, pDO, flag(0)),
ElemGravityOwner(uLabel, flag(0)),
needsAirProperties(false)
{
	NO_OP;
}

UserDefinedElem::~UserDefinedElem(void)
{
	NO_OP;
}

bool
UserDefinedElem::NeedsAirProperties(void) const
{
	return needsAirProperties;
}

void
UserDefinedElem::NeedsAirProperties(bool yesno)
{
	needsAirProperties = yesno;
}

Elem::Type 
UserDefinedElem::GetElemType(void) const
{
   	return Elem::LOADABLE;
}

AerodynamicElem::Type
UserDefinedElem::GetAerodynamicElemType(void) const
{
	return AerodynamicElem::AERODYNAMICLOADABLE;
}

