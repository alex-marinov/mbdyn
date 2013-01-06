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

/*
 * Registers a namespace that adds support for unit conversion
 * in math parser, based on UNIDATA's UDUNITS
 * <http://www.unidata.ucar.edu/software/udunits/>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "mathp.h"
#include "dataman.h"

extern "C" {
#include <udunits.h>
}

class UDUnitsNameSpace: public MathParser::NameSpace {
protected:
	MathParser::MathFunc_t f;

public:
	UDUnitsNameSpace(const char *path);
	~UDUnitsNameSpace(void);

	bool IsFunc(const std::string& fname) const;
	MathParser::MathFunc_t* GetFunc(const std::string& fname) const;
	TypedValue EvalFunc(MathParser::MathFunc_t *f,
		const MathParser::MathArgs& args) const;
};

static int
unit_convert(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_STRING);
	ASSERT(args[2]->Type() == MathParser::AT_STRING);
	ASSERT(args[3]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgString_t *arg1 = dynamic_cast<MathParser::MathArgString_t *>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgString_t *arg2 = dynamic_cast<MathParser::MathArgString_t *>(args[2]);
	ASSERT(arg2 != 0);

	MathParser::MathArgReal_t *arg3 = dynamic_cast<MathParser::MathArgReal_t *>(args[3]);
	ASSERT(arg3 != 0);

	int rc;

	utUnit	u_from;
	rc = utScan((*arg1)().c_str(), &u_from);
	switch (rc) {
	case 0:
		break;

	default:
		/* TODO: handle specific errors */
		silent_cerr("module-udunits: utScan could not interpret "
			"unit \"" << (*arg1)() << "\" in first arg"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	utUnit	u_to;
	rc = utScan((*arg2)().c_str(), &u_to);
	switch (rc) {
	case 0:
		break;

	default:
		silent_cerr("module-udunits: utScan could not interpret "
			"unit \"" << (*arg2)() << "\" in second arg"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	double slope, intercept;
	rc = utConvert(&u_from, &u_to, &slope, &intercept);
	switch (rc) {
	case 0:
		break;

	default:
		silent_cerr("module-udunits: utConvert failed (rc=" << rc << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (arg3->IsFlag(MathParser::AF_OPTIONAL_NON_PRESENT)) {
		if (intercept != 0) {
			silent_cerr("module-udunits: conversion between \"" << (*arg1)() << "\" "
				"and \"" << (*arg2)() << "\" "
				"has non-zero intercept" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		*out = slope;

	} else {
		*out = (*arg3)()*slope + intercept;
	}

	return 0;
}

UDUnitsNameSpace::UDUnitsNameSpace(const char *path)
: MathParser::NameSpace("units")
{
	int rc = utInit(path);

	if (path == 0) {
		path = "<unspecified>";
	}

	switch (rc) {
	case 0:
		// success
		break;

	case UT_ENOFILE:
		silent_cerr("utUnit could not find file \"" << path << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case UT_ESYNTAX:
		silent_cerr("utUnit found a syntax error "
			"in file \"" << path << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case UT_EUNKNOWN:
		silent_cerr("utUnit found an unknown specification "
			"in file \"" << path << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case UT_EIO:
		silent_cerr("utUnit hit an I/O error while reading "
			"file \"" << path << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case UT_EALLOC:
		silent_cerr("utUnit ran out of memory while reading "
			"file \"" << path << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	f.fname = "convert";
	f.args.resize(1 + 3);
	f.args[0] = new MathParser::MathArgReal_t;
	f.args[1] = new MathParser::MathArgString_t;
	f.args[2] = new MathParser::MathArgString_t;
	f.args[3] = new MathParser::MathArgReal_t(0., MathParser::AF_OPTIONAL);
	f.f = unit_convert;
	f.t = 0;
}

UDUnitsNameSpace::~UDUnitsNameSpace(void)
{
	utTerm();

	for (MathParser::MathArgs::iterator i = f.args.begin();
		i != f.args.end(); ++i)
	{
		delete *i;
	}
}

bool
UDUnitsNameSpace::IsFunc(const std::string& fname) const
{
	return GetFunc(fname) != 0;
}

MathParser::MathFunc_t*
UDUnitsNameSpace::GetFunc(const std::string& fname) const
{
	if (fname.compare(f.fname) == 0) {
		return const_cast<MathParser::MathFunc_t*>(&f);
	}

	return 0;
}

TypedValue 
UDUnitsNameSpace::EvalFunc(MathParser::MathFunc_t *f,
	const MathParser::MathArgs& args) const
{
	f->f(args);

	return TypedValue((*dynamic_cast<MathParser::MathArgReal_t*>(args[0]))());
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	MBDynParser	*pHP = (MBDynParser *)php;

	const char *path = 0;
	if (pHP->IsArg()) {
		path = pHP->GetFileName();
	}

	/* registers unit conversion namespace */
	MathParser::NameSpace *pNS = new UDUnitsNameSpace(path);
	int rc = pHP->GetMathParser().RegisterNameSpace(pNS);
	if (rc != 0) {
		delete pNS;
	}
	return rc;
}

