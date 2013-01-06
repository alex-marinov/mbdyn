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

//#include <mathp.h>
#include <tcl.h>

static Tcl_Interp *interp;
static int interp_cnt;

class TclPlugIn : public MathParser::PlugIn {
protected:
	TypedValue::Type type;
	Tcl_Obj *cmd;

public:
	TclPlugIn(MathParser& mp);
	~TclPlugIn(void);
	const char *sName(void) const;
	int Read(int argc, char *argv[]);
	TypedValue::Type GetType(void) const;
	TypedValue GetVal(void) const;
};

TclPlugIn::TclPlugIn(MathParser& mp)
: MathParser::PlugIn(mp), type(TypedValue::VAR_UNKNOWN),
cmd(0)
{
	if (!::interp) {
		::interp = Tcl_CreateInterp();
		if (!::interp) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	interp_cnt++;
}

TclPlugIn::~TclPlugIn(void)
{
	Tcl_DecrRefCount(cmd);

	if (--interp_cnt == 0) {
		if (::interp) {
			Tcl_DeleteInterp(interp);
		}
	}
}

const char *
TclPlugIn::sName(void) const
{
	return 0;
}

int
TclPlugIn::Read(int argc, char *argv[])
{
	char *s_type = argv[0];
	if (strcasecmp(s_type, "real") == 0) {
		type = TypedValue::VAR_REAL;

	} else if (strcasecmp(s_type, "integer") == 0) {
		type = TypedValue::VAR_INT;

	} else if (strcasecmp(s_type, "bool") == 0) {
		type = TypedValue::VAR_BOOL;

	} else {
		silent_cerr("unknown or unhandled type \"" << s_type << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	char *s_tcl = argv[1];
	if (strncasecmp(s_tcl, "file://", STRLENOF("file://")) == 0) {
		char *fname = &s_tcl[STRLENOF("file://")];
		FILE *fin;
		std::string s;
		char buf[1024];
		int cmdlen;

		fin = fopen(fname, "r");
		if (fin == 0) {
			silent_cerr("TclPlugIn::Read: unable to open file \"" << fname << "\"" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!fgets(buf, sizeof(buf), fin)) {
			silent_cerr("TclPlugIn::Read: unable to read from file \"" << fname << "\"" << std::endl);
			fclose(fin);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		s += buf;

		while (fgets(buf, sizeof(buf), fin)) {
			s += buf;
		}
		fclose(fin);
		
		cmd = Tcl_NewStringObj(s.c_str(), s.length());

	} else {

		/*
		 * check / escape string ?
		 */
		cmd = Tcl_NewStringObj(s_tcl, strlen(s_tcl));
	}

	if (cmd == 0) {
		silent_cerr("TclPlugIn::Read: Tcl_NewStringObj failed" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Tcl_IncrRefCount(cmd);

	return 0;
}

TypedValue::Type
TclPlugIn::GetType(void) const
{
	return type;
}

TypedValue 
TclPlugIn::GetVal(void) const
{
	Tcl_Obj *res;
	
	if (Tcl_EvalObjEx(interp, cmd, 0) != TCL_OK) {
		silent_cerr("TclPlugIn::GetVal: Tcl_EvalObjEx failed" 
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	res = Tcl_GetObjResult(interp);
	if (res == 0) {
		silent_cerr("TclPlugIn::GetVal: Tcl_GetObjResult failed" 
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	switch (type) {
	case TypedValue::VAR_INT: {
		int i;
		if (Tcl_GetIntFromObj(0, res, &i) != TCL_OK) {
			silent_cerr("TclPlugIn::GetVal: Tcl_GetIntFromObj failed"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		return TypedValue(i);
		}

	case TypedValue::VAR_REAL: {
		double d;
		if (Tcl_GetDoubleFromObj(0, res, &d) != TCL_OK) {
			silent_cerr("TclPlugIn::GetVal: Tcl_GetDoubleFromObj failed"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		return TypedValue(d);
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Tcl_ResetResult(interp);
}

static MathParser::PlugIn *
tcl_plugin(MathParser& mp, void *arg)
{
	MathParser::PlugIn *p = 0;
	
	SAFENEWWITHCONSTRUCTOR(p, TclPlugIn, TclPlugIn(mp));

	return p;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
#endif
	MBDynParser	*pHP = (MBDynParser *)php;

	pHP->GetMathParser().RegisterPlugIn("tcl", tcl_plugin, 0);

	return 0;
}

