/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_TCL

#include <tclpgin.h>

TclPlugIn::TclPlugIn(MathParser& mp, TypedValue::Type t)
: MathParser::PlugIn::PlugIn(mp), type(t), interp(NULL), cmd(NULL)
{
	interp = Tcl_CreateInterp();
	if (interp == NULL) {
		THROW(ErrGeneric());
	}
}

TclPlugIn::~TclPlugIn(void)
{
	/*
	 * FIXME: don't know how to destroy object!
	 */
	if (interp) {
		Tcl_DeleteInterp(interp);
	}
}

const char *
TclPlugIn::sName(void) const
{
	return NULL;
}

int
TclPlugIn::Read(int argc, char *argv[])
{
	if (strncasecmp(argv[0], "file://", 7) == 0) {
		FILE *fin;
		char buf[1024], *s;
		int cmdlen;

		fin = fopen(argv[0]+7, "r");
		if (fin == NULL) {
			cerr << "TclPlugIn::Read: error" << endl;
			THROW(ErrGeneric());
		}

		if (!fgets(buf, sizeof(buf), fin)) {
			cerr << "TclPlugIn::Read: error" << endl;
			THROW(ErrGeneric());
		}

		s = strdup(buf);
		cmdlen = strlen(buf);

		while (fgets(buf, sizeof(buf), fin)) {
			cmdlen += strlen(buf);
			s = (char *)realloc(s, cmdlen+1);
			if (s == NULL) {
				cerr << "TclPlugIn::Read: error" << endl;
				THROW(ErrGeneric());
			}
			strcat(s, buf);
		}
		
		cmd = Tcl_NewStringObj(s, cmdlen);

		free(s);
	} else {
		cmd = Tcl_NewStringObj(argv[0], strlen(argv[0]));
	}

	if (cmd == NULL) {
		cerr << "TclPlugIn::Read: error" << endl;
		THROW(ErrGeneric());
	}

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
	TypedValue rv(type);
	
	if (Tcl_GlobalEvalObj(interp, cmd) != TCL_OK) {
		cerr << "TclPlugIn::GetVal: Tcl_GlobalEvalObj: error" << endl;
		THROW(ErrGeneric());
	}

	if ((res = Tcl_GetObjResult(interp)) == NULL) {
		cerr << "TclPlugIn::GetVal: Tcl_GlobalEvalObj: error" << endl;
		THROW(ErrGeneric());
	}
	
	switch (type) {
	case TypedValue::VAR_INT: {
		Int d;
		if (Tcl_GetIntFromObj(NULL, res, &d) != TCL_OK) {
			cerr << "TclPlugIn::GetVal: Tcl_GetIntFromObj: error"
				<< endl;
			THROW(ErrGeneric());
		}
		rv = d;
		break;
	}
				  
	case TypedValue::VAR_REAL:
	default: {
		Real d;
		if (Tcl_GetDoubleFromObj(NULL, res, &d) != TCL_OK) {
			cerr << "TclPlugIn::GetVal: Tcl_GetDoubleFromObj: error"
				<< endl;
			THROW(ErrGeneric());
		}
		rv = d;
		break;
	}
	}

	Tcl_ResetResult(interp);

	return rv;
}

MathParser::PlugIn *
tcl_plugin(MathParser& mp, void *arg )
{
	MathParser::PlugIn *p = NULL;
	const char *type = (const char *)arg;
	
	if (type == NULL || strcasecmp(type, "real") == 0) {
		SAFENEWWITHCONSTRUCTOR(p, TclPlugIn::TclPlugIn,
				TclPlugIn::TclPlugIn(mp, TypedValue::VAR_REAL),
				MPmm);
	} else if (strcasecmp(type, "integer") == 0) {
		SAFENEWWITHCONSTRUCTOR(p, TclPlugIn::TclPlugIn,
				TclPlugIn::TclPlugIn(mp, TypedValue::VAR_INT), 
				MPmm);
	}
	return p;
}

#endif /* USE_TCL */

