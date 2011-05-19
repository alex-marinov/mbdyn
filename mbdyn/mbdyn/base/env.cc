/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include <sstream>

#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "ac/f2c.h"
#include "mathp.h"

/* environment */
extern char **environ;

static const char MBDYNPREFIX[] = "MBDYN";
static const int MBDYNPREFIXLEN = STRLENOF(MBDYNPREFIX);

extern void GetEnviron(MathParser& MP);

void
GetEnviron(MathParser& MP)
{
   	/* cerca la variabile MBDYNVARS */
   	char* p = NULL;
   	const char* s = "VARS";
	int l = MBDYNPREFIXLEN + strlen(s) + 1;
   	SAFENEWARR(p, char, l);
   	sprintf(p, "%s%s", MBDYNPREFIX, s);
   	char* e = getenv(p);
   	SAFEDELETEARR(p);

   	if (e != NULL) {
      		DEBUGCOUT("GetEnv: reading variable <" << e << ">" << std::endl);
		std::istringstream in(e);
      		InputStream In(in);
      		MP.GetLastStmt(In);
      		DEBUGCOUT("GetEnv: variable <" << e << "> read" << std::endl);
   	}

   	/* cerca le variabili definite singolarmente */
   	Table& T = MP.GetSymbolTable();
   	char** env = environ;
   	while (*env) {
      		if (strncmp(*env, MBDYNPREFIX, MBDYNPREFIXLEN) == 0) {
	 		DEBUGCOUT("reading var <" << *env << ">" << std::endl);
	 		long int i = 0;
	 		doublereal d = 0.;
	 		char* p = NULL;
	 		char* v = NULL;
	 		char* n = NULL;

	 		SAFESTRDUP(p, *env);
	 		v = std::strchr(p, '=');
	 		if (v == NULL) {
	    			silent_cerr("parse error in envvar <"
					<< p << ">" << std::endl);
	    			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	 		}

	 		*v = '\0';
	 		v++;

	 		if (strncmp(p, "MBDYNVARS", STRLENOF("MBDYNVARS")) == 0) {
				NO_OP;

			} else if (strncmp(p, "MBDYN_real_", STRLENOF("MBDYN_real_")) == 0) {
	    			n = p + STRLENOF("MBDYN_real_");
				char *endptr = NULL;
				errno = 0;
	    			d = strtod(v, &endptr);
				int save_errno = errno;
				if (endptr != NULL && endptr[0] != '\0') {
					silent_cerr("SetEnv: unable to parse "
							"real <" << v << "> "
							"for var <" << p << ">"
							<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);

				} else if (save_errno == ERANGE) {
					silent_cerr("SetEnv: real <" << v << "> "
						"for var <" << p << "> overflows"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
	    			DEBUGCOUT("setting real var <"
					<< n << "=" << d << ">" << std::endl);

	    			if ((T.Get(n)) == NULL) {
	       				if (T.Put(n, Real(d))  == NULL) {
		  				silent_cerr("SetEnv:"
							" error in insertion"
							" of real symbol <"
		    					<< n << ">" << std::endl);
		  				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       				}
	    			}

	 		} else if (strncmp(p, "MBDYN_integer_", STRLENOF("MBDYN_integer_")) == 0) {
	    			n = p + STRLENOF("MBDYN_integer_");
#ifdef HAVE_STRTOL
				char *endptr = NULL;
				errno = 0;
				i = strtol(v, &endptr, 10);
				int save_errno = errno;
				if (endptr != NULL && endptr[0] != '\0') {
					silent_cerr("SetEnv: unable to parse "
						"integer <" << v << "> "
						"for var <" << p << ">"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);

				} else if (save_errno == ERANGE) {
					silent_cerr("SetEnv: integer <" << v << "> "
						"for var <" << p << "> overflows"
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
#else /* !HAVE_STRTOL */
	    			i = atol(v);
#endif /* !HAVE_STRTOL */
	    			DEBUGCOUT("setting integer var <"
					<< n << "=" << i << ">" << std::endl);

	    			if ((T.Get(n)) == NULL) {
	       				if (T.Put(n, Int(i))  == NULL) {
		  				silent_cerr("SetEnv:"
							" error in insertion"
							" of integer symbol <"
		    					<< n << ">" << std::endl);
		  				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       				}
	    			}

	 		} else if (strncmp(p, "MBDYN_string_", STRLENOF("MBDYN_string_")) == 0) {
	    			n = p + STRLENOF("MBDYN_string_");
	    			if ((T.Get(n)) == NULL) {
	       				if (T.Put(n, TypedValue(std::string(v)))  == NULL) {
		  				silent_cerr("SetEnv:"
							" error in insertion"
							" of string symbol <"
		    					<< n << ">" << std::endl);
		  				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	       				}
	    			}

	 		} else {
	    			silent_cerr("unknown var type <"
					<< p << ">; skipping ..." << std::endl);
	 		}
      		}
      		env++;
   	}

   	/* altro ... */
}

