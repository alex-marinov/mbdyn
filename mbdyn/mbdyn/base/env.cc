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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sstream.h>

#include <mathp.h>

/* environment */
extern char **environ;

static const char MBDYNPREFIX[] = "MBDYN";
static const int MBDYNPREFIXLEN = strlen(MBDYNPREFIX);

extern void GetEnviron(MathParser& MP);

void
GetEnviron(MathParser& MP)
{
   	/* cerca la variabile MBDYNVARS */
   	char* p = NULL;
   	char* s = "VARS";
	int l = MBDYNPREFIXLEN+strlen(s)+1;
   	SAFENEWARR(p, char, l);
   	sprintf(p, "%s%s", MBDYNPREFIX, s);
   	char* e = getenv(p);
   	SAFEDELETEARR(p);     
   
   	if (e != NULL) {
      		DEBUGCOUT("GetEnv: reading variable <" << e << ">" << endl);
#if defined(HAVE_SSTREAM)
		std::istringstream in(e);
#else /* HAVE_STRSTREAM_H */
      		istrstream in(e);
#endif /* HAVE_STRSTREAM_H */
      		InputStream In(in);	 
      		MP.GetLastStmt(In);
      		DEBUGCOUT("GetEnv: variable <" << e << "> read" << endl);
   	}
   
   	/* cerca le variabili definite singolarmente */
   	Table& T = MP.GetSymbolTable();
   	char** env = environ;
   	while (*env) {
      		if (strncmp(*env, MBDYNPREFIX, MBDYNPREFIXLEN) == 0) {
	 		DEBUGCOUT("reading var <" << *env << ">" << endl);
	 		long int i = 0;
	 		double d = 0.;
	 		char* p = NULL;
	 		char* v = NULL;
	 		char* n = NULL;
	 
	 		SAFESTRDUP(p, *env);
	 		v = strchr(p, '=');
	 		if (v == NULL) {
	    			cerr << "parse error in envvar <" 
					<< p << ">" << endl;	  
	    			THROW(ErrGeneric());
	 		}

	 		*v = '\0';
	 		v++;
	 
	 		if (strncmp(p+5, "_real_", 6) == 0) {
	    			n = p+11;
	    			d = atof(v);
	    			DEBUGCOUT("setting real var <" 
					<< n << "=" << d << ">" << endl);
	    
	    			if ((T.Get(n)) == NULL) {	   
	       				if (T.Put(n, Real(d))  == NULL) {      
		  				cerr << "SetEnv:"
							" error in insertion"
							" of real symbol <"
		    					<< n << ">" << endl;
		  				THROW(ErrGeneric());
	       				}
	    			}
	    
	 		} else if (strncmp(p+5, "_integer_", 9) == 0) {
	    			n = p+14;
	    			i = atoi(v);
	    			DEBUGCOUT("setting integer var <" 
					<< n << "=" << i << ">" << endl);
	    
	    			if ((T.Get(n)) == NULL) {	   
	       				if (T.Put(n, Int(i))  == NULL) {
		  				cerr << "SetEnv:"
							" error in insertion"
							" of integer symbol <"
		    					<< n << ">" << endl;
		  				THROW(ErrGeneric());
	       				}
	    			}
	    
	 		} else {
	    			cerr << "unknown var type <" 
					<< p << ">; skipping ..." << endl;
	 		}
      		}
      		env++;
   	}
   
   	/* altro ... */
}

