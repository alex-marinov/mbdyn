/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/******************************************************************************

Macro di assert personalizzate

Uso: ASSERT(<expr>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore std::cerr il file e la riga
		solo se DEBUG e' definita.

Uso: ASSERTMSG(<expr>, <msg>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore std::cerr il file e la riga,
		seguiti dal messaggio <msg>, solo se DEBUG e' definita.

Entrambe chiamano la funzione _Assert(file, line, msg = NULL);
se msg e' definito, viene aggiunto in coda al messaggio di default

******************************************************************************/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <cassert>
#include <cstring>
#include <cstdlib>

#include "myassert.h"

/* flag di silent run (no output su stdout) */
int fSilent = 0;
int fPedantic = 0;

#ifdef DEBUG

long int debug_level = MYDEBUG_ANY;
long int DEFAULT_DEBUG_LEVEL = MYDEBUG_ANY;

void _Assert(const char* file, const int line, const char* msg)
{
   std::cout.flush();
   
   std::cerr << std::endl << "ASSERT fault in file " << file 
     << " at line " << line;
   if (msg) { 
      std::cerr << ':' << std::endl << msg; 
   }
   std::cerr << std::endl;
   
#ifdef DEBUG_STOP
   throw MyAssert::ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif   
   
   return;
}

std::ostream& _Out(std::ostream& out, const char* file, const int line)
{
   std::cout.flush();
   
   // out << "File <" << file << ">, line [" << line << "]: ";
   out << "[" << file << "," << line << "]: ";
   return out;
}

int get_debug_options(const char *const s, const debug_array da[])
{
   if (s == NULL || s[0] == '\0') {
      ::debug_level = DEFAULT_DEBUG_LEVEL;
      return 0;
   }
   
   const char* p = s;
   while (true) {
      const char* sep = std::strchr(p, ':');
      unsigned int l;
      if (sep != NULL) {
	 l = int(sep-p);
      } else {
	 l = strlen(p);
      }
      debug_array* w = (debug_array*)da;
      while (w->s != NULL) {
	 if (l == strlen(w->s) && strncmp(w->s, p, l) == 0) {
	    break;
	 }
	 w++;
      }
      if (w->s == NULL) {
	 if (l == 4 && strncmp("none", p, 4) == 0) {
	    ::debug_level = MYDEBUG_NONE;
	 } else if (l == 3 && strncmp("any", p, 3) == 0) {
	    ::debug_level = MYDEBUG_ANY;
	 } else {
	    silent_cerr("Unknown debug level \"");
	    for (unsigned int i = 0; i < l; i++) {
	       silent_cerr(p[i]);
	    }
	    silent_cerr("\"" << std::endl);
	 }
      } else {
	 ::debug_level |= w->l;
	 silent_cerr("debug level: " << w->s << std::endl);
      }
      if (sep == NULL) {
	 break;
      }
      p = sep+1;
   }
   
   return 0;
}

#endif /* DEBUG */

#if defined(__GNUC__) && (defined(_M_IX86) || defined(__x86_64)) && !defined(NDEBUG) && (defined(__CYGWIN__) || defined(_WIN32)) 
extern "C" void __assert_func (const char* file, int line, const char* func, const char* expr)
{
    std::cerr << "assertion " << expr << " failed: file " << file  << ":" << line  << ":" << func << std::endl;

    // debug break interrupt on x86 and x86_64 makes debugging easier on Windows
    __asm__ volatile ("int $3");

    abort();
}
#endif
