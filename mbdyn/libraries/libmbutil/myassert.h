/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

Uso:

		ASSERT(<expr>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore cerr il file e la riga
		solo se DEBUG e' definita.

		ASSERTMSG(<expr>, <msg>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore cerr il file e la riga,
		seguiti dal messaggio <msg>, solo se DEBUG e' definita.

		ASSERTBREAK(<expr>)
		ASSERTMSGBREAK(<expr>, <msg>)

		sono analoghe alle precedenti, ma in piu' chiamano exit(1) per
		interrompere il programma.


		Sono definite inoltre le macro COUT e CERR.
		Entrambe antepongono al messaggio che le segue una intestazione
		con il nome del file ed il numero di riga da cui proviene la chiamata.
		Prima della scrittura viene forzato il flush dello standard output
		nel	caso che i due flussi di uscita facessero capo alla stessa device.
                Infine le macro DEBUGCOUT(msg) e DEBUGCERR(msg) che scrivono
                il messaggio desiderato sui rispettivi flussi in uscita
                solo se e' definito #define DEBUG. Attenzione, se si concatenano
                piu' item con l'operatore << non occorre usare parentesi strane.

******************************************************************************/


#ifndef MYASSERT_H
#define MYASSERT_H

#include <stdlib.h>
#include <iostream>
#include <except.h>

#define NO_OP do {} while(0)

#ifdef DEBUG
/* predefined debug level flags (reserved from 0x0000001 to 0x00000080) */
const long int MYDEBUG_NONE                = 0x00000000;
const long int MYDEBUG_ANY                 = 0xFFFFFFFF;
const long int MYDEBUG_FNAMES              = 0x00000001;
#ifdef __GNUC__
const long int MYDEBUG_PRETTYFN            = 0x00000002;
#endif /* __GNUC__ */

/* debug level global var */
extern long int debug_level /* = MYDEBUG_ANY */ ;
extern long int DEFAULT_DEBUG_LEVEL;

struct debug_array {
   const char*  s;
   long int l;
};

extern int get_debug_options(const char *const s, const debug_array da[]);
   



class MyAssert {
 public:
   class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
};

extern void _Assert(const char* file, const int line, const char* msg = NULL);
extern std::ostream& _Out(std::ostream& out, const char* file, const int line);

#define ASSERT(expr) \
    do { \
        if (!(expr)) { \
            _Assert(__FILE__, __LINE__); \
        } \
    } while (0)


#define ASSERTBREAK(expr) \
    do { \
        if (!(expr)) { \
            _Assert(__FILE__, __LINE__); \
            throw MyAssert::ErrGeneric(MBDYN_EXCEPT_ARGS); \
        } \
    } while (0)

#define ASSERTMSG(expr, msg) \
    do { \
        if (!(expr)) { \
            _Assert(__FILE__, __LINE__, (msg)); \
        } \
    } while (0)

#define ASSERTMSGBREAK(expr, msg) \
    do { \
        if (!(expr)) { \
            _Assert(__FILE__, __LINE__, (msg)); \
            throw MyAssert::ErrGeneric(MBDYN_EXCEPT_ARGS); \
        } \
    } while (0)




#define COUT \
    _Out(std::cout, __FILE__, __LINE__)

#define CERR \
    _Out(std::cerr, __FILE__, __LINE__)



#define DEBUGCOUT(msg) \
    do { \
        _Out(std::cout, __FILE__, __LINE__) << msg; std::cout.flush(); \
    } while (0)

#define DEBUGCERR(msg) \
    do { \
        _Out(std::cerr, __FILE__, __LINE__) << msg; std::cerr.flush(); \
    } while (0)

#define DEBUG_LEVEL(level) \
    ((level) & ::debug_level)

#define DEBUG_LEVEL_MATCH(level) \
      (((level) & ::debug_level) == (level))

#define DEBUGLCOUT(level, msg) \
    do { \
        if (::debug_level & (level)) { \
            DEBUGCOUT(msg); \
        } \
    } while (0)

#define DEBUGLCERR(level, msg) \
    do { \
        if (::debug_level & (level)) { \
            DEBUGCERR(msg); \
        } \
    } while (0)

#define DEBUGLMCOUT(level, msg) \
    do { \
        if ((::debug_level & (level)) == (level)) { \
            DEBUGCOUT(msg); \
        } \
    } while (0)

#define DEBUGLMCERR(level, msg) \
    do { \
        if ((::debug_level & (level)) == (level)) { \
            DEBUGCERR(msg); \
        } \
    } while (0)

#if defined(__GNUC__)
#define DEBUGCOUTFNAME(fname) \
    do { \
        if (::debug_level & MYDEBUG_FNAMES) { \
            if (::debug_level & MYDEBUG_PRETTYFN) { \
                DEBUGCOUT("Entering <" << __PRETTY_FUNCTION__ << ">" << std::endl); \
            } else { \
                DEBUGCOUT("Entering <" << __FUNCTION__ << ">" << std::endl); \
            } \
        } \
    } while (0)
#else /* !__GNUC__ */
#define DEBUGCOUTFNAME(fname) \
    DEBUGLCOUT(MYDEBUG_FNAMES, "Entering <" << fname << ">" << std::endl)
#endif /* !__GNUC__ */

#else /* !DEBUG */
#define ASSERT(expr) \
    NO_OP

#define ASSERTBREAK(expr) \
    NO_OP

#define ASSERTMSG(expr, msg) \
    NO_OP

#define ASSERTMSGBREAK(expr, msg) \
    NO_OP

#define COUT \
    std::cout

#define CERR \
    std::cerr


#define DEBUGCOUT(msg) \
    do { } while (0)        

#define DEBUGCERR(msg) \
    do { } while (0)

#define DEBUG_LEVEL(level) \
    0

#define DEBUG_LEVEL_MATCH(level) \
    0

#define DEBUGLCOUT(level, msg) \
    do { } while (0)        

#define DEBUGLCERR(level, msg) \
    do { } while (0)

#define DEBUGLMCOUT(level, msg) \
    do { } while (0)        

#define DEBUGLMCERR(level, msg) \
    do { } while (0)

#define DEBUGCOUTFNAME(fname) \
    do { } while (0)

#endif /* !DEBUG */

#endif /* MYASSERT_H */
