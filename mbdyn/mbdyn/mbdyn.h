/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifndef MBDYN_H
#define MBDYN_H

/* Global typedefs (unused yet) */
typedef double mbReal;
typedef int mbInt;
typedef long int mbLong;
typedef long int mbFlag;

/* Math typedefs (deprecated; will be obsoleted) */
#ifndef HAVE_F2C_H	/* defined in "f2c.h" */
typedef long int flag;	/* boolean return value; will be obsoleted by "bool" */
#endif /* !HAVE_F2C_H */
#if 0
typedef unsigned int  Index;
#endif
typedef int           Int;
typedef long int      Lint;
typedef double        Real;

/* signal types */
#ifdef HAVE_SIGNAL
#ifndef HAVE___SIGHANDLER_T
#ifndef HAVE_SIGHANDLER_T
typedef void (*__sighandler_t)(int);
#else /* HAVE_SIGHANDLER_T */
typedef sighandler_t __sighandler_t;
#endif /* HAVE_SIGHANDLER_T */
#endif /* !HAVE___SIGHANDLER_T */
#endif /* HAVE_SIGNAL */

/* sig_atomic_t */
#ifndef HAVE_SIG_ATOMIC_T
typedef int sig_atomic_t;
#endif /* HAVE_SIG_ATOMIC_T */

/* replacement for bool */
#ifdef HAVE_BOOL
#ifdef NEED_BOOL_H
#include <bool.h>
#endif /* NEED_BOOL_H */
#else /* !HAVE_BOOL */
typedef int bool;
enum {
	false = 0,
	true = 1
};
#endif /* !HAVE_BOOL */

/* decides whether to build include capability in parser */
#if defined(HAVE_GETCWD) && defined(HAVE_CHDIR)
#define USE_INCLUDE_PARSER
#endif /* defined(HAVE_GETCWD) && defined(HAVE_CHDIR) */

/* Global variables */
extern int fSilent;
extern int fPedantic;

/* Global macros */
#ifdef __cplusplus

#define silent_output \
	(::fSilent > 0)

#define silent_cout(arg) \
    	do { \
        	if (::fSilent < 1) { \
            		std::cout << arg; \
        	} \
    	} while (0)

#define silent_cerr(arg) \
	do { \
		if (::fSilent < 2) { \
			std::cerr << arg; \
		} \
	} while (0)

#define pedantic_output \
	(::fPedantic > 0)

#define pedantic_cout(arg) \
    	do { \
        	if (::fPedantic > 2) { \
            		std::cout << arg; \
        	} \
    	} while (0)

#define pedantic_cerr(arg) \
	do { \
		if (::fPedantic > 1) { \
			std::cerr << arg; \
		} \
	} while (0)

#endif /* __cplusplus */

/* Debug levels (from 0x0001 to 0x0080 are reserved) */
static const long int MYDEBUG_INPUT               = 0x00000100;
static const long int MYDEBUG_ASSEMBLY            = 0x00000200;
static const long int MYDEBUG_DERIVATIVES         = 0x00000400;
static const long int MYDEBUG_FSTEPS              = 0x00000800;
static const long int MYDEBUG_MEM                 = 0x00001000;
static const long int MYDEBUG_MPI                 = 0x00002000;
static const long int MYDEBUG_PRED                = 0x00004000;
static const long int MYDEBUG_RESIDUAL            = 0x00008000;
static const long int MYDEBUG_SOL                 = 0x00010000;
static const long int MYDEBUG_INIT                = 0x00020000;
static const long int MYDEBUG_OUTPUT              = 0x00040000;
static const long int MYDEBUG_JAC                 = 0x00080000;

#ifdef USE_RTAI
/* visible to all */
extern void *mbdyn_rtai_task;
#endif /* USE_RTAI */

#endif /* MBDYN_H */

