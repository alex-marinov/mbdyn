/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifndef EXCEPT_H
#define EXCEPT_H

#include <stdio.h>

#include <ac/iostream>

#if 0
#ifdef USE_EXCEPTIONS
#define THROW(cl) \
    do {          \
        throw cl; \
    } while (0)
#else /* !USE_EXCEPTIONS */
#define THROW(cl)                                                   \
    do {                                                            \
        fprintf(stderr, "(%s,%d) aborting after call to exit()\n",  \
                __FILE__, __LINE__);                                \
        exit(EXIT_FAILURE);                                         \
    } while (0)
#endif /* !USE_EXCEPTIONS */
#endif

class NoErr {};
class ErrGeneric {};
class ErrInterrupted {};

class ErrOutOfRange {};
class ErrDivideByZero {};
class ErrMemory {
  public: 
    ErrMemory(void) {};
    ErrMemory(const char* const s) { silent_cerr(s << std::endl); };
    ErrMemory(std::ostream& out, const char* const s) { out << s << std::endl; };   
};

class EndOfFile {};
class ErrFile {};
class ErrFileSystem {};

class ErrNotAvailableYet {};
class ErrNotImplementedYet {};

#endif /* EXCEPT_H */

