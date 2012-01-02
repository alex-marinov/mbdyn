/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#ifndef MBSLEEP_H
#define MBSLEEP_H

#include "ac/f2c.h"
#include "math.h"

#ifdef __cplusplus
#include <ostream>
extern "C" {
#endif /* __cplusplus */

#if defined(HAVE_NANOSLEEP)
#include <time.h>
typedef struct timespec mbsleep_t;
#elif defined(HAVE_USLEEP)
#include <unistd.h>
typedef useconds_t mbsleep_t;
#else /* ! USLEEP */
typedef unsigned long mbsleep_t;
#endif /* ! SLEEP */

extern mbsleep_t mbsleep_init(long t);
extern int mbsleep_real2sleep(doublereal d, mbsleep_t *t);
extern int mbsleep_sleep2real(const mbsleep_t *t, doublereal *d);
extern int mbsleep(const mbsleep_t *t);

#ifdef __cplusplus
}

#ifdef HAVE_NANOSLEEP
extern std::ostream& operator << (std::ostream& out, const mbsleep_t& t);

extern bool operator < (const mbsleep_t& t1, const mbsleep_t& t2);
extern bool operator > (const mbsleep_t& t1, const mbsleep_t& t2);
extern bool operator <= (const mbsleep_t& t1, const mbsleep_t& t2);
extern bool operator >= (const mbsleep_t& t1, const mbsleep_t& t2);
extern bool operator == (const mbsleep_t& t1, const mbsleep_t& t2);
extern bool operator != (const mbsleep_t& t1, const mbsleep_t& t2);

extern bool operator < (const mbsleep_t& t1, const long& t2);
extern bool operator > (const mbsleep_t& t1, const long& t2);
extern bool operator <= (const mbsleep_t& t1, const long& t2);
extern bool operator >= (const mbsleep_t& t1, const long& t2);
extern bool operator == (const mbsleep_t& t1, const long& t2);
extern bool operator != (const mbsleep_t& t1, const long& t2);
/* otherwise operators are not needed */
#endif /* HAVE_NANOSLEEP */

#endif /* __cplusplus */

#endif /* MBSLEEP_H */
