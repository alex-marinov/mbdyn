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

#include "mbconfig.h"

#include <stdlib.h>

#include "mbsleep.h"

mbsleep_t
mbsleep_init(long t)
{
#if defined(HAVE_NANOSLEEP)
	mbsleep_t v = { t, 0 };
	return v;
#else /* ! HAVE_NANOSLEEP */
	return t;
#endif /* ! HAVE_NANOSLEEP */
}

int
mbsleep_real2sleep(doublereal d, mbsleep_t *t)
{
	if (d < 0) {
		return -1;
	}

#if defined(HAVE_NANOSLEEP)
	t->tv_sec = (time_t)floor(d);
	t->tv_nsec = (long)floor((d - t->tv_sec)*1000000000);
#elif defined(HAVE_USLEEP)
	*t = (unsigned long)floor(d*1000000);
#elif defined(HAVE_SLEEP)
	if (d > 1.) {
		*t = (unsigned int)floor(d);
	} else {
		*t = 1;
	}
#else /* ! SLEEP */
	*t = 0;
#endif /* ! SLEEP */

	return 0;
}

int
mbsleep_sleep2real(const mbsleep_t *t, doublereal *d)
{
#if defined(HAVE_NANOSLEEP)
	*d = ((doublereal)(t->tv_nsec))/1000000000 + t->tv_sec;
#elif defined(HAVE_USLEEP)
	*d = ((doublereal)*t)/1000000;
#elif defined(HAVE_SLEEP)
	*d = (doublereal)*t;
#else /* ! SLEEP */
	*d = 0.;
#endif /* ! SLEEP */

	return 0;
}

int
mbsleep(const mbsleep_t *t)
{
#if defined(HAVE_NANOSLEEP)
	return nanosleep(t, NULL);
#elif defined(HAVE_USLEEP)
	return usleep(*t);
#elif defined(HAVE_SLEEP)
	return sleep(*t);
#endif /* ! SLEEP */

	return -1;
}

