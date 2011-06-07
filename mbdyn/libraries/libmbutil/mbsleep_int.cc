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

#include "myassert.h"
#include "mbsleep.h"

#ifdef HAVE_NANOSLEEP
std::ostream&
operator << (std::ostream& out, const mbsleep_t& t)
{
	out << t.tv_sec;
	if (t.tv_nsec) {
		out << '.' << t.tv_nsec;
	}
	return out;
}

bool
operator < (const mbsleep_t& t1, const mbsleep_t& t2)
{
	if (t1.tv_sec != t2.tv_sec) {
		return t1.tv_sec < t2.tv_sec;
	}
	return t1.tv_nsec < t2.tv_nsec;
}

bool
operator > (const mbsleep_t& t1, const mbsleep_t& t2)
{
	if (t1.tv_sec != t2.tv_sec) {
		return t1.tv_sec > t2.tv_sec;
	}
	return t1.tv_nsec > t2.tv_nsec;
}

bool
operator <= (const mbsleep_t& t1, const mbsleep_t& t2)
{
	if (t1.tv_sec != t2.tv_sec) {
		return t1.tv_sec <= t2.tv_sec;
	}
	return t1.tv_nsec <= t2.tv_nsec;
}

bool
operator >= (const mbsleep_t& t1, const mbsleep_t& t2)
{
	if (t1.tv_sec != t2.tv_sec) {
		return t1.tv_sec >= t2.tv_sec;
	}
	return t1.tv_nsec >= t2.tv_nsec;
}

bool
operator == (const mbsleep_t& t1, const mbsleep_t& t2)
{
	return (t1.tv_sec == t2.tv_sec) && (t1.tv_nsec == t2.tv_nsec);
}

bool
operator != (const mbsleep_t& t1, const mbsleep_t& t2)
{
	return (t1.tv_sec != t2.tv_sec) || (t1.tv_nsec != t2.tv_nsec);
}

bool
operator < (const mbsleep_t& t1, const long& t2)
{
	if (t1.tv_sec != t2) {
		return t1.tv_sec < t2;
	}
	// NOTE: we assume that tv_nsec >= 0
	ASSERT(t1.tv_nsec >= 0);
	return false;
}

bool
operator > (const mbsleep_t& t1, const long& t2)
{
	if (t1.tv_sec != t2) {
		return t1.tv_sec > t2;
	}
	return t1.tv_sec > 0;
}

bool
operator <= (const mbsleep_t& t1, const long& t2)
{
	if (t1.tv_sec != t2) {
		return t1.tv_sec < t2;
	}
	// NOTE: we assume that tv_nsec >= 0
	ASSERT(t1.tv_nsec >= 0);
	return t1.tv_sec == 0;
}

bool
operator >= (const mbsleep_t& t1, const long& t2)
{
	if (t1.tv_sec != t2) {
		return t1.tv_sec > t2;
	}
	return t1.tv_sec > 0;
}

bool
operator == (const mbsleep_t& t1, const long& t2)
{
	return (t1.tv_sec == t2) && (t1.tv_nsec == 0);
}

bool
operator != (const mbsleep_t& t1, const long& t2)
{
	return (t1.tv_sec != t2) || (t1.tv_nsec != 0);
}

#endif // HAVE_NANOSLEEP
