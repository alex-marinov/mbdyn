/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include <cstdio>
#include <cstring>
#include "mbstrbuf.h"
//#include <stdio.h>

std::ostream&
mbstrbuf::stats(std::ostream& out)
{
	return out << "len=" << len << "; cursor=" << cursor << std::endl;
}


void
mbstrbuf::make_room(unsigned newlen)
{
	newlen = std::max(cursor + 2*newlen, 2*len);
	char *ptr = new char[newlen];
	memcpy(ptr, buf, cursor);
	delete[] buf;
	buf = ptr;
	len = newlen;
}

void 
mbstrbuf::return_cursor(unsigned newcursor)
{
	cursor = newcursor;
}


void
mbstrbuf::print_str(const char *str)
{
	int slen, buflen;

retry:;
	buflen = len - cursor;
	slen = strlen(str);
	if (slen >= buflen) {
		make_room(slen);
		goto retry;
	}

	memcpy(&buf[cursor], str, slen);
	cursor += slen;
	buf[cursor] = '\0';
}

void
mbstrbuf::print_double(const char *fmt, double d)
{
	int dlen, buflen;

retry:;
	buflen = len - cursor;
	dlen = snprintf(&buf[cursor], buflen, fmt, d);
	if (dlen >= buflen) {
		make_room(dlen);
		goto retry;
	}

	cursor += dlen;
}

const char *
mbstrbuf::get_buf(void) const
{
	return buf;
}

unsigned
mbstrbuf::get_len(void) const
{
	return len;
}

std::ostream&
operator << (std::ostream& out, const mbstrbuf& buf)
{
	return out << buf.buf;
}

#ifdef MAIN

int
main(void)
{
	mbstrbuf buf(10);
	buf.stats(std::cout) << std::endl;

	buf.print_double("%e", 10.2);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	return 0;
}

#endif // MAIN
