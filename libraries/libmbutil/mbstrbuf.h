/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef MBSTRBUF_H
#define MBSTRBUF_H

#include <iostream>

class mbstrbuf {
	friend std::ostream& operator << (std::ostream&, const mbstrbuf&);

private:
	unsigned len;
	unsigned cursor;
	char *buf;

public:
	mbstrbuf(unsigned l = 0) : len(l), cursor(0), buf(0) {
		if (len > 0) {
			buf = new char[len];
		}
	};

	std::ostream& stats(std::ostream& out);

	void make_room(unsigned newlen);
	void return_cursor(unsigned newcursor = 0);
	void print_str(const char *str);
	void print_double(const char *fmt, double d);
	const char *get_buf(void) const;
	unsigned get_len(void) const;
};

std::ostream& operator << (std::ostream& out, const mbstrbuf& buf);

#endif // MBSTRBUF_H
