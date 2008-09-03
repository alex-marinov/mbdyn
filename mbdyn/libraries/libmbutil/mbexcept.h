/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifndef MBEXCEPT_H
#define MBEXCEPT_H

#include <exception>
#include <stdexcept>
#include <sstream>

// Use this macro instead of the required set of args when declaring the constructor
// of classes derived from Err_base
#define MBDYN_EXCEPT_ARGS_DECL \
	const char *file, int line, const char *func, const char *cvs_header
// Use this macro to pass the required set of args thru to Err_base from the constructor
// of derived classes
#define MBDYN_EXCEPT_ARGS_PASSTHRU \
	file, line, func, cvs_header
// Use this macro to pass the required set of args to error classes derived from Err_base
#if __GNUC__ >= 2
#define MBDYN_EXCEPT_ARGS( cvs_header ) \
	__FILE__ , __LINE__ , __PRETTY_FUNCTION__, (cvs_header)
#else // ! __GNUC__
// FIXME: need to detect whether __func__ (C99) is available
#define MBDYN_EXCEPT_ARGS( cvs_header ) \
	__FILE__ , __LINE__ , "(unknown)", (cvs_header)
#endif // ! __GNUC__

// Example:

#if 0

class MyErr : public Err_base {
public:
	MyErr(MBDYN_EXCEPT_ARGS_DECL)
	: Err_base(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	virtual MyErr(void) throw() {};
	const char * what(void) const throw() {
		return (Err_base::str() + ": MyErr").c_str();
	};
};

class MyErrWithMessage : public Err_base {
private:
	std::string m_msg;

public:
	MyErrWithMessage(MBDYN_EXCEPT_ARGS_DECL, const std::string& msg)
	: Err_base(MBDYN_EXCEPT_ARGS_PASSTHRU), m_msg(msg) {};
	virtual MyErrWithMessage(void) throw() {};
	const char * what(void) const throw() {
		return (Err_base::str() + ": MyErrWithMessage(" + m_msg + ")").c_str();
	};
};

int
func(int i)
{
	if (i == 0) {
		throw MyErr(MBDYN_EXCEPT_ARGS("$Header$"), "i is zero!");
	}

	return i;
}

#endif

class Err_base : public std::exception {
private:
	std::string s;

public:
	Err_base(MBDYN_EXCEPT_ARGS_DECL)
	{
		std::stringstream ss;
		ss << "function=`" << func << "' "
			"[" << file << ":" << line;
		if (cvs_header && cvs_header[0]) {
			ss << " " << cvs_header;
		}
		ss << "]";
		s = ss.str();
	};

	virtual ~Err_base(void) throw() {};

	const std::string& str(void) const {
		return s;
	};

	const char * what(void) const throw() {
		return s.c_str();
	};
};

#endif // MBEXCEPT_H

