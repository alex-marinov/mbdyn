/* $Header$ */
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

#include <iostream>
#include <sstream>
#include <cstdlib>
#include "except.h"

class ErrPlain : public MBDynErrBase {
public:
	ErrPlain(MBDYN_EXCEPT_ARGS_DECL)
	: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	virtual ~ErrPlain(void) throw() {};
};

class ErrCode: public MBDynErrBase {
private:
	std::string s;
public:
	ErrCode(MBDYN_EXCEPT_ARGS_DECL, int code = -1)
	: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {
		std::stringstream ss;
		ss << MBDynErrBase::what() << " code=" << code;
		s = ss.str();
	};
	virtual ~ErrCode(void) throw() {};

	const char* what(void) const throw () {
		return s.c_str();
	};
};

int
main(int argc, const char *argv[])
{
	std::cout << "usage: testexcept [reason [code]]"
		<< std::endl << std::endl;
	switch (argc) {
	case 1:
		throw ErrPlain(MBDYN_EXCEPT_ARGS);
		break;

	case 2:
		throw ErrPlain(MBDYN_EXCEPT_ARGS, argv[1]);
		break;

	default:
		throw ErrCode(MBDYN_EXCEPT_ARGS, argv[1], std::atoi(argv[2]));
		break;
	}

	return 0;
}
