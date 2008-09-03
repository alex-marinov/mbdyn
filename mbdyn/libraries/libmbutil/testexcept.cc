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

#include <iostream>
#include "mbexcept.h"

class ErrPlain : public Err_base {
public:
	ErrPlain(MBDYN_EXCEPT_ARGS_DECL)
	: Err_base(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	virtual ~ErrPlain(void) throw() {};
};

class ErrReason : public Err_base {
private:
	std::string reason;

public:
	ErrReason(MBDYN_EXCEPT_ARGS_DECL, const std::string& r)
	: Err_base(MBDYN_EXCEPT_ARGS_PASSTHRU), reason(r) {};
	virtual ~ErrReason(void) throw() {};

	const char* what(void) const throw () {
		return (Err_base::str() + ": " + reason).c_str();
	};
};

int
main(int argc, const char *argv[])
{
	if (argc == 1) {
		throw ErrPlain(MBDYN_EXCEPT_ARGS("$Header$"));
	} else {
		throw ErrReason(MBDYN_EXCEPT_ARGS("$Header$"), argv[1]);
	}

	return 0;
}
