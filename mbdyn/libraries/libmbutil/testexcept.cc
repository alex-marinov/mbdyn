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

class Err : public std::exception, public Err_base {
protected:
	std::string reason;

public:
	Err(MBDYN_EXCEPT_ARGS_DECL, const std::string& r)
	: Err_base(MBDYN_EXCEPT_ARGS_PASSTHRU), reason(r) {};
	virtual ~Err(void) throw() {};

	const char* what(void) const throw () {
		return (Err_base::str() + " called Err (" + reason + ")").c_str();
	};
};

int
main(int argc, const char *argv[])
{
	if (argc > 1) {
		// passa la stringa CVS header...
		throw Err(MBDYN_EXCEPT_ARGS("$Header$"), argv[1]);

	} else {
		// ... non la passa
		throw Err(MBDYN_EXCEPT_ARGS(""), "for testing");
	}

	return 0;
}
