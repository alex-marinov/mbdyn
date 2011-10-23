/* $Header$ */
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

#include "demangle.h"

#include <malloc.h>
#ifdef HAVE_CXXABI_H
#include <cxxabi.h>
#endif

std::string
mbdyn_demangle(const char *name)
{
	std::string ret;

#ifdef HAVE_CXXABI_H
	const char* const demangled_name = name;
	int status = -1;
	char* res = abi::__cxa_demangle(name, NULL, NULL, &status);
	if (status == 0) {
		demangled_name = res;
	}
	ret = demangled_name;
#else
	ret = name;
#endif

	return ret;
}

std::string
mbdyn_demangle(const std::type_info& t)
{
	return mbdyn_demangle(t.name());
}

