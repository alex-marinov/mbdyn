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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#if defined(USE_RUNTIME_LOADING) && defined(HAVE_LTDL_H)
#include <ltdl.h>
#endif // USE_RUNTIME_LOADING && HAVE_LTDL_H

#include "dataman.h"
#include "mbdefs.h"

static bool done = false;

void
module_initialize(void)
{
	if (::done) {
		return;
	}

	::done = true;

#ifdef USE_RUNTIME_LOADING
	if (lt_dlinit()) {
		silent_cerr("unable to initialize run-time loading" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	 * NOTE: this macro is defined in mbdefs.h
	 */
	if (lt_dlsetsearchpath(MODULE_LOADPATH) != 0) {
		silent_cerr("unable to initialize load path" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // USE_RUNTIME_LOADING
}

void
module_finalize(void)
{
	if (!::done) {
		return;
	}

	::done = false;

#ifdef USE_RUNTIME_LOADING
	lt_dlexit();
#endif // USE_RUNTIME_LOADING
}
