/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003
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


#ifdef HAVE_CONFIG_H
#include <mbconfig.h> 		/* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/iostream>
#include <ac/f2c.h>

#include <dataman.h>
#include "loadable.h"

#include "module-aerodyn.h"

extern "C" int
__FC_DECL__(getrotorparams)(F_REAL *Omega, F_REAL *gamma, F_REAL *VHUB,
		F_REAL *tau)
{
	*Omega = 0.;
	*gamma = 0.;
	*VHUB = 0.;
	*tau = 0;

	return 0;
}

extern "C" int
__FC_DECL__(getbladeparams)(F_REAL *psi)
{
	*psi = 0.;

	return 0;
}

extern "C" int
__FC_DECL__(getelemparams)(F_INTEGER *MulTabLoc, F_REAL *phi, F_REAL *radius, 
		F_REAL *XGRND, F_REAL *YGRND, F_REAL *ZGRND)
{
	*MulTabLoc = 1;
	*phi = 0.;
	*radius = 0.;
	*XGRND = 0.;
	*YGRND = 0.;
	*ZGRND = 0.;
	
	return 0;
}

extern "C" int
__FC_DECL__(getvnvt)(F_REAL *VX, F_REAL *VY, F_REAL *VZ, F_REAL *VT,
		F_REAL *VNW, F_REAL *VNE)
{
	*VX = 0.;
	*VY = 0.;
	*VZ = 0.;
	*VT = 0.;
	*VNW = 0.;
	*VNE = 0.;

	return 0;
}

extern "C" int
__FC_DECL__(usrmes)(F_LOGICAL *Logical, F_CHAR msg[],
		F_INTEGER *code, F_CHAR level[])
{
	if (*Logical) {
		silent_cerr("AeroDyn [" << level << ":"
			<< *code << "] " << msg << std::endl);
	}

	return 0;
}

extern "C" int
__FC_DECL__(abort)(void)
{
	throw ErrGeneric();
}

