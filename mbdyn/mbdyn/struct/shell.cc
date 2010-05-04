/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2010
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Riccardo Vescovini	<vescovini@aero.polimi.it>
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

/*
 * Inspired by
 * Wojciech Witkowski
 * "4-Node combined shell element with semi-EAS-ANS strain interpolations
 * in 6-parameter shell theories with drilling degrees of freedom"
 * Comput Mech (2009) 43:307­319 DOI 10.1007/s00466-008-0307-x
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "constltp_impl.h"
#include "tpldrive_impl.h"
#include "shell.h"
#include "mynewmem.h"


Shell::Shell(unsigned uLabel, const DofOwner* pDO, flag fOut)
: Elem(uLabel, fOut),
ElemGravityOwner(uLabel, fOut),
ElemWithDofs(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut)
{
	NO_OP;
}

Shell::~Shell(void)
{
	NO_OP;
}

