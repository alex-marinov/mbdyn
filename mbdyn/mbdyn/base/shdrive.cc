/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Drive che usa  i gradi di liberta' */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "shdrive.h"

SHDriveCaller::SHDriveCaller(const DriveHandler* pDH,
	const DriveCaller *pFunc, const DriveCaller *pTrigger,
	const doublereal dVal0)
: DriveCaller(pDH), dVal0(dVal0)
{
	iSHDriveNumber = pDrvHdl->iSHInit(pFunc, pTrigger, dVal0);
}

SHDriveCaller::SHDriveCaller(const DriveHandler* pDH,
	integer iSHDriveNumber)
: DriveCaller(pDH), iSHDriveNumber(iSHDriveNumber)
{
	dVal0 = pDrvHdl->dGetSHVal0(iSHDriveNumber);
}

SHDriveCaller::~SHDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
SHDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		SHDriveCaller,
		SHDriveCaller(pDrvHdl, iSHDriveNumber));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
SHDriveCaller::Restart(std::ostream& out) const
{
	out
		<< " sample and hold, ",
			pDrvHdl->pGetSHFunc(iSHDriveNumber)->Restart(out)
		<< ", ",
			pDrvHdl->pGetSHTrigger(iSHDriveNumber)->Restart(out);
	if (dVal0) {
		out << ", initial value, " << dVal0;
	}

	return out;
}
