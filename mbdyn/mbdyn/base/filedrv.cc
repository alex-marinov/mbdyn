/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* file driver */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <fstream>

#include "dataman.h"
#include "mbdyn.h"
#include "filedrv.h"
#include "fixedstep.h"
#include "varstep.h"
#include "sockdrv.h"
#include "streamdrive.h"
#ifdef USE_SOCKET
#include "socketstreamdrive.h"
#endif // USE_SOCKET
#include "rtai_in_drive.h"

/* FileDrive - begin */

	

FileDrive::FileDrive(unsigned int uL, const DriveHandler* pDH,
	const std::string& s,
	integer nd, const std::vector<doublereal>& v0)
: Drive(uL, pDH), sFileName(s), iNumDrives(nd), pdVal(NULL)
{
   	SAFENEWARR(pdVal, doublereal, nd + 1);
	if (v0.empty()) {
   		for (int iCnt = 0; iCnt <= nd; iCnt++) {
      			pdVal[iCnt] = 0.;
   		}

	} else {
		ASSERT(v0.size() == unsigned(nd));
		pdVal[0] = 0.;
		for (int iCnt = 0; iCnt < nd; iCnt++) {
			pdVal[iCnt + 1] = v0[iCnt];
		}
	}
}


FileDrive::~FileDrive(void)
{
   	if (pdVal != NULL) {
      		SAFEDELETEARR(pdVal);
   	}
}


Drive::Type
FileDrive::GetDriveType(void) const
{
	return Drive::FILEDRIVE;
}

doublereal
FileDrive::dGet(const doublereal& /* t */ , int i) const
{
   	ASSERT(i > 0 && i <= iNumDrives);
   	return pdVal[i];
}

/* FileDrive - end */


/* FileDriveCaller - begin */

FileDriveCaller::FileDriveCaller(const DriveHandler* pDH,
		const FileDrive* p, integer i,
		const doublereal& da)
: DriveCaller(pDH), pFileDrive(p), iNumDrive(i), dAmplitude(da)
{
	ASSERT(pFileDrive != NULL);
	ASSERT(iNumDrive > 0 && iNumDrive <= pFileDrive->iGetNumDrives());
}


FileDriveCaller::~FileDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller* FileDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC,
			FileDriveCaller,
			FileDriveCaller(pDrvHdl,  pFileDrive,
				iNumDrive, dAmplitude));

	return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
FileDriveCaller::Restart(std::ostream& out) const
{
	out << " file, " << pFileDrive->GetLabel()
		<< ", " << iNumDrive;
	if (dAmplitude != 1.) {
		out << ", amplitude, " << dAmplitude;
	}
	return out;
}

/* FileDriveCaller - end */

