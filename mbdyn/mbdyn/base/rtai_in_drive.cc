/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_RTAI

#include <rtai_mbox_drive.h>

RTAIMailboxDrive::RTAIMailboxDrive(unsigned int uL, const DriveHandler* pDH,
			 integer nd)
: FileDrive(uL, pDH, "rtai_mailbox", nd), 
pdVal(NULL)
{
   	ASSERT(nd > 0);
   
   	SAFENEWARR(pdVal, doublereal, nd+1);

	/*
	 * initialize mailbox and so on
	 */
}

RTAIMailboxDrive::~RTAIMailboxDrive(void) 
{
   	if (pdVal != NULL) {
      		SAFEDELETEARR(pdVal);
   	}
}

void 
RTAIMailboxDrive::ServePending(void)
{
	/*
	 * store in pdVal the values of all the channels
	 * served by the mailbox
	 */
}

FileDrive::Type 
RTAIMailboxDrive::GetFileDriveType(void) const
{
   	return FileDrive::RTAI_MAILBOX;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
RTAIMailboxDrive::Restart(std::ostream& out) const
{
   	return out << "# RTAIMailboxDrive not implemented yet" << std::endl;
}
   
doublereal
RTAIMailboxDrive::dGet(const doublereal& /* t */ , int i) const
{
   	ASSERT(i > 0 && i <= iNumDrives);
   	return pdVal[i];
}

#endif /* USE_RTAI */

