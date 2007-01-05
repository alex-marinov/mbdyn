/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/* socket driver */

#ifndef RTAI_IN_DRIVE_H
#define RTAI_IN_DRIVE_H

/* include del programma */

#include "streamdrive.h"

#ifdef USE_RTAI

/* RTAIInDrive - begin */

class RTAIInDrive : public StreamDrive {
protected:
	
	/* FIXME: store restart info as well */
	const char *host;
	unsigned long node;
	int port;
	void *mbx;

public:
   	RTAIInDrive(unsigned int uL,
			 const DriveHandler* pDH,
			 const char* const sFileName,
			 const char *h,
			 integer nd, bool c, unsigned long /*int*/ n);
   
   	virtual ~RTAIInDrive(void);
   
   	virtual FileDrive::Type GetFileDriveType(void) const;

   	/* Scrive il contributo del DriveCaller al file di restart */
   	virtual std::ostream& Restart(std::ostream& out) const;
   
   	virtual void ServePending(const doublereal& t);
};

/* RTAIInDrive - end */

#endif /* USE_RTAI */

class DataManager;
class MBDynParser;

extern Drive *
ReadRTAIInDrive(DataManager *pDM, MBDynParser& HP, unsigned int uLabel);


#endif /* RTAI_IN_DRIVE_H */

