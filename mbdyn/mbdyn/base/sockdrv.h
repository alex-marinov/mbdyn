/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifndef SOCKDRV_H
#define SOCKDRV_H

#ifdef USE_SOCKET_DRIVES

const unsigned short int MBDynSocketDrivePort = 5555;


/* include del programma */

#include "auth.h"
#include "filedrv.h"
#include <fstream.h>


/* SocketDrive - begin */

class SocketDrive : public FileDrive {
protected:
   	enum { 
      		DEFAULT     = 0x0000,
      		INCREMENTAL = 0x0001,
      		IMPULSIVE   = 0x0002
   	};
   
   	unsigned short int Port;
   	AuthMethod* auth;
   
   	doublereal* pdVal;
   	int* pFlags;
   
   	int sock;
   
public:
   	SocketDrive(unsigned int uL, const DriveHandler* pDH,
	            unsigned short int p, AuthMethod* a, integer nd);
   
   	virtual ~SocketDrive(void);
   
   	virtual FileDriveType::Type GetFileDriveType(void) const;

   	/* Scrive il contributo del DriveCaller al file di restart */
   	virtual ostream& Restart(ostream& out) const;
   
   	virtual const doublereal& dGet(const doublereal& t, int i = 1) const;
   
   	virtual void ServePending(void);
};

/* SocketDrive - end */

#endif /* USE_SOCKET_DRIVES */

#endif /* SOCKDRV_H */

