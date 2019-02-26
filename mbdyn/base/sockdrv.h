/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef SOCKDRV_H
#define SOCKDRV_H

#ifdef USE_SOCKET

const unsigned short int MBDynSocketDrivePort = 5555;
extern const char *MBDynSocketDrivePath;

/* include del programma */

#include "auth.h"
#include "filedrv.h"

/* SocketDrive - begin */

class SocketDrive : public FileDrive {
protected:
   	enum {
      		DEFAULT     = 0x0000,
      		INCREMENTAL = 0x0001,
      		IMPULSIVE   = 0x0002
   	};

	int type;
	union {
 	  	unsigned short int Port;
		const char *Path;
	} data;
   	AuthMethod *auth;

   	int *pFlags;

   	int sock;

	void Init(void);

public:
   	SocketDrive(unsigned int uL, const DriveHandler* pDH,
	        unsigned short int p, AuthMethod* a,
		integer nd, const std::vector<doublereal>& v0);

   	SocketDrive(unsigned int uL, const DriveHandler* pDH,
	        const char *path,
		integer nd, const std::vector<doublereal>& v0);

   	virtual ~SocketDrive(void);

   	/* Scrive il contributo del DriveCaller al file di restart */
   	virtual std::ostream& Restart(std::ostream& out) const;

   	virtual void ServePending(const doublereal& t);
};

/* SocketDrive - end */

#endif /* USE_SOCKET */

class DataManager;
class MBDynParser;

struct SocketDR : public DriveRead {
public:
	virtual Drive *
	Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP);
};


#endif /* SOCKDRV_H */

