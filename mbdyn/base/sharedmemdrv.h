/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/sockdrv.h,v 1.29 2017/01/12 14:46:10 masarati Exp $ */
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

/* SharedMem driver */

#ifndef SHAREDMEMDRV_H
#define SHAREDMEMDRV_H

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

extern const char *MBDynSharedMemDrivePath;

#include "auth.h"
#include "filedrv.h"
#include "sock.h"

/* SharedMemDrive - begin */

class SharedMemDrive : public FileDrive {
protected:
   	enum {
      		DEFAULT     = 0x0000,
      		INCREMENTAL = 0x0001,
      		IMPULSIVE   = 0x0002
   	};

	const char *Name;
   	AuthMethod *auth;

   	int *pFlags;

   	boost::interprocess::shared_memory_object shm;
   	boost::interprocess::mapped_region region;
   	void * addr;
   	shared_memory_buffer *data;

	void Init(void);

public:
   	SharedMemDrive(unsigned int uL,
                const DriveHandler* pDH,
	            const char *path,
		        integer nd,
		        const std::vector<doublereal>& v0);

   	virtual ~SharedMemDrive(void);

   	/* Scrive il contributo del DriveCaller al file di restart */
   	virtual std::ostream& Restart(std::ostream& out) const;

   	virtual void ServePending(const doublereal& t);
};

/* SharedMemDrive - end */

#endif /* USE_SOCKET */

class DataManager;
class MBDynParser;

struct SharedMemDR : public DriveRead {
public:
	virtual Drive *
	Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP);
};


#endif /* SHAREDMEMDRV_H */

