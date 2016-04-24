/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#ifndef BUFFERSTREAMDRIVE_H
#define BUFFERSTREAMDRIVE_H

#include "streamdrive.h"

#ifdef USE_SOCKET

#include "usesock.h"

/* BufferStreamDrive - begin */

class BufferStreamDrive : public StreamDrive {
protected:
	unsigned int InputEvery;
	unsigned int InputCounter;

	std::vector<doublereal> buffer;

	StreamDriveEcho *pSDE;

public:
	BufferStreamDrive(unsigned int uL,
		const DriveHandler* pDH,
		integer nd, const std::vector<doublereal>& v0,
		StreamDrive::Modifier *pMod,
		unsigned int ie,
		StreamDriveEcho *pSDE);

	virtual ~BufferStreamDrive(void);

	const std::vector<doublereal>& GetBuf(void);

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

/* BufferStreamDrive - end */

#endif // USE_SOCKET

class DataManager;
class MBDynParser;

struct BufferStreamDR : public DriveRead {
public:
	BufferStreamDR(void) { NO_OP; };

	virtual Drive *
	Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP);
};

#endif /* BUFFERSTREAMDRIVE_H */

