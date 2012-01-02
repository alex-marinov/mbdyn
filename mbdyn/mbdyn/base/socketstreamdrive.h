/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/*
 * Author: Michele Attolico <attolico@aero.polimi.it>
 */

/* socket driver */

#ifndef SOCKETSTREAMDRIVE_H
#define SOCKETSTREAMDRIVE_H

#include "streamdrive.h"
#include "usesock.h"


/* SocketStreamDrive - begin */

class SocketStreamDrive : public StreamDrive {
protected:
	unsigned int InputEvery;
	bool bReceiveFirst;
	unsigned int InputCounter;

	UseSocket *pUS;
	int recv_flags;
	struct timeval SocketTimeout;

	std::string sOutFileName;
	std::ofstream outFile;
	int iPrecision;
	doublereal dShift;

public:
	SocketStreamDrive(unsigned int uL,
		const DriveHandler* pDH,
		UseSocket *pUS, bool c,
		const std::string& sFileName,
		integer nd, const std::vector<doublereal>& v0,
		unsigned int ie, bool bReceiveFirst,
		int flags,
		const struct timeval& st,
		const std::string& sOutFileName, int iPrecision,
		doublereal dShift);

	virtual ~SocketStreamDrive(void);

	virtual FileDrive::Type GetFileDriveType(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

/* SocketStreamDrive - end */

class DataManager;
class MBDynParser;

extern Drive* ReadSocketStreamDrive(DataManager* pDM,
		MBDynParser& HP, unsigned int uLabel);

#endif /* SOCKETSTREAMDRIVE_H */

