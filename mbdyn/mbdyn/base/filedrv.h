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

/* file driver */

#ifndef FILEDRV_H
#define FILEDRV_H

#include "mbdyn.h"
#include "drive.h"

/* FileDrive - begin */

class FileDrive : public Drive {
public:
	enum Type {
		UNKNOWN = -1,

		FIXEDSTEP = 0,
		VARIABLESTEP,
		SOCKET,
		STREAM,
		SOCKETSTREAM,
		RTAI_IN,

		LASTFILEDRIVE
	};

protected:
	std::string sFileName;
	integer iNumDrives;
   	doublereal* pdVal;

public:
	FileDrive(unsigned int uL, const DriveHandler* pDH,
		const std::string& s,
		integer nd, const std::vector<doublereal>& v0);
	virtual ~FileDrive(void);

	virtual Drive::Type GetDriveType(void) const;

	virtual FileDrive::Type GetFileDriveType(void) const = 0;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;

	virtual inline integer iGetNumDrives(void) const;

	virtual doublereal dGet(const doublereal& t, int i = 1) const;
};


inline integer
FileDrive::iGetNumDrives(void) const
{
	return iNumDrives;
}

/* FileDrive - end */


/* FileDriveCaller - begin */

class FileDriveCaller : public DriveCaller {
protected:
	FileDrive* pFileDrive;
	integer iNumDrive;
	doublereal dAmplitude;

public:
	FileDriveCaller(const DriveHandler* pDH, const FileDrive* p,
			integer i, const doublereal& da);
	virtual ~FileDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Restituisce il valore del driver */
	virtual inline doublereal dGet(const doublereal& dVal) const;
	/* virtual inline doublereal dGet(void) const; */
};


inline doublereal
FileDriveCaller::dGet(const doublereal& dVal) const
{
	return dAmplitude*pFileDrive->dGet(dVal, iNumDrive);
}

/* FileDriveCaller - end */

class DataManager;
class MBDynParser;

extern Drive* ReadFileDriver(DataManager* pDM,
		MBDynParser& HP, unsigned int uLabel);

#endif /* FILEDRV_H */

