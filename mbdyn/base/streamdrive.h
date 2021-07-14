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

#ifndef STREAMDRIVE_H
#define STREAMDRIVE_H

#include "filedrv.h"
#include "bufmod.h"

/* StreamDrive - begin */

class StreamDrive : public FileDrive {
public:
	class Modifier {
	public:
		Modifier(void);
		virtual ~Modifier(void);
		virtual size_t GetSize(void) const = 0;
		virtual void Modify(doublereal *out, const void *in) const = 0;
	};

	class Copy : public StreamDrive::Modifier {
	protected:
		integer m_iND;

	public:
		Copy(integer iND);
		size_t GetSize(void) const;
		void Modify(doublereal *out, const void *in) const;
	};

protected:
	/* MBox buffer */
	int size;
	std::vector<char> buf;

	bool create;

	const Modifier *pMod;

public:
   	StreamDrive(unsigned int uL,
		const DriveHandler* pDH,
		const std::string& sFileName,
		integer nd, const std::vector<doublereal>& v0,
		bool c, StreamDrive::Modifier *pmod);

   	virtual ~StreamDrive(void);
	void SetModifier(const Modifier *p);
	const StreamDrive::Modifier *pGetModifier(void) const;
};

extern StreamDrive::Modifier *
ReadStreamDriveModifier(MBDynParser& HP, integer nDrives);

/* StreamDrive - end */

/* StreamDriveEcho - begin */

class StreamDriveEcho {
private:
	const DriveHandler *pDrvHdl;
	std::string sOutFileName;
	std::vector<doublereal> echoBuf;
	std::ofstream outFile;
	int iPrecision;
	doublereal dShift;

public:
	StreamDriveEcho(const DriveHandler *pDrvHdl, std::string& sOutFileName, int iPrecision, doublereal dShift);
	~StreamDriveEcho(void);
	bool bIsEcho(void) const;
	bool Init(const std::string& msg, unsigned uLabel, unsigned nChannels);
	void EchoPrepare(const doublereal *pbuf, unsigned nChannels);
	void Echo(const doublereal *pbuf, unsigned nChannels);
};

extern StreamDriveEcho *
ReadStreamDriveEcho(const DataManager *pDM, MBDynParser& HP);

/* StreamDriveEcho - end */

class StreamDriveCopyCast : public StreamDrive::Modifier
{
protected:
	size_t m_size;
	std::vector<BufCast *> m_data;

public:
	StreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data);
	~StreamDriveCopyCast(void);

	size_t GetSize(void) const;
	void Modify(doublereal *out, const void *in) const;
};

#endif // STREAMDRIVE_H
