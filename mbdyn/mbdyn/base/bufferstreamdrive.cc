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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "filedrv.h"
#include "streamdrive.h"
#include "bufferstreamdrive.h"

BufferStreamDrive_base::BufferStreamDrive_base(unsigned int uL,
	const DriveHandler* pDH,
	integer nd, const std::vector<doublereal>& v0,
	StreamDrive::Modifier *pMod,
	unsigned int ie,
	StreamDriveEcho *pSDE)
: StreamDrive(uL, pDH, "buffer", nd, v0, true, pMod),
InputEvery(ie), InputCounter(ie - 1),
pSDE(pSDE)
{
	// NOTE: InputCounter is set to InputEvery - 1 so that input
	// is expected at initialization (initial time) and then every
	// InputEvery steps; for example, for InputEvery == 4, input
	// is expected at:
	//	initial time
	//	initial time + 4 * timestep
	//	initial time + 8 * timestep
	ASSERT(InputEvery > 0);

	InputCounter -= InputEvery;

	if (pSDE != 0) {
		pSDE->Init("BufferStreamDrive_base", uLabel, nd);
	}
}

BufferStreamDrive_base::~BufferStreamDrive_base(void)
{
	if (pSDE != 0) {
		delete pSDE;
	}
}

const integer
BufferStreamDrive_base::GetBufSize(void) const
{
	return size;
}

void
BufferStreamDrive_base::ServePending(const doublereal& t)
{
	/* read only every InputEvery steps */
	InputCounter++;
	if (InputCounter != InputEvery) {
		return;
	}
	InputCounter = 0;

	if (pSDE != 0) {
		pSDE->EchoPrepare(&pdVal[1], iNumDrives);
	}

	// copy values from buffer
	pMod->Modify(&pdVal[1], GetBufRaw());

	if (pSDE != 0) {
		pSDE->Echo(&pdVal[1], iNumDrives);
	}
}


BufferStreamDrive::BufferStreamDrive(unsigned int uL,
	const DriveHandler* pDH,
	integer nd, const std::vector<doublereal>& v0,
	StreamDrive::Modifier *pMod,
	unsigned int ie,
	StreamDriveEcho *pSDE)
: BufferStreamDrive_base(uL, pDH, nd, v0, pMod, ie, pSDE),
buffer(nd)
{
	NO_OP;
}

BufferStreamDrive::~BufferStreamDrive(void)
{
	NO_OP;
}

const doublereal *
BufferStreamDrive::GetBufRaw(void)
{
	// paranoid sanity check: callers of GetBuf() could have altered the size of the buffer...
	ASSERT(buffer.size() == iNumDrives);

	return &buffer[0];
}

std::vector<doublereal>&
BufferStreamDrive::GetBuf(void)
{
	// paranoid sanity check: callers of GetBuf() could have altered the size of the buffer...
	ASSERT(buffer.size() == iNumDrives);

	return buffer;
}

/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
BufferStreamDrive::Restart(std::ostream& out) const
{
	// input every, echo, ...
	out << "  file: " << uLabel << ", buffer stream, type, stl, " << iNumDrives << ";" << std::endl;
	return out;
}
   

BufferStreamDriveRaw::BufferStreamDriveRaw(unsigned int uL,
	const DriveHandler* pDH,
	integer nd, const std::vector<doublereal>& v0,
	StreamDrive::Modifier *pMod,
	unsigned int ie,
	StreamDriveEcho *pSDE,
	bool bOwnsMemory)
: BufferStreamDrive_base(uL, pDH, nd, v0, pMod, ie, pSDE),
m_bOwnsMemory(bOwnsMemory),
m_pBuffer(0)
{
	if (m_bOwnsMemory) {
		m_pBuffer = new doublereal[nd];
	}
}

BufferStreamDriveRaw::~BufferStreamDriveRaw(void)
{
	if (m_bOwnsMemory) {
		delete[] m_pBuffer;
	}
}

bool
BufferStreamDriveRaw::bOwnsMemory(void) const
{
	return m_bOwnsMemory;
}

void
BufferStreamDriveRaw::SetBufRaw(integer n, const doublereal *p)
{
	if (n != size) {
		// error
		std::ostringstream os;
		os << "setting buffer pointer in BufferStreamDriveRaw of wrong size (original=" << size << ", new=" << n << ")";
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, os.str());
	}

	if (m_bOwnsMemory) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "setting buffer pointer in BufferStreamDriveRaw that owns its memory");
	}

	if (m_pBuffer != 0) {
		// error; maybe we could simply replace it, couldn't we?
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "setting buffer pointer in BufferStreamDriveRaw that has already been set");
	}

	m_pBuffer = p;
}

const doublereal *
BufferStreamDriveRaw::GetBufRaw(void)
{
	ASSERT(m_pBuffer != 0);

	return m_pBuffer;
}

/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
BufferStreamDriveRaw::Restart(std::ostream& out) const
{
	// input every, echo, ...
	out << "  file: " << uLabel << ", buffer stream, type, raw, owns memory, " << (m_bOwnsMemory ? "yes" : "no" ) << iNumDrives << ";" << std::endl;
	return out;
}
   

/* legge i drivers tipo stream */

static Drive *
ReadBufferStreamDrive(const DataManager *pDM, MBDynParser& HP, unsigned uLabel)
{
	enum {
		STL,
		RAW
	} eType = STL;
	bool bOwnsMemory(true);

	if (HP.IsKeyWord("type")) {
		if (HP.IsKeyWord("raw")) {
			eType = STL;
			if (HP.IsKeyWord("owns" "memory")) {
				bOwnsMemory = HP.GetYesNoOrBool();
			}

		} else if (!HP.IsKeyWord("stl")) {
			silent_cerr("BufferStreamDrive"
				"(" << uLabel << "\"): "
				"invalid type at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	unsigned int InputEvery = 1;
	if (HP.IsKeyWord("input" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("BufferStreamDrive"
				"(" << uLabel << "\"): "
				"invalid \"input every\" value " << i
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		InputEvery = (unsigned int)i;
	}

	StreamDriveEcho *pSDE = ReadStreamDriveEcho(pDM, HP);

	int idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("BufferStreamDrive"
			"(" << uLabel << "\"): "
			"illegal number of channels " << idrives
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<doublereal> v0;
	if (HP.IsKeyWord("initial" "values")) {
		v0.resize(idrives);
		for (int i = 0; i < idrives; i++) {
			v0[i] = HP.GetReal();
		}
	}

	StreamDrive::Modifier *pMod(0);
	if (HP.IsKeyWord("modifier")) {
		pMod = ReadStreamDriveModifier(HP, idrives);
	}

	Drive* pDr = 0;
	switch (eType) {
	case STL:
		SAFENEWWITHCONSTRUCTOR(pDr, BufferStreamDrive,
			BufferStreamDrive(uLabel,
				pDM->pGetDrvHdl(),
				idrives, v0, pMod,
				InputEvery, pSDE));
		break;

	case RAW:
		SAFENEWWITHCONSTRUCTOR(pDr, BufferStreamDriveRaw,
			BufferStreamDriveRaw(uLabel,
				pDM->pGetDrvHdl(),
				idrives, v0, pMod,
				InputEvery, pSDE, bOwnsMemory));
		break;

	default:
		ASSERT(0);
		break;
	}

	return pDr;
}


Drive *
BufferStreamDR::Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP)
{
	return ReadBufferStreamDrive(pDM, HP, uLabel);
}
