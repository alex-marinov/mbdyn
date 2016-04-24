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

BufferStreamDrive::BufferStreamDrive(unsigned int uL,
	const DriveHandler* pDH,
	integer nd, const std::vector<doublereal>& v0,
	StreamDrive::Modifier *pMod,
	unsigned int ie,
	StreamDriveEcho *pSDE)
: StreamDrive(uL, pDH, 0, nd, v0, true, pMod),
InputEvery(ie), InputCounter(ie - 1),
buffer(nd),
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

	pSDE->Init("BufferStreamDrive", uLabel, nd);
}

BufferStreamDrive::~BufferStreamDrive(void)
{
	if (pSDE != 0) {
		delete pSDE;
	}
}

const std::vector<doublereal>&
BufferStreamDrive::GetBuf(void)
{
	return buffer;
}

/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
BufferStreamDrive::Restart(std::ostream& out) const
{
	out << "  file: " << uLabel << ", buffer stream, " << iNumDrives << ";" << std::endl;
	return out;
}
   
void
BufferStreamDrive::ServePending(const doublereal& t)
{
	/* read only every InputEvery steps */
	InputCounter++;
	if (InputCounter != InputEvery) {
		return;
	}
	InputCounter = 0;
	
	pSDE->EchoPrepare(&pdVal[1], iNumDrives);

	// copy values from buffer
	pMod->Modify(&pdVal[1], &buffer[0]);

	pSDE->Echo(&pdVal[1], iNumDrives);
}


/* legge i drivers tipo stream */

static Drive *
ReadBufferStreamDrive(const DataManager *pDM, MBDynParser& HP, unsigned uLabel)
{
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
	SAFENEWWITHCONSTRUCTOR(pDr, BufferStreamDrive,
		BufferStreamDrive(uLabel,
			pDM->pGetDrvHdl(),
			idrives, v0, pMod,
			InputEvery, pSDE));

	return pDr;
}


Drive *
BufferStreamDR::Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP)
{
	return ReadBufferStreamDrive(pDM, HP, uLabel);
}
