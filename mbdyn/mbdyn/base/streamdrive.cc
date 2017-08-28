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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "filedrv.h"
#include "streamdrive.h"

#include "bufmod.h"

/* StreamDrive - begin */

StreamDrive::StreamDrive(unsigned int uL,
	const DriveHandler* pDH,
	const std::string& sFileName,
	integer nd, const std::vector<doublereal>& v0,
	bool c, StreamDrive::Modifier *pmod)
: FileDrive(uL, pDH, sFileName, nd, v0),
create(c),
pMod(pmod)
{
   	ASSERT(nd > 0);

	if (pMod == 0) {
		pMod = new StreamDrive::Copy(nd);
	}

	// initialize mailbox and so on
	size = pMod->GetSize();
	buf.resize(size);
}

StreamDrive::~StreamDrive(void) 
{
	if (pMod) {
		delete pMod;
	}
}

void
StreamDrive::SetModifier(const Modifier *p)
{
	// FIXME: what if is not 0?
	if (pMod != 0) {
		delete pMod;
	}

	pMod = p;
}

const StreamDrive::Modifier *StreamDrive::pGetModifier(void) const
{
	return pMod;
}

StreamDrive::Modifier::Modifier(void)
{
	NO_OP;
}

StreamDrive::Modifier::~Modifier(void)
{
	NO_OP;
}

StreamDrive::Copy::Copy(integer iND)
: m_iND(iND)
{
	NO_OP;
}

size_t
StreamDrive::Copy::GetSize(void) const
{
	return (sizeof(doublereal)*m_iND);
}

void
StreamDrive::Copy::Modify(doublereal *out, const void *in) const
{
	const doublereal *pd_in = (const doublereal *)in;
	for (int i = 0; i < m_iND; i++) {
		out[i] = pd_in[i];
	}
}

/* StreamDrive - end */


/* StreamDriveCopyCast - begin */

/* moved to streamdrive.h FIXME: DELETE
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
};*/

StreamDriveCopyCast::StreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data)
: m_size(size), m_data(data)
{
#ifdef DEBUG
	size_t minsize = m_data[m_data.size() - 1]->offset() + m_data[m_data.size() - 1]->size();
	ASSERT(size >= minsize);
#endif
}

StreamDriveCopyCast::~StreamDriveCopyCast(void)
{
	for (std::vector<BufCast *>::iterator i = m_data.begin(); i != m_data.end(); ++i) {
		delete *i;
	}
}

size_t
StreamDriveCopyCast::GetSize(void) const
{
	return m_size;
}

void
StreamDriveCopyCast::Modify(doublereal *out, const void *in) const
{
	for (size_t i = 0; i != m_data.size(); ++i) {
		out[i] = m_data[i]->cast(in);
	}
}

/* StreamDriveCopyCast - end */

StreamDrive::Modifier *
ReadStreamDriveModifier(MBDynParser& HP, integer nDrives)
{
	StreamDrive::Modifier *pSDM(0);

	if (HP.IsKeyWord("copy" "cast")) {
		std::vector<BufCast *> data(nDrives);
		ReadBufCast(HP, data);
		size_t minsize = data[data.size() - 1]->offset() + data[data.size() - 1]->size();
		size_t size = minsize;
		if (HP.IsKeyWord("size")) {
			integer i = HP.GetInt();
			if (i <= 0) {
				silent_cerr("ReadStreamDriveModifier: invalid size " << i
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			size = size_t(i);
			if (size < minsize) {
				silent_cerr("ReadStreamDriveModifier: size " << size
					<< " is less than min size " << minsize
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		pSDM = new StreamDriveCopyCast(size, data);

	} else if (!HP.IsKeyWord("copy")) {
		// TODO: "copy" with byte swap for network-independent stuff...
		silent_cerr("ReadStreamDriveModifier: unknown modifier type at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pSDM;
}

/* StreamDriveEcho - begin */

StreamDriveEcho::StreamDriveEcho(const DriveHandler *pDrvHdl, std::string& sOutFileName, int iPrecision, doublereal dShift)
: pDrvHdl(pDrvHdl),
sOutFileName(sOutFileName),
iPrecision(iPrecision),
dShift(dShift)
{
	NO_OP;
}

StreamDriveEcho::~StreamDriveEcho(void)
{
	NO_OP;
}

bool
StreamDriveEcho::Init(const std::string& msg, unsigned uLabel, unsigned nChannels)
{
	outFile.open(sOutFileName.c_str());
	if (!outFile) {
		silent_cerr(msg << "(" << uLabel << "): "
			"unable to open echo file '" << sOutFileName << "'" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	echoBuf.resize(nChannels);

	if (iPrecision > 0) {
		outFile.precision(iPrecision);
	}
	outFile.setf(std::ios::scientific);

	outFile
		<< "# generated by SocketStreamDrive(" << uLabel << ")"
		<< std::endl;
	if (nChannels == 1) {
		outFile
			<< "# Time, Channel #1"
			<< std::endl;

	} else {
		outFile
			<< "# Time, Channels #1-" << nChannels
			<< std::endl;
	}

	return true;
}

void
StreamDriveEcho::EchoPrepare(const doublereal *pbuf, unsigned size)
{
	for (unsigned i = 1; i < size; i++) {
		echoBuf[i] = pbuf[i];
	}
}

void
StreamDriveEcho::Echo(const doublereal *pbuf, unsigned size)
{
	for (unsigned i = 0; i < size ; i++) {
		if (pbuf[i] != echoBuf[i]) {
			// changed; need to write
			outFile << (pDrvHdl->dGetTime() + dShift);
			for (unsigned j = 0; j < size ; j++) {
				outFile << " " << pbuf[j];
			}
			outFile << std::endl;
			break;
		}
	}
}

StreamDriveEcho *
ReadStreamDriveEcho(const DataManager *pDM, MBDynParser& HP)
{
	StreamDriveEcho *pSDE(0);

	std::string sOutFileName;
	int iPrecision = -1;
	doublereal dShift = 0.;

	if (HP.IsKeyWord("echo")) {
		const char *s = HP.GetFileName();
		if (s == NULL) {
			silent_cerr("ReadStreamDriveEcho: "
				"unable to parse echo file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		sOutFileName = s;

		if (HP.IsKeyWord("precision")) {
			iPrecision = HP.GetInt();
			if (iPrecision <= 0) {
				silent_cerr("ReadStreamDriveEcho: "
					"invalid echo precision " << iPrecision
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		if (HP.IsKeyWord("shift")) {
			dShift = HP.GetReal();
		}

		pSDE = new StreamDriveEcho(pDM->pGetDrvHdl(), sOutFileName, iPrecision, dShift);
	}

	return pSDE;
}

/* StreamDriveEcho - end */

