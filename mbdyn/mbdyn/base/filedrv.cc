/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
#include <ac/fstream>

#include <dataman.h>
#include <mbdyn.h>
#include <filedrv.h>
#include <fixedstep.h>
#include <sockdrv.h>
#include <streamdrive.h>
#ifdef USE_SOCKET
#include <socketstreamdrive.h>
#endif // USE_SOCKET
#include <rtai_in_drive.h>

/* FileDrive - begin */

	

FileDrive::FileDrive(unsigned int uL, const DriveHandler* pDH,
		const char* const s, integer nd)
: Drive(uL, pDH), sFileName(NULL), iNumDrives(nd), pdVal(NULL)
{
	ASSERT(s != NULL);
	SAFESTRDUP(sFileName, s);
   	SAFENEWARR(pdVal, doublereal, nd + 1);
   	for (int iCnt = 0; iCnt <= nd; iCnt++) {
      		pdVal[iCnt] = 0.;
   	}
}


FileDrive::~FileDrive(void)
{
	if (sFileName != NULL) {
		SAFEDELETEARR(sFileName);
	}

   	if (pdVal != NULL) {
      		SAFEDELETEARR(pdVal);
   	}
}


Drive::Type
FileDrive::GetDriveType(void) const
{
	return Drive::FILEDRIVE;
}

doublereal
FileDrive::dGet(const doublereal& /* t */ , int i) const
{
   	ASSERT(i > 0 && i <= iNumDrives);
   	return pdVal[i];
}

/* FileDrive - end */


/* FileDriveCaller - begin */

FileDriveCaller::FileDriveCaller(const DriveHandler* pDH,
		const FileDrive* p, integer i,
		const doublereal& da)
: DriveCaller(pDH), pFileDrive((FileDrive*)p), iNumDrive(i), dAmplitude(da)
{
	ASSERT(pFileDrive != NULL);
	ASSERT(iNumDrive > 0 && iNumDrive <= pFileDrive->iGetNumDrives());
}


FileDriveCaller::~FileDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller* FileDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC,
			FileDriveCaller,
			FileDriveCaller(pDrvHdl,  pFileDrive,
				iNumDrive, dAmplitude));

	return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
FileDriveCaller::Restart(std::ostream& out) const
{
	out << " file, " << pFileDrive->GetLabel()
		<< ", " << iNumDrive;
	if (dAmplitude != 1.) {
		out << ", amplitude, " << dAmplitude;
	}
	return out;
}

/* FileDriveCaller - end */


/* legge i drivers tipo file */

Drive*
ReadFileDriver(DataManager* pDM,
		MBDynParser& HP,
		unsigned int uLabel)
{
	Drive* pDr = NULL;

	const char* sKeyWords[] = {
		"fixed" "step",
		"socket",
		"stream",
		"socket" "stream",
		"rtai" "input",
		NULL
	};

	enum KeyWords {
		UNKNOWN = -1,

		FIXEDSTEP = 0,
		SOCKET,
		STREAM,
		SOCKETSTREAM,
		RTAIINPUT,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di drive */
	KeyWords CurrKeyWord = KeyWords(HP.GetWord());

	switch (CurrKeyWord) {
	case FIXEDSTEP:
		pDr = ReadFixedStepFileDrive(pDM, HP, uLabel);
		break;

	case SOCKET:
		pDr = ReadSocketDrive(pDM, HP, uLabel);
		break;

	case SOCKETSTREAM:
		pedantic_cout("\"socket stream\" is deprecated; "
				"use \"stream\" instead at line "
				<< HP.GetLineData() << std::endl);
		goto read_stream;

	case RTAIINPUT:
#ifndef USE_RTAI
	       	silent_cout("need --with-rtai to allow RTAI mailboxes; "
				"using stream..." << std::endl);
#endif /* ! USE_RTAI */
		/* fall thru... */

	case STREAM:
read_stream:;
#ifdef USE_RTAI
		if (mbdyn_rtai_task != NULL){
			silent_cout("starting RTAI input drive" << std::endl);
			pDr = ReadRTAIInDrive(pDM, HP, uLabel);
		} else 
#endif /* USE_RTAI */		
		{
#ifdef USE_SOCKET
			silent_cout("starting stream drive" << std::endl);
			pDr = ReadSocketStreamDrive(pDM, HP, uLabel);
#else // ! USE_SOCKET
			silent_cerr("stream drive " << uLabel
				<< " not allowed at line " << HP.GetLineData()
				<< " because apparently the current architecture "
				"does not support sockets" << std::endl);
			throw ErrGeneric();
#endif // ! USE_SOCKET
		}
		break;

	default:
		silent_cerr("unknown file drive at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	return pDr;
} /* End of ReadFileDriver */

