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

/*
 * Michele Attolico <attolico@aero.polimi.it>
*/


#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <netdb.h>

#include <dataman.h>
#include <filedrv.h>
#include <streamdrive.h>
#include "sock.h"

#include <socketstreamdrive.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>

#define DEFAULT_PORT	5500 /* FIXME: da definire meglio */
#define DEFAULT_HOST 	"127.0.0.1"

SocketStreamDrive::SocketStreamDrive(unsigned int uL,
		DataManager* pDM,
		const char* const sFileName,
		integer nd, unsigned int ie, bool c,
		unsigned short int p,
		const char* const h, int flags)
: StreamDrive(uL, pDM->pGetDrvHdl(), sFileName, nd, c),
InputEvery(ie), InputCounter(0), pUS(0), recv_flags(flags)
{
	ASSERT(InputEvery > 0);

	SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(h, p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamDrive::SocketStreamDrive(unsigned int uL,
		DataManager* pDM,
		const char* const sFileName,
		integer nd, unsigned int ie, bool c,
		const char* const p, int flags)
: StreamDrive(uL, pDM->pGetDrvHdl(), sFileName, nd, c),
InputEvery(ie), InputCounter(0), pUS(0), recv_flags(flags)
{
	ASSERT(InputEvery > 0);

	SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(p, c));
	if (c) {
		pDM->RegisterSocketUser(pUS);
	} else {
		pUS->Connect();
	}
}

SocketStreamDrive::~SocketStreamDrive(void)
{
	if (pUS != 0) {
		SAFEDELETE(pUS);
	}
}

FileDrive::Type
SocketStreamDrive::GetFileDriveType(void) const
{
	return FileDrive::SOCKETSTREAM;
}

/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
SocketStreamDrive::Restart(std::ostream& out) const
{
	out << "  file: " << uLabel << ", socket stream," 
		" stream drive name, \"" << sFileName << "\"";
	pUS->Restart(out);
	return out << ", " << iNumDrives << ";" << std::endl;
}
   
void
SocketStreamDrive::ServePending(const doublereal& t)
{
	/* by now, an abandoned drive is not read any more;
	 * should we retry or what? */
	if (pUS->Abandoned()) {
		return;
	}

	ASSERT(pUS->Connected());
	
	/* read only every InputEvery steps */
	InputCounter++;
	if (InputCounter != InputEvery) {
		return;
	}
	InputCounter = 0;

	/* FIXME: no receive at first step? */
	switch (recv(pUS->GetSock(), buf, size, recv_flags)) {
	case 0:
do_abandon:;
		silent_cout("SocketStreamDrive(" << sFileName << "): "
				<< "communication closed by host" << std::endl);
		pUS->Abandon();
		break;

	case -1: {
		int save_errno = errno;

		/* some errno values may be legal */
		switch (save_errno) {
		case ECONNRESET:
			goto do_abandon;
		}

		char *err_msg = strerror(save_errno);

		silent_cout("SocketStreamDrive(" << sFileName << ") "
				<< "failed (" << save_errno << ": " 
				<< err_msg << ")" << std::endl);
		throw ErrGeneric();
	}

	default: {	
		doublereal *rbuf = (doublereal *)buf;

		for (int i = 1; i <= iNumDrives; i++) {
			pdVal[i] = rbuf[i-1];
		}
		break;
	}
	}
}


/* legge i drivers tipo stream */

Drive*
ReadSocketStreamDrive(DataManager* pDM,
		MBDynParser& HP,
		unsigned int uLabel)
{
	
	bool create = false;
	unsigned short int port = DEFAULT_PORT;
	const char *name = 0;
	const char *host = 0;
	const char *path = 0;


	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "drive" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == 0) {
			silent_cerr("unable to read stream drive name "
				"for SocketStreamDrive(" << uLabel 
				<< ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();

		} 

		SAFESTRDUP(name, m);

	} else {
		silent_cerr("missing stream drive name "
			"for SocketStreamDrive(" << uLabel
			<< ") at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			create = true;
		} else if (HP.IsKeyWord("no")) {
			create = false;
		} else {
			silent_cerr("\"create\" must be \"yes\" or \"no\" "
				"for stream drive \"" << name << "\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
	}
		
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == 0) {
			silent_cerr("unable to read local path for "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
		
		SAFESTRDUP(path, m);
	}
	
	if (HP.IsKeyWord("port")) {
		if (path != 0) {
			silent_cerr("cannot specify a port "
					"for a local socket in "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();		
		}
		int p = HP.GetInt();
		/* Da sistemare da qui */
#ifdef IPPORT_USERRESERVED
		if (p <= IPPORT_USERRESERVED) {
			silent_cerr("SocketStreamDrive(" << uLabel << "): "
					"cannot listen on port " << port
					<< ": less than IPPORT_USERRESERVED=" 
					<< IPPORT_USERRESERVED
					<< " at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric();
		}
		/* if #undef'd, don't bother checking;
		 * the OS will do it for us */
#endif /* IPPORT_USERRESERVED */
		port = p;
	}

	
	if (HP.IsKeyWord("host")) {
		if (path != 0) {
			silent_cerr("cannot specify an allowed host "
					"for a local socket in "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("unable to read host for "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}

		SAFESTRDUP(host, h);

	} else if (!path && !create) {
		silent_cerr("host undefined for "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
			<< HP.GetLineData() << std::endl);
		silent_cerr("using default host: "
			<< DEFAULT_HOST << std::endl);
		SAFESTRDUP(host, DEFAULT_HOST);
	}

	/* we want to block until the whole chunk is received */
	int flags = 0;
#ifdef MSG_WAITAL
	flags |= MSG_WAITALL;
#endif /* MSG_WAITALL */
	while (HP.IsArg()) {
		if (HP.IsKeyWord("no" "signal")) {
			flags |= MSG_NOSIGNAL;

#if 0
		/* not honored by recv() */
		} else if (HP.IsKeyWord("non" "blocking")) {
			flags |= MSG_DONTWAIT;
#endif
		} else {
			break;
		}
	}

	unsigned int InputEvery = 1;
	if (HP.IsKeyWord("input" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("invalid input every value " << i
					<< " at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric();
		}
		InputEvery = (unsigned int)i;
	}

	int idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("illegal number of channels for "
				<< "SocketStreamDrive("
			 	<< uLabel << ") at line "
			<< std::endl);
		throw ErrGeneric();
	}

	Drive* pDr = 0;
	if (path == 0) {
		SAFENEWWITHCONSTRUCTOR(pDr, SocketStreamDrive,
				SocketStreamDrive(uLabel, 
				pDM,
				name, idrives, InputEvery, create, port, host,
				flags));

	} else {
		SAFENEWWITHCONSTRUCTOR(pDr, SocketStreamDrive,
				SocketStreamDrive(uLabel, 
				pDM,
				name, idrives, InputEvery, create, path,
				flags));
	}

	return pDr;
} /* End of ReadStreamDrive */

