/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef USE_SOCKET

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

#include "dataman.h"
#include "socketstream_out_elem.h"
#include "sock.h"

/* SocketStreamElem - begin */

SocketStreamElem::SocketStreamElem(unsigned int uL,
	const std::string& name,
	unsigned int oe,
	UseSocket *pUS,
	StreamContent *pSC,
	int flags, bool bSendFirst, bool bAbortIfBroken)
: Elem(uL, flag(0)),
StreamOutElem(uL, name, oe),
pUS(pUS), pSC(pSC), send_flags(flags),
bSendFirst(bSendFirst), bAbortIfBroken(bAbortIfBroken)
{
	NO_OP;
}

SocketStreamElem::~SocketStreamElem(void)
{
	if (pUS != 0) {
		SAFEDELETE(pUS);
	}

	if (pSC != 0) {
		SAFEDELETE(pSC);
	}
}

std::ostream&
SocketStreamElem::Restart(std::ostream& out) const
{   	
	return out << "# SocketStreamElem(" << GetLabel() << "): "
		"not implemented yet" << std::endl;
}	

void
SocketStreamElem::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (bSendFirst) {
		AfterConvergence(X, XP);
	}
}

void
SocketStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	/* by now, an abandoned element does not write any more;
	 * should we retry or what? */
	if (pUS->Abandoned()) {
		return;
	}

	ASSERT(pUS->Connected());

	/* output only every OutputEvery steps */
	OutputCounter++;
	if (OutputCounter != OutputEvery) {
		return;
	}
	OutputCounter = 0;

	// prepare the output buffer
	pSC->Prepare();

	if (send(pUS->GetSock(), pSC->GetBuf(), pSC->GetSize(), send_flags) != pSC->GetSize()) {
		int save_errno = errno;
		char *msg = strerror(save_errno);
		
		silent_cerr("SocketStreamElem(" << name << "): send() failed "
				"(" << save_errno << ": " << msg << ")"
				<< std::endl);

		if (bAbortIfBroken) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}

		pUS->Abandon();
	}
}

void
SocketStreamElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

Elem *
ReadSocketStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type)
{
	bool bCreate = false;
	unsigned short int port = -1;
	std::string name;
	std::string host;
	std::string path;

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == 0) {
			silent_cerr("unable to read stream name "
				"for SocketStreamElem(" << uLabel << ") "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} 
		
		name = m;

	} else {
		silent_cerr("missing stream name "
			"for SocketStreamElem(" << uLabel
			<< ") at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			bCreate = true;

		} else if (HP.IsKeyWord("no")) {
			bCreate = false;

		} else {
			silent_cerr("\"create\" must be either "
					"\"yes\" or \"no\" "
					"at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == 0) {
			silent_cerr("unable to read local path for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		path = m;
	}

	if (HP.IsKeyWord("port")) {
		if (!path.empty()) {
			silent_cerr("cannot specify a port "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}
		int p = HP.GetInt();
		/* Da sistemare da qui */

#ifdef IPPORT_USERRESERVED
		if (p <= IPPORT_USERRESERVED) {
			silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << "): "
				"cannot listen on port " << port
				<< ": less than IPPORT_USERRESERVED=" 
				<< IPPORT_USERRESERVED
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		/* if #undef'd, don't bother checking;
		 * the OS will do it for us */
#endif /* IPPORT_USERRESERVED */
		port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (!path.empty()) {
			silent_cerr("cannot specify an allowed host "
					"for a local socket in "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("unable to read host for "
				<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
				<< "(" << uLabel << ") at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		host = h;

	} else if (path.empty() && !bCreate) {
		/* INET sockets (!path) must be created if host is missing */
		silent_cerr("host undefined for "
			<< psElemNames[Elem::SOCKETSTREAM_OUTPUT]
			<< "(" << uLabel << ") at line "
			<< HP.GetLineData() << std::endl);
		host = DEFAULT_HOST;
		silent_cerr("using default host: "
			<< host << ":"
			<< (port == (unsigned short int)(-1) ? DEFAULT_PORT : port)
			<< std::endl);
	}

	int flags = 0;
	bool bSendFirst = true;
	bool bAbortIfBroken = false;
	while (HP.IsArg()) {
		if (HP.IsKeyWord("no" "signal")) {
			flags |= MSG_NOSIGNAL;

		} else if (HP.IsKeyWord("non" "blocking")) {
#ifdef MSG_DONTWAIT
			flags |= MSG_DONTWAIT;
#else /* !MSG_DONTWAIT */
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"MSG_DONTWAIT undefined; "
				"your mileage may vary" << std::endl);
#endif /* !MSG_DONTWAIT */

		} else if (HP.IsKeyWord("no" "send" "first")) {
			bSendFirst = false;

		} else if (HP.IsKeyWord("abort" "if" "broken")) {
			bAbortIfBroken = true;

		} else {
			break;
		}
	}

	unsigned int OutputEvery = 1;
	if (HP.IsKeyWord("output" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("invalid output every value " << i
					<< " at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		OutputEvery = (unsigned int)i;
	}

	StreamContent *pSC = ReadStreamContent(pDM, HP, type);

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* costruzione del nodo */
	UseSocket *pUS = 0;
	if (path.empty()) {
		if (port == (unsigned short int)(-1)) {
			port = DEFAULT_PORT;
			silent_cerr("port undefined; using default port "
				<< port << std::endl);
		}
      
		SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(host.c_str(), port, bCreate));

	} else {
		SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(path.c_str(), bCreate));
	}

	if (bCreate) {
		pDM->RegisterSocketUser(pUS);

	} else {
		pUS->Connect();
	}

	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
		SocketStreamElem(uLabel, name, OutputEvery,
			pUS, pSC, flags, bSendFirst, bAbortIfBroken));

	return pEl;
}

#endif // USE_SOCKET
