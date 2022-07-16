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

/*
 * Michele Attolico <attolico@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_SOCKET
#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
  #include <ws2tcpip.h>
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>
#endif /* _WIN32 */
#endif // USE_SOCKET

#include "dataman.h"
#include "socketstream_out_elem.h"
#include "socketstreammotionelem.h"
#include "sock.h"

#ifdef USE_RTAI
#include "rtai_out_elem.h"
#endif // USE_RTAI

#include "geomdata.h"

#ifdef USE_SOCKET

/* SocketStreamElem - begin */

SocketStreamElem::SocketStreamElem(unsigned int uL,
	const std::string& name,
	unsigned int oe,
	UseSocket *pUS,
	StreamContent *pSC,
	int flags, bool bSendFirst, bool bAbortIfBroken,
	StreamOutEcho *pSOE,
	bool bMsgDontWait)
: Elem(uL, flag(0)),
StreamOutElem(uL, name, oe),
pUS(pUS), pSC(pSC), send_flags(flags),
bSendFirst(bSendFirst), bAbortIfBroken(bAbortIfBroken),
bMsgDontWait(bMsgDontWait),
pSOE(pSOE)
{
	if (pSOE) {
		pSOE->Init("SocketStreamElem", uLabel, pSC->GetNumChannels());
	}
}

SocketStreamElem::~SocketStreamElem(void)
{
	if (pUS != 0) {
		SAFEDELETE(pUS);
	}

	if (pSC != 0) {
		SAFEDELETE(pSC);
	}

	if (pSOE != 0) {
		delete pSOE;
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
		// output imposed values (before "derivatives")
		OutputCounter = OutputEvery - 1;

		AfterConvergence(X, XP);
	}

	// do not send "derivatives"
	OutputCounter = -1;
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

	// check whether echo is needed
	if (pSOE) {
		pSOE->Echo((doublereal *)pSC->GetBuf(), pSC->GetNumChannels());
	}

	// int rc = sendn(pUS->GetSock(), pSC->GetOutBuf(), pSC->GetOutSize(), send_flags);
	ssize_t rc = pUS->send(pSC->GetOutBuf(), pSC->GetOutSize(), send_flags);
	if (rc == -1 || rc != pSC->GetOutSize()) {
		int save_errno = WSAGetLastError();
#ifdef _WIN32
            int test_errno = WSAEWOULDBLOCK;
#else
            int test_errno = EAGAIN;
#endif /* _WIN32 */
		if (save_errno == test_errno && bMsgDontWait) {
			// would block; continue (and discard...)
			return;
		}
		
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

#endif // USE_SOCKET

/*----------------------------------------------------------------------------
management of 'content type' for stream output element ('values','motion', etc)
------------------------------------------------------------------------------

Rearranged by Luca Conti (May 2017) on the basis of previous existing code
(fully working, just rearranged).

Edited in order to apply the same mechanism with 'readers' and 'maps' (std::map)
  already in use for constitutive laws and drives
*/

/*temporary data that allows to store all necessary information for creation
 of socket stream output element (used by functions in struct
SocketStreamElemRead)*/
struct SocketStreamOutputDataTmp {
	bool bIsRTAI;
	StreamOutEcho *pSOE;
	StreamContent *pSC;
	#if defined(HAVE_GETADDRINFO)
			struct addrinfo hints,*res;
			int rc;
	#elif defined(HAVE_GETHOSTBYNAME)
			struct hostent *he;
	#elif defined(HAVE_INET_ATON)
			struct in_addr addr;
	#endif // ! HAVE_GETADDRINFO && ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON
	std::string name;
	std::string path;
	std::string host;
	bool bCreate;
	bool bNonBlocking;
	bool bNoSignal;
	bool bSendFirst;
	bool bAbortIfBroken;
	unsigned short int port;
	int socket_type;
	int flags;
	unsigned int OutputEvery;
};

#ifdef USE_SOCKET
static const int sock_stream = SOCK_STREAM;
static const int sock_dgram = SOCK_DGRAM;
#else // ! USE_SOCKET
static const int sock_stream = 1;
static const int sock_dgram = 2;
#endif // ! USE_SOCKET


/*these functions are called by ReadSocketStreamElem in order to get and elaborate all data necessary to create the element*/
struct SocketStreamOutputElemCreator {
	virtual void getSocketStreamOutParam(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type, SocketStreamOutputDataTmp& socketStreamOutputDataTmp);
	virtual Elem* createSocketStreamOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type, SocketStreamOutputDataTmp& socketStreamOutputDataTmp);
};

Elem *
ReadSocketStreamElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type)
{
	SocketStreamOutputDataTmp socketStreamOutputDataTmp;
	SocketStreamOutputElemCreator socketStreamOutputElemCreator;
	socketStreamOutputElemCreator.getSocketStreamOutParam(pDM, HP, uLabel, type, socketStreamOutputDataTmp);
	socketStreamOutputDataTmp.pSOE = ReadStreamOutEcho(HP);
	socketStreamOutputDataTmp.pSC = ReadStreamContent(pDM, HP, type);
	Elem* pEl = socketStreamOutputElemCreator.createSocketStreamOutElem(pDM, HP, uLabel, type, socketStreamOutputDataTmp);
  return pEl;
}

void SocketStreamOutputElemCreator::getSocketStreamOutParam(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type, SocketStreamOutputDataTmp& socketStreamOutputDataTmp){
	socketStreamOutputDataTmp.bIsRTAI = false;//bool bIsRTAI(false);
#ifdef USE_RTAI
	if (::rtmbdyn_rtai_task != 0) {
		socketStreamOutputDataTmp.bIsRTAI = true;
	}
#endif // USE_RTAI
#ifndef USE_SOCKET
	if (!socketStreamOutputDataTmp.bIsRTAI) {
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"not allowed because apparently the current "
			"architecture does not support sockets "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // ! USE_SOCKET

	socketStreamOutputDataTmp.port = (unsigned short int)(-1);
	socketStreamOutputDataTmp.bCreate = false;
	bool bGotCreate(false);

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"unable to read stream name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (socketStreamOutputDataTmp.bIsRTAI && strlen(m) != 6) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"illegal stream name \"" << m << "\" "
				"(must be exactly 6 chars) "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		socketStreamOutputDataTmp.name = m;

	} else {
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"missing stream name "
			"at line " << HP.GetLineData() << std::endl);
		if (socketStreamOutputDataTmp.bIsRTAI) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("create")) {
		bGotCreate = true;
		if (!HP.GetYesNo(socketStreamOutputDataTmp.bCreate)) {
			silent_cerr("SocketStreamElem(" << uLabel << "):"
				"\"create\" must be either "
				"\"yes\" or \"no\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
#ifdef _WIN32
        silent_cerr("SocketStreamElem(" << uLabel << "): "
            "local sockets are not supported on Windows, you must use inet "
            "at line " << HP.GetLineData()
            << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else
		const char *m = HP.GetFileName();

		if (m == 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"unable to read local path "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		socketStreamOutputDataTmp.path = m;
#endif /* _WIN32 */
	}

	if (HP.IsKeyWord("port")) {
		if (!socketStreamOutputDataTmp.path.empty()) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"cannot specify a port for a local socket "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		int p = HP.GetInt();

		if (p <= 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"illegal port " << p << " "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		socketStreamOutputDataTmp.port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (!socketStreamOutputDataTmp.path.empty()) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"cannot specify an allowed host "
				"for a local socket "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const char *h;

		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		socketStreamOutputDataTmp.host = h;

	} else if (socketStreamOutputDataTmp.path.empty() && !socketStreamOutputDataTmp.bCreate) {
		/* INET sockets (!path) must be created if host is missing */
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"host undefined "
			"at line " << HP.GetLineData() << std::endl);
		socketStreamOutputDataTmp.host = DEFAULT_HOST;
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"using default host: "
			<< socketStreamOutputDataTmp.host << ":"
			<< (socketStreamOutputDataTmp.port == (unsigned short int)(-1) ? DEFAULT_PORT : socketStreamOutputDataTmp.port)
			<< std::endl);
	}

	socketStreamOutputDataTmp.socket_type = sock_stream;
	if (HP.IsKeyWord("socket" "type")) {
		if (HP.IsKeyWord("udp")) {
			socketStreamOutputDataTmp.socket_type = sock_dgram;

			if (!bGotCreate) {
				socketStreamOutputDataTmp.bCreate = true;
			}

		} else if (!HP.IsKeyWord("tcp")) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"invalid socket type "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (socketStreamOutputDataTmp.socket_type == sock_dgram && socketStreamOutputDataTmp.bCreate) {
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"socket type=udp incompatible with create=yes "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	socketStreamOutputDataTmp.bNonBlocking = false;
	socketStreamOutputDataTmp.bNoSignal = false;
	socketStreamOutputDataTmp.bSendFirst = true;
	socketStreamOutputDataTmp.bAbortIfBroken = false;
	while (HP.IsArg()) {
		if (HP.IsKeyWord("no" "signal")) {
			socketStreamOutputDataTmp.bNoSignal = true;

		} else if (HP.IsKeyWord("signal")) {
			socketStreamOutputDataTmp.bNoSignal = false;

		} else if (HP.IsKeyWord("blocking")) {
			socketStreamOutputDataTmp.bNonBlocking = false;

		} else if (HP.IsKeyWord("non" "blocking")) {
			socketStreamOutputDataTmp.bNonBlocking = true;

		} else if (HP.IsKeyWord("no" "send" "first")) {
			socketStreamOutputDataTmp.bSendFirst = false;

		} else if (HP.IsKeyWord("send" "first")) {
			socketStreamOutputDataTmp.bSendFirst = true;

		} else if (HP.IsKeyWord("abort" "if" "broken")) {
			socketStreamOutputDataTmp.bAbortIfBroken = true;

		} else if (HP.IsKeyWord("do" "not" "abort" "if" "broken")) {
			socketStreamOutputDataTmp.bAbortIfBroken = false;

		} else {
			break;
		}
	}

	socketStreamOutputDataTmp.OutputEvery = 1;
	if (HP.IsKeyWord("output" "every")) {
		int i = HP.GetInt();
		if (i <= 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"invalid output every value " << i << " "
		 		"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		socketStreamOutputDataTmp.OutputEvery = (unsigned int)i;
	}

}

Elem* SocketStreamOutputElemCreator::createSocketStreamOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type, SocketStreamOutputDataTmp& socketStreamOutputDataTmp) {
	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"semicolon expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Elem *pEl = 0;
	bool bMsgDontWait = false;
	
	// .log file output
	std::ostream& out = pDM->GetLogFile();
	out 
		<< "outputelement: " << uLabel 
		<< " stream";

	if (socketStreamOutputDataTmp.bIsRTAI) {
#ifdef USE_RTAI
		if (socketStreamOutputDataTmp.pSOE != 0) {
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"echo ignored in RTAI mode"
				<< std::endl);
		}

		unsigned long node = (unsigned long)-1;
#if defined(HAVE_GETADDRINFO)
		socketStreamOutputDataTmp.hints = {0};
		socketStreamOutputDataTmp.res = null;

		socketStreamOutputDataTmp.hints.ai_family = AF_INET;
		socketStreamOutputDataTmp.hints.ai_socktype = SOCK_STREAM; // FIXME: what about SOCK_DGRAM?
		socketStreamOutputDataTmp.rc = getaddrinfo(socketStreamOutputDataTmp.host.c_str(), NULL, &(socketStreamOutputDataTmp.hints), &(socketStreamOutputDataTmp.res));
		if (socketStreamOutputDataTmp.rc == 0) {
			node = ((struct sockaddr_in *)(socketStreamOutputDataTmp.res)->ai_addr)->sin_addr.s_addr;
			freeaddrinfo(socketStreamOutputDataTmp.res);
		}
#elif defined(HAVE_GETHOSTBYNAME)
		socketStreamOutputDataTmp.he = gethostbyname(socketStreamOutputDataTmp.host.c_str());
		if (socketStreamOutputDataTmp.he != NULL) {
			node = ((unsigned long *)(socketStreamOutputDataTmp.he)->h_addr_list[0])[0];
		}
#elif defined(HAVE_INET_ATON)
		if (inet_aton(socketStreamOutputDataTmp.host.c_str(), &(socketStreamOutputDataTmp.addr))) {
			node = socketStreamOutputDataTmp.addr.s_addr;
		}
#else // ! HAVE_GETADDRINFO && ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"host (RTAI RPC) not supported "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! HAVE_GETADDRINFO && ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON

		if (node == (unsigned long)-1) {
			silent_cerr("RTMBDynInDrive(" << uLabel << "): "
				"unable to convert host \"" << host << "\" to node" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		silent_cerr("starting RTMBDynOutputElement(" << uLabel << ")..."
			<< std::endl);
		SAFENEWWITHCONSTRUCTOR(pEl, RTMBDynOutElem,
			RTMBDynOutElem(uLabel,
				socketStreamOutputDataTmp.name, socketStreamOutputDataTmp.host, node,
				socketStreamOutputDataTmp.bCreate, socketStreamOutputDataTmp.pSC, socketStreamOutputDataTmp.bNonBlocking));

		out
			<< " " << "RTAI"
			<< " " << name
			<< " " << host
			<< " " << node
			<< " " << bCreate;

		WriteStreamContentLogOutput(pSC, out);

#endif // USE_RTAI

	} else {
		socketStreamOutputDataTmp.flags = 0;

#ifdef USE_SOCKET
		/* node creation */
		UseSocket *pUS = 0;
		if (socketStreamOutputDataTmp.path.empty()) {
			if (socketStreamOutputDataTmp.port == (unsigned short int)(-1)) {
				socketStreamOutputDataTmp.port = DEFAULT_PORT;
				silent_cerr("SocketStreamElem(" << uLabel << "): "
					"port undefined; using default port "
					<< socketStreamOutputDataTmp.port << " at line "
					<< HP.GetLineData() << std::endl);
			}

			SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(socketStreamOutputDataTmp.host.c_str(),
				socketStreamOutputDataTmp.port, socketStreamOutputDataTmp.socket_type, socketStreamOutputDataTmp.bCreate));
		
			// .log file output
			out 
				<< " " << "INET"
				<< " " << socketStreamOutputDataTmp.name
				<< " " << socketStreamOutputDataTmp.host
				<< " " << socketStreamOutputDataTmp.port;
			if (socketStreamOutputDataTmp.socket_type == sock_dgram) {
				out << " udp";
			} else {
				out << " tcp";
			}
			out << " " << socketStreamOutputDataTmp.bCreate;


		}
		else {
#ifdef _WIN32
			silent_cerr("SocketStreamElem(" << uLabel << "): "
				"Unix style local sockets are not supported on windows"
				<< std::endl);
#else
			SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(socketStreamOutputDataTmp.path.c_str(),
				socketStreamOutputDataTmp.socket_type, socketStreamOutputDataTmp.bCreate));
			
			// .log file output
			out 
				<< " " << "UNIX"
				<< " " << socketStreamOutputDataTmp.name
				<< " " << socketStreamOutputDataTmp.path;
			if (socketStreamOutputDataTmp.socket_type == sock_dgram) {
				out << " udp";
			} else {
				out << " tcp";
			}
			out << " " << socketStreamOutputDataTmp.bCreate;
#endif // _WIN32	
		}

		if (socketStreamOutputDataTmp.bCreate) {
			pDM->RegisterSocketUser(pUS);

		} else {
			pUS->Connect();
		}

#ifdef MSG_NOSIGNAL
		if (socketStreamOutputDataTmp.bNoSignal) {
			// NOTE: we assume MSG_NOSIGNAL is a macro...
			socketStreamOutputDataTmp.flags |= MSG_NOSIGNAL;

		} else {
			socketStreamOutputDataTmp.flags &= ~MSG_NOSIGNAL;
		}
#else // !MSG_NOSIGNAL
		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"MSG_NOSIGNAL undefined; "
			"your mileage may vary" << std::endl);
#endif // !MSG_NOSIGNAL

#ifdef MSG_DONTWAIT
		// NOTE: we assume MSG_DONTWAIT is a macro...
		if (socketStreamOutputDataTmp.bNonBlocking) {
			socketStreamOutputDataTmp.flags |= MSG_DONTWAIT;

		} else {
			socketStreamOutputDataTmp.flags &= ~MSG_DONTWAIT;
		}
#else // !MSG_DONTWAIT

                if (socketStreamOutputDataTmp.bNonBlocking) {
                    bMsgDontWait = true;
                }

		silent_cerr("SocketStreamElem(" << uLabel << "): "
			"MSG_DONTWAIT undefined; "
			"your mileage may vary" << std::endl);
#endif // !MSG_DONTWAIT

		silent_cerr("starting SocketStreamElem(" << uLabel << ")..."
			<< std::endl);
		SAFENEWWITHCONSTRUCTOR(pEl, SocketStreamElem,
			SocketStreamElem(uLabel, socketStreamOutputDataTmp.name, socketStreamOutputDataTmp.OutputEvery,
				pUS, socketStreamOutputDataTmp.pSC, socketStreamOutputDataTmp.flags,
				socketStreamOutputDataTmp.bSendFirst, socketStreamOutputDataTmp.bAbortIfBroken,
				socketStreamOutputDataTmp.pSOE, bMsgDontWait));

		out 
			<< " " << (!socketStreamOutputDataTmp.bNoSignal)
			<< " " << (!socketStreamOutputDataTmp.bNonBlocking)
			<< " " << socketStreamOutputDataTmp.bSendFirst
			<< " " << socketStreamOutputDataTmp.bAbortIfBroken
			<< " " << socketStreamOutputDataTmp.OutputEvery;
		
		WriteStreamContentLogOutput(socketStreamOutputDataTmp.pSC, out);

#endif // USE_SOCKET
	}

	return pEl;
}

void
WriteStreamContentLogOutput(const StreamContent* pSC, 
	std::ostream& out)
{
	const StreamContentMotion* pSCM = NULL;
	pSCM = dynamic_cast<const StreamContentMotion*>(pSC);

	if (pSCM == NULL) 
	{
		// pSC is of type StreamContent::VALUES
		out 
			<< " values"
			<< " " << pSC->GetNumChannels()
			<< std::endl;
	} else {
		// pSC is of type StreamContent::MOTION
		const unsigned uFlags = pSCM->uGetFlags();
		out 
			<< " motion"
			<< " " << bool(uFlags & GeometryData::X)
			<< " " << bool(uFlags & GeometryData::R)
			<< " " << bool(uFlags & GeometryData::RT)
			<< " " << bool(uFlags & GeometryData::V)
			<< " " << bool(uFlags & GeometryData::W);

		for (std::vector<const StructNode*>::const_iterator i = pSCM->nodes_begin(); 
				i != pSCM->nodes_end(); ++i) 
		{
			out << " " << (*i)->GetLabel();
		}

		out << std::endl;
	}
}
