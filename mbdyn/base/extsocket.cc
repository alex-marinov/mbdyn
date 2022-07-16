/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/extforce.cc,v 1.50 2017/01/12 14:46:09 masarati Exp $ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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

#include "extforce.h"


#ifdef USE_SOCKET

#include "dataman.h"
#include "extedge.h"
#include "except.h"
#include "solver.h"

#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <sys/stat.h>

#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
  #include <ws2tcpip.h>
#else
  #include <sys/types.h>
  #include <unistd.h>
#endif /* _WIN32 */

#include "extsocket.h"


ExtSocketHandler::ExtSocketHandler(UseSocket *pUS, mbsleep_t SleepTime,
	int recv_flags, int send_flags)
: ExtRemoteHandler(SleepTime, true, false),
pUS(pUS), recv_flags(recv_flags), send_flags(send_flags)
{
    pedantic_cout("ExtSocketHandler in constructor" << std::endl);
	NO_OP;
}

ExtSocketHandler::~ExtSocketHandler(void)
{
    pedantic_cout("ExtSocketHandler: in destructor" << std::endl);
	if (pUS->GetSock() != INVALID_SOCKET) {
		uint8_t u = ES_ABORT;

		// ignore result
		(void)sendn(pUS->GetSock(), (const char *)&u, sizeof(u), send_flags);
	}
	SAFEDELETE(pUS);
}

ExtFileHandlerBase::Negotiate
ExtSocketHandler::NegotiateRequest(void) const
{
    pedantic_cout("ExtSocketHandler: in NegotiateRequest" << std::endl);

	if (pUS->Create()) {
		return ExtFileHandlerBase::NEGOTIATE_SERVER;

	} else {
		return ExtFileHandlerBase::NEGOTIATE_CLIENT;
	}
}

bool
ExtSocketHandler::Prepare_pre(void)
{
	uint8_t u;
	ssize_t rc;

	pedantic_cout("ExtSocketHandler: in Prepare_pre" << std::endl);

	switch (NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_CLIENT:
	    /* This the server, send negotiation request */
		u = ES_NEGOTIATION;

		pedantic_cout("ExtSocketHandler::Prepare_pre -- NEGOTIATE_CLIENT -- sending negotiation request" << std::endl);

		rc = sendn(pUS->GetSock(), (const char *)&u, sizeof(u), send_flags);

		pedantic_cout("ExtSocketHandler::Prepare_pre -- NEGOTIATE_CLIENT -- negotiation request sent" << std::endl);

		if (rc == SOCKET_ERROR) {
			int save_errno = errno;
			silent_cerr("ExtSocketHandler::Prepare_pre -- NEGOTIATE_CLIENT -- negotiation request send() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			return (bOK = false);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler::Prepare_pre -- NEGOTIATE_CLIENT -- negotiation request send() failed "
				"(sent " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			return (bOK = false);
		}
		break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER:
        /* This the client, check if negotiation request has been sent */
	    pedantic_cout("ExtSocketHandler::Prepare_pre -- NEGOTIATE_SERVER -- trying to receive negotiation request" << std::endl);

		rc = recvn(pUS->GetSock(), (char *)&u, sizeof(u), recv_flags);

		if (rc == SOCKET_ERROR) {
			return (bOK = false);

		} else if (rc != sizeof(u)) {
			return (bOK = false);
		}

		if (u != ES_NEGOTIATION) {
			silent_cerr("ExtSocketHandler::Prepare_pre -- NEGOTIATE_SERVER -- negotiation request recv() failed"
				<< std::endl);
			return (bOK = false);
		}
		break;

	default:
		ASSERT(0);
	}

	pedantic_cout("ExtSocketHandler::Prepare_pre -- NEGOTIATE_SERVER -- negotiation started" << std::endl);

	return true;
}

void
ExtSocketHandler::Prepare_post(bool ok)
{
    pedantic_cout("ExtSocketHandler: in Prepare_post, ok is " << ok << std::endl);
	if (NegotiateRequest()) {
		uint8_t u = ok ? ES_OK : ES_ABORT;
        pedantic_cout("ExtSocketHandler: sending negotiation response" << std::endl);
		ssize_t rc = sendn(pUS->GetSock(), (const char *)&u, sizeof(u),
			send_flags);
		if (rc == SOCKET_ERROR) {
			int save_errno = WSAGetLastError();
			silent_cerr("ExtSocketHandler: negotiation response send() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: negotiation response send() failed "
				"(sent " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pedantic_cout("ExtSocketHandler: negotiation response sent" << std::endl);

	} else {
		uint8_t u;
		pedantic_cout("ExtSocketHandler: receiving negotiation response" << std::endl);
		ssize_t rc = recvn(pUS->GetSock(), (char *)&u, sizeof(u),
			recv_flags);
		if (rc == SOCKET_ERROR) {
			int save_errno = WSAGetLastError();
			silent_cerr("ExtSocketHandler: negotiation response recv() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: negotiation response recv() failed "
				"(received " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		}
        pedantic_cout("ExtSocketHandler: negotiation response received" << std::endl);
		if (u != ES_OK) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

void
ExtSocketHandler::AfterPredict(void)
{
    pedantic_cout("ExtSocketHandler: in AfterPredict" << std::endl);
	bLastReadForce = false;
	bReadForces = true;
}

bool
ExtSocketHandler::Send_pre(SendWhen when)
{
    pedantic_cout("ExtSocketHandler: in Send_pre" << std::endl);
	if (!bReadForces || !bOK) {
		return false;
	}

	uint8_t u;
	if (when == SEND_AFTER_CONVERGENCE) {
		u = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;
	} else {
		u = ES_REGULAR_DATA;
	}
	pedantic_cout("ExtSocketHandler: sending when to send data" << std::endl);
	ssize_t rc = sendn(pUS->GetSock(), (const char *)&u, sizeof(u), send_flags);
	if (rc == SOCKET_ERROR) {
		int save_errno = WSAGetLastError();
		silent_cerr("ExtSocketHandler: send() failed "
			"(" << save_errno << ": " << strerror(save_errno) << ")"
			<< std::endl);
		mbdyn_set_stop_at_end_of_iteration();
		return (bOK = false);

	} else if (rc != sizeof(u)) {
		silent_cerr("ExtSocketHandler: send() failed "
			"(sent " << rc << " bytes "
			"instead of " << sizeof(u) << ")"
			<< std::endl);
		mbdyn_set_stop_at_end_of_iteration();
		return (bOK = false);
	}

	pedantic_cout("ExtSocketHandler: sent when to send data" << std::endl);

	return true;
}

void
ExtSocketHandler::Send_post(SendWhen when)
{
    pedantic_cout("ExtSocketHandler: in Send_post" << std::endl);
#if 0
	if (when == SEND_AFTER_CONVERGENCE) {
		bReadForces = true;
	}
#endif
	NO_OP;
}

bool
ExtSocketHandler::Recv_pre(void)
{
    pedantic_cout("ExtSocketHandler: in Recv_pre" << std::endl);
	if (!bReadForces) {
		return false;
	}

	uint8_t u = 0;

#ifdef _WIN32
    int test_errno = WSAEWOULDBLOCK;
#else
    int test_errno = EAGAIN;
#endif /* _WIN32 */

    int save_errno = 0;

	if (SleepTime != 0) {
		for ( ; ; ) {
			ssize_t rc;
#ifdef _WIN32
            // ensure no blocking for duration of recv call
            if (pUS->IsBlocking()){
                // make non-blocking
                unsigned long mode = 1;
                int result = ioctlsocket(pUS->GetSock(), FIONBIO, &mode);
            }
            rc = recvn(pUS->GetSock(), (char *)&u, sizeof(u),
				recv_flags );
            save_errno = WSAGetLastError();
            if (pUS->IsBlocking()){
                // make blocking again
                unsigned long mode = 0;
                int result = ioctlsocket(pUS->GetSock(), FIONBIO, &mode);
            }
#else
			rc = recvn(pUS->GetSock(), (char *)&u, sizeof(u),
				recv_flags | MSG_DONTWAIT);
            save_errno = WSAGetLastError();
#endif /* _WIN32 */
			if (rc != SOCKET_ERROR) {
				break;
			}

			if (save_errno != test_errno) {
				silent_cerr("ExtSocketHandler: "
					"recv() failed (" << save_errno << ": "
					<< sock_err_string(save_errno) << ")"
					<< std::endl);
				mbdyn_set_stop_at_end_of_iteration();
				return (bOK = false);
			}

			if (mbdyn_stop_at_end_of_iteration()) {
				return (bOK = false);
			}

			mbsleep(&SleepTime);
		}

	} else {
	    // wait (block) until the status code is returned
		ssize_t rc = recvn(pUS->GetSock(), (char *)&u, sizeof(u), recv_flags);

		if (rc == SOCKET_ERROR) {
			int save_errno = WSAGetLastError();

			if (WSAGetLastError() != test_errno) {
				silent_cerr("ExtSocketHandler: "
					"recv() failed (" << save_errno << ": "
					<< sock_err_string(save_errno) << ")"
					<< std::endl);
				mbdyn_set_stop_at_end_of_iteration();
				return (bOK = false);
			}

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: "
				"recv()=" << rc << " (expected " << sizeof(u) << ")"
				<< std::endl);
			mbdyn_set_stop_at_end_of_iteration();
			return (bOK = false);
		}
	}

	return ActOnCmd(u);
}

int
ExtSocketHandler::GetOutFileDes(void)
{
    pedantic_cout("ExtSocketHandler: in GetOutFileDes" << std::endl);
	return pUS->GetSock();
}

int
ExtSocketHandler::GetSendFlags(void) const
{
	return send_flags;
}

int
ExtSocketHandler::GetInFileDes(void)
{
	return pUS->GetSock();
}

int
ExtSocketHandler::GetRecvFlags(void) const
{
	return recv_flags;
}

#endif /* USE_SOCKET */


ExtFileHandlerBase *
ReadExtSocketHandler(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel)
{
    pedantic_cout("In ReadExtSocketHandler" << std::endl);
#ifdef USE_SOCKET
	ExtFileHandlerBase *pEFH = 0;

	bool bCreate = false;
	unsigned short int port = (unsigned short int)-1;
	std::string host;
	std::string path;

	if (HP.IsKeyWord("create")) {
		if (!HP.GetYesNo(bCreate)) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"\"create\" must be either \"yes\" or \"no\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
#ifdef _WIN32
        silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"local sockets are not supported on Windows, you must use inet "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#else /* _WIN32 */
		const char *m = HP.GetFileName();

		if (m == 0) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"unable to read local path "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		path = m;
#endif /* _WIN32 */
	}

	if (HP.IsKeyWord("port")) {
		if (!path.empty()) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot specify port "
				"for a local socket "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		int p = HP.GetInt();
		/* Da sistemare da qui */
#ifdef IPPORT_USERRESERVED
		if (p <= IPPORT_USERRESERVED) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot listen on reserved port "
				<< port << ": less than "
				"IPPORT_USERRESERVED=" << IPPORT_USERRESERVED
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
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot specify host for a local socket "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const char *h;

		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		host = h;

	} else if (path.empty() && !bCreate) {
		silent_cerr("ExtSocketHandler"
			"(" << uLabel << "): "
			"host undefined "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int socket_type = SOCK_STREAM;
	if (HP.IsKeyWord("socket" "type")) {
		if (HP.IsKeyWord("udp")) {
			socket_type = SOCK_DGRAM;

		} else if (!HP.IsKeyWord("tcp")) {
			silent_cerr("ExtSocketHandler(" << uLabel << "\"): "
				"invalid socket type "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if ((socket_type == SOCK_DGRAM) && !bCreate) {
		silent_cerr("ExtSocketHandler(" << uLabel << "\"): "
			"socket type=upd incompatible with create=no "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// we want to block until the whole chunk is received
	int recv_flags = 0;
	int send_flags = 0;
#ifdef MSG_WAITALL
	recv_flags |= MSG_WAITALL;
#endif /* MSG_WAITALL */

	while (HP.IsArg()) {
		if (HP.IsKeyWord("signal")) {
#ifdef MSG_NOSIGNAL
			recv_flags &= ~MSG_NOSIGNAL;
			send_flags &= ~MSG_NOSIGNAL;
#else // ! MSG_NOSIGNAL
			silent_cout("ExtSocketHandler"
				"(" << uLabel << "): "
				"MSG_NOSIGNAL not defined (ignored) "
				"at line " << HP.GetLineData()
				<< std::endl);
#endif // ! MSG_NOSIGNAL

		// not honored by recv(2)
		} else if (HP.IsKeyWord("no" "signal")) {
#ifdef MSG_NOSIGNAL
			recv_flags |= MSG_NOSIGNAL;
			send_flags |= MSG_NOSIGNAL;
#else // ! MSG_NOSIGNAL
			silent_cout("ExtSocketHandler"
				"(" << uLabel << "): "
				"MSG_NOSIGNAL not defined (ignored) "
				"at line " << HP.GetLineData()
				<< std::endl);
#endif // ! MSG_NOSIGNAL

		} else {
			break;
		}
	}

	UseSocket *pUS = 0;
	if (path.empty()) {
		if (port == (unsigned short int)(-1)) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"port missing"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(host, port, socket_type, bCreate));

	} else {
#ifdef _WIN32
        // error should be thrown before we reach here
#else /* _WIN32 */
		SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(path, socket_type, bCreate));
#endif /* _WIN32 */
	}
    pedantic_cout("In ReadExtSocketHandler before if (bCreate)" << std::endl);
	if (bCreate) {
        pedantic_cout("In ReadExtSocketHandler bCreate true so RegisterSocketUser" << std::endl);
		pDM->RegisterSocketUser(pUS);

	} else {
	    pedantic_cout("In ReadExtSocketHandler bCreate false so Connect" << std::endl);
		pUS->Connect();
	}

	mbsleep_t SleepTime = mbsleep_init(0);
	std::streamsize Precision = 0;
	ReadExtFileParams(pDM, HP, uLabel, SleepTime, Precision);
	// NOTE: so far, precision is ignored

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtSocketHandler,
		ExtSocketHandler(pUS, SleepTime, recv_flags, send_flags));
    pedantic_cout("In ReadExtSocketHandler at end" << std::endl);
	return pEFH;
#else // ! USE_SOCKET
	silent_cerr("ExtSocketHandler not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
}
