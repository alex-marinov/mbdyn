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

#ifdef USE_SOCKET

#include "myassert.h"
#include "mynewmem.h"
#include "mbsleep.h"

#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#include <stdlib.h>
#include <errno.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <time.h>

#include "usesock.h"
#include "sock.h"

#define DEFAULT_PORT	5500 /* FIXME: da definire meglio */
#define DEFAULT_HOST 	"127.0.0.1"

UseSocket::UseSocket(bool c)
: sock(-1),
socket_type(MBDYN_DEFAULT_SOCKET_TYPE),
create(c),
connected(false),
abandoned(false),
socklen(0)
{
	NO_OP;
}

UseSocket::UseSocket(int t, bool c)
: sock(-1),
socket_type(t),
create(c),
connected(false), // t == SOCK_DGRAM),
abandoned(false),
socklen(0)
{
	NO_OP;
}

UseSocket::~UseSocket(void)
{
	if (sock != -1) {
		int status;
		time_t t = time(NULL);

		if (socket_type == SOCK_STREAM) {
			status = shutdown(sock, SHUT_RDWR);
			int save_errno = errno;
			if (status == 0) {
				pedantic_cout("UseSocket::~UseSocket: shutdown: "
					"socket=" << sock
					<< " status=" << status
					<< " time=" << asctime(localtime(&t))
					<< std::endl);

			} else {
				char *msg = strerror(save_errno);
				silent_cerr("UseSocket::~UseSocket: shutdown "
					"socket=" << sock
					<< " status=" << status
					<< " time=" << asctime(localtime(&t))
					<< " error=" << save_errno << " (" << msg << ")"
					<< std::endl);
			}
		}

		status = close(sock);
		pedantic_cout("UseSocket::~UseSocket: close: "
			"socket=" << sock
			<< " status=" << status
			<< " time=" << asctime(localtime(&t))
			<< std::endl);
	}
}

std::ostream&
UseSocket::Restart(std::ostream& out) const
{
	const char *sType;
	switch (socket_type) {
	case SOCK_STREAM:
		sType = "tcp";
		break;
	case SOCK_DGRAM:
		sType = "udp";
		break;
	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	out << ", socket type, " << sType;

	if (create) {
		out << ", create, yes";
	}

	return out;
}

int
UseSocket::GetSock(void) const
{
	return sock;
}

void
UseSocket::SetSock(int s)
{
	sock = s;
}

void
UseSocket::PostConnect(void)
{
	struct linger lin;
	lin.l_onoff = 1;
	lin.l_linger = 0;
	
	if (setsockopt(GetSock(), SOL_SOCKET, SO_LINGER, &lin, sizeof(lin))) {
		int save_errno = errno;
		char *msg = strerror(save_errno);

      		silent_cerr("UseSocket::PostConnect: setsockopt() failed "
			"(" << save_errno << ": " << msg << ")"
			<< std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
UseSocket::Connect(void)
{
	if (socket_type != SOCK_STREAM) {
		return;
	}

	// FIXME: retry strategy should be configurable
	int count = 600;
	mbsleep_t timeout;
	mbsleep_real2sleep(0.1, &timeout);

	for ( ; count > 0; count--) {
		if (connect(sock, GetSockaddr(), GetSocklen()) < 0) {
			int save_errno = errno;
			switch (save_errno) {
			case ECONNREFUSED:	// inet
			case ENOENT:		// unix
				/* Socket does not exist yet; retry */
				mbsleep(&timeout);
				continue;
			}

			/* Connect failed */
			char *msg = strerror(save_errno);
			silent_cerr("UseSocket::Connect: connect() failed "
				"(" << save_errno << ": " << msg << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* Success */
		connected = true;
		PostConnect();
		return;
	}

	silent_cerr("UseSocket(): connection timed out"
		<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
UseSocket::ConnectSock(int s)
{
	SetSock(s);

	connected = true;

	PostConnect();
}

bool
UseSocket::Create(void) const
{
	return create;
}

bool
UseSocket::Connected(void) const
{
	return connected;
}

void
UseSocket::Abandon(void)
{
	abandoned = true;
}

bool
UseSocket::Abandoned(void) const
{
	return abandoned;
}

socklen_t&
UseSocket::GetSocklen(void) const
{
	return socklen;
}

ssize_t
UseSocket::send(const void *buf, size_t len, int flags)
{
	switch (socket_type) {
	case SOCK_STREAM:
		return ::send(sock, buf, len, flags);

	case SOCK_DGRAM:
		return ::sendto(sock, buf, len, flags, GetSockaddr(), GetSocklen());

	default:
		ASSERT(0);
	}

	return -1;
}

ssize_t
UseSocket::recv(void *buf, size_t len, int flags)
{
	switch (socket_type) {
	case SOCK_STREAM:
		return ::recv(sock, buf, len, flags);

	case SOCK_DGRAM:
		// struct sockaddr src_addr = { 0 };
		// socklen_t addrlen = 0;
		// return ::recvfrom(sock, buf, len, flags, &src_addr, &addrlen);
		return ::recvfrom(sock, buf, len, flags, 0, 0);

	default:
		ASSERT(0);
	}

	return -1;
}

UseInetSocket::UseInetSocket(const std::string& h, unsigned short p, bool c)
: UseSocket(c),
host(h),
port(p)
{
	UseInetSocket_int();
}

UseInetSocket::UseInetSocket(const std::string& h, unsigned short p, int t, bool c)
: UseSocket(t, c),
host(h),
port(p)
{
	UseInetSocket_int();
}

void
UseInetSocket::UseInetSocket_int(void)
{
   	ASSERT(port > 0);

	socklen = sizeof(addr);

	/* if 0, means INET on localhost */
	if (host.empty()) {
		host = DEFAULT_HOST;
	}

	if (create) {
		int			save_errno;
		
   		sock = mbdyn_make_inet_socket_type(0, host.c_str(), port, socket_type,
			1, &save_errno);
		
		if (sock == -1) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("UseInetSocket(" << host << ":" << port << "): "
				"socket() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

   		} else if (sock == -2) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("UseInetSocket(" << host << ":" << port << "): "
				"bind() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

   		} else if (sock == -3) {
      			silent_cerr("UseInetSocket(" << host << ":" << port << "): "
				"illegal host name \"" << host << "\" "
				"(" << save_errno << ")"
				<< std::endl);
      			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   		}

		if (socket_type == SOCK_STREAM) {
   			if (listen(sock, 1) < 0) {
				save_errno = errno;
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("UseInetSocket(" << host << ":" << port << "): "
					"listen() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
      				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
   		}
	}
}

UseInetSocket::~UseInetSocket(void)
{
	NO_OP;
}

std::ostream&
UseInetSocket::Restart(std::ostream& out) const
{
	UseSocket::Restart(out);
	out << ", port, " << port;
	if (!host.empty()) {
		out << ", host, " << "\"" << host << "\"";
	}
	return out;
}

void
UseInetSocket::Connect(void)
{
	if (connected) {
		return;
	}

#if 0
	sock = socket(PF_INET, SOCK_STREAM, 0);
	if (sock < 0) {
		silent_cerr("UseSocket(): socket() failed "
				<< "\"" << host << ":" << port << "\""
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	addr.sin_family = AF_INET;
	addr.sin_port = htons(port);
	if (inet_aton(host.c_str(), &addr.sin_addr) == 0) {
		silent_cerr("UseSocket(): unknown host "
				"\"" << host << ":" << port << "\""
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);	
	}
#endif

	int save_errno;
	sock = mbdyn_make_inet_socket_type(&addr, host.c_str(), port, socket_type, 0, &save_errno);
	if (sock < 0) {
		const char	*err_msg = strerror(save_errno);

		silent_cerr("UseInetSocket(" << host << ":" << port << "): socket() failed "
			"(" << save_errno << ": " << err_msg << ")"
			<< std::endl);
      		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pedantic_cout("connecting to inet socket "
			"\"" << inet_ntoa(addr.sin_addr) 
			<< ":" << ntohs(addr.sin_port)
			<< "\" ..." << std::endl);

	UseSocket::Connect();
}
	
void
UseInetSocket::ConnectSock(int s)
{
	UseSocket::ConnectSock(s);
	
	silent_cout("INET connection to port=" << port << " by client "
			<< "\"" << inet_ntoa(addr.sin_addr)
			<< ":" << ntohs(addr.sin_port) << "\""
			<< std::endl);
}

struct sockaddr *
UseInetSocket::GetSockaddr(void) const
{
	return (struct sockaddr *)&addr;
}
	
UseLocalSocket::UseLocalSocket(const std::string& p, bool c)
: UseSocket(c), path(p)
{
	UseLocalSocket_int();
}

UseLocalSocket::UseLocalSocket(const std::string& p, int t, bool c)
: UseSocket(t, c), path(p)
{
	UseLocalSocket_int();
}

void
UseLocalSocket::UseLocalSocket_int(void)
{
	ASSERT(!path.empty());

	socklen = sizeof(addr);
	
	if (create) {
		int		save_errno;

		sock = mbdyn_make_named_socket_type(0, path.c_str(), socket_type, 1, &save_errno);
		
		if (sock == -1) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("UseLocalSocket(" << path << "): "
				"socket() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

   		} else if (sock == -2) {
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("UseLocalSocket(" << path << "): "
				"bind() failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   		}

		if (socket_type == SOCK_STREAM) {
   			if (listen(sock, 1) < 0) {
				save_errno = errno;
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("UseLocalSocket(" << path << "): "
					"listen() failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
      				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   			}
		}
	}	 
}

UseLocalSocket::~UseLocalSocket(void)
{
	if (!path.empty() && create) {
		unlink(path.c_str());
	}
}

std::ostream&
UseLocalSocket::Restart(std::ostream& out) const
{
	return UseSocket::Restart(out) << ", path, " << "\"" << path << "\"";
}

void
UseLocalSocket::Connect(void)
{
	if (connected) {
		return;
	}

	sock = socket(PF_LOCAL, socket_type, 0);
	if (sock < 0) {
		int save_errno = errno;
		char *msg = strerror(save_errno);

		silent_cerr("UseSocket(): socket() failed "
				"(" << save_errno << ": " << msg << ")"
				<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	addr.sun_family = AF_UNIX;
	if (path.size() >= sizeof(addr.sun_path)) {
		silent_cerr("UseSocket(): path=\"" << path << "\" exceeds allowable size "
			"of LOCAL socket path (" << sizeof(addr.sun_path) << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	memcpy(addr.sun_path, path.c_str(), path.size() + 1); // terminating '\0'

	pedantic_cout("connecting to local socket \"" << path << "\" ..."
		<< std::endl);

	UseSocket::Connect();
}

void
UseLocalSocket::ConnectSock(int s)
{
	UseSocket::ConnectSock(s);
	
	silent_cout("LOCAL connection to path=\"" << path << "\""
		<< std::endl);
}

struct sockaddr *
UseLocalSocket::GetSockaddr(void) const
{
	return (struct sockaddr *)&addr;
}

#endif // USE_SOCKET
