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

#ifdef USE_SOCKET

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
  #include <ws2tcpip.h>
#else
  /* Assume that any non-Windows platform uses POSIX-style sockets instead. */
  // #include <fcntl.h> // TODO: nothing appears to be used from fcntl.h?
  #include <unistd.h>
  #include <fcntl.h>
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <netinet/in.h>
  #include <netdb.h>
  #include <sys/un.h>
  #include <arpa/inet.h>
  #include <errno.h>
  #include <sys/poll.h>
  #include <netinet/tcp.h>
#endif

#include "sock.h"

int
mbdyn_make_inet_socket(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int dobind, int *perrno)
{
	return mbdyn_make_inet_socket_type(sock, name, hostname, port, MBDYN_DEFAULT_SOCKET_TYPE, dobind, perrno);
}

int
mbdyn_host2inet_addr(struct sockaddr_in *name, const char *hostname, unsigned short int port, int socket_type, int *perrno)
{
	if (hostname) {
#if defined(HAVE_GETADDRINFO)
		char portbuf[sizeof("65535") + 1];
		struct addrinfo hints = { 0 }, *res = NULL;
		int rc;

		rc = snprintf(portbuf, sizeof(portbuf), "%d", (int)port);
		if (rc > STRLENOF("65535")) {
			return -4;
		}

		hints.ai_family = AF_INET;
		hints.ai_socktype = socket_type;
		rc = getaddrinfo(hostname, portbuf, &hints, &res);
		if (rc != 0) {
			*perrno = WSAGetLastError();
			return -3;
		}

		name->sin_addr = ((struct sockaddr_in *)res->ai_addr)->sin_addr;

		freeaddrinfo(res);

#elif defined(HAVE_GETHOSTBYNAME)
		struct hostent *hostinfo;

		/* TODO: use getnameinfo() if available */
		hostinfo = gethostbyname(hostname);
		if (hostinfo == NULL) {
			*perrno = h_errno;
			return -3;
		}

		name->sin_addr = *(struct in_addr *)hostinfo->h_addr;
#elif defined(HAVE_INET_ATON)
		struct in_addr addr;
		if (inet_aton(hostname, &addr) == 0) {
			*perrno = WSAGetLastError();
			return -3;
		}
		name->sin_addr = addr.s_addr;
#else
		return -3;
#endif
	} else {
		name->sin_addr.s_addr = htonl(INADDR_ANY);
	}

	return 0;
}

int
mbdyn_make_inet_socket_type(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int socket_type, int dobind, int *perrno)
{
	struct sockaddr_in tmpname = { 0 };

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Give the socket a name. */
   	name->sin_family = AF_INET;
   	name->sin_port = htons(port);

	int rc = mbdyn_host2inet_addr(name, hostname, port, socket_type, perrno);
	if (rc != 0) {
		return rc;
	}

   	/* Create the socket. */
	*sock = socket(PF_INET, socket_type, 0);
   	//bool isblocking = true;
   	if (*sock == INVALID_SOCKET) {
		if (perrno) {
			*perrno = WSAGetLastError();
		}
		return -1;
   	}

   	/* disable Nagle's algorithm */
	int flag = 1;
	int result = setsockopt(*sock,            /* socket affected */
				IPPROTO_TCP,     /* set option at TCP level */
				TCP_NODELAY,     /* name of option */
				(char *) &flag,  /* the cast is historical cruft */
				sizeof(int));    /* length of option value */

	if (dobind) {
		rc = bind(*sock, (struct sockaddr *) name, sizeof(struct sockaddr_in));
		if (rc == SOCKET_ERROR) {
			if (perrno) {
				*perrno = WSAGetLastError();
			}
			return -2;
		}
   	}

   	return 0;
}



char* sock_err_string (int err)
{
#ifdef _win32
    char msg[100];
    winsock_err_string (err, msg);
#else
    char* msg = strerror(err);
#endif /* _win32 */

    return msg;
}

#ifdef _WIN32
void winsock_err_string (int err, char* msg)
{
    switch (err) {
    case 0:                  strcpy (msg, "No error"); break;
    case WSAEINTR:           strcpy (msg, "Interrupted system call"); break;
    case WSAEBADF:           strcpy (msg, "Bad file number"); break;
    case WSAEACCES:          strcpy (msg, "Permission denied"); break;
    case WSAEFAULT:          strcpy (msg, "Bad address"); break;
    case WSAEINVAL:          strcpy (msg, "Invalid argument"); break;
    case WSAEMFILE:          strcpy (msg, "Too many open sockets"); break;
    case WSAEWOULDBLOCK:     strcpy (msg, "Operation would block"); break;
    case WSAEINPROGRESS:     strcpy (msg, "Operation now in progress"); break;
    case WSAEALREADY:        strcpy (msg, "Operation already in progress"); break;
    case WSAENOTSOCK:        strcpy (msg, "Socket operation on non-socket"); break;
    case WSAEDESTADDRREQ:    strcpy (msg, "Destination address required"); break;
    case WSAEMSGSIZE:        strcpy (msg, "Message too long"); break;
    case WSAEPROTOTYPE:      strcpy (msg, "Protocol wrong type for socket"); break;
    case WSAENOPROTOOPT:     strcpy (msg, "Bad protocol option"); break;
    case WSAEPROTONOSUPPORT: strcpy (msg, "Protocol not supported"); break;
    case WSAESOCKTNOSUPPORT: strcpy (msg, "Socket type not supported"); break;
    case WSAEOPNOTSUPP:      strcpy (msg, "Operation not supported on socket"); break;
    case WSAEPFNOSUPPORT:    strcpy (msg, "Protocol family not supported"); break;
    case WSAEAFNOSUPPORT:    strcpy (msg, "Address family not supported"); break;
    case WSAEADDRINUSE:      strcpy (msg, "Address already in use"); break;
    case WSAEADDRNOTAVAIL:   strcpy (msg, "Can't assign requested address"); break;
    case WSAENETDOWN:        strcpy (msg, "Network is down"); break;
    case WSAENETUNREACH:     strcpy (msg, "Network is unreachable"); break;
    case WSAENETRESET:       strcpy (msg, "Net connection reset"); break;
    case WSAECONNABORTED:    strcpy (msg, "Software caused connection abort"); break;
    case WSAECONNRESET:      strcpy (msg, "Connection reset by peer"); break;
    case WSAENOBUFS:         strcpy (msg, "No buffer space available"); break;
    case WSAEISCONN:         strcpy (msg, "Socket is already connected"); break;
    case WSAENOTCONN:        strcpy (msg, "Socket is not connected"); break;
    case WSAESHUTDOWN:       strcpy (msg, "Can't send after socket shutdown"); break;
    case WSAETOOMANYREFS:    strcpy (msg, "Too many references, can't splice"); break;
    case WSAETIMEDOUT:       strcpy (msg, "Connection timed out"); break;
    case WSAECONNREFUSED:    strcpy (msg, "Connection refused"); break;
    case WSAELOOP:           strcpy (msg, "Too many levels of symbolic links"); break;
    case WSAENAMETOOLONG:    strcpy (msg, "File name too long"); break;
    case WSAEHOSTDOWN:       strcpy (msg, "Host is down"); break;
    case WSAEHOSTUNREACH:    strcpy (msg, "No route to host"); break;
    case WSAENOTEMPTY:       strcpy (msg, "Directory not empty"); break;
    case WSAEPROCLIM:        strcpy (msg, "Too many processes"); break;
    case WSAEUSERS:          strcpy (msg, "Too many users"); break;
    case WSAEDQUOT:          strcpy (msg, "Disc quota exceeded"); break;
    case WSAESTALE:          strcpy (msg, "Stale NFS file handle"); break;
    case WSAEREMOTE:         strcpy (msg, "Too many levels of remote in path"); break;
    case WSASYSNOTREADY:     strcpy (msg, "Network system is unavailable"); break;
    case WSAVERNOTSUPPORTED: strcpy (msg, "Winsock version out of range"); break;
    case WSANOTINITIALISED:  strcpy (msg, "WSAStartup not yet called"); break;
    case WSAEDISCON:         strcpy (msg, "Graceful shutdown in progress"); break;
    case WSAHOST_NOT_FOUND:  strcpy (msg, "Host not found"); break;
    case WSANO_DATA:         strcpy (msg, "No host data of that type was found"); break;

    }
}

#else
/* errno cannot be used on windows for sockets, so we use
   WSAGetLastError but on non-windows define it to just return
   errno */
int WSAGetLastError()
{
    return errno;
}

int
mbdyn_make_named_socket(struct sockaddr_un *name, const char *path,
	int dobind, int *perrno)
{
	return mbdyn_make_named_socket_type(name, path, MBDYN_DEFAULT_SOCKET_TYPE, dobind, perrno);
}

int
mbdyn_make_named_socket_type(struct sockaddr_un *name, const char *path,
	int socket_type, int dobind, int *perrno)
{
   	int sock = -1;

   	struct sockaddr_un tmpname = { 0 };
	socklen_t size;

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Create the socket. */
   	sock = socket(PF_LOCAL, socket_type, 0);
   	if (sock < 0) {
		if (perrno) {
			*perrno = errno;
		}
      		return -1;
   	}

   	/* Give the socket a name. */
   	name->sun_family = AF_LOCAL;
   	strncpy(name->sun_path, path, sizeof(name->sun_path));
#ifdef HAVE_OFFSETOF
	size = (offsetof(struct sockaddr_un, sun_path)
			+ strlen(name->sun_path) + 1);
#else /* HAVE_OFFSETOF */
	size = sizeof(struct sockaddr_un);
#endif /* !HAVE_OFFSETOF */

   	if (dobind) {
		int rc = bind(sock, (struct sockaddr *)name, size);
		if (rc < 0) {
			if (perrno) {
				*perrno = errno;
			}
      			return -2;
		}
   	}

   	return sock;
}
#endif /* _WIN32 */

#endif /* USE_SOCKET */
