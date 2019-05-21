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

#ifndef SOCK_H
#define SOCK_H

#include "mbconfig.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef HAVE_SOCKLEN_T
typedef int socklen_t;
#endif /* HAVE_SOCKLEN_T */

#ifndef _WIN32
/* We define the SOCKET type if not defined already */
#ifndef SOCKET_TYPEDEF
#define SOCKET_TYPEDEF
typedef int SOCKET;
#endif /* SOCKET_TYPEDEF */
#endif /* ! _WIN32 */

/** Creates an inet socket of default type (SOCK_STREAM)
 *
 *  param name Output sockaddr_in structure which will have a name of the socket added
 *    a new sockaddr_in structure will be created if this is NULL
 *  param hostname Input A pointer to a NULL-terminated ANSI string that contains a host
 *    (node) name or a numeric host address string. The numeric host address string is an
 *    IPv4 address or an IPv6 hex address.
 *  param port Input port number on which which to connect
 *  param dobind Input flag determining whether to bind to the socket
 *  param perror Output error returned when the socket was created
 *
 */
extern int
mbdyn_make_inet_socket(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int dobind, int *perror);

/** Creates an inet socket
 *
 *  param name Output sockaddr_in structure which will have a name of the socket added
 *    a new sockaddr_in structure will be created if this is NULL
 *  param hostname Input A pointer to a NULL-terminated ANSI string that contains a host
 *    (node) name or a numeric host address string. The numeric host address string is an
 *    IPv4 address or an IPv6 hex address.
 *  param port Input port number on which which to connect
 *  param socket_type Input type of socket to be created, can be SOCK_STREAM or
 *    SOCK_DGRAM
 *  param dobind Input flag determining whether to bind to the socket
 *  param perror Output error returned when the socket was created
 *
 */
extern int
mbdyn_make_inet_socket_type(SOCKET* sock, struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int socket_type, int dobind, int *perror);

char* sock_err_string (int err);

#ifdef _WIN32

/* h_errno is not used on windows for sockets, we can't redefine errno, but
   it should be safe enough to redefine h_errno which is
   networking specific */
#define h_errno WSAGetLastError()


void winsock_err_string (int err, char* msg);

#else

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif /* ! SOCKET_ERROR */

/* winsock sockets return INVALID_SOCKET rather that -1
   so we define it here for BSD sockets */
#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif /* INVALID_SOCKET */

/* errno does not work on windows with winsock, you must call
   WSAGetLastError instead. On other platforms we define this function
   here, which simply returns errno */
extern int
WSAGetLastError(void);

/** Creates a local (unix) socket of default type (SOCK_STREAM)
 *
 *  param name Output sockaddr_un structure which will have a name of the socket added
 *    a new sockaddr_in structure will be created if this is NULL
 *  param path Input A pointer to a NULL-terminated ANSI string that contains a path
 *    for the socket connection
 *  param port Input port number on which which to connect
 *  param dobind Input flag determining whether to bind to the socket
 *  param perror Output error returned when the socket was created
 *
 */
extern int
mbdyn_make_named_socket(struct sockaddr_un *name, const char *path,
	int dobind, int *perror);

/** Creates a local (unix) socket
 *
 *  param name Output sockaddr_un structure which will have a name of the socket added
 *    a new sockaddr_un structure will be created if this is NULL
 *  param path Input A pointer to a NULL-terminated ANSI string that contains a path
 *    for the socket connection
 *  param port Input port number on which which to connect
 *  param socket_type Input type of socket to be created, can be SOCK_STREAM or
 *    SOCK_DGRAM
 *  param dobind Input flag determining whether to bind to the socket
 *  param perror Output error returned when the socket was created
 *
 */
extern int
mbdyn_make_named_socket_type(struct sockaddr_un *name, const char *path,
	int socket_type, int dobind, int *perror);
#endif /* _WIN32 */

/* These have been erased from FreeBSD's headers, so we have to 
 * work around with this ugly hack. This is not exactly the same,
 * but it seems to work.
 * I chose this place for the definitions, since there are different
 * files using these macros and so we just need to define them once.
 */
#ifndef IPPORT_USERRESERVED
#define IPPORT_USERRESERVED IPPORT_RESERVED
#endif /* ! IPPORT_USERRESERVED */

/* MBDyn's default socket type is tcp (SOCK_STREAM) (not udp, SOCK_DGRAM) */
#define MBDYN_DEFAULT_SOCKET_TYPE SOCK_STREAM

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* SOCK_H */
