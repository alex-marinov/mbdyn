/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
#include <errno.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
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

int
mbdyn_make_inet_socket(struct sockaddr_in *name, const char *hostname,
	unsigned short int port, int dobind, int *perrno)
{
   	int sock = -1;
	struct sockaddr_in tmpname;

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Give the socket a name. */
   	name->sin_family = AF_INET;
   	name->sin_port = htons(port);
	
	if (hostname) {
#if defined(HAVE_GETADDRINFO)
		char portbuf[sizeof("65535")];
		struct addrinfo hints = { 0 }, *res = NULL;
		int rc;

		rc = snprintf(portbuf, sizeof(portbuf), "%d", port);
		if (rc > STRLENOF("65535")) {
			return -4;
		}

		hints.ai_family = AF_INET;
		hints.ai_socktype = SOCK_STREAM;
		rc = getaddrinfo(hostname, portbuf, &hints, &res);
		if (rc != 0) {
			*perrno = errno;
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
			*perrno = errno;
			return -3;
		}
		name->sin_addr = addr.s_addr;
#else
		return -3;
#endif
	} else {
		name->sin_addr.s_addr = htonl(INADDR_ANY);
	}

   	/* Create the socket. */
   	sock = socket(PF_INET, SOCK_STREAM, 0);
   	if (sock < 0) {
		if (perrno) {
			*perrno = errno;
		}
      		return -1;
   	}

	if (dobind) {
		int rc = -1;
		rc = bind(sock, (struct sockaddr *) name, sizeof(struct sockaddr_in));
		if (rc < 0) {
			if (perrno) {
				*perrno = errno;
			}
			return -2;
		}
   	}

   	return sock;
}

int
mbdyn_make_named_socket(struct sockaddr_un *name, const char *path,
	int dobind, int *perrno)
{
   	int sock = -1;

   	struct sockaddr_un tmpname;
	socklen_t size;

	if (name == NULL) {
		name = &tmpname;
	}

	if (perrno) {
		*perrno = 0;
	}

   	/* Create the socket. */
   	sock = socket(PF_LOCAL, SOCK_STREAM, 0);
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
   	if (dobind && bind(sock, (struct sockaddr *) name, size) < 0) {
		if (perrno) {
			*perrno = errno;
		}
      		return -2;
   	}

   	return sock;
}

#endif /* USE_SOCKET */
