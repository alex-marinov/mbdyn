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
#include "dataman.h"
#include "dataman_.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
  /* See http://stackoverflow.com/questions/12765743/getaddrinfo-on-win32 */
  #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0501  /* Windows XP. */
  #endif
  #include <winsock2.h>
  #include <ws2tcpip.h>
#else
  #include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>
#endif /* _WIN32 */

#include "sock.h"

void
DataManager::RegisterSocketUser(UseSocket *pUS)
{
	ASSERT(pUS->GetSock() != INVALID_SOCKET);
	SocketUsers[SocketUsers.size()] = pUS;
}

void
DataManager::DeleteSocketUsers(void)
{
	for (std::map<int, UseSocket *>::iterator i = SocketUsers.begin();
			i != SocketUsers.end(); ++i)
	{
		delete i->second;
	}
	SocketUsers.clear();
}

void
DataManager::WaitSocketUsers(void)
{
	if (SocketUsers.empty()) {
		return;
	}

	time_t finalTime = 0;
	if (SocketUsersTimeout != 0) {
		finalTime = time(NULL) + SocketUsersTimeout;
	}

	/* declare two sets of sockets, zeroing them out with FD_ZERO */
	fd_set	active_set, read_set;
	FD_ZERO(&active_set);

	/* insert the already registered sockets into active set */
	std::map<int, UseSocket *>::iterator ri;
	std::map<int, UseSocket *>::const_iterator re = SocketUsers.end();
	int nactive; /* count how many are added */
	for (nactive = 0, ri = SocketUsers.begin(); ri != re; ++nactive, ++ri) {
		FD_SET(ri->second->GetSock(), &active_set);
	}

	/* wait for all registered */
	while (nactive > 0) {

        pedantic_cout ("DataManager::WaitSocketUsers(): "
						"In connection loop with " << nactive
						<< " active sockets to process (nactive)" << std::endl);

		struct timeval	timeout, *timeoutp = NULL;
        /* check if user specified timeout for connection to be made
           successfully has passed. If finalTime is zero, there is
           no timeout and we will keep trying forever. */
		if (finalTime != 0) {
			timeout.tv_sec = time(NULL);
			if (timeout.tv_sec >= finalTime) {
do_timeout:;
				silent_cerr("DataManager::WaitSocketUsers(): "
						"timeout " << SocketUsersTimeout
						<< " s exceeded" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			timeout.tv_sec = finalTime - timeout.tv_sec;
			timeout.tv_usec = 0;

			timeoutp = &timeout;
		}

		read_set = active_set;

		/* FD_SETSIZE determines the maximum number of descriptors in a set. */
		int a = select(FD_SETSIZE, &read_set, 0, 0, timeoutp);

        pedantic_cout ("DataManager::WaitSocketUsers(): "
						"Called select on all active sockets. 'a' is " << a << std::endl);

		switch (a) {
		case SOCKET_ERROR: {
			int save_errno = WSAGetLastError();
			char *msg = sock_err_string(save_errno);

			silent_cerr("select() failed: " << save_errno << " "
				"(" << msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		case 0:
			 goto do_timeout;
		}

		pedantic_cout ("DataManager::WaitSocketUsers(): "
						   "Looping on active to see if connection is possible" << std::endl);

		/* loop on active to see what is being connected */
		for (int i = 0; i < nactive && a != SOCKET_ERROR; i++) {

			UseSocket *pUS = SocketUsers[i];

			int isset = FD_ISSET(pUS->GetSock(), &read_set);

            pedantic_cout ("DataManager::WaitSocketUsers(): "
						   "Socket: " << i
						   << " Info: " << pUS->GetSockaddrStr()
						   << " isset: " << isset
						   << std::endl);

			if (isset) {

                pedantic_cout ("DataManager::WaitSocketUsers(): "
						   "Attempting to accept socket: " << i
						   << " Info: " << pUS->GetSockaddrStr()
						   << std::endl);

				SOCKET sock;

				sock = accept(pUS->GetSock(), pUS->GetSockaddr(),
						&pUS->GetSocklen());

				if (sock == INVALID_SOCKET) {
					int save_errno = WSAGetLastError();
					char *msg = sock_err_string(save_errno);

					silent_cerr("accept() failed: "
							<< save_errno << " "
							"(" << msg << ")"
							<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* remove accepted from set */
				FD_CLR(pUS->GetSock(), &active_set);

#ifdef _WIN32
				closesocket (pUS->GetSock());
#else
				close(pUS->GetSock());
#endif /* _WIN32 */

				SocketUsers.erase(i);

				/* register as connected */
				pUS->ConnectSock(sock);
				nactive--;
				a--;
			}
		}
	}

    pedantic_cout ("DataManager::WaitSocketUsers(): "
					"Leaving processing loop with " << nactive << " active sockets" << std::endl);
}

#endif // USE_SOCKET
