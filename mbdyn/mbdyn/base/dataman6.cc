/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>

void
DataManager::RegisterSocketUser(UseSocket *pUS)
{
	ASSERT(pUS->GetSock() != -1);
	SocketUsers[pUS->GetSock()] = pUS;
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
	
	fd_set	active_set, read_set;
	FD_ZERO(&active_set);

	/* insert registered sockets in set */
	std::map<int, UseSocket *>::iterator ri;
	std::map<int, UseSocket *>::const_iterator re = SocketUsers.end();
	int n;
	for (n = 0, ri = SocketUsers.begin(); ri != re; ++n, ++ri) {
		FD_SET(ri->first, &active_set);
	}

	/* wait for all registered */
	while (n > 0) {
		struct timeval	timeout, *timeoutp = 0;

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
		int a = select(FD_SETSIZE, &read_set, 0, 0, timeoutp);
		switch (a) {
		case -1: {
			int save_errno = errno;
			char *msg = strerror(save_errno);

			silent_cerr("select() failed: " << save_errno << " "
				"(" << msg << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		case 0:
			 goto do_timeout;
		}

		/* loop on active to see what is being connected */
		for (int i = 0; i < FD_SETSIZE && a > 0; i++) {
			if (FD_ISSET(i, &read_set)) {
				UseSocket *pUS = SocketUsers[i];
				int sock;
				
				sock = accept(i, pUS->GetSockaddr(),
						&pUS->GetSocklen());
				if (sock < 0) {
					int save_errno = errno;
					char *msg = strerror(save_errno);

					silent_cerr("accept() failed: "
							<< save_errno << " "
							"(" << msg << ")"
							<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* remove accepted from set */
				FD_CLR(i, &active_set);
				close(i);
				SocketUsers.erase(i);

				/* register as connected */
				pUS->ConnectSock(sock);
				n--;
				a--;
			}
		}
	}
}

#endif // USE_SOCKET
