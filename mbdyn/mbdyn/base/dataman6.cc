/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

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
DataManager::WaitSocketUsers(void)
{
	int n;
	
	fd_set	active_set, read_set;
	FD_ZERO(&active_set);

	/* insert registered sockets in set */
	std::map<int, UseSocket *>::iterator ri;
	std::map<int, UseSocket *>::const_iterator re = SocketUsers.end();
	for (n = 0, ri = SocketUsers.begin(); ri != re; n++, ri++) {
		FD_SET(ri->first, &active_set);
	}

	/* wait for all registered */
	while (n > 0) {
		read_set = active_set;
		if (select(FD_SETSIZE, &read_set, NULL, NULL, NULL) < 0) {
			int save_errno = errno;
			char *msg = strerror(save_errno);

			silent_cerr("select() failed: " << save_errno << " "
				"(" << msg << ")" << std::endl);
			throw ErrGeneric();
		}

		/* loop on active to see what is being connected */
		for (int i = 0; i < FD_SETSIZE; i++) {
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
					throw ErrGeneric();
				}

				/* remove accepted from set */
				FD_CLR(i, &active_set);
				close(i);

				/* register as connected */
				pUS->ConnectSock(sock);
				n--;
			}
		}
	}
}

