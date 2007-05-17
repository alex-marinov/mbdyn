/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "sock.h"
#include "s2s.h"

int
main(int argc, char *argv[])
{
	struct s2s_t	s2s;

	try {
		s2s.parse(argc, argv);
		s2s.prepare();

	} catch (...) {
		s2s.shutdown();
		exit(EXIT_FAILURE);
	}

	if (s2s.nChannels == 0) {
		std::cin.getline(s2s.buf, sizeof(s2s.buf));

		std::istringstream	str(s2s.buf);

		for (;; s2s.nChannels++) {
			double	d;

			str >> d;

			if (!str) {
				break;
			}
			
			s2s.dbuf.insert(s2s.dbuf.end(), d);
		}

		send(s2s.sock, (char *)&s2s.dbuf[0], sizeof(double)*s2s.nChannels, 0);

	} else {
		s2s.dbuf.resize(s2s.nChannels);
	}

	while (true) {
		for (int i = 0; i < s2s.nChannels; i++) {
			std::cin >> s2s.dbuf[i];
			if (!std::cin) {
				goto done;
			}
		}
		
		send(s2s.sock, (char *)&s2s.dbuf[0], sizeof(double)*s2s.nChannels, 0);
	}

done:;
	s2s.shutdown();
	exit(EXIT_SUCCESS);
}

#else // ! USE_SOCKET

#include <iostream>

int
main(void)
{
	std::cerr << "sockets not available" << std::endl;

	return EXIT_FAILURE;
}

#endif // ! USE_SOCKET
