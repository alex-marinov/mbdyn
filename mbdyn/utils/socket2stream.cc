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

#include <cstring>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
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
	s2s_t	s2s;

	try {
		s2s.parse(argc, argv);
		s2s.prepare();

	} catch (...) {
		s2s.shutdown();
		exit(EXIT_FAILURE);
	}

	if (s2s.nChannels == 0) {
		s2s.nChannels = 1;
	}

	s2s.dbuf.resize(s2s.nChannels);
	const char	*sep = " ";
	while (true) {
		int len = s2s.recv(0);

		switch (len) {
		case -1: {
			int		save_errno = errno;

			if (save_errno == EAGAIN && !s2s.is_blocking()) {
				continue;
			}
				
			const char	*err_msg = strerror(save_errno);
			silent_cerr("recv(" << s2s.sock << ",\"" << s2s.buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
		}

		case 0:
			goto done;

		default:
			break;
		}

		for (int i = 0; i < s2s.nChannels - 1; i++) {
			std::cout << s2s.dbuf[i] << sep;
		}
		std::cout << s2s.dbuf[s2s.nChannels - 1] << std::endl;
	}

done:;
	s2s.shutdown();
	exit(EXIT_SUCCESS);
}

