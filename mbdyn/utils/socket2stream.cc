/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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
	try {
		s2s_prepare(argc, argv);

	} catch (...) {
		s2s_shutdown(EXIT_FAILURE);
	}

	if (::s2s.nChannels == 0) {
		::s2s.nChannels = 1;
	}

	double	dbuf[::s2s.nChannels];
	char	*sep = " ";
	while (true) {
		int len = recv(::s2s.sock, (char *)dbuf, sizeof(double)*::s2s.nChannels, 0);

		switch (len) {
		case -1: {
			int		save_errno = errno;
			const char	*err_msg = strerror(save_errno);
			
			silent_cerr("recv(" << ::s2s.sock << ",\"" << ::s2s.buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
		}

		case 0:
			goto done;

		default:
			break;
		}

		for (int i = 0; i < ::s2s.nChannels - 1; i++) {
			std::cout << dbuf[i] << sep;
		}
		std::cout << dbuf[::s2s.nChannels - 1] << std::endl;
	}

done:;
	s2s_shutdown(EXIT_SUCCESS);
}

