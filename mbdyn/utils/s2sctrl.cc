/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
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
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
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
	s2s_t	s2s_measures;
	s2s_t	s2s_controls;

	try {
		s2s_measures.prepare();
		s2s_controls.prepare();

	} catch (...) {
		s2s_measures.shutdown();
		s2s_controls.shutdown();
		exit(EXIT_FAILURE);
	}

	if (s2s_measures.nChannels == 0) {
		// error?
	}

	if (s2s_controls.nChannels == 0) {
		// error?
	}

	s2s_measures.dbuf.resize(s2s_measures.nChannels);
	s2s_controls.dbuf.resize(s2s_controls.nChannels);

	while (true) {
		// read new measures
		int len = s2s_measures.recv(0);

		// check sanity
		switch (len) {
		case -1: {
			int		save_errno = errno;
			const char	*err_msg = strerror(save_errno);
			
			silent_cerr("recv(" << s2s_measures.sock << ",\"" << s2s_measures.buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
		}

		case 0:
			goto done;

		default:
			break;
		}

#if 0
		// now s2s_measures.dbuf contains s2s_measures.nChannels measures 
		for (int i = 0; i < s2s_measures.nChannels - 1; i++) {
			measures[i] = s2s_measures.dbuf[i];
		}
#endif

		// use them

#if 0
		// write s2s_controls.nChannels controls to s2s_controls.dbuf
		for (int i = 0; i < s2s_controls.nChannels - 1; i++) {
			s2s_controls.dbuf[i] = controls[i];
		}
#endif

		// send new controls
		s2s_controls.send(0);
	}

done:;
	s2s_measures.shutdown();
	s2s_controls.shutdown();
	exit(EXIT_SUCCESS);
}

