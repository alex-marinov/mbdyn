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

#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <error.h>
#include <sys/ioctl.h>

#include <iostream>

int
mbdyn_task2cpu(int cpu)
{
	static bool	disabled = false;
	int		fd = -1;

	if (!disabled) {
#ifdef HAVE_TASK2CPU
		fd = open("/dev/TASK2CPU", O_RDWR);
		if (fd != -1) {
			ioctl(fd, 0, cpu);
			close(fd);

		} else {
			int save_errno = errno;
			char *err_msg = strerror(save_errno);

			silent_cerr("Error opening /dev/TASK2CPU ("
					<< save_errno << ": " << err_msg << ";"
				       " ignored)" << std::endl);
			disabled = true;
		}
#else /* ! HAVE_TASK2CPU */
		silent_cerr("/dev/TASK2CPU is not available" << std::endl);
		disabled = true;
#endif /* ! HAVE_TASK2CPU */
	}

	return (fd == -1);
}

