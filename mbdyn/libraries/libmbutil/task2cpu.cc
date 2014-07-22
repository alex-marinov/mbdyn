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

#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#ifdef HAVE_TASK2CPU
#include <sys/ioctl.h>
#endif // HAVE_TASK2CPU

#include "ac/pthread.h"

#ifdef USE_PTHREAD_SETAFFINITY_NP
#include <sched.h>
#endif

#include <iostream>

static bool		mbdyn_task2cpu_disabled = false;
#ifdef HAVE_THREADS
static pthread_mutex_t	mbdyn_task2cpu_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif /* HAVE_THREADS */

int
mbdyn_task2cpu(int cpu)
{
	int		fd = -1;

#ifdef HAVE_THREADS
	pthread_mutex_lock(&::mbdyn_task2cpu_mutex);
#endif /* HAVE_THREADS */
	if (!::mbdyn_task2cpu_disabled) {
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
			::mbdyn_task2cpu_disabled = true;
		}
#elif defined(USE_PTHREAD_SETAFFINITY_NP)
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(cpu, &cpuset);

		const pthread_t thread = pthread_self();

		int s = pthread_setaffinity_np(thread, sizeof(cpuset), &cpuset);

		if (0 != s) {
			silent_cerr("warning: pthread_setaffinity_np failed with status " << s << std::endl);
			goto failed;
		}

#ifdef DEBUG
		s = pthread_getaffinity_np(thread, sizeof(cpuset), &cpuset);

		if (0 != s) {
			DEBUGCERR("warning: pthread_getaffinity_np failed with status " << s << std::endl);
			goto failed;
		}

		for (int i = 0; i < CPU_SETSIZE; ++i) {
			if (CPU_ISSET(i, &cpuset)) {
				DEBUGCERR("thread " << cpu + 1 << " running at CPU " << i << std::endl);
			}
		}
#endif
		fd = 1; // return 0

	failed:
		if (fd == -1) {
			::mbdyn_task2cpu_disabled = true;
		}
#else  /* ! HAVE_TASK2CPU */
		silent_cerr("/dev/TASK2CPU or pthread_getaffinity_np are not available" << std::endl);
		::mbdyn_task2cpu_disabled = true;
#endif /* ! HAVE_TASK2CPU */
	}
#ifdef HAVE_THREADS
	pthread_mutex_unlock(&::mbdyn_task2cpu_mutex);
#endif /* HAVE_THREADS */

	return (fd == -1);
}

