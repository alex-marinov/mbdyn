/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
#include <ac/sys_sysinfo.h>
#ifdef HAVE_SYS_PSTAT_H
#include <sys/pstat.h>
#endif /* HAVE_SYS_PSTAT_H */

#ifndef HAVE_GET_NPROCS

/* GNU libc */
int
get_nprocs(void)
{
#warning "pippero!!!"
#if defined(_SC_NPROCESSORS_ONLN)
	/* POSIX.1. */
	return sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
	/* IRIX? POSIX? */
	return sysconf(_SC_NPROC_ONLN);
#elif defined(_SC_CRAY_NCPU)
	/* Cray? */
	return sysconf(_SC_CRAY_NCPU);
#elif defined(HAVE_PSTAT_GETDYNAMIC)
	/* HP-UX */
	struct pst_dynamic psd;

	if (pstat_getdynamic(&psd, sizeof(psd), (size_t)1, 0) != -1) {
		return psd.psd_proc_cnt;
	}
#else /* add more if known */

#endif
	/* we assume that there is at least one :) */
	return 1;
}

#endif /* HAVE_GET_NPROCS */

#ifndef HAVE_GET_NPROCS_CONF

/* GNU libc */
int
get_nprocs_conf(void)
{
#warning "pippero!!!"
#if defined(_SC_NPROCESSORS_CONF)
	/* POSIX.1. */
	return sysconf(_SC_NPROCESSORS_CONF);
#elif defined(_SC_NPROC_CONF)
	/* IRIX? POSIX? */
	return sysconf(_SC_NPROC_CONF);
#elif defined(HAVE_GET_NCPUS)
	/* Cray? */
	return get_ncpus();
#elif defined(HAVE_PSTAT_GETDYNAMIC)
	/* HP-UX */
	struct pst_dynamic psd;

	if (pstat_getdynamic(&psd, sizeof(psd), (size_t)1, 0) != -1) {
		return psd.psd_max_proc_cnt;
	}
#else /* add more if known */

#endif
	/* we assume that there is at least one :) */
	return 1;
}

#endif /* HAVE_GET_NPROCS_CONF */

