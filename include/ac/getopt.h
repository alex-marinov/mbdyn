/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/*
 * This header file comes from the GNU implementation of the libc 
 * `getopt' function. The body of the function comes from the OpenLDAP
 * 
 * 		http://www.openldap.org
 * 
 * implementation for those systems that do not have a native one.
 * The intellectual property of the OpenLDAP implementation remains 
 * with the original Authors, whose contribution is kindly appreciated.
 *
 * The original copyright statement follows.
 */

/* Declarations for getopt.
 * Copyright (C) 1989,90,91,92,93,94,96,97,98 Free Software Foundation, Inc.
 * This file is part of the GNU C Library.
 *
 * The GNU C Library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * The GNU C Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with the GNU C Library; see the file COPYING.LIB.  If not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.  */

#ifndef GETOPT_H
#define GETOPT_H

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else /* !HAVE_GETOPT_H */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <unistd.h>

extern int getopt(int argc, char * const argv[], const char *optstring);

extern char *optarg;
extern int optind, opterr, optopt;

#if 0
extern int getopt_long(int argc, 
                       char * const argv[],
		       const char *optstring,
		       const struct option *longopts, int *longindex);
#endif /* 0 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* !HAVE_GETOPT_H */

#endif /* GETOPT_H */

