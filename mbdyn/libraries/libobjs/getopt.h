/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 */

#ifndef GETOPT_H
#define GETOPT_H

#ifndef HAVE_GETOPT_H

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

#endif /* !HAVE_GETOPT_H */

#endif /* GETOPT_H */

