/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/*
 * The body of the `getopt' function comes from the OpenLDAP
 * 
 * 		http://www.openldap.org
 * 
 * implementation for those systems that do not have a native one.
 * The intellectual property of the OpenLDAP implementation remains 
 * with the original Authors, whose contribution is kindly appreciated.
 */

/*
 * Copyright 1998-2015 The OpenLDAP Foundation, All Rights Reserved.
 * COPYING RESTRICTIONS APPLY, see COPYRIGHT file
 */
/*
	getopt.c

	modified public-domain AT&T getopt(3)
	modified by Kurt Zeilenga for inclusion into OpenLDAP
	modified by Pierangelo Masarati for inclusion with MBDyn
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifndef HAVE_GETOPT

#include <stdio.h>

#include <string.h>
#include <unistd.h>

#ifdef HAVE_IO_H
#include <io.h>
#endif

#ifndef STDERR_FILENO
#define STDERR_FILENO 2
#endif

int opterr = 1;
int optind = 1;
int optopt;
char * optarg;

static void ERR (char * const argv[], const char * s, char c)
{
	char errbuf[2];

#ifdef DF_TRACE_DEBUG
printf("DF_TRACE_DEBUG: 	static void ERR () in getopt.c\n");
#endif
	if (opterr)
	{
		errbuf[0] = c;
		errbuf[1] = '\n';
		(void) write(STDERR_FILENO,argv[0],strlen(argv[0]));
		(void) write(STDERR_FILENO,s,strlen(s));
		(void) write(STDERR_FILENO,errbuf,sizeof errbuf);
	}
}

int getopt (int argc, char * const argv [], const char * opts)
{
	static int sp = 1, error = (int) '?';
	static char sw = '-', eos = '\0', arg = ':';
	register char c, * cp;

#ifdef DF_TRACE_DEBUG
printf("DF_TRACE_DEBUG: 	int getopt () in getopt.c\n");
#endif
	if (sp == 1)
	{
		if (optind >= argc || argv[optind][0] != sw
		|| argv[optind][1] == eos)
			return EOF;
		else if (strcmp(argv[optind],"--") == 0)
		{
			optind++;
			return EOF;
		}
	}
	c = argv[optind][sp];
	optopt = (int) c;
	if (c == arg || (cp = strchr(opts,c)) == NULL)
	{
		ERR(argv,": illegal option--",c);
		if (argv[optind][++sp] == eos)
		{
			optind++;
			sp = 1;
		}
		return error;
	}
	else if (*++cp == arg)
	{
		if (argv[optind][sp + 1] != eos)
			optarg = &argv[optind++][sp + 1];
		else if (++optind >= argc)
		{
			ERR(argv,": option requires an argument--",c);
			sp = 1;
			return error;
		}
		else
			optarg = argv[optind++];
		sp = 1;
	}
	else
	{
		if (argv[optind][++sp] == eos)
		{
			sp = 1;
			optind++;
		}
		optarg = NULL;
	}
	return (int) c;
}
#endif /* HAVE_GETOPT */

