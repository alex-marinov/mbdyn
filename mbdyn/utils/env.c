/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern char **environ;

int 
main(void)
{
   	char** s = environ;
   	int c = 0;

	if (s == NULL || *s == NULL) {
		printf("environment is empty!\n");
		return EXIT_SUCCESS;
	}
	
  	do {
	 	long int i = 0;
	 	double d = 0.;
	 	char* p = NULL;
	 	char* v = NULL;
	 	char* n = NULL;

		if (strncmp(*s, "MBDYN", 5) != 0) {
			continue;
		}

		c++;
	 
	 	printf("%s\n", *s);

	 	p = strdup(*s);
	 	if (p == NULL) {
	    		fprintf(stderr, "error 1\n");
	    		exit(EXIT_FAILURE);
	 	}
	 	v = strchr(p, '=');
	 	if (v == NULL) {
	    		fprintf(stderr, "error 2\n");
	    		exit(EXIT_FAILURE);
	 	}

	 	*v = '\0';
	 	v++;
	 
	 	if (strncmp(p+5, "_real_", 6) == 0) {
	    		n = p+11;
	    		d = atof(v);
	    		printf("env=%s, var=%s, val=%s(%e, real)\n", 
				*s, n, v, d);
	 	} else if (strncmp(p+5, "_integer_", 9) == 0) {
	    		n = p+14;
	    		i = atoi(v);
	    		printf("env=%s, var=%s, val=%s(%ld, integer)\n", 
				*s, n, v, i);
	 	} else {
	    		fprintf(stderr, "error 3\n");
	    		exit(EXIT_FAILURE);
	 	}      	
   	} while (*++s);

	printf("%d entries\n", c);
	
   	return EXIT_SUCCESS;
}

