/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

int get_line(FILE* fin, char* buf, int buf_size)
{
   char* s = buf;
   int c = 0;
   
   while ((c = fgetc(fin)) != '\n') {
      if (c == EOF) {
	 return 1;
      }      
      *s++ = c;
      if (s >= buf+buf_size-1) {
	 return -1;
      }
   }
   *s = '\0';
   
   return 0;
}

char* get_word(char* buf, char** next)
{
   char* s = buf;
   char* beg = NULL;
   
   if (s == NULL) {
      *next = NULL;
      return NULL;
   }
   
   while (isspace(*s)) {
      s++;
   }
   
   if (*s == '\0') {
      *next = NULL;
      return NULL;
   }
   
   beg = s;
   while (!isspace(*s)) {
      s++;
   }
   
   if (*s == '\0') {
      *next = NULL;
   } else {
      *next = s;
   }
   
   return beg;
}


int main(int argn, const char* const argv[])
{
   FILE* fin = NULL;
   int isf = 0;
   double* pd = NULL;
   double** pdv = NULL;
   int* pi = NULL;
   int ifirst = 0;
   int icurr = 0;
   int nrows = 0;
   int ncols = 0;
#define BUFSIZE 1024
   static char buf[BUFSIZE];
   char* s = NULL;
   char* next = NULL;
   int i = 0;
   int j = 0;
   
   int i1 = 0;
   int i2 = i1+1;
   int i3 = i2+1;
   
   double dt = 0.;
   double dt2 = 0.;
   
   if (argn == 1 || argn > 3) {
      fprintf(stderr, "usage: %s <file> <dt>\n", argv[0]);
      exit(EXIT_FAILURE);
   } else if (argn == 2) {
      fin = stdin;
      dt = atof(argv[1]);
   } else {
      if (!strcmp(argv[1], "-")) {
	 fin = stdin;
      } else {	 
	 fin = fopen(argv[1], "r");
	 if (fin == NULL) {
	    fprintf(stderr, "%s: file <%s> doesn't exist\n", argv[0], argv[1]);
	    exit(EXIT_FAILURE);
	 }
	 isf = 1;
      }	  
      dt = atof(argv[2]);
   }
   
   if (dt <= 0.) {
      fprintf(stderr, "%s: illegal dt = %f\n", argv[0], dt);
      exit(EXIT_FAILURE);
   }
   dt2 = 2.*dt;

   if (get_line(fin, buf, BUFSIZE) == -1) {
      fprintf(stderr, "%s: line 1 is longer that %d in file <%s>\n", 
	      argv[0], BUFSIZE, argv[1]);
      exit(EXIT_FAILURE);
   }

   s = get_word(buf, &next);
   if (s == NULL) {
      fprintf(stderr, "%s: unable to read first label in file <%s>\n",
	      argv[0], argv[1]);
      exit(EXIT_FAILURE);      
   }
   ifirst = atoi(s);
   
   while ((s = get_word(next, &next))) {
      ncols++;
   }
   
   nrows++;
   while (1) {
      if (get_line(fin, buf, BUFSIZE) == -1) {
	 fprintf(stderr, "%s: line %d is longer that %d in file <%s>\n", 
		 argv[0], nrows+1, BUFSIZE, argv[1]);
	 exit(EXIT_FAILURE);
      }
      
      if (sscanf(buf, "%d", &icurr) < 1) {
	 fprintf(stderr, "%s: unable to read label %d in file <%s>\n",
		 argv[0], nrows+1, argv[1]);
	 exit(EXIT_FAILURE);
      }
      if (icurr == ifirst) {
	 break;
      }
      nrows++;      
   }
   
   pd = (double*)malloc(sizeof(double)*(ncols*nrows*3));
   pdv = (double**)malloc(sizeof(double*)*(nrows*3));
   pi = (int*)malloc(sizeof(int)*nrows);
   
   if (pd == NULL || pdv == NULL) {
      fprintf(stderr, "%s: out of memory?\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   
   for (i = 3*nrows; i-- > 0; ) {
      pdv[i] = pd+i*ncols;
   }

   
   
   rewind(fin);

   for (i = 0; i < nrows; i++) {
      get_line(fin, buf, BUFSIZE);
      next = buf;
      s = get_word(next, &next);
      pi[i] = atoi(s);
      for (j = 0; j < ncols; j++) {
	 s = get_word(next, &next);
	 pdv[i1*nrows+i][j] = atof(s);
      }
   }
   
   for (i = 0; i < nrows; i++) {      
      get_line(fin, buf, BUFSIZE);
      next = buf;
      s = get_word(next, &next);
      printf("%8d", pi[i]);
      for (j = 0; j < ncols; j++) {
	 s = get_word(next, &next);
	 pdv[i2*nrows+i][j] = atof(s);
	 printf("%16.8e", (pdv[i2*nrows+i][j]-pdv[i1*nrows+i][j])/dt);
      }
      printf("\n");
   }
   
   while (1) {
      for (i = 0; i < nrows; i++) {	 
	 if (get_line(fin, buf, BUFSIZE) == 1) {
	    goto last_line;
	 }      
	 next = buf;
	 s = get_word(next, &next);
	 printf("%8d", pi[i]);
	 for (j = 0; j < ncols; j++) {
	    s = get_word(next, &next);
	    pdv[i3*nrows+i][j] = atof(s);
	    printf("%16.8e", (pdv[i3*nrows+i][j]-pdv[i1*nrows+i][j])/dt2);
	 }
	 printf("\n");
      }      
      i1 = (i1+1)%3;
      i2 = (i2+1)%3;
      i3 = (i3+1)%3;
   }
   
last_line:
  
   for (i = 0; i < nrows; i++) {
      printf("%8d", pi[i]);
      for (j = 0; j < ncols; j++) {	
	 printf("%16.8e", (pdv[i2*nrows+i][j]-pdv[i1*nrows+i][j])/dt);
      }
      printf("\n");
   }

   if (isf) {
      fclose(fin);
   }

   return 0;
}
   
