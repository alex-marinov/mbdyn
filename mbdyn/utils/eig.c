#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <malloc.h>
#include <myf2c.h>

extern "C" int 
dgegv_(char* jobvl,
       char* jobvr, 
       integer* n, 
       doublereal* a,
       integer* lda, 
       doublereal* b, 
       integer* ldb, 
       doublereal* alphar, 
       doublereal* alphai, 
       doublereal* beta, 
       doublereal* vl, 
       integer* ldvl, 
       doublereal* vr, 
       integer* ldvr, 
       doublereal* work, 
       integer* lwork, 
       integer* info);

static int 
read_matrix_header(char* buf, integer* nrows, integer* ncols)
{
   char* p = NULL;
   char* q = NULL;
   
   p = strchr(buf, ':');
   if (p == NULL) {
      return -1;
   }
   p++;
   
   q = strchr(p, 'x');      
   if (q == NULL) {
      return -1;
   }
   *q = '\0';
   *ncols = atoi(q+1);  
   *nrows = atoi(p);
   
   return 0;
}

int
get_matrix(FILE* in, integer nrows, integer ncols, doublereal* m)
{
   char buf[1024];
   int i;
   int j;
   
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
	 char* p = buf;
	 
	 while (isspace(*p = fgetc(in))) {
	 }
	 p++;
	 while (!isspace(*p = fgetc(in))) {
	    p++;
	 }
	 *p = '\0';
	 
	 m[i+nrows*j] = strtod(buf, NULL);
      }
   }
   
   return 0;
}

int
main(int argc, char* argv[])
{
   FILE* in = NULL;
   char buf[1024];
   int i;
   int j;
   
   FILE* out = NULL;
   
   integer nrows = 0;
   integer ncols = 0;
   
   doublereal* pd = NULL;
   doublereal* a = NULL;
   doublereal* b = NULL;
   doublereal* alphar = NULL;
   doublereal* alphai = NULL;
   doublereal* beta = NULL;
   doublereal* vl = NULL;
   doublereal* vr = NULL;
   doublereal* work = NULL;   
   
   integer n = 0;
   integer lwork = 0;
   integer info = 0;
   
   char sl[4] = "V";
   char sr[4] = "V";
      
   int rc = 0;
   
   if (argc < 3) {
      fprintf(stderr, "eig: missing file name\ncommand: %s", argv[0]);
      for (i = 1; i < argc; i++) {
	 fprintf(stderr, " %s", argv[i]);
      }
      fprintf(stderr, "\n");
      return -1;
   }
   
   in = fopen(argv[1], "r");
   if (in == NULL) {
      fprintf(stderr, "eig: unable to open file '%s'\n", argv[1]);
      return -1;
   }

   /* matrix A */
   do {
      if (fgets(buf, sizeof(buf), in) == NULL) {
	 fprintf(stderr, "eig: unable to read matrix A header\n");
	 return -1;
      }
   } while (buf[0] != '#');
   
   if (read_matrix_header(buf, &nrows, &ncols) != 0) {
      fprintf(stderr, "eig: corrupted header for matrix A\n");    
   }
   
   if (nrows <= 0 || ncols <= 0 || nrows != ncols) {
      fprintf(stderr, "eig: illegal size %dx%d for matrix A\n", nrows, ncols);
      return -1;
   }
   
   n = nrows;
   lwork = 8*n;
   
   pd = (doublereal*)malloc(sizeof(doublereal)*((4*n+3)*n+lwork));
   if (pd == NULL) {
      fprintf(stderr, "eig: memory exausted?\n");
      return -1;
   }
   
   a = pd;
   b = a+n*n;
   vl = b+n*n;
   vr = vl+n*n;
   alphar = vr+n*n;
   alphai = alphar+n;
   beta = alphai+n;
   work = beta+n;

   get_matrix(in, n, n, a);
   
   /* matrix B */
   do {
      if (fgets(buf, sizeof(buf), in) == NULL) {
	 fprintf(stderr, "eig: unable to read matrix B header\n");
	 return -1;
      }
   } while (buf[0] != '#');
   
   if (read_matrix_header(buf, &nrows, &ncols) != 0) {
      fprintf(stderr, "eig: corrupted header for matrix B\n");
   }
   
   if (nrows <= 0 || ncols <= 0 || nrows != ncols || nrows != n) {
      fprintf(stderr, "eig: illegal size %dx%d for matrix B\n", nrows, ncols);
      return -1;      
   }

   get_matrix(in, n, n, b);
   
   dgegv_(sl,
	  sr, 
	  &n, 
	  a,
	  &n,
	  b, 
	  &n,
	  alphar, 
	  alphai, 
	  beta, 
	  vl, 
	  &n, 
	  vr, 
	  &n, 
	  work, 
	  &lwork, 
	  &info);
   
   if (info != 0) {
      fprintf(stderr, "error return from dgegv\n");
      rc = -1;
      goto error_return;
   }
   
   out = fopen(argv[2], "w");
   if (out == NULL) {
      fprintf(stderr, "eig: unable to open file '%s'\n", argv[2]);
      return -1;
   }

   fprintf(out, "# vector AlphaR: %d\n", n);
   for (i = 0; i < n; i++) {
      fprintf(out, "%20.8e\n", alphar[i]);
   }
   
   fprintf(out, "# vector AlphaI: %d\n", n);
   for (i = 0; i < n; i++) {
      fprintf(out, "%20.8e\n", alphai[i]);
   }
   
   fprintf(out, "# vector Beta: %d\n", n);
   for (i = 0; i < n; i++) {
      fprintf(out, "%20.8e\n", beta[i]);
   }
   
   fprintf(out, "# matrix VL: %dx%d\n", n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	 fprintf(out, "%20.8e", vl[i+n*j]);
      }
      fprintf(out, "\n");
   }
   
   fprintf(out, "# matrix VR: %dx%d\n", n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	 fprintf(out, "%20.8e", vr[i+n*j]);
      }
      fprintf(out, "\n");
   }
   
error_return:
   if (pd != NULL) {
      free(pd);
   }
   
   return rc;
}
