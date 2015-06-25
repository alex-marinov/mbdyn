/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#ifdef USE_HARWELL

#include <stdlib.h>
#include <stdio.h>
#include "ac/f2c.h"

#include "harwout.h"

/* ma28ad - begin */

int __FC_DECL__(w28ad9)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE N OUT OF RANGE = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28AD_99999, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE N OUT OF RANGE = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad8)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE NZ NON POSITIVE = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28AD_99998, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE NZ NON POSITIVE = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad7)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE LICN TOO SMALL = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28AD_99997, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE LICN TOO SMALL = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad6)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE LIRN TOO SMALL = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28AD_99996, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE LIRN TOO SMALL = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad5)(integer* lp)
{
#ifdef __cplusplus
   char* v[3];
   
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE INDICES FOUND OUT OF RANGE";
   v[1] = "\n";
   v[2] = NULL;
		   
   harwell_error(MA28AD_99995, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE INDICES FOUND OUT OF RANGE\n");
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad4)(integer* lp, 
	    integer* i1, doublereal* d, integer* i2, integer* i3)
{
#ifdef __cplusplus
   char* v[11];
   char buf_i1[256];
   char buf_d[256];
   char buf_i2[256];
   char buf_i3[256];
   
   sprintf(buf_i1, "%6d", *i1);
   sprintf(buf_d, "%12.4e", *d);
   sprintf(buf_i2, "%8d", *i2);
   sprintf(buf_i3, "%8d", *i3);
         
   v[0] = " ";
   v[1] = buf_i1;
   v[2] = "TH ELEMENT WITH VALUE ";
   v[3] = buf_d;
   v[4] = "\n";
   v[5] = " IS OUT OF RANGE WITH INDICES ";
   v[6] = buf_i2;
   v[7] = ", ";
   v[8] = buf_i3;
   v[9] = "\n";
   v[10] = NULL;
		   
   harwell_error(MA28AD_99994, v);
#else /* !__cplusplus */   
   printf(" %6dTH ELEMENT WITH VALUE %12.4e\n"
	  " IS OUT OF RANGE WITH INDICES %8d, %8d\n",
	  *i1, *d, *i2, *i3);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad3)(integer* lp, integer* i1, integer* i2, doublereal* d)
{
#ifdef __cplusplus
   char* v[8];
   char buf_i1[256];
   char buf_i2[256];
   char buf_d[256];  
   
   sprintf(buf_i1, "%8d", *i1);
   sprintf(buf_i2, "%8d", *i2);
   sprintf(buf_d, "%12.4e", *d);
         
   v[0] = " DUPLICATE ELEMENT IN POSITION ";
   v[1] = buf_i1;
   v[2] = ",";
   v[3] = buf_i1;
   v[4] = " WITH VALUE ";
   v[5] = buf_d;
   v[6] = "\n";
   v[7] = NULL;

   harwell_error(MA28AD_99993, v);
#else /* !__cplusplus */   
   printf(" DUPLICATE ELEMENT IN POSITION %8d,%8d WITH VALUE %12.4e\n",
	  *i1, *i2, *d);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad2)(integer* lp)
{
#ifdef __cplusplus
   char* v[3];
   
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE ERROR RETURN FROM MC23A/AD";
   v[1] = "\n";
   v[2] = NULL;
		   
   harwell_error(MA28AD_99992, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE ERROR RETURN FROM MC23A/AD\n");
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28ad1)(integer* lp)
{
#ifdef __cplusplus
   char* v[3];
   
   v[0] = " ERROR RETURN FROM MA28A/AD BECAUSE ERROR RETURN FROM MA30A/AD";
   v[1] = "\n";
   v[2] = NULL;
		   
   harwell_error(MA28AD_99991, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28A/AD BECAUSE ERROR RETURN FROM MA30A/AD\n");
#endif /* !__cplusplus */
   return 0;
}

/* ma28ad - end */


/* ma28bd - begin */

int __FC_DECL__(w28bd9)(integer* lp, integer* i1, integer* i2)
{
#ifdef __cplusplus
   char* v[7];
   char buf_i1[256];
   char buf_i2[256];
   
   sprintf(buf_i1, "%4d", *i1);
   sprintf(buf_i2, "%7d", *i2);
         
   v[0] = " ERROR RETURN FROM MA28B/BD WITH IFLAG=";
   v[1] = buf_i1;
   v[2] = "\n";
   v[3] = buf_i1;
   v[4] = " ENTRIES DROPPED FROM STRUCTURE BY MA28A/AD";
   v[5] = "\n";
   v[6] = NULL;

   harwell_error(MA28BD_99999, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD WITH IFLAG=%4d\n"
	  "%7d ENTRIES DROPPED FROM STRUCTURE BY MA28A/AD\n",
	  *i1, *i2);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28bd8)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28B/BD BECAUSE N OUT OF RANGE = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28BD_99998, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD BECAUSE N OUT OF RANGE = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28bd7)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28B/BD BECAUSE NZ NON POSITIVE = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28BD_99997, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD BECAUSE NZ NON POSITIVE = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28bd6)(integer* lp, integer* i)
{
#ifdef __cplusplus
   char* v[4];
   char buf_i[256];
   
   sprintf(buf_i, "%10d", *i);
   v[0] = " ERROR RETURN FROM MA28B/BD BECAUSE LICN TOO SMALL = ";
   v[1] = buf_i;
   v[2] = "\n";
   v[3] = NULL;
		   
   harwell_error(MA28BD_99996, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD BECAUSE LICN TOO SMALL = %10d\n", *i);
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28bd5)(integer* lp)
{
#ifdef __cplusplus
   char* v[3];
   
   v[0] = " ERROR RETURN FROM MA28B/BD BECAUSE ERROR RETURN FROM MA30B/BD";
   v[1] = "\n";
   v[2] = NULL;
		   
   harwell_error(MA28BD_99995, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD BECAUSE ERROR RETURN FROM MA30B/BD\n");
#endif /* !__cplusplus */
   return 0;
}

int __FC_DECL__(w28bd4)(integer* lp)
{
#ifdef __cplusplus
   char* v[3];
   
   v[0] = " ERROR RETURN FROM MA28B/BD";
   v[1] = "\n";
   v[2] = NULL;
		   
   harwell_error(MA28BD_99994, v);
#else /* !__cplusplus */   
   printf(" ERROR RETURN FROM MA28B/BD\n");
#endif /* !__cplusplus */
   return 0;
}

/* ma28bd - end */

#endif /* USE_HARWELL */

