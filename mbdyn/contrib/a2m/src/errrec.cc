/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


#include <errrec.h>
#include <errors.h>

void error_recovery (int ncount, FILE* logfile)
{
   /* Nel caso una card non venga riconosciuta del tutto */
   fprintf (logfile,"<UNRECOGNIZED CARD BEHIND>\n");
}


void out_error (int err_num, const char* message)
{
   fprintf (errfile,"\nLINE %5d :", ncount);
   fprintf (errfile,"%s [%s]",Error_code[err_num],message);
   nerr++;
   return;
}

void out_warning (const char* message)
{
   fprintf (errfile,"\nLINE %5d :", ncount-1);
   fprintf (errfile,"WARNING ! %s",message);
   nwarnings++;
   return;
}

void out_warning (int err_num, const char* message)
{
   fprintf (errfile,"\nLINE %5d :[WARNING]:", ncount);
   fprintf (errfile,"%s [%s]",Error_code[err_num],message);
   nwarnings++;
   return;
}

void out_table (const char* message)
{
   fprintf (tablefile,message);
   return;
}

int yyerror(char *s)
{
   /* Questo blocco scrive nel file di log del lexer il token che ha
    * generato l'errore nella sua esatta posizione */
   fprintf (logfile,"(BAD TOKEN)"); 	
   out_error (2,yytext);
   return -1;
}
