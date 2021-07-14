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



// HELPERS.CC - FUNZIONI DI VARIA UTILITA'

#include <stdlib.h>
#include <helpers.h>
#include <string.h>
#include <getopt.h>

extern Boolean VERBOSE_MODE;
extern Boolean OVERWRITE_MODE;
extern Boolean THROW_MODE;
extern Boolean SKIP_MODE;
extern Boolean EXTENDED_MATRIX_DISPLAY;
extern Boolean SPECIFY_FILE;
extern Boolean REMOVE_REMARK;

void strupper (char* s)
{
   for (int i=0;i<strlen(s);i++)
     s[i]=toupper(s[i]);
}


char * CheckInputFile (int argn, char *argv[], FILE*& erf, FILE*& msg, 
		       FILE*& tbl)
{
	FILE *test = NULL;
        const char* fname=NULL;
   
	while (1) {
		int opt;

		opt = getopt(argn, argv, "voksefhr");
		
		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'v' : 
			VERBOSE_MODE=Y; 
			break;
		case 'o' : 
			OVERWRITE_MODE=Y; 
			break;
		case 'k' : 
			THROW_MODE=N; 
			break;
		case 's' : 
			SKIP_MODE=Y; 
			break;
	        case 'r' :
		        REMOVE_REMARK=Y;
		        break;
		case 'e' : 
			EXTENDED_MATRIX_DISPLAY=Y; 
			break;
		case 'f' : 
			SPECIFY_FILE=Y; 
		        fname=argv[optind];
		        break;
		case 'h' :
			/* stampa l'help di sistema */
			cout << "Options for a++:" << endl
	<< "-v  : verbose mode (display all adams card with all entries)"
	<< endl << "-o  : overwrite mode (overwrites card with same label)"
	<< endl << "-k  : throw mode (erase card with syntax form error)"
	<< endl << "-e  : extended matrix display"
	<< endl << "-f  : specify file (prompt for user output filename)"
	<< endl << "-r  : remove remarks from the output file"
	<< endl;
			exit(EXIT_SUCCESS);

		}
	}

	if (argv[optind] == NULL) {
		return NULL;
	}
        int l = strlen(argv[optind]);
	char *base_input = new char[l+1];
	strcpy(base_input, argv[optind]);
	
	if (l > 4 && strcasecmp(base_input+l-4, ".adm") == 0) {
		base_input[l-4] = '\0';
	}
	char *base_error = new char[l+4+1];
	char *base_messg = new char[l+4+1];
        char *base_table = new char[l+4+1];
        strcpy(base_error, base_input);
	strcat(base_error, ".err");
	strcpy(base_messg, base_input);
	strcat(base_messg, ".msg");
        strcpy(base_table, base_input);
        strcat(base_table,".ref");

	erf=fopen(base_error,"w");
	msg=fopen(base_messg,"w");
        tbl=fopen(base_table,"w");
	
        cout << "Base input is:" << base_input << endl;
        return base_input;
}
