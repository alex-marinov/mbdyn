
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
	const char *base = argv[optind];
	char *base_error = new char[strlen(base)+4+1];
	char *base_messg = new char[strlen(base)+4+1];
	char *base_input = new char[strlen(base)+1];
        char *base_table = new char[strlen(base)+4+1];
        strcpy(base_error, base);
	strcat(base_error, ".err");
	strcpy(base_messg, base);
	strcat(base_messg, ".msg");
        strcpy(base_table, base);
        strcat(base_table,".ref");
	strcpy(base_input, base);

	erf=fopen(base_error,"w");
	msg=fopen(base_messg,"w");
        tbl=fopen(base_table,"w");
	
        cout << "Base input is:" << base_input << endl;
        return base_input;
}
