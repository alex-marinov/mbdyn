
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
