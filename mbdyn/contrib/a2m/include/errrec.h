#ifndef ERRREC_H
#define ERRREC_H

// MODULO DI GESTIONE DEGLI ERRORI
// ERR(OR)REC(OGNIZE) 

#include <stdio.h>
#include <defs.h>

// Variabili relative ai file
extern FILE* errfile;
extern FILE* logfile;
extern FILE* messagefile;
extern FILE* tablefile;
//// Variabili relative a lex
extern int ncount;
extern const char* yytext;
// Variabili relative al bison
extern const char** Token_Table;
extern inline const char* Find_Token(int);
// Variabili relative alla gestione errori
extern int nerr;
extern int nwarnings;
extern Boolean ES;

void error_recovery (int, FILE*);
void out_error (int,const char*);
int  yyerror (char*);
void out_warning (const char*);
void out_warning (int, const char*);
void out_table (const char*);

#endif
