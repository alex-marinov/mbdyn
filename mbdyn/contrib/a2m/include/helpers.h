#ifndef HELPERS_H
#define HELPERS_H

#include <defs.h>
#include <ctype.h>
#include <stdio.h>

void strupper (char *s);
//Boolean CheckInputFile (char*,FILE*&,FILE*&);
char *CheckInputFile (int,char* argv[],FILE*&,FILE*&,FILE*&);
#endif

