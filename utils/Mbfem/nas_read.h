#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <map>
#include <vector>
#include <math.h>

#include <myassert.h>
#include <matvec3.h>
#include <fullmh.h>

#include "block.h"

int readBgpdt(FILE*, int*, int, char*,std::map<int,int>,std::map<int,int>
		,std::map<int,int>,std::map<int,int>,std::map<int,Block>*,std::vector<int>*);
int readop2(char*,std::map<int,int>,std::map<int,int>,
		std::map<int,int>,std::map<int,int>,std::map<int,Block>*,std::vector<int>*);
int iopen(FILE*, char*);
int iheadr(FILE*, char*, int*);
int inasread(FILE*, int*, int, int, int*);
int namecomp(char*, char**, int);
FILE* open_nasout(char*,int,std::vector<int>);
void nas_cards(FILE*,std::map<int,Block>,int,int,int);
void text_dump(char*,int,int,std::map<int,Block>);
void readmat(char*,std::vector<int>,std::map<int,Block>*);
void rdir_dump(char*,int,int,std::map<int,Block>);
