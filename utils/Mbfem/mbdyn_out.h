#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <map>
#include <vector>

#include <myassert.h>
#include <matvec3.h>

#include "block.h"

std::map<int,int> read_mov(char*,std::map<int,int>,std::map<int,int>,std::vector<int>,std::map<int,Block>*,bool,double);
void read_frc(char*,std::map<int,int>,std::vector<int>,std::map<int,int>,std::map<int,Block>*,double);
