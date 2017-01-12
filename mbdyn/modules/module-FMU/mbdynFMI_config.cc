/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
        Copyright (C) 2016(-2017) all rights reserved.
        The copyright of this patch is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described 
        in the GNU Public License version 2.1

*/

#include "mbdynFMI_config.h"
#include <stdio.h>
#include <string.h>

void importlogger(jm_callbacks* c, jm_string module, jm_log_level_enu_t log_level, jm_string message)
{

       printf("module = %s, log level = %d: %s\n", module, log_level, message);

}

void setup_callbacks(jm_callbacks* callbacks){
	callbacks->malloc = malloc;
	callbacks->calloc = calloc;
	callbacks->realloc = realloc;
	callbacks->free = free;
	callbacks->logger = importlogger;
	callbacks->log_level = jm_log_level_debug;
	callbacks->context = 0;
	printf("Callback Setup Done! \n");
}

std::string UncompressLocation(const char* location){
        int length = strlen(location);
        int i;

        for (i=length; i>0; i--){
                if(location[i]==47){
                        break;
                }
        }

	std::string destination(location);
//        char* destination = (char*) malloc(length*sizeof(char) );
	destination.resize(i+1);
        return destination.c_str();
}

	
