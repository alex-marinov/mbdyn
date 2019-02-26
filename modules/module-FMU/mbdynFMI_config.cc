/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
        Copyright (C) 2016-2017 all rights reserved.
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

	
