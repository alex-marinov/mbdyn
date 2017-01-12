/*
        AUTHOR: Devyesh Tandon <devyeshtandon+mbdyn@gmail.com>
        Copyright (C) 2016(-2017) all rights reserved.
        The copyright of this patch is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described 
        in the GNU Public License version 2.1

*/

#include "mbdynFMI.h"

#define BUFFER 1000

void importlogger(jm_callbacks* c, jm_string module, jm_log_level_enu_t log_level, jm_string message);

void setup_callbacks(jm_callbacks* callbacks);

int error(const char* test, int value);

int sumde(int x, int y);

std::string UncompressLocation(const char* location);
