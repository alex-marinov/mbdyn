# $Header$
###############################################################################
#
# Add module dependencies here; e.g., if in dir 'module-X' there's 
# a file 'myfunc.cc' write:
#
# MODULE_DEPENDENCIES=myfunc.o
#
# Then put this file in directory 'module-X' as 'Makefile.inc'
#
###############################################################################

# FIXME: this flag is needed by g++ 4.8 to include <cstdint>
MODULE_INCLUDE="-std=c++0x"

