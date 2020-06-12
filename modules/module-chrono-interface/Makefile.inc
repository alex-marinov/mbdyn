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

#where I generate the shared lib including C_E models, and some funcs. 			
LIBPATH=/home/runsen/Documents/GSoC/C::E_models

#the header for the funcs “mbdyn_ce.h”.
MODULE_INCLUDE = -I .
#link the shared library “libCE_models.so”
MODULE_LINK= -L$(LIBPATH) -lCE_models
#MODULE_LINK = -L$(LIBPATH) -lCE_models -Wl, -rpath=$(LIBPATH) -lCE_models 
#MODULE_LINK = -lCE_models


