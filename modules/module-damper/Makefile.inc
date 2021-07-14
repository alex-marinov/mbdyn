## $Header$
MODULE_DEPENDENCIES = nlrheo_damper.lo
# Needs Gnu Scientific Library
#MODULE_LINK = -lgsl
# Needs C blas; comment if C blas is already linked otherwise
MODULE_LINK = -lgsl -lgslcblas

# This is required to build both the module and the main
all: libmodule-damper.la nlrheo_main

# Rule for the main (FIXME: should export the CXXCOMPILE...)
nlrheo_main: nlrheo_main.cc nlrheo_damper.c $(MBDYN_MODULE_SRCDIR)/nlrheo_damper.h
	-$(CXX) -g -O0 -I$(MBDYN_MODULE_SRCDIR) -o $@ nlrheo_main.cc nlrheo_damper.c $(MODULE_LINK)
