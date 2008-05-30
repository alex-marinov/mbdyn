## $Header$
MODULE_DEPENDENCIES = nlrheo_damper.o
# Needs Gnu Scientific Library
#MODULE_LINK = -lgsl
# Needs C blas; comment if C blas is already linked otherwise
MODULE_LINK = -lgsl -lgslcblas

nlrheo_main: nlrheo_main.cc nlrheo_damper.cc $(MBDYN_MODULE_SRCDIR)/nlrheo_damper.h
	-$(CXX) -g -O0 -I$(MBDYN_MODULE_SRCDIR) -o $@ nlrheo_main.cc nlrheo_damper.cc $(MODULE_LINK)
