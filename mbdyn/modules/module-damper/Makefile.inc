## $Header$
MODULE_DEPENDENCIES = nlrheo_damper.o
# Needs Gnu Scientific Library
#MODULE_LINK = -lgsl
# Needs C blas; comment if C blas is already linked otherwise
MODULE_LINK = -lgsl -lgslcblas
