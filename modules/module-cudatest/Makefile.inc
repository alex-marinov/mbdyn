## $Header$
MODULE_DEPENDENCIES=cudatest.lo
MODULE_LINK= -L/usr/local/cuda/lib -lcudart -lcuda
#NVCC=nvcc
#
#cudatest.lo:
#	ln -sf $(srcdir)/module-cudatest/cudatest.cu .
#	../libtool --mode=compile $(NVCC) cudatest.cu

