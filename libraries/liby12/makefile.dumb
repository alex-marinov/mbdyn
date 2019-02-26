# dumb makefile to compile the bare library

# replace with your fortran compiler and ar/ranlib
F77=g77
AR=ar
RANLIB=ranlib

# do not change below this line

SRCS=y12mae.f y12maf.f y12mbe.f y12mbf.f y12mce.f y12mcf.f y12mde.f \
y12mdf.f y12mfe.f y12mge.f y12mhe.f
OBJS=y12mae.o y12maf.o y12mbe.o y12mbf.o y12mce.o y12mcf.o y12mde.o \
y12mdf.o y12mfe.o y12mge.o y12mhe.o
LIB=liby12.a

$(LIB): $(SRCS)
	$(FC) -c $(SRCS)
	$(AR) ru $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

