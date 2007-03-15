include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libann.a library
lib_LTLIBRARIES = libann.la

libann_la_SOURCES = \
ActivationFunction.c \
ActivationFunction.h \
ann.c \
ann.h \
matrix.c \
matrix.h \
scale_factor.c

libann_la_LIBADD = @LIBS@
libann_la_LDFLAGS =

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
