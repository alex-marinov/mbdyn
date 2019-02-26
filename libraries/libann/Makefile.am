include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libmbann.a library
noinst_LTLIBRARIES = libmbann.la

libmbann_la_SOURCES = \
ActivationFunction.c \
ActivationFunction.h \
ann.c \
ann.h \
matrix.c \
matrix.h

libmbann_la_LIBADD = @LIBS@
libmbann_la_LDFLAGS =

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
