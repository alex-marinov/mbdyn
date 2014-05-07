include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libcolamd.a library
noinst_LTLIBRARIES = libmbdyncolamd.la

libmbdyncolamd_la_SOURCES = \
colamd.c \
colamd.h

libmbdyncolamd_la_LIBADD = @LIBS@
libmbdyncolamd_la_LDFLAGS =

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
