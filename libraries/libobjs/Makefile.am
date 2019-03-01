include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libobjs.la

libobjs_la_SOURCES = \
get_nprocs.c \
getopt.c

libobjs_la_LIBADD = @LIBS@
libobjs_la_LDFLAGS =

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath

include $(top_srcdir)/build/bot.mk
