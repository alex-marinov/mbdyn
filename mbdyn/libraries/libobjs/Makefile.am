include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libobjs.la

libobjs_la_SOURCES = \
dgegv.c \
getopt.c \
getopt.h

libobjs_la_LIBADD = @LIBS@
libobjs_la_LDFLAGS =

INCLUDES = -I../../mbdyn \
	-I../../libraries/libmbutil \
	-I../../libraries/libmbmath

include $(top_srcdir)/build/bot.mk
