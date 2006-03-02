include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libcolamd.a library
lib_LTLIBRARIES = libcolamd.la

libcolamd_la_SOURCES = \
colamd.c \
colamd.h

libcolamd_la_LIBADD = @LIBS@
libcolamd_la_LDFLAGS =

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
