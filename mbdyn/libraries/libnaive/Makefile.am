include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libcolamd.a library
lib_LTLIBRARIES = libnaive.la

libnaive_la_SOURCES = \
mthrdslv.c \
mthrdslv.h \
pmthrdslv.c

libnaive_la_LIBADD = @LIBS@
libnaive_la_LDFLAGS =

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
