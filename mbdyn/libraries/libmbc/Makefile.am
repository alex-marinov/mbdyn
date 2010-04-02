include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libmbc.a library
noinst_LTLIBRARIES = libmbc.la
#lib_LTLIBRARIES = libmbc.la

libmbc_la_SOURCES = \
mbc.c \
mbc.h \
sock.c \
sock.h

libmbc_la_LIBADD = @LIBS@
libmbc_la_LDFLAGS =

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
