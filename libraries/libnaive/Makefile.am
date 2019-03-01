include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libnaive.a library
noinst_LTLIBRARIES = libnaive.la

libnaive_la_SOURCES = \
mthrdslv.c \
mthrdslv.h \
pmthrdslv.c \
pmthrdslv.h

libnaive_la_LIBADD = @LIBS@ \
@ATOMIC_OPS_LIBS@
libnaive_la_LDFLAGS =

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
