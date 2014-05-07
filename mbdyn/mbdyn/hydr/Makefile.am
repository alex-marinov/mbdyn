include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libhydr.a library
noinst_LTLIBRARIES = libhydr.la
libhydr_la_SOURCES = \
actuator.cc \
actuator.h \
hfluid.cc \
hfluid.h \
hfluid_.h \
hminor.cc \
hminor.h \
hutils.cc \
hutils.h \
pipe.cc \
pipe.h \
preselem.cc \
preselem.h \
presnode.h \
valve.cc \
valve.h

libhydr_la_LIBADD = @LIBS@
libhydr_la_LDFLAGS = -static

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base \
-I$(srcdir)/../../mbdyn/struct \
-I$(srcdir)/../../mbdyn/aero \
-I$(srcdir)/../../mbdyn/elec \
-I$(srcdir)/../../mbdyn/hydr

include $(top_srcdir)/build/bot.mk
