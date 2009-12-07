include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libaero.a library
noinst_LTLIBRARIES = libaero.la
libaero_la_SOURCES = \
aerod2.f \
aerod2.h \
aerodata.cc \
aerodata.h \
aerodata_impl.cc \
aerodata_impl.h \
aerodc81.c \
aerodc81.h \
aerodyn.cc \
aerodyn.h \
aerodyn_.h \
aeroelem.cc \
aeroelem.h \
aeroext.h \
aeroext.cc \
aeromodal.cc \
aeromodal.h \
c81data.cc \
c81data.h \
cyclocopter.cc \
cyclocopter.h \
genfm.cc \
genfm.h \
indvel.cc \
indvel.h \
instruments.cc \
instruments.h \
rotor.cc \
rotor.h \
windprof.cc \
windprof.h

libaero_la_LIBADD = @LIBS@
libaero_la_LDFLAGS = -static

INCLUDES = \
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
