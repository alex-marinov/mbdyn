include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libaero.a library
noinst_LTLIBRARIES = libaero.la
libaero_la_SOURCES = \
aerod2.f \
aerod2.h \
aerodata.cc \
aerodata.h \
aerodc81.c \
aerodc81.h \
aerodyn.cc \
aerodyn.h \
aeroelem.cc \
aeroelem.h \
aeromodal.cc \
aeromodal.h \
c81data.cc \
c81data.h \
rotor.cc \
rotor.h	\
aeroext.h \
aeroext.cc

libaero_la_LIBADD = @LIBS@
libaero_la_LDFLAGS = -static

noinst_PROGRAMS = testc81
testc81_SOURCES = testc81.cc
testc81_LDADD = libaero.la ../../libraries/libmbutil/libmbutil.la

INCLUDES = \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../libraries/libmbwrap \
-I../../mbdyn \
-I../../mbdyn/base \
-I../../mbdyn/struct \
-I../../mbdyn/aero \
-I../../mbdyn/elec \
-I../../mbdyn/hydr

include $(top_srcdir)/build/bot.mk
