## Process this file with automake to produce Makefile.in

# Build libhydr.a library
lib_LTLIBRARIES = libhydr.la
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

