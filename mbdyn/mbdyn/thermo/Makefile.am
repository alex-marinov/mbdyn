include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libelec.a library
noinst_LTLIBRARIES = libthermo.la
libthermo_la_SOURCES = \
therm.cc \
therm.h \
thermalcapacitance.h \
thermalcapacitance.cc \
thermalnode.cc \
thermalnode.h \
thermalresitance.h \
thermalresistance.cc \
thermalsource.h \
thermalsource.cc

libthermo_la_LIBADD = @LIBS@
libthermo_la_LDFLAGS = -static

INCLUDES = \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../libraries/libmbwrap \
-I../../libraries/libcppad \
-I../../mbdyn \
-I../../mbdyn/base \
-I../../mbdyn/struct \
-I../../mbdyn/aero \
-I../../mbdyn/elec \
-I../../mbdyn/thermo \
-I../../mbdyn/hydr

include $(top_srcdir)/build/bot.mk
