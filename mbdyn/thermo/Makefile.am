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
thermalresistance.h \
thermalresistance.cc \
thermalsource.h \
thermalsource.cc

libthermo_la_LIBADD = @LIBS@
libthermo_la_LDFLAGS = -static

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../libraries/libcppad \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base \
-I$(srcdir)/../../mbdyn/struct \
-I$(srcdir)/../../mbdyn/aero \
-I$(srcdir)/../../mbdyn/elec \
-I$(srcdir)/../../mbdyn/hydr \
-I$(srcdir)/../../mbdyn/thermo

include $(top_srcdir)/build/bot.mk
