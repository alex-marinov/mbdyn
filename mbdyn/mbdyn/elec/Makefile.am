include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libelec.a library
noinst_LTLIBRARIES = libelec.la
libelec_la_SOURCES = \
discctrl.cc \
discctrl.h \
elec.cc \
elec.h \
elecnode.cc \
elecnode.h \
forgfact.cc \
forgfact.h \
genel.cc \
genel.h \
genel_.h \
genfilt.cc \
genfilt.h \
gpc.cc \
gpc.h \
id.cc \
id.h \
px.cc \
px.h \
rottrim.cc \
rottrim.h \
swashpl.cc \
swashpl.h

libelec_la_LIBADD = @LIBS@
libelec_la_LDFLAGS = -static

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

# Program that tests Deadbeat Predictive Control
noinst_PROGRAMS = dpc
dpc_SOURCES =
dpc_LDADD = dpc.o \
../base/libbase.la \
../../libraries/libmbutil/libmbutil.la \
../../libraries/libmbmath/libmbmath.la \
@MESCHACH_LIBS@ @LAPACK_LIBS@ @F2C_LIBS@ @LIBS@

# Explicit compile rules to enforce the use of CXX!
dpc.o : gpc.cc
	-$(CXXCOMPILE) $(CXXFLAGS) -DTEST_DPC -c gpc.cc && mv gpc.o dpc.o

dpc : dpc.o
	-$(CXXLINK) -o dpc $(dpc_LDADD)

include $(top_srcdir)/build/bot.mk
