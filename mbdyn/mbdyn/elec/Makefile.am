include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libelec.a library
noinst_LTLIBRARIES = libelec.la
libelec_la_SOURCES = \
accelerometer.cc \
accelerometer.h \
displacement.cc \
displacement.h \
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
motor.cc \
motor.h \
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

include $(top_srcdir)/build/bot.mk
