## Process this file with automake to produce Makefile.in

# Build libelec.a library
lib_LTLIBRARIES = libelec.la
libelec_la_SOURCES = \
discctrl.cc  discctrl.h  elec.cc    elec.h   elecnode.cc \
elecnode.h   forgfact.cc  forgfact.h  genel.cc   genel.h  genel_.h \
genfilt.cc   genfilt.h    gpc.cc      gpc.h      id.cc    id.h \
px.cc        px.h         swashpl.cc  swashpl.h

libelec_la_LIBADD = @LIBS@
libelec_la_LDFLAGS = -static

INCLUDES = -I../../include -I../../libraries/libmbutil -I../../libraries/libmbmath -I../../mbdyn -I../../mbdyn/base -I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec -I../../mbdyn/hydr

