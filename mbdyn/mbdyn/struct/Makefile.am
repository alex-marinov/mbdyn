## Process this file with automake to produce Makefile.in

# Build libstruct.a library
lib_LTLIBRARIES = libstruct.la
libstruct_la_SOURCES = \
accj.cc      accj.h      autostr.cc   autostr.h \
beam.cc      beam.h       body.cc      body.h      drvhinge.cc  drvhinge.h \
drvj.cc      drvj.h       genj.cc      genj.h      gravity.cc   gravity.h \
inline.cc    inline.h     inplanej.cc  inplanej.h  inplcntj.h   joint.cc \
joint.h      joint_.h     kin.cc      kin.h        kinj.cc \
kinj.h       modal.cc     modal.h      planedj.cc  planedj.h    planej.cc \
planej.h     prismj.cc    prismj.h     pzbeam.cc   pzbeam.h     rodj.cc \
rodj.h       spherj.cc    spherj.h     strnode.cc  strnode.h    univj.cc \
univj.h      vehj.cc      vehj.h       vehj2.cc    vehj2.h \
beam2.h      beam2.cc     pzbeam2.h    pzbeam2.cc  shapefnc.h   shapefnc.cc \
hbeam.h      hbeam.cc     hbeam_interp.cc

libstruct_la_LIBADD = @LIBS@
libstruct_la_LDFLAGS = -static

INCLUDES = -I../../include -I../../libraries/libmbutil -I../../libraries/libmbmath -I../../mbdyn -I../../mbdyn/base -I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec -I../../mbdyn/hydr

