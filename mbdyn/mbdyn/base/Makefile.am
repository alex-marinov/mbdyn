## Process this file with automake to produce Makefile.in

# Build libbase.a library
lib_LTLIBRARIES = libbase.la

libbase_la_SOURCES = \
dataman5.cc  elman.cc     loadable.cc  reffrm.h \
auth.cc      dataman_.h   loadable.h   restart.h \
auth.h       demangled.h  mangled.h \
bool.h       dofdrive.cc  enums.cc     mbpar.cc \
bulk.cc      dofdrive.h   mbpar.h \
bulk.h       dofman.cc    memmans.h \
constltp.h   dofown.h     env.cc       multistp.cc \
constltp_.h  dofpgin.cc   filedrv.cc   multistp.h \
contactj.h   dofpgin.h    filedrv.h    mysleep.c    mysleep.h    shape.h \
damper.f     drive.cc     force.cc     node.cc      sockdrv.cc \
damper.h     drive.h      force.h      node.h       sockdrv.h \
dataman.cc   drive_.cc    integr.h     nodeman.cc   tpldrive.cc \
dataman.h    drive_.h     intmeth.cc   output.cc    tpldrive.h   tpldrive_.h\
dataman2.cc  driven.h     intmeth.h    output.h     vec3drv.h \
dataman3.cc  elem.cc      j2p.cc       posrel.h \
dataman4.cc  elem.h       j2p.h        readclaw.h \
schur.cc     schur.h      schurdataman.cc           schurdataman.h \
shockabsorber.h

libbase_la_LIBADD = @LIBS@
libbase_la_LDFLAGS = -static

INCLUDES = -I../../libraries/libmbutil -I../../libraries/libmbmath -I../../libraries/libobjs -I../../mbdyn -I../../mbdyn/base -I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec -I../../mbdyn/hydr

