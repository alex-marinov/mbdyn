include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libbase.a library
noinst_LTLIBRARIES = libbase.la

libbase_la_SOURCES = \
auth.cc \
auth.h \
bulk.cc \
bulk.h \
constltp.h \
constltp_.h \
contactj.h \
damper.f \
damper.h \
dataman.cc \
dataman.h \
dataman2.cc \
dataman3.cc \
dataman4.cc \
dataman5.cc \
dataman_.h \
dofdrive.cc \
dofdrive.h \
dofman.cc \
dofown.cc \
dofown.h \
dofpgin.cc \
dofpgin.h \
drive.cc \
drive.h \
driven.cc \
driven.h \
drive_.cc \
drive_.h \
elem.cc \
elem.h \
elman.cc \
enums.cc \
env.cc \
external.h \
external.cc \
filedrv.cc \
filedrv.h \
force.cc \
force.h \
integr.cc \
integr.h \
j2p.cc \
j2p.h \
loadable.cc \
loadable.h \
mbpar.cc \
mbpar.h \
node.cc \
node.h \
nodeman.cc \
nonlin.cc \
nonlin.h \
nonlinpb.h \
output.cc \
output.h \
posrel.h \
readclaw.h \
reffrm.h \
resforces.cc \
resforces.h \
restart.h \
sah.cc \
sah.h \
schurdataman.cc \
schurdataman.h \
shape.h \
shockabsorber.h \
simentity.cc \
simentity.h \
sockdrv.cc \
sockdrv.h \
solver.cc \
solver.h \
stepsol.cc \
stepsol.h \
symcltp.h \
tpldrive.cc \
tpldrive.h \
tpldrive_.h \
vec3drv.h

libbase_la_LIBADD = @LIBS@ $(GINACLIB_LIBS)
libbase_la_LDFLAGS = -static

INCLUDES = \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../libraries/libmbwrap \
-I../../mbdyn \
-I../../mbdyn/base \
-I../../mbdyn/struct \
-I../../mbdyn/aero \
-I../../mbdyn/elec \
-I../../mbdyn/hydr \
$(GINACLIB_CPPFLAGS)

include $(top_srcdir)/build/bot.mk
