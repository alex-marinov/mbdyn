include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libbase.a library
noinst_LTLIBRARIES = libbase.la

libbase_la_SOURCES = \
adams_res.cc \
auth.cc \
auth.h \
bicg.cc \
bicg.h \
bulk.cc \
bulk.h \
chaco_interface.cc \
chaco_interface.h \
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
fixedstep.cc \
fixedstep.h \
force.cc \
force.h \
gmres.cc \
gmres.h \
integr.cc \
integr.h \
j2p.cc \
j2p.h \
loadable.cc \
loadable.h \
mbpar.cc \
mbpar.h \
mbrtai_utils.c \
mbrtai_utils.h \
metiswrap.cc \
metiswrap.h \
mfree.cc \
mfree.h \
motionview_res.cc \
node.cc \
node.h \
nodeman.cc \
nonlin.cc \
nonlin.h \
nonlinpb.h \
nr.cc \
nr.h \
output.cc \
output.h \
posrel.h \
precond.cc \
precond.h \
precond_.h \
readclaw.h \
reffrm.cc \
reffrm.h \
resforces.cc \
resforces.h \
restart.h \
rtai_in_drive.cc \
rtai_in_drive.h \
rtai_out_elem.cc \
rtai_out_elem.h \
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
socketstreamdrive.cc \
socketstreamdrive.h \
streamdrive.cc \
streamdrive.h \
solver.cc \
solver.h \
solverdiagnostics.cc \
solverdiagnostics.h \
stepsol.cc \
stepsol.h \
stepsol.hc \
symcltp.h \
tpldrive.cc \
tpldrive.h \
tpldrive_.h \
vec3drv.h \
thirdorderstepsol.h \
thirdorderstepsol.cc \
ScalarFunctions.h \
ScalarFunctionsImpl.h \
ScalarFunctionsImpl.cc

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
