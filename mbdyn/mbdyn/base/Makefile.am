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
demangled.h \
dofdrive.cc \
dofdrive.h \
dofman.cc \
dofown.h \
dofpgin.cc \
dofpgin.h \
drive.cc \
drive.h \
driven.h \
drive_.cc \
drive_.h \
elem.cc \
elem.h \
elman.cc \
enums.cc \
env.cc \
filedrv.cc \
filedrv.h \
force.cc \
force.h \
integr.h \
intmeth.cc \
intmeth.h \
j2p.cc \
j2p.h \
loadable.cc \
loadable.h \
mangled.h \
mbpar.cc \
mbpar.h \
multistp.cc \
multistp.h \
node.cc \
node.h \
nodeman.cc \
output.cc \
output.h \
posrel.h \
readclaw.h \
reffrm.h \
resforces.cc \
resforces.h \
restart.h \
schurdataman.cc \
schurdataman.h \
shape.h \
shockabsorber.h \
sockdrv.cc \
sockdrv.h \
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
