include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libbase.a library
noinst_LTLIBRARIES = libbase.la

libbase_la_SOURCES = \
auth.cc \
auth.h \
bicg.cc \
bicg.h \
binaryoutput.cc \
binary_output_file.cc \
bufferstream_out_elem.cc \
bufferstream_out_elem.h \
bufferstreamdrive.cc \
bufferstreamdrive.h \
bulk.cc \
bulk.h \
constltp.h \
constltp_ann.h \
constltp_axw.h \
constltp_dcw.h \
constltp_impl.cc \
constltp_impl.h \
constltp_nlp.cc \
constltp_nlp.h \
constltp_nlsf.cc \
constltp_nlsf.h \
contactj.h \
converged.cc \
converged.h \
dataman.cc \
dataman.h \
dataman2.cc \
dataman3.cc \
dataman4.cc \
dataman6.cc \
dataman_.h \
datamanforward.h \
ddrive.cc \
ddrive.h \
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
extedge.cc \
extedge.h \
external.cc \
external.h \
extforce.cc \
extforce.h \
extsocket.cc \
extsocket.h \
filedrv.cc \
filedrv.h \
fixedstep.cc \
fixedstep.h \
force.cc \
force.h \
gmres.cc \
gmres.h \
hint.h \
hint_impl.cc \
hint_impl.h \
invdataman.cc \
invdyn.h \
invdyn.cc \
invsolver.cc \
invsolver.h \
invsolwrap.cc \
invsolwrap.h \
j2p.cc \
j2p.h \
linesearch.h \
linesearch.cc \
loadable.cc \
loadable.h \
mbpar.cc \
mbpar.h \
mfree.cc \
mfree.h \
modelns.cc \
modelns.h \
modules.cc \
modules.h \
mtdataman.cc \
mtdataman.h \
nestedelem.cc \
nestedelem.h \
netcdf_output_file.cc \
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
privdrive.cc \
privdrive.h \
privpgin.cc \
privpgin.h \
rbk.cc \
rbk.h \
rbk_impl.cc \
rbk_impl.h \
readlinsol.cc \
readlinsol.h \
reffrm.cc \
reffrm.h \
resforces.cc \
resforces.h \
restart.h \
rtposixsolver.cc \
rtposixsolver.h \
rtsolver.cc \
rtsolver.h \
sah.cc \
sah.h \
scalarvalue.cc \
scalarvalue.h \
shape.cc \
shape.h \
shape_impl.h \
shdrive.cc \
shdrive.h \
shockabsorber.h \
simentity.cc \
simentity.h \
sockdrv.cc \
sockdrv.h \
socketstreamdrive.cc \
socketstreamdrive.h \
socketstream_out_elem.cc \
socketstream_out_elem.h \
streamdrive.cc \
streamdrive.h \
streamoutelem.cc \
streamoutelem.h \
solver.cc \
solver.h \
solver_impl.h \
solverdiagnostics.cc \
solverdiagnostics.h \
stepsol.cc \
stepsol.h \
stepsol.hc \
symcltp.h \
tpldrive.cc \
tpldrive.h \
tpldrive_impl.cc \
tpldrive_impl.h \
vec3drv.h \
thirdorderstepsol.h \
thirdorderstepsol.cc \
TimeStepControl.h \
TimeStepControl.cc \
units.cc \
userelem.cc \
userelem.h \
usesock.cc \
usesock.h \
varstep.cc \
varstep.h \
ScalarFunctions.h \
ScalarFunctionsImpl.h \
ScalarFunctionsImpl.cc

nodist_libbase_la_SOURCES =

if BUILD_STATIC_MODULES
nodist_libbase_la_SOURCES += \
$(srcdir)/../../modules/module-minmaxdrive/module-minmaxdrive.cc \
$(srcdir)/../../modules/module-minmaxdrive/module-minmaxdrive.h \
$(srcdir)/../../modules/module-nodedistdrive/module-nodedistdrive.cc \
$(srcdir)/../../modules/module-nodedistdrive/module-nodedistdrive.h \
$(srcdir)/../../modules/module-hfelem/module-hfelem.cc \
$(srcdir)/../../modules/module-hfelem/module-hfelem.h
## add more...

if USE_OCTAVE
nodist_libbase_la_SOURCES += \
$(srcdir)/../../modules/module-octave/module-octave.cc \
$(srcdir)/../../modules/module-octave/module-octave.h \
$(srcdir)/../../modules/module-octave/octave_object.cc \
$(srcdir)/../../modules/module-octave/octave_object.h
endif
endif

if USE_GINAC
libbase_la_SOURCES += \
ginaccltp.h \
ginacdrive.cc \
ginacdrive.h
endif

if USE_RTAI
libbase_la_SOURCES += \
mbrtai_utils.c \
mbrtai_utils.h \
rtai_in_drive.cc \
rtai_in_drive.h \
rtai_out_elem.cc \
rtai_out_elem.h \
rtaisolver.cc \
rtaisolver.h
endif

if USE_JDQZ
libbase_la_SOURCES += \
eigjdqz.cc \
eigjdqz.h
endif

if USE_SCHUR
libbase_la_SOURCES += \
schurdataman.cc \
schurdataman.h
endif

libbase_la_LIBADD = @LIBS@ @SECURITY_LIBS@ $(GINACLIB_LIBS)
libbase_la_LDFLAGS = -static

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../libraries/libann \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base \
-I$(srcdir)/../../mbdyn/struct \
-I$(srcdir)/../../mbdyn/aero \
-I$(srcdir)/../../mbdyn/elec \
-I$(srcdir)/../../mbdyn/thermo \
-I$(srcdir)/../../mbdyn/hydr \
-I$(srcdir)/../../modules \
$(GINACLIB_CPPFLAGS) \
@OCTAVE_INCLUDE@

noinst_PROGRAMS = inusetest
inusetest_SOURCES = inusetest.cc
inusetest_LDADD = @ATOMIC_OPS_LIBS@ \
../../libraries/libmbutil/libmbutil.la

include $(top_srcdir)/build/bot.mk
