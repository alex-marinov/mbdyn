include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libstruct.a library
noinst_LTLIBRARIES = libstruct.la
libstruct_la_SOURCES = \
accj.cc \
accj.h \
autostr.cc \
autostr.h \
beam.cc \
beam.h \
beam2.cc \
beam2.h \
beamslider.h \
beamslider.cc \
body.cc \
body.h \
body_vm.cc \
body_vm.h \
brake.cc \
brake.h \
distance.h \
distance.cc \
drvdisp.cc \
drvdisp.h \
drvhinge.cc \
drvhinge.h \
drvj.cc \
drvj.h \
friction.cc \
friction.h \
genj.cc \
genj.h \
gimbal.cc \
gimbal.h \
gravity.cc \
gravity.h \
hbeam.cc \
hbeam.h \
hbeam_interp.cc \
hbeam_interp.h \
inertia.cc \
inertia.h \
inline.cc \
inline.h \
inplanej.cc \
inplanej.h \
inplcntj.h \
joint.cc \
joint.h \
joint_.h \
jointreg.cc \
jointreg.h \
membrane.h \
membrane.cc \
membraneeas.h \
membraneeas.cc \
modal.cc \
modal.h \
modaledge.cc \
modaledge.h \
modalext.cc \
modalext.h \
modalmappingext.cc \
modalmappingext.h \
modalforce.cc \
modalforce.h \
planej.cc \
planej.h \
point_contact.cc \
point_contact.h \
prismj.cc \
prismj.h \
pzbeam.cc \
pzbeam.h \
pzbeam2.cc \
pzbeam2.h \
rodbezj.cc \
rodbezj.h \
rodj.cc \
rodj.h \
shapefnc.cc \
shapefnc.h \
shell.cc \
shell.h \
shell.hc \
shelleasans.cc \
shelleasans.h \
shelleas.cc \
shelleas.h \
socketstreammotionelem.cc \
socketstreammotionelem.h \
spherj.cc \
spherj.h \
stredge.cc \
stredge.h \
strext.cc \
strext.h \
strmappingext.cc \
strmappingext.h \
strforce.cc \
strforce.h \
strforce_impl.h \
strnode.cc \
strnode.h \
tdclw.cc \
tdclw.h \
totalequation.cc \
totalequation.h \
totalj.cc \
totalj.h \
univj.cc \
univj.h \
vb.cc \
vb.h \
vehj.cc \
vehj.h \
vehj2.cc \
vehj2.h \
vehj3.cc \
vehj3.h \
vehj4.cc \
vehj4.h

nodist_libstruct_la_SOURCES =

if BUILD_STATIC_MODULES
# wheel2 and more
nodist_libstruct_la_SOURCES += \
$(srcdir)/../../modules/module-wheel2/module-wheel2.cc \
$(srcdir)/../../modules/module-wheel2/module-wheel2.h \
\
$(srcdir)/../../modules/module-asynchronous_machine/module-asynchronous_machine.cc \
$(srcdir)/../../modules/module-asynchronous_machine/module-asynchronous_machine.h \
\
$(srcdir)/../../modules/module-inline_friction/module-inline_friction.cc \
$(srcdir)/../../modules/module-inline_friction/module-inline_friction.h \
\
$(srcdir)/../../modules/module-nodedistdrive/module-nodedistdrive.cc \
$(srcdir)/../../modules/module-nodedistdrive/module-nodedistdrive.h \
\
$(srcdir)/../../modules/module-multi_step_drive/module-multi_step_drive.cc \
$(srcdir)/../../modules/module-multi_step_drive/module-multi_step_drive.h \
\
$(srcdir)/../../modules/module-switch_drive/module-switch_drive.cc \
$(srcdir)/../../modules/module-switch_drive/module-switch_drive.h

if USE_AUTODIFF
nodist_libstruct_la_SOURCES += \
$(srcdir)/../../modules/module-uni_in_plane/module-uni_in_plane.h \
$(srcdir)/../../modules/module-uni_in_plane/module-uni_in_plane.cc \
\
$(srcdir)/../../modules/module-ballbearing_contact/module-ballbearing_contact.h \
$(srcdir)/../../modules/module-ballbearing_contact/module-ballbearing_contact.cc
endif

if USE_SPARSE_AUTODIFF
nodist_libstruct_la_SOURCES += \
$(srcdir)/../../modules/module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.cc \
$(srcdir)/../../modules/module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.h \
\
$(srcdir)/../../modules/module-hydrodynamic_plain_bearing2/module-hydrodynamic_plain_bearing2.cc \
$(srcdir)/../../modules/module-hydrodynamic_plain_bearing2/module-hydrodynamic_plain_bearing2.h \
\
$(srcdir)/../../modules/module-journal_bearing/module-journal_bearing.h \
$(srcdir)/../../modules/module-journal_bearing/module-journal_bearing.cc \
\
$(srcdir)/../../modules/module-triangular_contact/module-triangular_contact.h \
$(srcdir)/../../modules/module-triangular_contact/module-triangular_contact.cc
endif

# add more...
endif

# screw
if MBDYN_DEVEL
nodist_libstruct_la_SOURCES += \
screwjoint.h \
screwjoint.cc
endif

libstruct_la_LIBADD = @LIBS@
libstruct_la_LDFLAGS = -static

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base \
-I$(srcdir)/../../mbdyn/struct \
-I$(srcdir)/../../mbdyn/aero \
-I$(srcdir)/../../mbdyn/elec \
-I$(srcdir)/../../mbdyn/hydr

include $(top_srcdir)/build/bot.mk

