include $(top_srcdir)/build/top.mk
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
distance.h \
distance.cc \
drvhinge.cc \
drvhinge.h \
drvj.cc \
drvj.h \
friction.cc \
genj.cc \
genj.h \
gravity.cc \
gravity.h \
hbeam.cc \
hbeam.h \
hbeam_interp.cc \
hbeam_interp.h \
inline.cc \
inline.h \
inplanej.cc \
inplanej.h \
inplcntj.h \
joint.cc \
joint.h \
joint_.h \
kin.cc \
kin.h \
kinj.cc \
kinj.h \
modal.cc \
modal.h \
module-wheel2.cc \
modules.h \
planej.cc \
planej.h \
prismj.cc \
prismj.h \
pzbeam.cc \
pzbeam.h \
pzbeam2.cc \
pzbeam2.h \
rodj.cc \
rodj.h \
shapefnc.cc \
shapefnc.h \
spherj.cc \
spherj.h \
strforce.cc \
strforce.h \
strnode.cc \
strnode.h \
univj.cc \
univj.h \
vehj.cc \
vehj.h \
vehj2.cc \
vehj2.h


libstruct_la_LIBADD = @LIBS@
libstruct_la_LDFLAGS = -static

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

module-wheel2.cc:
	if test -f module-wheel2.cc ; then \
		rm -f module-wheel2.cc ; \
	else \
		ln -s ../../modules/module-wheel2/module-wheel2.cc . ; \
	fi

	
