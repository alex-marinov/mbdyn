include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libaero.a library
noinst_LTLIBRARIES = libaero.la
libaero_la_SOURCES = \
aerod2.f \
aerod2.h \
aerodata.cc \
aerodata.h \
aerodata_impl.cc \
aerodata_impl.h \
aerodc81.c \
aerodc81.h \
aerodyn.cc \
aerodyn.h \
aerodyn_.h \
aeroelem.cc \
aeroelem.h \
aeroext.h \
aeroext.cc \
aeromodal.cc \
aeromodal.h \
c81data.cc \
c81data.h \
genfm.cc \
genfm.h \
gust.cc \
gust.h \
indvel.cc \
indvel.h \
instruments.cc \
instruments.h \
rotor.cc \
rotor.h \
windprof.cc \
windprof.h

nodist_libaero_la_SOURCES =

#if MBDYN_DEVEL
#nodist_libaero_la_SOURCES += \
#endif

libaero_la_LIBADD = @LIBS@
libaero_la_LDFLAGS = -static

if BUILD_STATIC_MODULES
nodist_libaero_la_SOURCES += \
$(srcdir)/../../modules/module-cyclocopter/module-cyclocopter.cc \
$(srcdir)/../../modules/module-cyclocopter/module-cyclocopter.h

if BUILD_CHARM
# charm
nodist_libaero_la_SOURCES += \
$(srcdir)/../../modules/module-charm/module-charm.cc \
$(srcdir)/../../modules/module-charm/mbcharm.h \
$(srcdir)/../../modules/module-charm/mbcharm_f.f \
$(srcdir)/../../modules/module-charm/charmdebug.f \
$(srcdir)/../../modules/module-charm/charm_common.c \
$(srcdir)/../../modules/module-charm/CharmWP.c

# NOTE: charm_common.c and CharmWP.c belong to CHARM's developers' kit
# and thus are not distributed with MBDyn

# NOTE: libWPModule.a must be in the load path; use LDFLAGS
libaero_la_LIBADD += -lWPModule

# add more...
endif
endif

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base \
-I$(srcdir)/../../mbdyn/struct \
-I$(srcdir)/../../mbdyn/aero \
-I$(srcdir)/../../mbdyn/elec \
-I$(srcdir)/../../mbdyn/hydr

include $(top_srcdir)/build/bot.mk
