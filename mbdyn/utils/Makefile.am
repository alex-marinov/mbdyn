include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build mbdyn utils
bin_PROGRAMS = \
autopilot \
cl \
dae-intg \
deriv \
env \
eu2rot \
femgen \
intg \
mbplot \
posrel \
rot2eu \
rot2eup \
rot2phi \
rotrel \
sockdrive \
vec2rot

LIBS = \
../libraries/libobjs/libobjs.la \
../libraries/libmbwrap/libmbwrap.la \
@Y12_LIBS@ \
../libraries/libmbmath/libmbmath.la \
../libraries/libmbutil/libmbutil.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@ \
@LIBS@

autopilot_SOURCES = autopilot.cc
cl_SOURCES = cl.cc
dae_intg_SOURCES = dae-intg.cc dae-intg.h
deriv_SOURCES = deriv.c
env_SOURCES = env.c
eu2rot_SOURCES = eu2rot.cc
femgen_SOURCES = femgen.f
intg_SOURCES = intg.cc
mbplot_SOURCES = mbplot.c
posrel_SOURCES = posrel.cc
rot2eup_SOURCES = rot2eup.cc
rot2eu_SOURCES = rot2eu.cc
rot2phi_SOURCES = rot2phi.cc
rotrel_SOURCES = rotrel.cc
sockdrive_SOURCES = sockdrive.cc
vec2rot_SOURCES = vec2rot.cc

intg_LDADD = -ldl
dae_intg_LDADD = -ldl

INCLUDES = \
-I../libraries/libmbutil \
-I../libraries/libmbmath \
-I../libraries/libmbwrap \
-I../mbdyn

# Rule to compile modules; run "make <module-name>.so"
%.so: %.o
	$(CXXLINK) -shared -o $@ $<

include $(top_srcdir)/build/bot.mk
