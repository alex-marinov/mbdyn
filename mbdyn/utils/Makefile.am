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
loadable \
posrel \
rot2eu \
rot2eup \
rot2phi \
rotrel \
sockdrive

LIBS = \
../libraries/libmbwrap/libmbwrap.la \
../libraries/libmbmath/libmbmath.la \
../libraries/libmbutil/libmbutil.la \
@Y12_LIBS@ \
@HARWELL_LIBS@ \
../libraries/libmbwrap/libmbwrap.la \
@MESCHACH_LIBS@ \
@F2C_LIBS@

autopilot_SOURCES = autopilot.cc
cl_SOURCES = cl.cc
dae_intg_SOURCES = dae-intg.cc
deriv_SOURCES = deriv.c
env_SOURCES = env.c
eu2rot_SOURCES = eu2rot.cc
femgen_SOURCES = femgen.f
intg_SOURCES = intg.cc
loadable_SOURCES = loadable.cc
posrel_SOURCES = posrel.cc
rot2eup_SOURCES = rot2eup.cc
rot2eu_SOURCES = rot2eu.cc
rot2phi_SOURCES = rot2phi.cc
rotrel_SOURCES = rotrel.cc
sockdrive_SOURCES = sockdrive.cc

intg_LDADD = -ldl
dae_intg_LDADD = -ldl

INCLUDES = \
-I../libraries/libmbutil \
-I../libraries/libmbmath \
-I../libraries/libmbwrap \
-I../libraries/libharwell \
-I../mbdyn

