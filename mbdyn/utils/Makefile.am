## Process this file with automake to produce Makefile.in

# Build mbdyn utils
bin_PROGRAMS = rotrel eu2rot rot2eu rot2eup rot2phi \
deriv cl loadable posrel sockdrive autopilot intg dae-intg env

LIBS = ../libraries/libmbmath/libmbmath.la ../libraries/libmbutil/libmbutil.la

rotrel_SOURCES = rotrel.cc
eu2rot_SOURCES = eu2rot.cc
rot2eu_SOURCES = rot2eu.cc
rot2eup_SOURCES = rot2eup.cc
rot2phi_SOURCES = rot2phi.cc
deriv_SOURCES = deriv.c
cl_SOURCES = cl.cc
loadable_SOURCES = loadable.cc
posrel_SOURCES = posrel.cc
sockdrive_SOURCES = sockdrive.cc
autopilot_SOURCES = autopilot.cc
intg_SOURCES = intg.cc
dae_intg_SOURCES = dae-intg.cc
env_SOURCES = env.c

LIBS = ../libraries/libmbutil/libmbutil.la ../libraries/libmbmath/libmbmath.la @F2C_LIBS@
intg_LDADD = -ldl
dae_intg_LDADD = -ldl

INCLUDES = -I../libraries/libmbutil -I../libraries/libmbmath -I../mbdyn -I../mbdyn/base -I../mbdyn/struct -I../mbdyn/aero -I../mbdyn/elec -I../mbdyn/hydr

