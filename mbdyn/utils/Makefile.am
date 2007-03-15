include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build mbdyn utils
bin_PROGRAMS = \
ann_sf \
ann_sim \
ann_tr \
autopilot \
c81merge \
c81test \
cl \
crypt \
dae-intg \
deriv \
env \
eu2rot \
femgen \
intg \
logproc \
mbplot \
posrel \
rot2eu \
rot2eup \
rot2phi \
rotrel \
sockdrive \
socket2stream \
stream2socket \
trim \
vec2rot

# do not install these
noinst_PROGRAMS = \
playground

MYLIBS = \
../libraries/libobjs/libobjs.la \
../libraries/libmbwrap/libmbwrap.la \
@Y12_LIBS@ \
../libraries/libcolamd/libcolamd.la \
../libraries/libann/libann.la \
../libraries/libmbmath/libmbmath.la \
../libraries/libmbutil/libmbutil.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@BLAS_LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@ \
@LIBS@

ann_sf_SOURCES = ann_sf.c
ann_sim_SOURCES = ann_sim.c
ann_tr_SOURCES = ann_tr.c
autopilot_SOURCES = autopilot.c
c81merge_SOURCES = c81merge.cc
c81test_SOURCES = c81test.cc
crypt_SOURCES = crypt.c
cl_SOURCES = cl.cc
dae_intg_SOURCES = dae-intg.cc dae-intg.h
deriv_SOURCES = deriv.c
env_SOURCES = env.c
eu2rot_SOURCES = eu2rot.cc
femgen_SOURCES = femgen.f
intg_SOURCES = intg.cc
playground_SOURCES = playground.cc
logproc_SOURCES = logproc.c
mbplot_SOURCES = mbplot.c
posrel_SOURCES = posrel.cc
rot2eup_SOURCES = rot2eup.cc
rot2eu_SOURCES = rot2eu.cc
rot2phi_SOURCES = rot2phi.cc
rotrel_SOURCES = rotrel.cc
sockdrive_SOURCES = sockdrive.c
socket2stream_SOURCES = socket2stream.cc s2s.cc s2s.h
stream2socket_SOURCES = stream2socket.cc s2s.cc s2s.h
trim_SOURCES = trim.cc s2s.cc s2s.h
vec2rot_SOURCES = vec2rot.cc

ann_sf_LDADD = $(MYLIBS)
ann_sim_LDADD = $(MYLIBS)
ann_tr_LDADD = $(MYLIBS)
autopilot_LDADD = $(MYLIBS) @SECURITY_LIBS@
# need to fool libtool in getting libbase linked twice ;)
# hope this is not fixed for a while...
c81merge_LDADD = \
../mbdyn/base/libbase.la \
../mbdyn/aero/libaero.la \
../mbdyn/../mbdyn/base/libbase.la \
$(MYLIBS)
c81test_LDADD = ../mbdyn/aero/libaero.la \
../libraries/libmbutil/libmbutil.la \
$(MYLIBS)
crypt_LDADD = $(MYLIBS) @SECURITY_LIBS@
cl_LDADD = $(MYLIBS)
dae_intg_LDADD = $(MYLIBS)
eu2rot_LDADD = $(MYLIBS)
intg_LDADD = $(MYLIBS)
playground_LDADD = $(MYLIBS)
posrel_LDADD = $(MYLIBS)
rot2eup_LDADD = $(MYLIBS)
rot2eu_LDADD = $(MYLIBS)
rot2phi_LDADD = $(MYLIBS)
rotrel_LDADD = $(MYLIBS)
sockdrive_LDADD = $(MYLIBS) @SECURITY_LIBS@
socket2stream_LDADD = $(MYLIBS)
stream2socket_LDADD = $(MYLIBS)
trim_LDADD = $(MYLIBS)
vec2rot_LDADD = $(MYLIBS)

logproc_LDFLAGS = -static

INCLUDES = \
-I../include \
-I$(srcdir)/../include \
-I$(srcdir)/../libraries/libmbutil \
-I$(srcdir)/../libraries/libmbmath \
-I$(srcdir)/../libraries/libmbwrap \
-I$(srcdir)/../libraries/libann \
-I$(srcdir)/../mbdyn \
-I$(srcdir)/../mbdyn/aero \
-I$(srcdir)/../mbdyn/base

# Rule to compile modules; run "make <module-name>.so"
%.so: %.o
	$(CXXLINK) -shared -o $@ $<

include $(top_srcdir)/build/bot.mk

