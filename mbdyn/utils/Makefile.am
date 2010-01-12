include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build mbdyn utils
bin_PROGRAMS = \
ann_in \
ann_sf \
ann_sim \
ann_tr \
c81merge \
c81test \
cl \
crypt \
dae-intg \
deriv \
env \
eu2rot \
eu2phi \
femgen \
intg \
posrel \
rot2eu \
rot2eup \
rot2phi \
rotrel \
test_modalext \
test_modalext_edge \
vec2rot

if USE_RTAI
bin_PROGRAMS += \
logproc
endif

if USE_X_ANN
if MBDYN_DEVEL
bin_PROGRAMS += \
mbfem
endif
endif

if USE_SOCKET
bin_PROGRAMS += \
autopilot \
sockdrive \
socket2stream \
stream2socket \
test_modalext_socket \
trim

if MBDYN_DEVEL
bin_PROGRAMS += \
s2sctrl
endif
endif

# do not install these
noinst_PROGRAMS = \
playground \
test_modalext \
test_modalext_edge

if USE_SOCKET
noinst_PROGRAMS += \
test_modalext_socket
endif

MYLIBS = \
../libraries/libmbwrap/libmbwrap.la \
@Y12_LIBS@ \
../libraries/libann/libann.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@TCL_LIBS@ \
@LIBS@

ann_in_SOURCES = ann_in.c
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
eu2phi_SOURCES = eu2phi.cc
femgen_SOURCES = femgen.f
intg_SOURCES = intg.cc intg.h
playground_SOURCES = playground.cc
logproc_SOURCES = logproc.c
posrel_SOURCES = posrel.cc
rot2eup_SOURCES = rot2eup.cc
rot2eu_SOURCES = rot2eu.cc
rot2phi_SOURCES = rot2phi.cc
rotrel_SOURCES = rotrel.cc
nodist_s2sctrl_SOURCES = s2sctrl.cc s2s.cc s2s.h
sockdrive_SOURCES = sockdrive.c
socket2stream_SOURCES = socket2stream.cc s2s.cc s2s.h
stream2socket_SOURCES = stream2socket.cc s2s.cc s2s.h
test_modalext_SOURCES = test_modalext.cc
test_modalext_edge_SOURCES = test_modalext_edge.c
test_modalext_socket_SOURCES = test_modalext_socket.c
trim_SOURCES = trim.cc s2s.cc s2s.h
vec2rot_SOURCES = vec2rot.cc
nodist_mbfem_SOURCES = \
Mbfem/block.cc \
Mbfem/block.h \
Mbfem/force3.h \
Mbfem/mbdyn_out.cc \
Mbfem/mbdyn_out.h \
Mbfem/mbfem.cc \
Mbfem/mls.cc \
Mbfem/mls.h \
Mbfem/nas_read.cc \
Mbfem/nas_read.h \
Mbfem/node3.h

ann_in_LDADD = $(MYLIBS)
ann_sf_LDADD = $(MYLIBS)
ann_sim_LDADD = $(MYLIBS)
ann_tr_LDADD = $(MYLIBS)
autopilot_LDADD = $(MYLIBS) @SECURITY_LIBS@

# need to fool libtool in getting libbase linked twice ;)
# hope this is not fixed for a while...
DO_FOOL_LIBTOOL =
if FOOL_LIBTOOL
DO_FOOL_LIBTOOL += ../mbdyn/base/libbase.la
endif

c81merge_LDADD = \
$(DO_FOOL_LIBTOOL) \
../mbdyn/aero/libaero.la \
../mbdyn/base/libbase.la \
$(MYLIBS)
c81test_LDADD = ../mbdyn/aero/libaero.la \
../libraries/libmbutil/libmbutil.la \
$(MYLIBS)
crypt_LDADD = $(MYLIBS) @SECURITY_LIBS@
cl_LDADD = $(MYLIBS)
dae_intg_LDADD = $(MYLIBS)
eu2rot_LDADD = $(MYLIBS)
eu2phi_LDADD = $(MYLIBS)
intg_LDADD = $(MYLIBS)
playground_LDADD = $(MYLIBS)
posrel_LDADD = $(MYLIBS)
rot2eup_LDADD = $(MYLIBS)
rot2eu_LDADD = $(MYLIBS)
rot2phi_LDADD = $(MYLIBS)
rotrel_LDADD = $(MYLIBS)
s2sctrl_LDADD = $(MYLIBS)
sockdrive_LDADD = $(MYLIBS) @SECURITY_LIBS@
socket2stream_LDADD = $(MYLIBS)
stream2socket_LDADD = $(MYLIBS)
test_modalext_socket_LDADD = $(MYLIBS)
trim_LDADD = $(MYLIBS)
vec2rot_LDADD = $(MYLIBS)
mbfem_LDADD = $(MYLIBS) @X_ANN_LIBS@

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
# NOTE: this is a GNU make extension
%.so: %.o
	$(CXXLINK) -shared -o $@ $<

include $(top_srcdir)/build/bot.mk

