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
eu2rot \
eu2phi \
femgen \
intg \
posrel \
printenv \
rot2eu \
rot2eup \
rot2phi \
rotrel \
vec2rot

if USE_RTAI
bin_PROGRAMS += \
logproc
endif

if USE_ANN
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
trim

if MBDYN_DEVEL
bin_PROGRAMS += \
s2s_gpc \
s2s_dss
endif
endif

# do not install these
noinst_PROGRAMS = \
playground \
test_modalext \
test_modalext_edge \
test_strext_edge

if USE_OCTAVE
noinst_PROGRAMS += \
mbdyn_util_oct.oct
endif

if USE_SOCKET
noinst_PROGRAMS += \
test_modalext_socket \
test_strext_socket \
test_strext_socket_f \
test_strext_socket_cxx \
test_strext_socket_f_cxx
endif

MYLIBS = \
../libraries/libmbwrap/libmbwrap.la \
@Y12_LIBS@ \
../libraries/libmbc/libmbc.la \
../libraries/libann/libmbann.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
@KLU_LIBS@ \
@MESCHACH_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
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
printenv_SOURCES = env.c
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
s2s_dss_SOURCES = s2s.cc s2s.h
nodist_s2s_dss_SOURCES = s2s_dss.cc
s2s_gpc_SOURCES = s2s.cc s2s.h
nodist_s2s_gpc_SOURCES = s2s_gpc.cc GPC.c GPC.h udu.c udu.h
sockdrive_SOURCES = sockdrive.c
socket2stream_SOURCES = socket2stream.cc s2s.cc s2s.h
stream2socket_SOURCES = stream2socket.cc s2s.cc s2s.h
test_modalext_SOURCES = test_modalext.cc
test_modalext_edge_SOURCES = test_modalext_edge.c
test_modalext_socket_SOURCES = test_modalext_socket.c
test_strext_edge_SOURCES = test_strext_edge.c
test_strext_socket_SOURCES = \
test_strext_socket.c \
test_strext_socket_lib.c \
test_strext_socket_lib.h
test_strext_socket_f_SOURCES = \
test_strext_socket_f.f \
test_strext_socket_f_main.c \
test_strext_socket_lib.c \
test_strext_socket_lib.h
test_strext_socket_cxx_SOURCES = \
test_strext_socket.c \
test_strext_socket_lib_cxx.cc \
test_strext_socket_lib.h
test_strext_socket_f_cxx_SOURCES = \
test_strext_socket_f.f \
test_strext_socket_f_main.c \
test_strext_socket_lib_cxx.cc \
test_strext_socket_lib.h
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

mbdyn_util_oct_oct_SOURCES = mbdyn_util_oct.cc
if USE_OCTAVE
OCTAVE_INCLUDE = -I$(top_srcdir)/include -I$(top_srcdir)/libraries/libmbutil -I$(top_srcdir)/libraries/libmbmath -I$(top_srcdir)/mbdyn -I../include
OCTAVE_LIBS = -L../libraries/libmbmath/.libs -L../libraries/libmbutil/.libs -lmbmath -lmbutil
OCTAVE_FUNCS = \
rotation_matrix_to_rotation_vector.oct \
rotation_vector_to_rotation_matrix.oct \
rotation_matrix_to_euler123.oct \
euler123_to_rotation_matrix.oct \
rotation_matrix_to_euler313.oct \
euler313_to_rotation_matrix.oct \
rotation_matrix_to_euler321.oct \
euler321_to_rotation_matrix.oct

mbdyn_util_oct.o: mbdyn_util_oct.cc
	mkoctfile $(OCTAVE_INCLUDE) -c -o $@ $<

mbdyn_util_oct.oct: mbdyn_util_oct.o
	mkoctfile -o $@ $< $(OCTAVE_LIBS)
endif

ann_in_LDADD = $(MYLIBS)
ann_sf_LDADD = $(MYLIBS)
ann_sim_LDADD = $(MYLIBS)
ann_tr_LDADD = $(MYLIBS)
autopilot_LDADD = $(MYLIBS) @SECURITY_LIBS@
mbdyn_util_oct_oct_LDADD = $(MYLIBS)

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
s2s_dss_LDADD = $(MYLIBS)
s2s_gpc_LDADD = $(MYLIBS)
sockdrive_LDADD = $(MYLIBS) @SECURITY_LIBS@
socket2stream_LDADD = $(MYLIBS)
stream2socket_LDADD = $(MYLIBS)
test_modalext_socket_LDADD = $(MYLIBS)
test_strext_socket_LDADD = $(MYLIBS)
test_strext_socket_f_LDADD = $(MYLIBS)
test_strext_socket_cxx_LDADD = $(MYLIBS)
test_strext_socket_f_cxx_LDADD = $(MYLIBS)
trim_LDADD = $(MYLIBS)
vec2rot_LDADD = $(MYLIBS)
mbfem_LDADD = $(MYLIBS) @ANN_LIBS@

logproc_LDFLAGS = -static

INCLUDES = \
-I../include \
-I$(srcdir)/../include \
-I$(srcdir)/../libraries/libmbc \
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

install-exec-local-oct:
	$(mkinstalldirs) $(DESTDIR)$(datarootdir) $(DESTDIR)$(octavedir)
	$(INSTALL_SCRIPT) -m 755 mbdyn_util_oct.oct $(DESTDIR)$(octavedir)
	@echo "creating symlinks..."
	@for f in $(OCTAVE_FUNCS) ; do \
		ln -s $(DESTDIR)$(octavedir)/mbdyn_util_oct.oct $(DESTDIR)$(octavedir)/$$f ; \
	done

install-exec-local: install-exec-local-oct

include $(top_srcdir)/build/bot.mk

