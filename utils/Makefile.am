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

if HAVE_GETPASS
bin_PROGRAMS += \
sockdrive \
socket2stream \
stream2socket \
trim
endif

if HAVE_TERMIOS
bin_PROGRAMS += \
autopilot
endif

if MBDYN_DEVEL
if !WINDOWS
bin_PROGRAMS += \
s2s_gpc \
s2s_dss
endif
endif

endif

if INSTALL_TEST_PROGRAMS
bin_PROGRAMS += \
playground \
test_modalext \
test_modalext_edge \
test_strext_edge
noinst_PROGRAMS =
else
# do not install these
noinst_PROGRAMS = \
playground \
test_modalext \
test_modalext_edge \
test_strext_edge
endif


if USE_NETCDF4
noinst_PROGRAMS += \
ncprint
endif

if USE_OCTAVE
noinst_PROGRAMS += \
mbdyn_util_oct.oct
endif

noinst_LTLIBRARIES = libfemgen.la

if USE_SOCKET
if INSTALL_TEST_PROGRAMS
bin_PROGRAMS += \
test_modalext_socket \
test_strext_socket \
test_strext_socket_f \
test_strext_socket_cxx \
test_strext_socket_f_cxx
else
noinst_PROGRAMS += \
test_modalext_socket \
test_strext_socket \
test_strext_socket_f \
test_strext_socket_cxx \
test_strext_socket_f_cxx
endif

# fake library to tweak the selected linker
noinst_LTLIBRARIES += libstrext_socket_f.la
endif

MYLIBS = \
../libraries/libmbutil/libmbutil.la \
../libraries/libmbwrap/libmbwrap.la \
@Y12_LIBS@

if BUILD_MBC
MYLIBS += \
../libraries/libmbc/libmbc.la
endif

MYLIBS += \
../libraries/libann/libmbann.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
@KLU_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
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
femgen_SOURCES = femgen.c
libfemgen_la_SOURCES = femgen_f.f
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
test_strext_socket_cxx_SOURCES = \
test_strext_socket.c \
test_strext_socket_lib_cxx.cc \
test_strext_socket_lib.h
# most of the code is in libstrext_socket_f.la to force linking with C linker
test_strext_socket_f_SOURCES = \
test_strext_socket_lib.c \
test_strext_socket_lib.h
# most of the code is in libstrext_socket_f.la to force linking with C++ linker
test_strext_socket_f_cxx_SOURCES = \
test_strext_socket_lib_cxx.cc \
test_strext_socket_lib.h
libstrext_socket_f_la_SOURCES = \
test_strext_socket_f.f \
test_strext_socket_f_main.c \
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

ncprint_SOURCES = ncprint.cc

mbdyn_util_oct_oct_SOURCES = mbdyn_util_oct.cc
if USE_OCTAVE
OCTAVE_XINCLUDE = -I$(top_srcdir)/include -I$(top_srcdir)/libraries/libmbutil -I$(top_srcdir)/libraries/libmbmath -I$(top_srcdir)/mbdyn -I../include
OCTAVE_XLIBS = -L../libraries/libmbmath/.libs -L../libraries/libmbutil/.libs -lmbmath -lmbutil
OCTAVE_XFUNCS = \
rotation_matrix_to_rotation_vector.oct \
rotation_vector_to_rotation_matrix.oct \
rotation_matrix_to_euler123.oct \
euler123_to_rotation_matrix.oct \
rotation_matrix_to_euler313.oct \
euler313_to_rotation_matrix.oct \
rotation_matrix_to_euler321.oct \
euler321_to_rotation_matrix.oct
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
femgen_LDADD = libfemgen.la
intg_LDADD = $(MYLIBS)
logproc_LDADD = @RTAI_LDFLAGS@
playground_LDADD = $(MYLIBS)
ncprint_LDADD = $(MYLIBS) @NETCDF_LIBS@
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
test_modalext_LDADD = $(MYLIBS)
test_modalext_edge_LDADD = $(MYLIBS)
test_modalext_socket_LDADD = $(MYLIBS)
test_strext_edge_LDADD = $(MYLIBS)
test_strext_socket_LDADD = $(MYLIBS)
test_strext_socket_f_LDADD = $(MYLIBS) libstrext_socket_f.la
test_strext_socket_cxx_LDADD = $(MYLIBS)
test_strext_socket_f_cxx_LDADD = $(MYLIBS) libstrext_socket_f.la
trim_LDADD = $(MYLIBS)
vec2rot_LDADD = $(MYLIBS)
mbfem_LDADD = $(MYLIBS) @ANN_LIBS@

# fool linking on MinGW/MSYS (WTF?!?)
# test_strext_edge_LINK = $(CXXLD) -o $@

logproc_LDFLAGS = -static

AM_CPPFLAGS = \
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

if USE_OCTAVE
PKG_NAME:=mbdyn_util_oct
COPYRIGHT:=pierangelo.masarati@polimi.it
PKG_AUTHOR:=r.resch@a1.net
PKG_DIR:=$(PKG_NAME)
PKG_FILE:=$(PKG_NAME).tar.gz

mbdyn_util_oct.o: mbdyn_util_oct.cc
	$(OCTAVE_MKOCTFILE) $(OCTAVE_XINCLUDE) -c -o $@ $<

mbdyn_util_oct.oct: mbdyn_util_oct.o
	$(OCTAVE_MKOCTFILE) -o $@ $< $(OCTAVE_XLIBS)

# workaround for Cygwin and MinGW
mbdyn_util_oct.oct.exe: mbdyn_util_oct.oct

$(PKG_FILE): $(top_srcdir)/utils/PKG_ADD mbdyn_util_oct.oct
	@rm -rf $(PKG_DIR) $(PKG_FILE) 
	@$(MKDIR_P) $(PKG_DIR)/package/inst
	@echo "Copyright $(COPYRIGHT)" > $(PKG_DIR)/package/COPYING
	@echo "Name: $(PKG_NAME)" > $(PKG_DIR)/package/DESCRIPTION
	@echo "Version: $(VERSION)" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Date: 2019-04-07" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Author: $(PKG_AUTHOR)" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Maintainer: $(PKG_AUTHOR)" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Title: MBDyn Octave package" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Description: MBDyn rotation handling functions" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "License: GPLv2" >> $(PKG_DIR)/package/DESCRIPTION
	@echo "Categories: Multibody dynamics" >> $(PKG_DIR)/package/DESCRIPTION
	@$(INSTALL_SCRIPT) mbdyn_util_oct.oct $(PKG_DIR)/package/inst/
	@$(INSTALL_SCRIPT) $(top_srcdir)/utils/PKG_ADD $(PKG_DIR)/package/PKG_ADD
	@cd $(PKG_DIR); tar -zcf $(PKG_FILE) package
	@cp -f $(PKG_DIR)/$(PKG_FILE) $(PKG_FILE)
	@rm -rf $(PKG_DIR)
	@printf "Octave package file \"%s\" created\n" $(PKG_FILE)

install-exec-local-oct: $(PKG_FILE)
	@octave-cli -qfHW --eval "pkg install $<"
	@printf "package \"%s\" installed\n" $(PKG_NAME)
else
install-exec-local-oct:
endif

install-exec-local: install-exec-local-oct

include $(top_srcdir)/build/bot.mk

