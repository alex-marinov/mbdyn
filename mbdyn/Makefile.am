include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build static libraries in subdirs
##SUBDIRS = @XSUBDIRS@
SUBDIRS = base struct aero elec hydr thermo

# Build mbdyn executable
bin_PROGRAMS = mbdyn
mbdyn_SOURCES = \
mbdyn.cc  \
mbdyn.h

# need to fool libtool in getting libbase linked twice ;)
# hope this is not fixed for a while...

DO_FOOL_LIBTOOL =
if FOOL_LIBTOOL
DO_FOOL_LIBTOOL += ../mbdyn/base/libbase.la
endif

mbdyn_LDADD = \
$(DO_FOOL_LIBTOOL) \
aero/libaero.la \
struct/libstruct.la \
elec/libelec.la \
thermo/libthermo.la \
hydr/libhydr.la \
base/libbase.la

if BUILD_MBC
mbdyn_LDADD += \
../libraries/libmbc/libmbc_static.la
endif

mbdyn_LDADD += \
../libraries/libmbwrap/libmbwrap.la \
../libraries/libann/libmbann.la \
$(GINACLIB_LIBS) \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@HARWELL_LIBS@ \
@ARPACK_LIBS@ \
@JDQZ_LIBS@ \
@LAPACK_LIBS@ \
@UMFPACK_LIBS@ \
@KLU_LIBS@ \
@MESCHACH_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@MPI_LIBS@ \
@CHACO_LIBS@ \
@METIS_LIBS@ \
@SECURITY_LIBS@ \
@NETCDF_LIBS@ \
@UDUNITS_LIBS@ \
@RTAI_LDFLAGS@ \
@RT_LIBS@ \
@ANN_LIBS@ \
@OCTAVE_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
$(GINACLIB_LIBS)

### note: to avoid dynamic linking, in mbdyn_LDADD use
### ../libraries/libmbc/.libs/libmbc.a

AM_CPPFLAGS = \
-I../include \
-I$(srcdir)/../include \
-I$(srcdir)/../libraries/libmbutil \
-I$(srcdir)/../libraries/libmbmath \
-I$(srcdir)/../libraries/libmbwrap \
-I$(srcdir)/../mbdyn/base \
-I$(srcdir)/../mbdyn/struct \
-I../mbdyn/thermo \
-I$(srcdir)/../mbdyn/aero \
-I$(srcdir)/../mbdyn/elec \
-I$(srcdir)/../mbdyn/hydr \
$(GINACLIB_CPPFLAGS)

include $(top_srcdir)/build/bot.mk
