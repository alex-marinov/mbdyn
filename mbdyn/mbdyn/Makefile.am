include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build static libraries in subdirs
##SUBDIRS = @XSUBDIRS@
SUBDIRS = base struct aero elec hydr thermo

# Build mbdyn executable
bin_PROGRAMS = mbdyn
mbdyn_SOURCES = \
dummy_c.c \
dummy_f.f  \
mbdyn.cc  \
mbdyn.h

# need to fool libtool in getting libbase linked twice ;)
# hope this is not fixed for a while...
mbdyn_LDADD = \
../mbdyn/base/libbase.la \
aero/libaero.la \
../libraries/libmbwrap/libmbwrap.la \
struct/libstruct.la \
elec/libelec.la \
thermo/libthermo.la \
hydr/libhydr.la \
base/libbase.la  \
$(GINACLIB_LIBS)\
../libraries/libmbmath/libmbmath.la \
../libraries/libmbutil/libmbutil.la \
@Y12_LIBS@ \
../libraries/libcolamd/libmbdyncolamd.la \
../libraries/libnaive/libnaive.la \
../libraries/libann/libann.la \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
../libraries/libobjs/libobjs.la \
@MESCHACH_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@MPI_LIBS@ \
@BLAS_LIBS@ \
@F2C_LIBS@ \
@CHACO_LIBS@ \
@METIS_LIBS@ \
@SECURITY_LIBS@ \
@THREAD_LIBS@ \
@TCL_LIBS@ \
@NETCDF_LIBS@ \
@UDUNITS_LIBS@ \
@X_ANN_LIBS@ \
@LIBS@ \
$(GINACLIB_LIBS)

INCLUDES = \
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
