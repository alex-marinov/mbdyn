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

DO_FOOL_LIBTOOL =
if FOOL_LIBTOOL
DO_FOOL_LIBTOOL += ../mbdyn/base/libbase.la
endif

mbdyn_LDADD = \
$(DO_FOOL_LIBTOOL) \
../mbdyn/base/libbase.la  \
aero/libaero.la \
struct/libstruct.la \
elec/libelec.la \
thermo/libthermo.la \
hydr/libhydr.la \
base/libbase.la  \
../libraries/libmbwrap/libmbwrap.la \
../libraries/libann/libann.la \
$(GINACLIB_LIBS) \
@Y12_LIBS@ \
@HARWELL_LIBS@ \
@ARPACK_LIBS@ \
@LAPACK_LIBS@ \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
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
@RT_LIBS@ \
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
