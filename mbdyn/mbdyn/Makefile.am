include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build static libraries in subdirs
##SUBDIRS = @XSUBDIRS@
SUBDIRS = base struct aero elec hydr

# Build mbdyn executable
bin_PROGRAMS = mbdyn
mbdyn_SOURCES = \
dummy_c.c \
dummy_f.f  \
mbdyn.cc  \
mbdyn.h

mbdyn_LDADD = \
../mbdyn/base/libbase.la \
aero/libaero.la \
../libraries/libmbwrap/libmbwrap.la \
struct/libstruct.la \
elec/libelec.la \
hydr/libhydr.la \
base/libbase.la \
../libraries/libmbmath/libmbmath.la \
../libraries/libmbutil/libmbutil.la \
@Y12_LIBS@ \
@HARWELL_LIBS@ \
@UMFPACK_LIBS@ \
../libraries/libobjs/libobjs.la \
@MESCHACH_LIBS@ @LAPACK_LIBS@ @MPI_LIBS@ @F2C_LIBS@ \
@METIS_LIBS@ @SECURITY_LIBS@ @TCL_LIBS@ @LIBS@ \
$(GINACLIB_LIBS)

INCLUDES = \
-I../libraries/libmbutil \
-I../libraries/libmbmath \
-I../libraries/libmbwrap \
-I../mbdyn/base \
-I../mbdyn/struct \
-I../mbdyn/aero \
-I../mbdyn/elec \
-I../mbdyn/hydr \
$(GINACLIB_CPPFLAGS)
include $(top_srcdir)/build/bot.mk
