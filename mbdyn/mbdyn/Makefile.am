## Process this file with automake to produce Makefile.in

# Build static libraries in subdirs
##SUBDIRS = @XSUBDIRS@
SUBDIRS = base struct aero elec hydr

# Build mbdyn executable
bin_PROGRAMS = mbdyn
mbdyn_SOURCES = mbdyn.cc mbdyn.h dummy_f.f dummy_c.c
## mbdyn_LDADD = @XLIBOBJS@
mbdyn_LDADD = \
	base/libbase.la \
	struct/libstruct.la \
	aero/libaero.la \
	elec/libelec.la \
	hydr/libhydr.la \
	../libraries/libmbmath/libmbmath.la \
	../libraries/libmbutil/libmbutil.la \
	../libraries/liby12/liby12.la \
	../libraries/libobjs/libobjs.la \
	@MESCHACH_LIBS@ @LAPACK_LIBS@ @MPI_LIBS@ @F2C_LIBS@ \
	@METIS_LIBS@ @PAM_LIBS@ @TCL_LIBS@ @LIBS@
mbdyn_LDFLAGS = -static

INCLUDES = -I../libraries/libmbutil -I../libraries/libmbmath -I../mbdyn/base -I../mbdyn/struct -I../mbdyn/aero -I../mbdyn/elec -I../mbdyn/hydr

