include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libmbmath.a library
noinst_LTLIBRARIES = libmbmath.la

libmbmath_la_SOURCES = \
bisec.cc \
bisec.h \
ccmh.cc \
ccmh.h \
cscmhtpl.h \
dgeequ.cc \
dgeequ.h \
dirccmh.cc \
dirccmh.h \
fullmh.cc \
fullmh.h \
geomdata.h \
interp.h \
interp.cc \
JacSubMatrix.cc \
JacSubMatrix.h \
ldl.h \
ls.cc \
ls.h \
matvec.h \
matvec3.cc \
matvec3.h \
matvec3n.cc \
matvec3n.h \
matvec6.cc \
matvec6.h \
matvecass.h \
matvecexp.cc \
matvecexp.h \
mbcomm.h \
mh.cc \
mh.h \
MLS.cc \
MLS.h \
naivemh.cc \
naivemh.h \
Rot.cc \
Rot.hh \
RotCoeff.hc \
RotCoeff.hh \
solman.cc \
solman.h \
spmapmh.cc \
spmapmh.h \
spmh.cc \
spmh.h \
stlvh.cc \
stlvh.h \
submat.cc \
submat.h \
tpls.h \
vh.cc \
vh.h \
epetravh.h \
epetravh.cc \
epetraspmh.h \
epetraspmh.cc \
sp_exp_dof_map.h \
sp_gradient_base.h \
sp_gradient_expr.h \
sp_gradient_func.h \
sp_gradient_fwd.h \
sp_gradient.h \
sp_gradient_op.h \
sp_gradient_util.h \
sp_matrix_base_fwd.h \
sp_matrix_base.h \
sp_gradient.cc \
sp_gradient_spmh.h \
sp_gradient_spmh.cc \
siconosmh.h \
siconosmh.cc

if USE_SCHUR
libmbmath_la_SOURCES += \
schsolman.cc \
schsolman.h \
schurmh.h \
schurmh.cc
endif

if USE_FORTRAN
libmbmath_la_SOURCES += \
ldl.f \
spdata.h \
spdata_.f
else
libmbmath_la_SOURCES += \
ldl_.c
endif

libmbmath_la_LIBADD = @LIBS@ @FCLIBS@ @ANN_LIBS@ @TRILINOS_LIBS@
libmbmath_la_LDFLAGS =

noinst_PROGRAMS = matmultest itertest dgeequtest subtest

noinst_PROGRAMS += \
sp_gradient_test

matmultest_SOURCES = matmultest.cc
matmultest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@TRILINOS_LIBS@ \
@LIBS@

itertest_SOURCES = itertest.cc
itertest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@TRILINOS_LIBS@ \
@LIBS@

dgeequtest_SOURCES = dgeequtest.cc
dgeequtest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@TRILINOS_LIBS@ \
@LIBS@

subtest_SOURCES = subtest.cc
subtest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@TRILINOS_LIBS@ \
@LIBS@

sp_gradient_test_SOURCES = \
sp_gradient_test.cc \
sp_gradient_test_func.h \
sp_gradient_test_func.cc

sp_gradient_test_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@PASTIX_LIBS@ \
@QRUPDATE_LIBS@ \
@SUITESPARSEQR_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@TRILINOS_LIBS@ \
@RT_LIBS@

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/liby12 \
-I$(srcdir)/../../libraries/libcolamd \
-I$(srcdir)/../../libraries/libmbutil \
-I$(srcdir)/../../libraries/libmbmath \
-I$(srcdir)/../../libraries/libnaive \
-I$(srcdir)/../../libraries/libmbwrap \
-I$(srcdir)/../../mbdyn \
-I$(srcdir)/../../mbdyn/base

include $(top_srcdir)/build/bot.mk
