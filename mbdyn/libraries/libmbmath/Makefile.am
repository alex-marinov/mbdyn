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
dgeequ.cc \
dgeequ.h \
dirccmh.cc \
dirccmh.h \
fullmh.cc \
fullmh.h \
geomdata.h \
gradient.h \
interp.h \
interp.cc \
JacSubMatrix.cc \
JacSubMatrix.h \
ldl.f \
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
spdata.h \
spdata_.f \
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
gradient.h \
matvec.h

if USE_SCHUR
libmbmath_la_SOURCES += \
schsolman.cc \
schsolman.h \
schurmh.h \
schurmh.cc
endif

libmbmath_la_LIBADD = @LIBS@ @FCLIBS@ @ANN_LIBS@
libmbmath_la_LDFLAGS =

noinst_PROGRAMS = matmultest itertest dgeequtest subtest matvectest

if USE_AUTODIFF
noinst_PROGRAMS += \
gradienttest
endif

matmultest_SOURCES = matmultest.cc
matmultest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@

itertest_SOURCES = itertest.cc
itertest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@

dgeequtest_SOURCES = dgeequtest.cc
dgeequtest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@

subtest_SOURCES = subtest.cc
subtest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@

gradienttest_SOURCES = gradienttest.cc
gradienttest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@ \
@RT_LIBS@

matvectest_SOURCES = matvectest.cc matvectestfor.f
matvectest_LDADD = \
libmbmath.la \
../libmbutil/libmbutil.la \
../libcolamd/libmbdyncolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@Y12_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@FCLIBS@ \
@LIBS@ \
@THREAD_LIBS@ \
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
