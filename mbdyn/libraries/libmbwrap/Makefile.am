include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libmbwrap.a library
lib_LTLIBRARIES = libmbwrap.la

libmbwrap_la_SOURCES = \
harwlib.h \
harwout.c \
harwout.h \
harwout_.cc \
harwrap.cc \
harwrap.h \
lapackwrap.cc \
lapackwrap.h \
mschwrap.cc \
mschwrap.h \
naivewrap.cc \
naivewrap.h \
parnaivewrap.cc \
parnaivewrap.h \
superluwrap.cc \
superluwrap.h \
umfpackwrap.cc \
umfpackwrap.h \
taucswrap.cc \
taucswrap.h \
y12wrap.cc \
y12wrap.h \
wsmpwrap.cc \
wsmpwrap.h

libmbwrap_la_LIBADD = \
../libmbmath/libmbmath.la \
../libmbutil/libmbutil.la \
../libobjs/libobjs.la \
../liby12/liby12.la \
../libcolamd/libcolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@WSMP_LIBS@ \
@LAPACK_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@ \
@THREAD_LIBS@

libmbwrap_la_LDFLAGS =

noinst_PROGRAMS = wraptest cctest

wraptest_SOURCES = wraptest.cc
wraptest_LDADD = \
./libmbwrap.la \
../libmbmath/libmbmath.la \
../libmbutil/libmbutil.la \
../libobjs/libobjs.la \
../liby12/liby12.la \
../libcolamd/libcolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@TAUCS_LIBS@ \
@LAPACK_LIBS@ \
@METIS_LIBS@ \
@BLAS_LIBS@ \
@F2C_LIBS@ \
@LIBS@ \
@TCL_LIBS@ \
@THREAD_LIBS@

cctest_SOURCES = cctest.cc
cctest_LDADD = \
./libmbwrap.la \
../libmbmath/libmbmath.la \
../libmbutil/libmbutil.la \
../libobjs/libobjs.la \
../liby12/liby12.la \
../libcolamd/libcolamd.la \
../libnaive/libnaive.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@SUPERLU_LIBS@ \
@BLAS_LIBS@ \
@LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@ \
@THREAD_LIBS@

INCLUDES = \
-I../../mbdyn \
-I../libmbutil \
-I../libmbmath \
-I../liby12 \
-I../libcolamd \
-I../libnaive

include $(top_srcdir)/build/bot.mk
