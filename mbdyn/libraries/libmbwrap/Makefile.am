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
mschwrap.cc \
mschwrap.h \
umfpackwrap.cc \
umfpackwrap.h \
y12wrap.cc \
y12wrap.h

libmbwrap_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbwrap_la_LDFLAGS =

noinst_PROGRAMS = wraptest
wraptest_SOURCES = test.cc
wraptest_LDADD = \
./libmbwrap.la \
../libmbmath/libmbmath.la \
../libmbutil/libmbutil.la \
../libobjs/libobjs.la \
../liby12/liby12.la \
@UMFPACK_LIBS@ \
@MESCHACH_LIBS@ \
@HARWELL_LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@

INCLUDES = \
-I../../mbdyn \
-I../libmbutil \
-I../libmbmath \
-I../liby12

include $(top_srcdir)/build/bot.mk
