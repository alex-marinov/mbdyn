## Process this file with automake to produce Makefile.in

# Build libmbwrap.a library
lib_LTLIBRARIES = libmbwrap.la

libmbwrap_la_SOURCES = \
harwout.c \
harwout.h \
harwout_.cc \
harwrap.cc \
harwrap.h \
mschwrap.cc \
mschwrap.h \
y12wrap.cc \
y12wrap.h

libmbwrap_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbwrap_la_LDFLAGS =

INCLUDES = \
-I../../libraries/liby12 \
-I../../libraries/libharwell \
-I../../libraries/libmbmath \
-I../../libraries/libmbutil \
-I../../mbdyn

bin_PROGRAMS = t
t_SOURCES = test.cc
t_LDADD = \
./libmbwrap.la \
../libmbmath/libmbmath.la \
../libmbutil/libmbutil.la \
../libobjs/libobjs.la \
../liby12/liby12.la \
../libharwell/libharwell.la \
./libmbwrap.la \
@MESCHACH_LIBS@ \
@F2C_LIBS@ \
@TCL_LIBS@

INCLUDES = \
-I../../mbdyn \
-I../libmbutil \
-I../libmbmath \
-I../liby12 \
-I../libharwell

