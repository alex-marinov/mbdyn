## Process this file with automake to produce Makefile.in

# Build libmbmath.a library
lib_LTLIBRARIES = libmbmath.la

libmbmath_la_SOURCES = \
fullmh.cc    fullmh.h     harwlib.h    harwlib_.f      spdata_.f \
harwout.c    harwout.h    harwout_.cc  harwrap.cc      harwrap.h   ldl.f \
ldl.h        matvec3.cc   matvec3.h    matvec3n.cc     matvec3n.h  matvec6.cc \
matvec6.h    mschwrap.cc  mschwrap.h   myf2c.h         mymath.h    solman.cc \
solman.h     submat.cc    submat.h     schsolman.cc    schsolman.h

libmbmath_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbmath_la_LDFLAGS =

INCLUDES = -I../../libraries/libmbutil -I../../libraries/libmbmath -I../../mbdyn/base -I../../mbdyn -I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec -I../../mbdyn/hydr

