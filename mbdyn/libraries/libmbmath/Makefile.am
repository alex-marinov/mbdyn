## Process this file with automake to produce Makefile.in

# Build libmbmath.a library
lib_LTLIBRARIES = libmbmath.la

libmbmath_la_SOURCES = \
fullmh.cc \
fullmh.h \
ldl.f \
ldl.h \
matvec3.cc \
matvec3.h \
matvec3n.cc \
matvec3n.h \
matvec6.cc \
matvec6.h \
matvecexp.cc \
matvecexp.h \
Rot.cc \
Rot.hh \
RotCoeff.hc \
RotCoeff.hh \
schsolman.cc \
schsolman.h \
solman.cc \
solman.h \
sparsemh.cc \
sparsemh.h \
spdata_.f \
spmapmh.h \
submat.cc \
submat.h

libmbmath_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbmath_la_LDFLAGS =

INCLUDES = \
-I../../libraries/liby12 \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../libraries/libmbwrap \
-I../../mbdyn

