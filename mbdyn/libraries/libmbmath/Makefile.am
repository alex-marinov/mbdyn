## Process this file with automake to produce Makefile.in

# Build libmbmath.a library
lib_LTLIBRARIES = libmbmath.la

libmbmath_la_SOURCES = \
fullmh.cc \
fullmh.h \
harwout.c \
harwout.h \
harwout_.cc \
harwrap.cc \
harwrap.h \
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
mschwrap.cc \
mschwrap.h \
myf2c.h \
mymath.h \
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
submat.cc \
submat.h \
y12wrap.cc \
y12wrap.h

libmbmath_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbmath_la_LDFLAGS =

INCLUDES = \
-I../../libraries/liby12 \
-I../../libraries/libharwell \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../mbdyn
