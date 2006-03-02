include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libmbmath.a library
lib_LTLIBRARIES = libmbmath.la

libmbmath_la_SOURCES = \
ccmh.cc \
ccmh.h \
dirccmh.cc \
dirccmh.h \
fullmh.cc \
fullmh.h \
interp.h \
interp.cc \
JacSubMatrix.cc \
JacSubMatrix.h \
ldl.f \
ldl.h \
ls.cc \
ls.h \
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
naivemh.cc \
naivemh.h \
Rot.cc \
Rot.hh \
RotCoeff.hc \
RotCoeff.hh \
schsolman.cc \
schsolman.h \
schurmh.h \
schurmh.cc \
solman.cc \
solman.h \
spdata.h \
spdata_.f \
spmapmh.cc \
spmapmh.h \
spmh.cc \
spmh.h \
submat.cc \
submat.h \
vh.cc \
vh.h

libmbmath_la_LIBADD = @LIBS@ @F2C_LIBS@
libmbmath_la_LDFLAGS =

INCLUDES = \
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
