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

