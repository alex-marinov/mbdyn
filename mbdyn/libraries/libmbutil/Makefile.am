## Process this file with automake to produce Makefile.in

# Build libmbutil.a library
lib_LTLIBRARIES = libmbutil.la
libmbutil_la_SOURCES = \
destr.h     except.h   filename.h  gauss.cc     gauss.h \
input.cc    input.h    llist.cc    llist.h      llstack.cc  llstack.h \
mathp.cc    mathp.h    mathtyp.h   myassert.cc  myassert.h  mynewmem.cc \
mynewmem.h  mystddef.h parser.cc   parser.h     parsinc.cc  parsinc.h \
search.h    stack.h    table.cc    table.h \
veciter.h   withlab.h  fn_UNIX.cc

libmbutil_la_LIBADD = @LIBS@
libmbutil_la_LDFLAGS =

INCLUDES = -I../../include -I../../libraries/libmbutil -I../../libraries/libmbmath -I../../mbdyn -I../../mbdyn/base -I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec -I../../mbdyn/hydr

