include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build libmbutil.a library
lib_LTLIBRARIES = libmbutil.la

libmbutil_la_SOURCES = \
destr.h \
dummypgin.cc \
dummypgin.h \
except.h \
filename.h \
fn_UNIX.cc \
gauss.cc \
gauss.h \
input.cc \
input.h \
llist.cc \
llist.h \
llstack.cc \
llstack.h \
mathp.cc \
mathp.h \
mathtyp.h \
mbsasl.c \
mbsasl.h \
myassert.cc \
myassert.h \
mynewmem.cc \
mynewmem.h \
mysleep.c\
mysleep.h\
parser.cc \
parser.h \
parsinc.cc \
parsinc.h \
search.h \
stack.h \
table.cc \
table.h \
tclpgin.cc \
tclpgin.h \
veciter.h \
withlab.h

libmbutil_la_LIBADD = @LIBS@
libmbutil_la_LDFLAGS =

INCLUDES = \
-I../../include \
-I../../libraries/libmbutil \
-I../../libraries/libmbmath \
-I../../mbdyn

noinst_PROGRAMS = mbsasltest
mbsasltest_SOURCES = mbsasl.c mbsasl.h mbsasltest.c
mbsasltest_LDADD = @THREAD_LIBS@ @SECURITY_LIBS@

include $(top_srcdir)/build/bot.mk
