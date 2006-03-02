include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in

# Build liby12.a library
lib_LTLIBRARIES = liby12.la

liby12_la_SOURCES = \
y12lib.h \
y12mae.f \
y12maf.f \
y12mbe.f \
y12mbf.f \
y12mce.f \
y12mcf.f \
y12mde.f \
y12mdf.f \
y12mfe.f \
y12mge.f \
y12mhe.f

liby12_la_LIBADD = @LIBS@ @F2C_LIBS@
liby12_la_LDFLAGS =

noinst_PROGRAMS = y12test
y12test_SOURCES = test.c
y12test_LDADD = ./liby12.la
INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../libmbmath \
-I$(srcdir)/../../mbdyn

include $(top_srcdir)/build/bot.mk
