## Process this file with automake to produce Makefile.in

# Build liby12.a library
lib_LTLIBRARIES = liby12.la

liby12_la_SOURCES = \
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

bin_PROGRAMS = t
t_SOURCES = test.c
t_LDADD = ./liby12.la
INCLUDES = \
-I../libmbmath \
-I../../mbdyn


