include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = \
f2c-int.h \
f2c.h \
float.h \
fstream \
getopt.h \
iomanip \
iostream \
lapack.h \
math.h \
pthread.h \
spinlock.h \
sstream \
sys_sysinfo.h

noinst_PROGRAMS = integertest
integertest_SOURCES = integertest.cc
INCLUDES = -I../../mbdyn

export SH_CXXCOMPILE=$(CXXCOMPILE)
export SH_CXXLD=$(CXXLD)

all-local: integertest.cc FORCE
	-@$(SHELL) integertest.sh

FORCE:

include $(top_srcdir)/build/bot.mk
