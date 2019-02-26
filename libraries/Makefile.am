include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

SUBDIRS = \
libmbutil \
libobjs \
libcolamd \
libnaive

if BUILD_Y12
SUBDIRS += liby12
endif

SUBDIRS += \
libmbmath \
libmbwrap \
libann

if BUILD_MBC
SUBDIRS += libmbc
endif

include $(top_srcdir)/build/bot.mk
