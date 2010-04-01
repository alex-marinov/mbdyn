include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

SUBDIRS = \
libmbc \
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

include $(top_srcdir)/build/bot.mk
