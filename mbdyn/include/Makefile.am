include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = mbconfig.h.in mbdefs.h.in

SUBDIRS = ac

all-local: mbdefs.h FORCE

mbdefs.h: mbdefs.h.in Makefile
	sed \
		-e "s;%MODULE_LOADPATH%;@MODULE_LOADPATH@;" \
		-e "s;%BINPATH%;${bindir};" \
		-e "s;%OCTAVEPATH%;${octavedir};" \
		$(srcdir)/$@.in > $@

FORCE:

include $(top_srcdir)/build/bot.mk
