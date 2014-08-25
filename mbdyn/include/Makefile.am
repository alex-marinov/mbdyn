include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = mbconfig.h.in mbdefs.h.in

SUBDIRS = ac

all-local: mbdefs.h FORCE

if USE_OCTAVE
octbindir="${bindir}/octave-`octave-config --version`"
else
octbindir="${bindir}"
endif

mbdefs.h: mbdefs.h.in Makefile
	sed \
		-e "s;%MODULE_LOADPATH%;@MODULE_LOADPATH@;" \
		-e "s;%BINPATH%;${bindir};" \
		-e "s;%OCTAVEPATH%;${octavedir};" \
		-e "s;%OCTAVEBINPATH%;${octbindir};" \
		$(srcdir)/$@.in > $@

uninstall-local: mbdefs.h
	rm -f mbdefs.h

FORCE:

include $(top_srcdir)/build/bot.mk
