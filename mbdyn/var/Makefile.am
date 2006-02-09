include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = \
naca0012.c81 \
abs2rel.awk \
mbdyn2pod.awk \
mbdyn2van.awk \
mbdyn2vol.awk \
periodize.awk \
pod.m

sharedir = "$(prefix)/share"
octavedir = "$(sharedir)/octave"
c81dir = "$(sharedir)/c81"
awkdir = "$(sharedir)/awk"

install-data-local: $(EXTRA_DIST)
	$(mkinstalldirs) $(DESTDIR)$(sharedir) $(DESTDIR)$(awkdir) $(DESTDIR)$(c81dir) $(DESTDIR)$(octavedir) $(DESTDIR)$(bindir)
	@for f in $(EXTRA_DIST) ; do \
		case $$f in \
		*.awk) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(awkdir) ;; \
		*.c81) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(c81dir) ;; \
		*.m) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(octavedir) ;; \
		*.sh) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(bindir) ;; \
		*.in) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

include $(top_srcdir)/build/bot.mk

