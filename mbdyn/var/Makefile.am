include $(top_srcdir)/build/top.mk
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = \
abs2rel.awk \
mbdyn2pod.awk \
mbdyn2van.awk \
mbdyn2vol.awk \
periodize.awk \
naca0012.c81 \
mbdyn2easyanim.sh.in \
pod.m

sharedir = "$(prefix)/share"
octavedir = "$(sharedir)/octave"
c81dir = "$(sharedir)/c81"
awkdir = "$(sharedir)/awk"

generated = mbdyn2easyanim.sh

all-local: $(generated)

$(generated):
	@sed -e "s;@AWKDIR@;$(awkdir);" $@.in > $@

install-data-local: $(EXTRA_DIST) $(generated)
	$(mkinstalldirs) $(DESTDIR)$(sharedir) $(DESTDIR)$(awkdir) $(DESTDIR)$(c81dir) $(DESTDIR)$(octavedir) $(DESTDIR)$(bindir)
	@for f in $^ ; do \
		case $$f in \
		*.awk) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(awkdir) ;; \
		*.c81) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(c81dir) ;; \
		*.m) $(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(octavedir) ;; \
		*.sh) $(INSTALL_SCRIPT) -m 555 $$f $(DESTDIR)$(bindir) ;; \
		*.in) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

distclean-local:
	rm -f $(generated)

include $(top_srcdir)/build/bot.mk

