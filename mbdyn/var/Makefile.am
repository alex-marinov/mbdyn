include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

EXTRA_DIST = \
abs2rel.awk \
mbdyn2easyanim.awk \
mbdyn2pod.awk \
periodize.awk \
 \
naca0012.c81 \
 \
mbdyn2easyanim.sh.in \
 \
gradnorm.m \
pod.m \
rotation_increment_gibbs_rodriguez.m \
rotation_perturbation_matrix_gibbs_rodriguez.m \
mbdyn_derivative.m \
 \
mbd.vim

generated = mbdyn2easyanim.sh

all-local: $(generated)

$(generated): FORCE
	@echo "generating $(generated)..."
	@for g in $(generated) ; do \
		sed \
			-e "s;@AWKDIR@;$(awkdir);" \
			$(top_srcdir)/var/$$g.in > $$g ; \
	done

install-data-local: $(EXTRA_DIST) $(generated)
	$(mkinstalldirs) $(DESTDIR)$(datarootdir) $(DESTDIR)$(awkdir) $(DESTDIR)$(c81dir) $(DESTDIR)$(octavedir) $(DESTDIR)$(bindir)
	@for f in $^ ; do \
		case $$f in \
		*.awk) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(awkdir)" ; $(INSTALL_SCRIPT) -m 444 $(top_srcdir)/var/$$f $(DESTDIR)$(awkdir) ;; \
		*.c81) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(c81dir)" ; $(INSTALL_SCRIPT) -m 444 $(top_srcdir)/var/$$f $(DESTDIR)$(c81dir) ;; \
		*.m) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(octavedir)" ; $(INSTALL_SCRIPT) -m 444 $(top_srcdir)/var/$$f $(DESTDIR)$(octavedir) ;; \
		*.sh) ;; \
		*.in) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

install-exec-local: $(EXTRA_DIST) $(generated)
	$(mkinstalldirs) $(DESTDIR)$(datarootdir) $(DESTDIR)$(awkdir) $(DESTDIR)$(c81dir) $(DESTDIR)$(octavedir) $(DESTDIR)$(bindir)
	@for f in $^ ; do \
		case $$f in \
		*.awk) ;; \
		*.c81) ;; \
		*.m) ;; \
		*.sh) echo "$(INSTALL_SCRIPT) -m 555 $$f $(DESTDIR)$(bindir)" ; $(INSTALL_SCRIPT) -m 555 $$f $(DESTDIR)$(bindir) ;; \
		*.in) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

# uninstall-local: $(EXTRA_DIST) $(generated)
uninstall-local:
	@for f in $(EXTRA_DIST) $(generated) ; do \
		case $$f in \
		*.awk) rm -f $(DESTDIR)$(awkdir)/$$f ;; \
		*.c81) rm -f $(DESTDIR)$(c81dir)/$$f ;; \
		*.m) rm -f $(DESTDIR)$(octavedir)/$$f ;; \
		*.sh) rm -f $(DESTDIR)$(bindir)/$$f ;; \
		*.in) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

distclean-local:
	rm -f $(generated)

.PHONY: FORCE

include $(top_srcdir)/build/bot.mk

