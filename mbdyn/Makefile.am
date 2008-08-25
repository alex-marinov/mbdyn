include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

codeasterdir = "$(sharedir)/Aster"
codeasterfiles = \
contrib/CodeAster/cms/cms.py

EXTRA_DIST = CHANGES \
contrib/SimulinkInterface/Makefile \
contrib/SimulinkInterface/MBDynCom.mdl \
contrib/SimulinkInterface/README \
contrib/SimulinkInterface/setup.m \
contrib/SimulinkInterface/sfun_mbdyn_com_read.c \
contrib/SimulinkInterface/sfun_mbdyn_com_write.c \
contrib/SimulinkInterface/sfun_mbdyn_start.c \
	contrib/SimulinkInterface/examples/pendulum \
	contrib/SimulinkInterface/examples/pendulum.mdl \
	contrib/SimulinkInterface/examples/README \
	contrib/SimulinkInterface/tmp/README \
contrib/Patches \
$(codeasterfiles)

install-data-local-codeaster: $(codeasterfiles)
	$(mkinstalldirs) $(DESTDIR)$(sharedir) $(DESTDIR)$(codeasterdir)
	@for f in $^ ; do \
		case $$f in \
		*.py) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(codeasterdir)" ; $(INSTALL_SCRIPT) -m 444 $(srcdir)/$$f $(DESTDIR)$(codeasterdir) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

install-data-local: install-data-local-codeaster

dist-hook:
	rm -rf `find $(distdir) -name CVS`



SUBDIRS = include build etc libraries manual mbdyn modules utils var

include $(top_srcdir)/build/bot.mk

