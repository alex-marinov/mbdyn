include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

octavedir = "$(datarootdir)/octave"
c81dir = "$(datarootdir)/c81"
awkdir = "$(datarootdir)/awk"

export octavedir c81dir awkdir

codeaster_dir = "$(datarootdir)/Aster"
codeaster_srcdir = contrib/CodeAster
codeaster_files = \
$(codeaster_srcdir)/cms/cms.py

scicos_dir = "$(datarootdir)/Scicos"
scicos_srcdir = contrib/ScicosInterface
scicos_files = \
$(scicos_srcdir)/buildblocks.sci \
$(scicos_srcdir)/loadblocks.sci \
$(scicos_srcdir)/sockread_intf.sci \
$(scicos_srcdir)/sockwrite_intf.sci \
$(scicos_srcdir)/sockread_comp.c \
$(scicos_srcdir)/sockwrite_comp.c \
$(scicos_srcdir)/README \
$(scicos_srcdir)/TODO \
$(scicos_srcdir)/examples/pendulum \
$(scicos_srcdir)/examples/Pendulum.cosf

simulink_dir = "$(datarootdir)/Simulink"
simulink_srcdir = contrib/SimulinkInterface
simulink_files = \
$(simulink_srcdir)/Makefile \
$(simulink_srcdir)/MBDynCom.mdl \
$(simulink_srcdir)/mbdyn.sh \
$(simulink_srcdir)/README \
$(simulink_srcdir)/setup.m \
$(simulink_srcdir)/sfun_mbdyn_com_read.c \
$(simulink_srcdir)/sfun_mbdyn_com_write.c \
$(simulink_srcdir)/sfun_mbdyn_start.c \
	$(simulink_srcdir)/examples/pendulum \
	$(simulink_srcdir)/examples/pendulum.mdl \
	$(simulink_srcdir)/examples/README \
	$(simulink_srcdir)/tmp/README

EXTRA_DIST = CHANGES \
contrib/Patches \
$(codeaster_files) \
$(scicos_files) \
$(simulink_files)

install-data-local-codeaster: $(codeaster_files)
	$(mkinstalldirs) $(DESTDIR)$(datarootdir) $(DESTDIR)$(codeaster_dir)
	@for f in $^ ; do \
		case $$f in \
		*.py) echo "$(INSTALL_SCRIPT) -m 444 $$f $(DESTDIR)$(codeaster_dir)" ; $(INSTALL_SCRIPT) -m 444 $(srcdir)/$$f $(DESTDIR)$(codeaster_dir) ;; \
		*) echo "($$f ignored)" ;; \
		esac ; \
	done

install-data-local: install-data-local-codeaster

dist-hook:
	rm -rf `find $(distdir) -name CVS`



SUBDIRS = include build etc libraries manual mbdyn modules utils var

include $(top_srcdir)/build/bot.mk

