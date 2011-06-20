include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

export LTCXXCOMPILE = @LIBTOOL@ --mode=compile @CXX@ @DEFS@ \
	-I. -I@srcdir@ \
	-I../../include \
	-I../$(top_srcdir)/include \
	-I../$(top_srcdir)/libraries/libmbutil \
	-I../$(top_srcdir)/libraries/libmbmath \
	-I../$(top_srcdir)/libraries/libmbwrap \
	-I../$(top_srcdir)/libraries/libobjs \
	-I../$(top_srcdir)/mbdyn \
	-I../$(top_srcdir)/mbdyn/base \
	-I../$(top_srcdir)/mbdyn/struct \
	-I../$(top_srcdir)/mbdyn/aero \
	-I../$(top_srcdir)/mbdyn/elec \
	-I../$(top_srcdir)/mbdyn/hydr \
	-prefer-pic \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CXXFLAGS) @CXXFLAGS@
CXXLD = @CXX@
export CXXLINK = @LIBTOOL@ --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	@CXXFLAGS@ @LDFLAGS@ -shared @FCLIBS@
export LTCOMPILE = @LIBTOOL@ --mode=compile @CC@ @DEFS@ \
	-I. -I@srcdir@ \
	-I../../include \
	-I../$(top_srcdir)/include \
	-I../$(top_srcdir)/libraries/libmbutil \
	-I../$(top_srcdir)/libraries/libmbmath \
	-I../$(top_srcdir)/libraries/libmbwrap \
	-I../$(top_srcdir)/libraries/libobjs \
	-I../$(top_srcdir)/mbdyn \
	-I../$(top_srcdir)/mbdyn/base \
	-I../$(top_srcdir)/mbdyn/struct \
	-I../$(top_srcdir)/mbdyn/aero \
	-I../$(top_srcdir)/mbdyn/elec \
	-I../$(top_srcdir)/mbdyn/hydr \
	-prefer-pic \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CFLAGS) @CFLAGS@
export LTF77COMPILE = @LIBTOOL@ --mode=compile @F77@ $(AM_FFLAGS) @FFLAGS@ \
	-prefer-pic
export LTFCCOMPILE = @LIBTOOL@ --mode=compile @FC@ $(AM_FFLAGS) @FFLAGS@ \
	-prefer-pic
export LIBTOOL = @LIBTOOL@
export MODULES_DIRS=@MODULES_DIRS@
export INSTALL_PROGRAM=${INSTALL}
export LIBEXEC_DIR=$(libexecdir)

.PHONY: all all__ clean distclean

# Note: this hack is required because macro "LIBTOOL" may expand
# to sh ../libtool, so, when we change directory to the modules,
# the effective libtool is not found.

all: libtool all__

libtool: 
	-@ln -sf ../libtool .
	
all__:
	@echo "Configured modules: $(MODULES_DIRS)"
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Playing with $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(srcdir)/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f,f90,f95} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f|\*.f90|\*.f95) ;; *) if test ! -f $$i/$$f ; then ( cd $$i && ln -sf ../$$isrc/$$f $$f ) ; fi ; ;; esac ; \
			done ; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			( cd $$i && $(MAKE) MBDYN_MODULE_SRCDIR=../$$isrc ) ; \
			rm -f $$i/Makefile ; \
		done \
	fi

install: all__
	$(mkinstalldirs) $(libexecdir)
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Installing $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(srcdir)/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f,f90,f95} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f|\*.f90|\*.f95) ;; *) if test ! -f $$i/$$f ; then ( cd $$i && ln -sf ../$$isrc/$$f $$f ) ; fi ; ;; esac ; \
			done ; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			( cd $$i && $(MAKE) install ) ; \
			rm -f $$i/Makefile ; \
		done \
	fi

clean clean-so:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Cleaning $$i ..."; \
			cp $(srcdir)/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) $@ ; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi

distclean:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Cleaning $$i ..."; \
			rm -rf $$i/*.so $$i/*.lo $$i/*.o ; \
		done \
	fi
	rm -f Makefile libtool

# Add distributed modules below
EXTRA_DIST = Makefile.module Makefile.template \
 \
module-aerodyn/NREL_AeroDyn.h \
module-aerodyn/Makefile.inc \
module-aerodyn/mbdyn_ad.f90 \
module-aerodyn/module-aerodyn.cc \
module-aerodyn/README \
 \
module-asynchronous_machine/asynchronous_machine.mbdyn \
module-asynchronous_machine/module-asynchronous_machine.cc \
module-asynchronous_machine/module-asynchronous_machine.h \
module-asynchronous_machine/doc/Doxyfile \
module-asynchronous_machine/doc/Makefile \
 \
module-constlaw/module-constlaw.cc \
module-constlaw/README \
 \
module-constlaw-f90/Makefile.inc \
module-constlaw-f90/module-constlaw-f90.cc \
module-constlaw-f90/README \
module-constlaw-f90/usrsub.f90 \
module-constlaw-f90/usrsub.h \
 \
module-constlaw-f95/Makefile.inc \
module-constlaw-f95/module-constlaw-f95.cc \
module-constlaw-f95/README \
module-constlaw-f95/usrsub.f95 \
module-constlaw-f95/usrsub.h \
 \
module-drive/module-drive.cc \
module-drive/README \
\
module-hydrodynamic_plain_bearing/diffsizes.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force_dv.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force.f95 \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing_force.h \
module-hydrodynamic_plain_bearing/Makefile.inc \
module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.cc \
module-hydrodynamic_plain_bearing/module-hydrodynamic_plain_bearing.h \
module-hydrodynamic_plain_bearing/hydrodynamic_plain_bearing.mbdyn \
 \
module-mds/mds \
module-mds/module-mds.cc \
 \
module-scalarfunc/module-scalarfunc.cc \
module-scalarfunc/README \
 \
module-tclpgin/Makefile.inc \
module-tclpgin/module-tclpgin.cc \
 \
module-template/module-template.cc \
 \
module-template2/module-template2.cc \
 \
module-udunits/Makefile.inc \
module-udunits/module-udunits.cc \
module-udunits/README \
 \
module-wheel2/module-wheel2.cc \
module-wheel2/module-wheel2.h

include $(top_srcdir)/build/bot.mk
