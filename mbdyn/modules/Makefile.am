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
	-fPIC \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CXXFLAGS) @CXXFLAGS@
CXXLD = @CXX@
export CXXLINK = @LIBTOOL@ --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	@CXXFLAGS@ @LDFLAGS@ -shared @F2C_LIBS@
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
	-fPIC \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CFLAGS) @CFLAGS@
export LTF77COMPILE = @LIBTOOL@ --mode=compile @F77@ $(AM_FFLAGS) @FFLAGS@ \
	-fPIC
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
			isrc=$(top_srcdir)/modules/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f) ;; *) if test ! -f $$i/$$f ; then ln -sf ../$$isrc/$$f $$i/$$f ; fi ; ;; esac ; \
			done ; \
			cp $(top_srcdir)/modules/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) ; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi

install:
	$(mkinstalldirs) $(libexecdir)
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Installing $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(top_srcdir)/modules/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f) ;; *) if test ! -f $$i/$$f ; then ln -sf ../$$isrc/$$f $$i/$$f ; fi ; ;; esac ; \
			done ; \
			cp $(top_srcdir)/modules/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) install; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi

clean:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Cleaning $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(top_srcdir)/modules/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f) ;; *) if test ! -f $$i/$$f ; then ln -sf ../$$isrc/$$f $$i/$$f ; fi ; ;; esac ; \
			done ; \
			cp $(top_srcdir)/modules/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) clean ; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi

distclean:
	-@if test x"$(MODULES_DIRS)" != x ; then \
		for i in $(MODULES_DIRS) ; do \
			echo "Dist-cleaning $$i ..."; \
			export MBDYN_MODULE=$$i ; \
			isrc=$(top_srcdir)/modules/$$i ; \
			if test -d $$isrc ; then mkdir -p $$i ; fi ; \
			for fsrc in $$isrc/*.{c,cc,f} ; do \
				f=`basename $$fsrc` ; \
				case $$f in \*.c|\*.cc|\*.f) ;; *) if test ! -f $$i/$$f ; then ln -sf ../$$isrc/$$f $$i/$$f ; fi ; ;; esac ; \
			done ; \
			cp $(top_srcdir)/modules/Makefile.module $$i/Makefile ; \
			cd $$i && $(MAKE) distclean ; \
			cd .. ; \
			rm -f $$i/Makefile ; \
		done \
	fi
	rm -f Makefile libtool

# Add dsitributed modules below
EXTRA_DIST = Makefile.module Makefile.template \
module-template/module-template.cc \
module-wheel2/module-wheel2.cc

include $(top_srcdir)/build/bot.mk
