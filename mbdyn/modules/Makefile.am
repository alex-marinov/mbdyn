## Process this file with automake to produce Makefile.in
## Do not forget to run automake ( --add-missing ) when you modify Makefile.am

export LTCXXCOMPILE = @LIBTOOL@ --mode=compile @CXX@ @DEFS@ \
	-I. -I@srcdir@ -I../../include \
	-I../../libraries/libmbutil -I../../libraries/libmbmath \
	-I../../libraries/libobjs -I../../mbdyn -I../../mbdyn/base \
	-I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec \
	-I../../mbdyn/hydr -fPIC \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CXXFLAGS) @CXXFLAGS@
CXXLD = @CXX@
export CXXLINK = @LIBTOOL@ --mode=link $(CXXLD) $(AM_CXXFLAGS) \
	@CXXFLAGS@ @LDFLAGS@ -shared
export LTF77COMPILE = @LIBTOOL@ --mode=compile @F77@ $(AM_FFLAGS) @FFLAGS@ \
	-fPIC


.PHONY: all clean distclean

all:
	-for i in @MODULES_DIRS@ ; do \
		export MBDYN_MODULE=$$i ; \
		cp Makefile.module $$i/Makefile ; \
		cd $$i && $(MAKE) ; \
		cd .. ; \
		rm -f $$i/Makefile ; \
	done

clean:
	-for i in @MODULES_DIRS@ ; do \
		export MBDYN_MODULE=$$i ; \
		cp Makefile.module $$i/Makefile ; \
		cd $$i && $(MAKE) clean ; \
		cd .. ; \
		rm -f $$i/Makefile ; \
	done

distclean: clean
	rm -f Makefile

