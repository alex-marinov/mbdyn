include $(top_srcdir)/build/top.mk
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
	@CXXFLAGS@ @LDFLAGS@ -shared -lg2c
export LTCOMPILE = @LIBTOOL@ --mode=compile @CC@ @DEFS@ \
	-I. -I@srcdir@ -I../../include \
	-I../../libraries/libmbutil -I../../libraries/libmbmath \
	-I../../libraries/libobjs -I../../mbdyn -I../../mbdyn/base \
	-I../../mbdyn/struct -I../../mbdyn/aero -I../../mbdyn/elec \
	-I../../mbdyn/hydr -fPIC \
	$(AM_CPPFLAGS) @CPPFLAGS@ $(AM_CFLAGS) @CFLAGS@
export LTF77COMPILE = @LIBTOOL@ --mode=compile @F77@ $(AM_FFLAGS) @FFLAGS@ \
	-fPIC
export MODULES_DIRS=@MODULES_DIRS@

.PHONY: all all__ clean distclean

# Note: this hack is required because macro "LIBTOOL" may expand
# to sh ../libtool, so, when we change directory to the modules,
# the effective libtool is not found.

all: all__ libtool

libtool: 
	-@ln -sf ../libtool .
	
all__:
	-@for i in $(MODULES_DIRS) ; do \
		echo "Playing with $$i ..."; \
		export MBDYN_MODULE=$$i ; \
		cp Makefile.module $$i/Makefile ; \
		cd $$i && $(MAKE) ; \
		cd .. ; \
		rm -f $$i/Makefile ; \
	done

clean:
	-@for i in $(MODULES_DIRS) ; do \
		echo "Cleaning $$i ..."; \
		export MBDYN_MODULE=$$i ; \
		cp Makefile.module $$i/Makefile ; \
		cd $$i && $(MAKE) clean ; \
		cd .. ; \
		rm -f $$i/Makefile ; \
	done

distclean:
	-@for i in $(MODULES_DIRS) ; do \
		echo "Dist-cleaning $$i ..."; \
		export MBDYN_MODULE=$$i ; \
		cp Makefile.module $$i/Makefile ; \
		cd $$i && $(MAKE) distclean ; \
		cd .. ; \
		rm -f $$i/Makefile ; \
	done
	rm -f Makefile libtool

include $(top_srcdir)/build/bot.mk
