include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libmbc.a library
#noinst_LTLIBRARIES = libmbc.la
lib_LTLIBRARIES = libmbc.la

libmbc_la_SOURCES = \
mbc.c \
mbc.h \
mbcxx.cc \
mbcxx.h \
sock.c \
sock.h

libmbc_la_LIBADD = @LIBS@
libmbc_la_LDFLAGS =

include_HEADERS = \
mbc.h \
mbcxx.h

INCLUDES = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../mbdyn

if USE_PYTHON
lib_LTLIBRARIES += libmbc_py.la
libmbc_py_la_SOURCES = \
mbc_py.i \
mbc_py.h \
mbc_py.cc \
mbc_py_global.i \
mbc_py_global.h \
mbc_py_wrap.cc

SWIG_FLAGS = @SWIG_LIB@ $(AX_SWIG_PYTHON_OPT) $(AX_SWIG_PYTHON_CPPFLAGS)
INCLUDES += $(PYTHON_CPPFLAGS) -I$(PYTHON_SITE_PKG)/numpy/core/include/numpy

mbc_py_wrap_dep = \
$(srcdir)/mbc_py.i \
$(srcdir)/mbc_py_global.i \
$(srcdir)/mbc_py_global.h
mbc_py_wrap.cc: $(mbc_py_wrap_dep)
	$(SWIG) $(SWIG_FLAGS) -o $@ $<

libmbc_py_la_LIBADD = @LIBS@ libmbc.la
libmbc_py_la_LDFLAGS =

install-exec-local:
	$(mkinstalldirs) $(libexecdir)/mbpy
	if test -f .libs/libmbc_py.so ; then \
		cp .libs/libmbc_py.so $(libexecdir)/mbpy/_mbc_py.so ; \
	fi
	$(install_sh_DATA) mbc_py.py $(srcdir)/mbc_py_interface.py $(libexecdir)/mbpy/

clean-local:
	rm -f mbc_py_wrap.cc mbc_py.py
endif

include $(top_srcdir)/build/bot.mk
