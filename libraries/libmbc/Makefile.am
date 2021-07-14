include $(top_srcdir)/build/top.mk
## $Header$
## Process this file with automake to produce Makefile.in

# Build libmbc.a library
noinst_LTLIBRARIES = libmbc_static.la

### libmbc_static.la is linked by mbdyn; libmbc.la is installed
### NOTE: all mbdyn needs are the "sock" stuff
libmbc_static_la_SOURCES = \
mbc_dummy.c \
sock.c \
sock.h

libmbc_static_la_LIBADD = @LIBS@
libmbc_static_la_LDFLAGS =


lib_LTLIBRARIES = libmbc.la

libmbc_la_SOURCES = \
mbc.c \
mbc.h \
mbcxx.cc \
mbcxx.h

libmbc_la_LIBADD = libmbc_static.la @LIBS@
libmbc_la_LDFLAGS =

include_HEADERS = \
mbc.h \
mbcxx.h

AM_CPPFLAGS = \
-I../../include \
-I$(srcdir)/../../include \
-I$(srcdir)/../../libraries/libmbc \
-I$(srcdir)/../../mbdyn

EXTRA_DIST = mbc_py_interface.py

if USE_PYTHON
lib_LTLIBRARIES += _mbc_py.la
_mbc_py_la_SOURCES = \
mbc_py.i \
mbc_py.h \
mbc_py.cc \
mbc_py_global.i \
mbc_py_global.h

nodist__mbc_py_la_SOURCES = \
mbc_py_wrap.cc

SWIG_FLAGS = @SWIG_LIB@ $(AX_SWIG_PYTHON_OPT) $(AX_SWIG_PYTHON_CPPFLAGS)
AM_CPPFLAGS += $(PYTHON_CPPFLAGS) $(PYTHON_NUMPY_HEADERS)

mbc_py_wrap_dep = \
$(srcdir)/mbc_py.i \
$(srcdir)/mbc_py_global.i \
$(srcdir)/mbc_py_global.h
mbc_py_wrap.cc: $(mbc_py_wrap_dep)
	$(SWIG) $(SWIG_FLAGS) -o $@ $<

_mbc_py_la_LIBADD = @LIBS@ libmbc.la
_mbc_py_la_LDFLAGS = -module

install-exec-local:
	$(mkinstalldirs) $(libexecdir)/mbpy
	$(install_sh_PROGRAM) .libs/_mbc_py.so $(libexecdir)/mbpy/
	$(install_sh_DATA) mbc_py.py $(srcdir)/mbc_py_interface.py $(libexecdir)/mbpy/

# remove _mbc_py.* because not directly usable; _mbc_py.so already in $(libexecdir)/mbpy/
install-exec-hook:
	rm -f $(libdir)/_mbc_py.*

clean-local:
	rm -f mbc_py_wrap.cc mbc_py.py
endif

include $(top_srcdir)/build/bot.mk
