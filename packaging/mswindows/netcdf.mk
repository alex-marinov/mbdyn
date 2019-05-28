# This file is part of MXE. See LICENSE.md for licensing information.

PKG             := netcdf
$(PKG)_WEBSITE  := https://www.unidata.ucar.edu/software/netcdf/
$(PKG)_DESCR    := NetCDF
$(PKG)_IGNORE   :=
# $(PKG)_VERSION  := 4.3.0
# 4.3.0 $(PKG)_CHECKSUM := e796413d27da6b053e07a18f567a1d0c23d2a317cef905faa2a05fe4f725fc63
#$(PKG)_VERSION  := 4.1.2
# 4.1.2 $(PKG)_CHECKSUM := 0c8df55f5967be2cd6ede1a1e8a4a96a69d6b716ec649b3eba9640dbbcda16b9
$(PKG)_VERSION  := 4.1.3
$(PKG)_CHECKSUM := 315bc385b0750dd85b9a122194382db52f432dea1390de9f7afa11cf91869213
$(PKG)_SUBDIR   := $(PKG)-$($(PKG)_VERSION)
$(PKG)_FILE     := $(PKG)-$($(PKG)_VERSION).tar.gz
#$(PKG)_URL      := ftp://ftp.unidata.ucar.edu/pub/netcdf/old/$($(PKG)_FILE)
$(PKG)_URL      := ftp://ftp.unidata.ucar.edu/pub/netcdf/$($(PKG)_FILE)
$(PKG)_DEPS     := cc curl hdf4 hdf5 portablexdr zlib

#define $(PKG)_UPDATE
#    $(WGET) -q -O- 'https://www.unidata.ucar.edu/downloads/netcdf/index.jsp' | \
#    grep netcdf | \
#    $(SED) -n 's,.*href="netcdf-\([0-9_]*\)">.*,\1,p' | \
#    head -1 | \
#    tr '_' '.'
#endef

# NetCDF uses '#ifdef IGNORE' as a synonym to '#if 0' in several places.
# IGNORE is assumed to never be defined, but winbase.h defines it...
# We just replace '#ifdef IGNORE' with '#if 0' to work around this.

# $(SED) -i -e 's/#ifdef IGNORE/#if 0/' libsrc4/nc4hdf.c libsrc4/ncfunc.c libsrc/attr.c ncgen/cvt.c && \

define $(PKG)_BUILD
    cd '$(1)' && \
        $(SED) -i -e 's/#ifdef IGNORE/#if 0/' nc_test4/tst_files.c \
                                                libcdmr/ast_runtime.c \
                                                libcdmr/ast_internal.c \
                                                libcdmr/cceparselex.h \
                                                libcdmr/ncaux.c \
                                                libcdmr/crce.c \
                                                libcdmr/ast_internal.h \
                                                libcdmr/nccrgetvara.c \
                                                ncgen/genbin.c \
                                                ncgen/genchar.c \
                                                ncgen/genj.c \
                                                ncgen/semantics.c \
                                                libsrc/attr.c \
                                                libsrc/attr.m4 \
                                                oc/curlfunctions.c \
                                                oc/ocuri.c \
                                                oc/oc.c \
                                                oc/ocxdr_stdio.c \
                                                oc/occlientparams.c \
                                                libdispatch/nclog.c \
                                                libdispatch/nc_uri.c \
                                                libsrc4/nc4hdf.c \
                                                libsrc4/nc4file.c \
                                                libsrc4/ncfunc.c \
                                                libdap2/dapattr3.c \
                                                libdap2/cdf4.c \
                                                libdap2/getvara4.c \
                                                libdap2/dapalign.c \
                                                libdap2/dapodom.c \
                                                libdap2/common34.c \
                                                libdap2/ncdap3.h \
                                                libdap2/daputil.h \
                                                libdap2/daputil.c \
                                                libdap2/getvara3.c \
                                                libdap2/ncdap3a.c \
                                                libdap2/constraints4.c \
                                                libdap2/constraints3.c \
                                                libdap2/ncdap4.c \
                                                libdap2/cdf3.c \
        && \
        ./configure \
            $(MXE_CONFIGURE_OPTS) \
            --enable-netcdf-4 \
            --enable-hdf4 \
            --with-hdf5 \
            --with-zlib \
            --with-cxx-compat \
            --disable-testsets \
            --disable-examples \
            CPPFLAGS="-D_DLGS_H -DWIN32_LEAN_AND_MEAN" \
            LIBS="-lmfhdf -ldf -lportablexdr -lws2_32"

    $(MAKE) -C '$(1)' -j '$(JOBS)' LDFLAGS=-no-undefined

    $(MAKE) -C '$(1)' -j 1 install
endef
