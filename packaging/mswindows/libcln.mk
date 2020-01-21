PKG             := libcln
$(PKG)_WEBSITE  := https://www.ginac.de/CLN
$(PKG)_DESCR    := CLN - Class Library for Numbers
$(PKG)_IGNORE   :=

$(PKG)_VERSION  := 0023432
$(PKG)_GH_CONF  := jrheinlaender/cln/branches/w64
$(PKG)_CHECKSUM := d814cd9cb9ff14f97e09d2766c5b25a2c680b01d0dee5da0c7bc51bf54641aea

#$(PKG)_VERSION  := 0e72f1f
#$(PKG)_CHECKSUM := e94fc1fd4a05be12e3fb25241260a7226bd5f2c423c80e4faa70ab80e9591431
#$(PKG)_GH_CONF  := jrheinlaender/cln/branches/w64

#$(PKG)_VERSION  := 1.3.6
#$(PKG)_CHECKSUM := f492530e8879bda529009b6033e1923c8f4aae843149fc28c667c20b094d984a
#$(PKG)_SUBDIR   := cln-$($(PKG)_VERSION)
#$(PKG)_FILE     := cln-$($(PKG)_VERSION).tar.bz2
#$(PKG)_URL      := $($(PKG)_WEBSITE)/$($(PKG)_FILE)

$(PKG)_DEPS     := cc gmp libiconv

#define $(PKG)_UPDATE
#    echo 'TODO: write update script for $(PKG).' >&2;
#    echo $($(PKG)_VERSION)
#endef

define $(PKG)_BUILD_STATIC
    cd '$(SOURCE_DIR)' && autoreconf -fi
#CXXFLAGS='$$CXXFLAGS -std=c++11'
    cd '$(BUILD_DIR)' && '$(SOURCE_DIR)/configure' \
        $(MXE_CONFIGURE_OPTS) ac_cv_c_bigendian=no
    cd '$(BUILD_DIR)' &&  $(MAKE) -j '$(JOBS)'
    cd '$(BUILD_DIR)' && $(MAKE) -j 1 install
endef

define $(PKG)_BUILD_SHARED
    cd '$(SOURCE_DIR)' && autoreconf -fi
#CXXFLAGS='$$CXXFLAGS -std=c++11'
    cd '$(BUILD_DIR)' && '$(SOURCE_DIR)/configure' \
        $(MXE_CONFIGURE_OPTS) ac_cv_c_bigendian=no
    cd '$(BUILD_DIR)' &&  $(MAKE) -j '$(JOBS)' LDFLAGS='-no-undefined'
    cd '$(BUILD_DIR)' && $(MAKE) -j 1 install
endef
