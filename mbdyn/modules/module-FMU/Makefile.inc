FMIL = /home/marco/Mbdyn/FMILibrary-2.0.1/Build
#FMIL = /home/marco/Mbdyn/FMILibrary-2.0.1-install
FMII = /home/marco/Mbdyn/FMILibrary-2.0.1-install

MODULE_INCLUDE = -I$(FMII)/include \
        -I$(FMIL)/../ThirdParty/Minizip/minizip \
        -I$(FMIL)/../src/Util/include \
        -I$(FMIL) \
        -I$(FMIL)/../src/Import/include/  \
	-I$(FMIL)/../src/Import/src/FMI1/ \
	-I$(FMIL)/../src/Import/src/FMI2/ \
	-I$(FMIL)/../src/XML/include/ \
	-I$(FMIL)/../src/CAPI/src/ \
	-I$(FMIL)/../src/CAPI/include

MODULE_LINK =  -L$(FMIL) -lfmizip -ljmutils -lfmilib\
          -L$(FMIL)/zlib/ -lz \
          -L$(FMIL)/Config.cmake/Minizip/ -lminizip  \
          -L$(FMIL)/ThirdParty/c99_snprintf/ -lc99snprintf \
	   $(FMIL)/libfmilib_shared.so 

MODULE_DEPENDENCIES = mbdynFMI.lo \
       mbdynFMI_config.lo\


