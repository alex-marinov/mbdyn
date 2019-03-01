# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2017
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
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


