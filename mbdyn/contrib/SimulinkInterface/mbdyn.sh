#! /bin/sh
# $Header$

# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2011
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

# This shell script is used by the simulink interface to execute MBDyn
# with an appropriately fed environment.  Customize as required by your
# installation of MBDyn and Simulink.  Please report errors and issues
# to mbdyn-users@mbdyn.org

# TODO: get values from configure...
INSTALLDIR=/usr/local/mbdyn
#INSTALLDIR=/usr/local
BINDIR="${INSTALLDIR}/bin"

# This is the location of MBDyn when run from within the build tree
#MBDYN=../../mbdyn/mbdyn

# This is the location of MBDyn when installed in the default location
MBDYN="${BINDIR}/mbdyn"

# Add further variable declarations
# for example, if MBDyn is linked with a version of libstdc++ that differs
# from the one distributed with matlab, put its location into LD_LIBRARY_PATH
# in most cases:
LD_LIBRARY_PATH=/usr/lib
# in most cases on x86_64:
#LD_LIBRARY_PATH=/usr/lib64
# but, usually, it's enough to let the system do its job!
#LD_LIBRARY_PATH=

exec "${MBDYN}" $@
#exec -a mbdyn "${MBDYN}" $@
# log output in 'log.txt' (enable "verbose" from interface)
#exec "${MBDYN}" $@ > log.txt 2>&1

