#!/bin/sh

# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2004
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

# on good compilers, compilation should fail
${SH_CXXCOMPILE} -DINTEGER_TEST -o integertest.o integertest.cc >/dev/null 2>&1
if test "$?" != 0 ; then
	# echo "integertest: sizeof(integer) == sizeof(int)"
	exit 0
fi

# if the compilation succeeded, linking should be straightforward
${SH_CXXLD} -o integertest integertest.o >/dev/null 2>&1
if test "$?" != 0 ; then
	echo "         ###################################"
	echo "         ### unable to link integertest! ###"
	echo "         ###################################"
	exit 0
fi

# in any case, if execution fails they differ
exec ./integertest >/dev/null 2>&1
if test "$?" != 0 ; then
	echo "         ###############################################"
	echo "         ### WARNING: sizeof(integer) != sizeof(int) ###"
	echo "         ###############################################"
	exit 0
fi

exit 0
