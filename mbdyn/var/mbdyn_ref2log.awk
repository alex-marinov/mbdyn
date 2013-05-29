# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2013
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

# usage: awk -f mbdyn_ref2log.awk [-v rot={euler123|euler313|euler321|phi|mat}] INPUT_FILE
#
# INPUT_FILE is a .rfm file

BEGIN {
	if (rot == "") {
		rot = "euler123";

	} else if (rot != "euler123" && rot != "euler313" && rot != "euler321" && rot != "phi" && rot != "mat") {
		exit 1;
	}
}

{
	printf "structural node: %s %s %s %s %s %s %s %s", $1, $2, $3, $4, rot, $5, $6, $7;
	if (rot == "mat") {
		printf " %s %s %s %s %s %s", $8, $9, $10, $11, $12, $13;
	}
	print "";
}

