# MBDyn (C) is a multibody analysis code. 
# $Header$
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2012
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
# Author:	Pierangelo Masarati	<masarati@aero.polimi.it>
# 
# Skip undesired labels
function skip_label(l)
{
	return 0;
}

BEGIN {
	FirstStep = 1;
	FirstLabel = -1;
} 
{
	# if first label, end line
	if ($1 == FirstLabel) {
		printf("\n");
	}

	# set first label
	if (FirstStep) {
		FirstLabel = $1
		FirstStep = 0;
	}

	# skip undesired labels (modify the body of the function at will)
	if (skip_label($1)) {
		next;
	}

	# collect values
	if (last == 0) {
		l = NF;
	} else {
		l = last;
	}

	for (j = 2; j <= l; j++) {
		printf(" %13.6e", $(j));
	}
}
END {
	# end last line
	printf("\n");
}
