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
# Prepare MBDyn's *.log file for loading within EasyAnim
#
# EasyAnim is a simple package for the visualization of mechanical systems
# it can be downloaded from http://mecara.fpms.ac.be/EasyDyn/
# For information, contact
#
# Prof. Olivier VERLINDEN
# Faculte Polytechnique de Mons
# Department of Theoretical Mechanics, Dynamics and Vibrations
# 31, Bd Dolez; 7000 MONS (BELGIUM)
# E-mail: Olivier.Verlinden@fpms.ac.be
#
# Usage:
#
# awk -f mbdyn2vol.awk <file>.log > <file>.vol

BEGIN {
	strnode_num = 0;
	beam2_num = 0;
	beam3_num = 0;
}

/structural node:/ {
	strnode_label[strnode_num] = $3;
	strnode[$3] = strnode_num;
	strnode[$3, 1] = $4;
	strnode[$3, 2] = $5;
	strnode[$3, 3] = $6;
	strnode_num++;
}

/beam2:/ {
	if (!($3 in strnode)) {
		print "structural node("$3") requested by beam2("$2") as node 1 not found" > "/dev/stderr";
		exit;
	}
	if (!($7 in strnode)) {
		print "structural node("$7") requested by beam2("$2") as node 2 not found" > "/dev/stderr";
		exit;
	}
	beam2_label[beam2_num] = $2;
	beam2[$2, 1] = $3;
	beam2[$2, 1, 1] = $4;
	beam2[$2, 1, 2] = $5;
	beam2[$2, 1, 3] = $6;
	beam2[$2, 2] = $7;
	beam2[$2, 2, 1] = $8;
	beam2[$2, 2, 2] = $9;
	beam2[$2, 2, 3] = $10;
	beam2_num++;
}

/beam3:/ {
	if (!($3 in strnode)) {
		print "structural node("$3") requested by beam3("$2") as node 1 not found" > "/dev/stderr";
		exit;
	}
	if (!($7 in strnode)) {
		print "structural node("$7") requested by beam3("$2") as node 2 not found" > "/dev/stderr";
		exit;
	}
	if (!($11 in strnode)) {
		print "structural node("$11") requested by beam3("$2") as node 3 not found" > "/dev/stderr";
		exit;
	}
	beam3_label[beam3_num] = $2;
	beam3[$2, 1] = $3;
	beam3[$2, 1, 1] = $4;
	beam3[$2, 1, 2] = $5;
	beam3[$2, 1, 3] = $6;
	beam3[$2, 2] = $7;
	beam3[$2, 2, 1] = $8;
	beam3[$2, 2, 2] = $9;
	beam3[$2, 2, 3] = $10;
	beam3[$2, 3] = $11;
	beam3[$2, 3, 1] = $12;
	beam3[$2, 3, 2] = $13;
	beam3[$2, 3, 3] = $14;
	beam3_num++;
}

END {
	printf("# this is a comment\n");
	printf("%d\n", strnode_num);
	for (i = 0; i < strnode_num; i++) {
		label = strnode_label[i];
		printf("%d %e %e %e\n", label, strnode[label, 1], strnode[label, 2], strnode[label, 3]);
	}
	printf("%d\n", beam2_num + 2*beam3_num);
	for (i = 0; i < beam2_num; i++) {
		label = beam2_label[i];
		printf("%d %d %d %d\n", label, beam2[label, 1], beam2[label, 2], 14);
	}
	for (i = 0; i < beam3_num; i++) {
		label = beam3_label[i];
		printf("%d_1 %d %d %d\n", label, beam3[label, 1], beam3[label, 2], 14);
		printf("%d_2 %d %d %d\n", label, beam3[label, 2], beam3[label, 3], 14);
	}
	printf("%d\n", 0);
}
