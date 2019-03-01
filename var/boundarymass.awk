# $Header$
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
# bm.awk: helper to extract data for boundary mass calculations
# requires "node" to be set to the label of the node that's being
# modified.

BEGIN {
	if (node == 0) {
		print "\"node\" must be set";
		exit;
	}
}
/ RECORD GROUP 1, HEADER/ {
	getline;
	getline;
	nnodes = $2;
	nmodes = $3;
}
/ RECORD GROUP 2, FINITE ELEMENT NODE LIST/ {
	cnt = 0;
	gotit = 0;
	for (; gotit == 0; ) {
		getline;
		if (substr($0, 1, 2) == "**") {
			break;
		}
		for (col = 1; col <= NF; col++) {
			if ($col == node) {
				gotit = 1;
				cnt += col - NF;
				break;
			}
		}
		cnt += NF;
	}

	if (gotit == 0) {
		printf "node %d not found\n", node;
		exit
	}
	printf "% node=%d index=%d\n", node, cnt;
}
/ RECORD GROUP 5, NODAL X COORDINATES/ {
	for (row = 1; row <= cnt; row++) {
		getline;
	}
	X = $1;
}
/ RECORD GROUP 6, NODAL Y COORDINATES/ {
	for (row = 1; row <= cnt; row++) {
		getline;
	}
	Y = $1;
}
/ RECORD GROUP 7, NODAL Z COORDINATES/ {
	for (row = 1; row <= cnt; row++) {
		getline;
	}
	Y = $1;
	printf "% node=%d position=[%e, %e, %e];\n", node, X, Y, Z;
}
/    NORMAL MODE SHAPE #/ {
	mode = $6;
	printf "X%d = [ ...\n", mode;
	for (row = 1; row <= nnodes; row++) {
		getline;
		if (row == cnt) {
			shape[mode] = $0;
		}
		printf "        %s ...\n", $0;
	}
	printf "];\n";
}
/ RECORD GROUP 9, MODAL MASS MATRIX/ {
	printf "M = [ ...\n";
	for (row = 1; row <= nmodes; row++) {
		getline;
		printf "        %s; ...\n", $0;
	}
	printf "];\n";
}
/ RECORD GROUP 10, MODAL STIFFNESS MATRIX/ {
	printf "K = [ ...\n";
	for (row = 1; row <= nmodes; row++) {
		getline;
		printf "        %s; ...\n", $0;
	}
	printf "];\n";
}
/ RECORD GROUP 11, DIAGONAL OF LUMPED MASS MATRIX/ {
	printf "X = [ ...\n";
	for (row = 1; row <= nmodes; row++) {
		printf "        %s; ...\n", shape[row];
	}
	printf "];\n";
	for (row = 1; row <= cnt; row++) {
		getline;
	}
	printf "BM = diag([%s]);\n", $0;
	printf "M = M - X*BM*X';\n";
	printf "[v, l] = eig(M\\K);\n";
	printf "[l, I] = sort(diag(l));\n";
	printf "v = v(:, I);\n";
	printf "v = v*diag(1./sqrt(diag(v'*M*v)));\n";
	printf "m = diag(diag(v'*M*v));\n";
	printf "k = diag(diag(v'*K*v));\n";
	printf "x = v'*X;\n";
	printf "XX = [ ...\n";
	for (col = 1; col <= nmodes; col++) {
		printf "        X%d' ...\n", col;
	}
	printf "];\n";
	printf "xx = XX*v;\n";

	printf "disp('** RECORD GROUP 8, MODE SHAPES')\n";
	printf "for cnt = 1:%d,\n", nmodes;
	printf "        disp(sprintf('**    NORMAL MODE SHAPE #%%3d', cnt));\n";
	printf "        for row = 1:%d,\n", nnodes;
	printf "                disp(sprintf('%% 18.10E', xx(6*row-6+1:6*row, cnt)));\n";
	printf "        end;\n";
	printf "end\n";
	printf "disp('**')\n";

	printf "disp('** RECORD GROUP 9, MODAL MASS MATRIX')\n";
	printf "for row = 1:%d,\n", nmodes;
	printf "        disp(sprintf('%% 18.10E', m(row, :)));\n";
	printf "end\n";
	printf "disp('**')\n";

	printf "disp('** RECORD GROUP 10, MODAL STIFFNESS MATRIX')\n";
	printf "for row = 1:%d,\n", nmodes;
	printf "        disp(sprintf('%% 18.10E', k(row, :)));\n";
	printf "end\n";
	printf "disp('**')\n";
}
