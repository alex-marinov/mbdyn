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
# Prepare MBDyn ASCII output (.log & .mov files) for loading within EasyAnim
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
# Used by mbdyn2easyanim.sh; don't use directly.
#
# Recognized vars (set with -v <var>=<value>):
#
# showNone		suppresses all additional entities generation
# showAll		enables all additional entities generation
# showHinge		generates node at hinge and offsets for hinges
# noShowDistance	suppresses generation of distance joints
# noShowRod		suppresses generation of rod joints
# showDeformableJoint	generates nodes at joints and offsets for def. joint
# noShowBeam		suppresses generation of beam elements
# noShowAero		suppresses generation of aerodynamic elements
# every			outputs data every many time steps
# triad_length		if > 0, triads are drawn for each node with that length

# multiplies a 3x3 matrix times a 3 vector, adds a second 3 vector
# and stores result in another vector
function mat3_mul_vec3_add_vec3(m, v1, v2, r)
{
	r[1] = m[1, 1]*v1[1] + m[1, 2]*v1[2] + m[1, 3]*v1[3] + v2[1];
	r[2] = m[2, 1]*v1[1] + m[2, 2]*v1[2] + m[2, 3]*v1[3] + v2[2];
	r[3] = m[3, 1]*v1[1] + m[3, 2]*v1[2] + m[3, 3]*v1[3] + v2[3];
}

# trasforms 3 angles as output by MBDyn (1,2,3 sequence) into R matrix
function euler123_2R(Alpha, Beta, Gamma, R,   dCosAlpha, dSinAlpha, dCosBeta, dSinBeta, dCosGamma, dSinGamma) {
	dCosAlpha = cos(Alpha*AngleScale);
	dSinAlpha = sin(Alpha*AngleScale);
	dCosBeta = cos(Beta*AngleScale);
	dSinBeta = sin(Beta*AngleScale);
	dCosGamma = cos(Gamma*AngleScale);
	dSinGamma = sin(Gamma*AngleScale);

	R[1, 1] = dCosBeta*dCosGamma;
	R[2, 1] = dCosAlpha*dSinGamma+dSinAlpha*dSinBeta*dCosGamma;
	R[3, 1] = dSinAlpha*dSinGamma-dCosAlpha*dSinBeta*dCosGamma;
	R[1, 2] = -dCosBeta*dSinGamma;
	R[2, 2] = dCosAlpha*dCosGamma-dSinAlpha*dSinBeta*dSinGamma;
	R[3, 2] = dSinAlpha*dCosGamma+dCosAlpha*dSinBeta*dSinGamma;
	R[1, 3] = dSinBeta;
	R[2, 3] = -dSinAlpha*dCosBeta;
	R[3, 3] = dCosAlpha*dCosBeta;
}

# trasforms 3 angles as output by MBDyn (3,1,3 sequence) into R matrix
function euler313_2R(Alpha, Beta, Gamma, R,   dCosAlpha, dSinAlpha, dCosBeta, dSinBeta, dCosGamma, dSinGamma) {
	dCosAlpha = cos(Alpha*AngleScale);
	dSinAlpha = sin(Alpha*AngleScale);
	dCosBeta = cos(Beta*AngleScale);
	dSinBeta = sin(Beta*AngleScale);
	dCosGamma = cos(Gamma*AngleScale);
	dSinGamma = sin(Gamma*AngleScale);

	R[1, 1] = dCosAlpha*dCosGamma - dSinAlpha*dCosBeta*dSinGamma;
	R[2, 1] = dSinAlpha*dCosGamma + dCosAlpha*dCosBeta*dSinGamma;
	R[3, 1] = dSinBeta*dSinGamma;
	R[1, 2] = -dCosAlpha*dSinGamma - dSinAlpha*dCosBeta*dCosGamma;
	R[2, 2] = -dSinAlpha*dSinGamma + dCosAlpha*dCosBeta*dCosGamma;
	R[3, 2] = dSinBeta*dCosGamma;
	R[1, 3] = dSinAlpha*dSinBeta;
	R[2, 3] = -dCosAlpha*dSinBeta;
	R[3, 3] = dCosBeta;
}

# trasforms 3 angles as output by MBDyn (3,2,1 sequence) into R matrix
function euler321_2R(Alpha, Beta, Gamma, R,   dCosAlpha, dSinAlpha, dCosBeta, dSinBeta, dCosGamma, dSinGamma) {
	dCosAlpha = cos(Alpha*AngleScale);
	dSinAlpha = sin(Alpha*AngleScale);
	dCosBeta = cos(Beta*AngleScale);
	dSinBeta = sin(Beta*AngleScale);
	dCosGamma = cos(Gamma*AngleScale);
	dSinGamma = sin(Gamma*AngleScale);

	R[1, 1] = dCosAlpha*dCosBeta;
	R[2, 1] = dSinAlpha*dCosBeta;
	R[3, 1] = -dSinBeta;
	R[1, 2] = -dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma;
	R[2, 2] = dCosAlpha*dCosGamma + dSinAlpha*dSinBeta*dSinGamma;
	R[3, 2] = dCosBeta*dSinGamma;
	R[1, 3] = dSinAlpha*dSinGamma + dCosAlpha*dSinBeta*dCosGamma;
	R[2, 3] = -dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma;
	R[3, 3] = dCosBeta*dCosGamma;
}

# trasforms 3 angles as output by MBDyn (1,2,3 sequence) into R matrix
function phi2R(phi1, phi2, phi3, R,   dPhi, da, db) {
	dPhi = sqrt(phi1*phi1 + phi2*phi2 + phi3*phi3);
	if (dPhi > 1e-6) {
		da = sin(dPhi)/dPhi;
		db = (1 - cos(dPhi))/(dPhi*dPhi);

		R[1, 1] = 1. - db*(phi2*phi2 + phi3*phi3);
		R[2, 1] = da*phi3 + db*phi2*phi1;
		R[3, 1] = -da*phi2 + db*phi3*phi1;
		R[1, 2] = -da*phi3 + db*phi1*phi2;
		R[2, 2] = 1. - db*(phi3*phi3 + phi1*phi1);;
		R[3, 2] = da*phi1 + db*phi3*phi2;
		R[1, 3] = da*phi2 + db*phi1*phi3;
		R[2, 3] = -da*phi1 + db*phi2*phi3;
		R[3, 3] = 1. - db*(phi1*phi1 + phi2*phi2);;

	} else if (dPhi > 1e-12) {
		R[1, 1] = 1.;
		R[2, 1] = phi3;
		R[3, 1] = -phi2;
		R[1, 2] = -phi3;
		R[2, 2] = 1.;
		R[3, 2] = phi1;
		R[1, 3] = phi2;
		R[2, 3] = -phi1;
		R[3, 3] = 1.;

	} else {
		R[1, 1] = 1.;
		R[2, 1] = 0.;
		R[3, 1] = 0.;
		R[1, 2] = 0.;
		R[2, 2] = 1.;
		R[3, 2] = 0.;
		R[1, 3] = 0.;
		R[2, 3] = 0.;
		R[3, 3] = 1.;
	}
}

function nodeprop_add(name, color, radius) {
	nodeprop[nodeprop_num, "name"] = name;
	nodeprop[nodeprop_num, "color"] = color;
	nodeprop[nodeprop_num, "radius"] = radius;
	nodeprop_num++;
}

function node_absolute_add(name, x, y, z, prop) {
	node[node_num] = name;
	node[node_num, "relative"] = -1;
	node[node_num, 1] = x;
	node[node_num, 2] = y;
	node[node_num, 3] = z;
	strnode[name, 1] = x;
	strnode[name, 2] = y;
	strnode[name, 3] = z;
	node[node_num, "prop"] = prop;
	node_num++;
}

function node_add(name, relative, x, y, z, prop) {
	node[node_num] = name;
	node[node_num, "relative"] = relative;
	node[node_num, 1] = x;
	node[node_num, 2] = y;
	node[node_num, 3] = z;
	node[node_num, "prop"] = prop;
	node_num++;
}

function edgeprop_add(name, color, radius) {
	edgeprop[edgeprop_num, "name"] = name;
	edgeprop[edgeprop_num, "color"] = color;
	edgeprop[edgeprop_num, "radius"] = radius;
	edgeprop_num++;
}

function edge_add(name, node1label, node2label, prop) {
	edge[edge_num] = name;
	edge[edge_num, 1] = node1label;
	edge[edge_num, 2] = node2label;
	edge[edge_num, "prop"] = prop;
	edge_num++;
}

function sideprop_add(name, color) {
	sideprop[sideprop_num, "name"] = name;
	sideprop[sideprop_num, "color"] = color;
	sideprop_num++;
}

function side_add(name, numnodes, nodelabels, prop) {
	side[side_num] = name;
	side[side_num, "N"] = numnodes;
	for (i = 1; i <= numnodes; i++) {
		# check if exists?
		if (!(i in nodelabels)) {
			print "node #" i " of " numnodes " missing while adding side\"" name "\"" > /dev/stderr;
			exit 1;
		}
		side[side_num, i] = nodelabels[i];
	}
	side[side_num, "prop"] = prop;
	side_num++;
}

BEGIN {
	isvan = 0;
	steps = 0;

	if (every) {
		if (every < 0) {
			every = 0;

		} else {
			every--;
		}
	}

	# Note: by default, distance, rod, beam and aero elements are output
	# generic hinges and deformable joints can be optionally output
	# by setting "-v showHinge=1" and "-v showDeformableJoint=1"
	# to disable all output but that of structural nodes, create a .awk
	# that sets the related members of "show" to 0 and include it with -f.

	# safe defaults?
	show["hinge"] = 0;
	show["deformablejoint"] = show["hinge"];
	show["totaljoint"] = show["hinge"];
	show["totalpinjoint"] = show["hinge"];
	show["distance"] = 1;
	show["rod"] = 1;
	show["beam"] = 1;
	show["shell"] = 1;
	show["aero"] = 1;

	if (showNone) {
		for (i in show) {
			show[i] = 0;
		}

	} else if (showAll) {
		for (i in show) {
			show[i] = 1;
		}

	} else {
		if (showHinge) {
			show["hinge"] = 1;
		}

		if (noShowDistance) {
			show["distance"] = 0;
		}

		if (noShowRod) {
			show["rod"] = 0;
		}

		if (showDeformableJoint) {
			show["deformablejoint"] = 1;
		}

		if (showTotalJoint) {
			show["totaljoint"] = 1;
		}

		if (showTotalPinJoint) {
			show["totalpinjoint"] = 1;
		}

		if (noShowBeam) {
			show["beam"] = 0;
		}

		if (noShowShell) {
			show["shell"] = 0;
		}

		if (noShowAero) {
			show["aero"] = 0;
		}
	}

	node_num = 0;
	edge_num = 0;
	side_num = 0;

	strnode_num = 0;

	j_rod_num = 0;
	j_distance_num = 0;
	# more joints...

	beam2_num = 0;
	beam3_num = 0;

	shell4_num = 0;

	aero0_num = 0;
	aero2_num = 0;
	aero3_num = 0;

	deg2rad = 0.017453293;
	rad2deg = 57.29578;
	AngleScale = deg2rad;

	volfile = file ".vol";
	vanfile = file ".van";

	nodeprop_num = 0;
	edgeprop_num = 0;
	sideprop_num = 0;

	# node props
	nodeprop_add("hide", 0, 0);
	nodeprop_add("hinge_node", 2, 1.5);
	nodeprop_add("distance_node", 1, 1.);
	nodeprop_add("rod_node", 1, 1.);
	nodeprop_add("deformablejoint_node", 1, 1.);
	nodeprop_add("totaljoint_node", 2, 1.5);
	nodeprop_add("beam_node", 1, 1.);
	nodeprop_add("shell_node", 0, 0.);
	nodeprop_add("aero_node", 0, 0.);

	# edge props
	edgeprop_add("hinge_offset", 10, .5);
	edgeprop_add("distance_edge", 1, 1.);
	edgeprop_add("distance_offset", 12, .5);
	edgeprop_add("rod_edge", 1, 1.);
	edgeprop_add("rod_offset", 12, .5);
	edgeprop_add("deformablejoint_edge", 5, 2.5);
	edgeprop_add("deformablejoint_offset", 12, .5);
	edgeprop_add("totaljoint_edge", 1, 1.);
	edgeprop_add("totaljoint_offset", 10, .5);
	edgeprop_add("beam_edge", 14, 1.);
	edgeprop_add("beam_offset", 12, .5);
	edgeprop_add("triad_x", 4, 1.);
	edgeprop_add("triad_y", 2, 1.);
	edgeprop_add("triad_z", 1, 1.);

	# side props
	sideprop_add("shell", 12);
	sideprop_add("aero", 14);
}

isvan == 0 && /^structural node:/ {
	strnode_label[strnode_num] = $3;
	strnode[$3] = strnode_num;
	strnode[$3, 1] = $4;
	strnode[$3, 2] = $5;
	strnode[$3, 3] = $6;
	if ($7 == "phi") {
		strnode[$3, "orientation"] = "phi";
		strnode[$3, 4] = $8;
		strnode[$3, 5] = $9;
		strnode[$3, 6] = $10;

	} else if ($7 == "euler123") {
		strnode[$3, "orientation"] = "euler123";
		strnode[$3, 4] = $8;
		strnode[$3, 5] = $9;
		strnode[$3, 6] = $10;

	} else if ($7 == "euler313") {
		strnode[$3, "orientation"] = "euler313";
		strnode[$3, 4] = $8;
		strnode[$3, 5] = $9;
		strnode[$3, 6] = $10;

	} else if ($7 == "euler321") {
		strnode[$3, "orientation"] = "euler321";
		strnode[$3, 4] = $8;
		strnode[$3, 5] = $9;
		strnode[$3, 6] = $10;

	} else if ($7 == "mat") {
		strnode[$3, "orientation"] = "mat";
		strnode[$3, 4] = $8;
		strnode[$3, 5] = $9;
		strnode[$3, 6] = $10;
		strnode[$3, 7] = $11;
		strnode[$3, 8] = $12;
		strnode[$3, 9] = $13;
		strnode[$3, 10] = $14;
		strnode[$3, 11] = $15;
		strnode[$3, 12] = $16;

	} else {
		# error!
	}
	strnode_num++;

	if (!exclude["structural node", $3]) {
		node_absolute_add($3, $4, $5, $6, "default");

		if (triad_length > 0.) {
			node_add("triad_x_" $3, $3, triad_length, 0., 0., "hide");
			node_add("triad_y_" $3, $3, 0., triad_length, 0., "hide");
			node_add("triad_z_" $3, $3, 0., 0., triad_length, "hide");

			edge_add("triad_x_" $3, $3, "triad_x_" $3, "triad_x");
			edge_add("triad_y_" $3, $3, "triad_y_" $3, "triad_y");
			edge_add("triad_z_" $3, $3, "triad_z_" $3, "triad_z");
		}
	}
}

isvan == 0 && (/sphericalhinge:/ || /revolutehinge/ || /universalhinge/) && show["hinge"] {
	if (!exclude["joint", $2]) {
		j_hinge_label[j_hinge_num] = $2;
		j_hinge[$2] = j_hinge_num;
		j_hinge[$2, 1] = $3;
		j_hinge[$2, 1, 1] = $4;
		j_hinge[$2, 1, 2] = $5;
		j_hinge[$2, 1, 3] = $6;
		j_hinge[$2, 2] = $16;
		j_hinge[$2, 2, 1] = $17;
		j_hinge[$2, 2, 2] = $18;
		j_hinge[$2, 2, 3] = $19;
		j_hinge_num++;

		label1 = j_hinge[$2, 1];
		label2 = j_hinge[$2, 2];

		if (j_hinge[$2, 1, 1] != 0. || j_hinge[$2, 1, 2] != 0. || j_hinge[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "hinge_" $2 "_point1";
			node_add(label, j_hinge[$2, 1],
				j_hinge[$2, 1, 1],
				j_hinge[$2, 1, 2],
				j_hinge[$2, 1, 3],
				"hinge_node");

			edge_add("hinge_" $2 "_offset1", label1, label, "hinge_offset");

			label1 = label;
		}

		if (j_hinge[$2, 2, 1] != 0. || j_hinge[$2, 2, 2] != 0. || j_hinge[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "hinge_" $2 "_point2";
			node_add(label, j_hinge[$2, 2],
				j_hinge[$2, 2, 1],
				j_hinge[$2, 2, 2],
				j_hinge[$2, 2, 3],
				"hinge_node");

			edge_add("hinge_" $2 "_offset2", label2, label, "hinge_offset");

			label2 = label;
		}
	}
}

isvan == 0 && /^distance:/ && show["distance"] {
	if (!exclude["joint", $2]) {
		j_distance_label[j_distance_num] = $2;
		j_distance[$2] = j_distance_num;
		j_distance[$2, 1] = $3;
		j_distance[$2, 1, 1] = $4;
		j_distance[$2, 1, 2] = $5;
		j_distance[$2, 1, 3] = $6;
		j_distance[$2, 2] = $7;
		j_distance[$2, 2, 1] = $8;
		j_distance[$2, 2, 2] = $9;
		j_distance[$2, 2, 3] = $10;
		j_distance_num++;

		label1 = j_distance[$2, 1];
		label2 = j_distance[$2, 2];

		if (j_distance[$2, 1, 1] != 0. || j_distance[$2, 1, 2] != 0. || j_distance[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "distance_" $2 "_point1";
			node_add(label, j_distance[$2, 1],
				j_distance[$2, 1, 1],
				j_distance[$2, 1, 2],
				j_distance[$2, 1, 3],
				"distance_node");

			edge_add("distance_" $2 "_offset1", label1, label, "distance_offset");

			label1 = label;
		}

		if (j_distance[$2, 2, 1] != 0. || j_distance[$2, 2, 2] != 0. || j_distance[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "distance_" $2 "_point2";
			node_add(label, j_distance[$2, 2],
				j_distance[$2, 2, 1],
				j_distance[$2, 2, 2],
				j_distance[$2, 2, 3],
				"distance_node");

			edge_add("distance_" $2 "_offset2", label2, label, "distance_offset");

			label2 = label;
		}

		edge_add("distance_" $2, label1, label2, "distance_edge");
	}
}

isvan == 0 && /^rod:/ && show["rod"]{
	if (!exclude["joint", $2]) {
		j_rod_label[j_rod_num] = $2;
		j_rod[$2] = j_rod_num;
		j_rod[$2, 1] = $3;
		j_rod[$2, 1, 1] = $4;
		j_rod[$2, 1, 2] = $5;
		j_rod[$2, 1, 3] = $6;
		j_rod[$2, 2] = $7;
		j_rod[$2, 2, 1] = $8;
		j_rod[$2, 2, 2] = $9;
		j_rod[$2, 2, 3] = $10;
		j_rod_num++;

		label1 = j_rod[$2, 1];
		label2 = j_rod[$2, 2];

		if (j_rod[$2, 1, 1] != 0. || j_rod[$2, 1, 2] != 0. || j_rod[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "rod_" $2 "_point1";
			node_add(label, j_rod[$2, 1],
				j_rod[$2, 1, 1],
				j_rod[$2, 1, 2],
				j_rod[$2, 1, 3],
				"rod_node");

			edge_add("rod_" $2 "_offset1", label1, label, "rod_offset");

			label1 = label;
		}

		if (j_rod[$2, 2, 1] != 0. || j_rod[$2, 2, 2] != 0. || j_rod[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "rod_" $2 "_point2";
			node_add(label, j_rod[$2, 2],
				j_rod[$2, 2, 1],
				j_rod[$2, 2, 2],
				j_rod[$2, 2, 3],
				"rod_node");

			edge_add("rod_" $2 "_offset2", label2, label, "rod_offset");

			label2 = label;
		}

		edge_add("rod_" $2, label1, label2, "rod_edge");
	}
}

isvan == 0 && (/deformablejoint:/ || /deformabledisplacementjoint:/) && show["deformablejoint"] {
	if (!exclude["joint", $2]) {
		j_deformablejoint_label[j_deformablejoint_num] = $2;
		j_deformablejoint[$2] = j_deformablejoint_num;
		j_deformablejoint[$2, 1] = $3;
		j_deformablejoint[$2, 1, 1] = $4;
		j_deformablejoint[$2, 1, 2] = $5;
		j_deformablejoint[$2, 1, 3] = $6;
		j_deformablejoint[$2, 2] = $16;
		j_deformablejoint[$2, 2, 1] = $17;
		j_deformablejoint[$2, 2, 2] = $18;
		j_deformablejoint[$2, 2, 3] = $19;
		j_deformablejoint_num++;

		label1 = j_deformablejoint[$2, 1];
		label2 = j_deformablejoint[$2, 2];

		if (j_deformablejoint[$2, 1, 1] != 0. || j_deformablejoint[$2, 1, 2] != 0. || j_deformablejoint[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "deformablejoint_" $2 "_point1";
			node_add(label, j_deformablejoint[$2, 1],
				j_deformablejoint[$2, 1, 1],
				j_deformablejoint[$2, 1, 2],
				j_deformablejoint[$2, 1, 3],
				"deformablejoint_node");

			edge_add("deformablejoint_" $2 "_offset1", label1, label, "deformablejoint_offset");

			label1 = label;
		}

		if (j_deformablejoint[$2, 2, 1] != 0. || j_deformablejoint[$2, 2, 2] != 0. || j_deformablejoint[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "deformablejoint_" $2 "_point2";
			node_add(label, j_deformablejoint[$2, 2],
				j_deformablejoint[$2, 2, 1],
				j_deformablejoint[$2, 2, 2],
				j_deformablejoint[$2, 2, 3],
				"deformablejoint_node");

			edge_add("deformablejoint_" $2 "_offset2", label2, label, "deformablejoint_offset");

			label2 = label;
		}

		edge_add("deformablejoint_" $2, label1, label2, "deformablejoint_edge");
	}
}

isvan == 0 && (/^totaljoint:/ || /^totalequation:/ || /^totalreaction/) && show["totaljoint"] {
	if (!exclude["joint", $2]) {
		j_totaljoint_label[j_totaljoint_num] = $2;
		j_totaljoint[$2] = j_totaljoint_num;
		j_totaljoint[$2, 1] = $3;
		j_totaljoint[$2, 1, 1] = $4;
		j_totaljoint[$2, 1, 2] = $5;
		j_totaljoint[$2, 1, 3] = $6;
		j_totaljoint[$2, 2] = $25;
		j_totaljoint[$2, 2, 1] = $26;
		j_totaljoint[$2, 2, 2] = $27;
		j_totaljoint[$2, 2, 3] = $28;
		j_totaljoint_num++;

		label1 = j_totaljoint[$2, 1];
		label2 = j_totaljoint[$2, 2];

		if (j_totaljoint[$2, 1, 1] != 0. || j_totaljoint[$2, 1, 2] != 0. || j_totaljoint[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "totaljoint_" $2 "_point1";
			node_add(label, j_totaljoint[$2, 1],
				j_totaljoint[$2, 1, 1],
				j_totaljoint[$2, 1, 2],
				j_totaljoint[$2, 1, 3],
				"totaljoint_node");

			edge_add("totaljoint_" $2 "_offset1", label1, label, "totaljoint_offset");

			label1 = label;
		}

		if (j_totaljoint[$2, 2, 1] != 0. || j_totaljoint[$2, 2, 2] != 0. || j_totaljoint[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "totaljoint_" $2 "_point2";
			node_add(label, j_totaljoint[$2, 2],
				j_totaljoint[$2, 2, 1],
				j_totaljoint[$2, 2, 2],
				j_totaljoint[$2, 2, 3],
				"totaljoint_node");

			edge_add("totaljoint_" $2 "_offset2", label2, label, "totaljoint_offset");

			label2 = label;
		}

		edge_add("totaljoint_" $2, label1, label2, "totaljoint_edge");
	}
}

isvan == 0 && /^totalpinjoint:/ && show["totalpinjoint"] {
	if (!exclude["joint", $2]) {
		j_totalpinjoint_label[j_totalpinjoint_num] = $2;
		j_totalpinjoint[$2] = j_totalpinjoint_num;
		j_totalpinjoint[$2, 1] = $3;
		j_totalpinjoint[$2, 1, 1] = $4;
		j_totalpinjoint[$2, 1, 2] = $5;
		j_totalpinjoint[$2, 1, 3] = $6;
		j_totalpinjoint[$2, 2, 1] = $25;
		j_totalpinjoint[$2, 2, 2] = $26;
		j_totalpinjoint[$2, 2, 3] = $27;
		j_totalpinjoint_num++;

		label1 = j_totalpinjoint[$2, 1];

		if (j_totalpinjoint[$2, 1, 1] != 0. || j_totalpinjoint[$2, 1, 2] != 0. || j_totalpinjoint[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "totalpinjoint_" $2 "_point1";
			node_add(label, label1,
				j_totalpinjoint[$2, 1, 1],
				j_totalpinjoint[$2, 1, 2],
				j_totalpinjoint[$2, 1, 3],
				"totaljoint_node");

			edge_add("totalpinjoint_" $2 "_offset1", label1, label, "totaljoint_offset");

			label1 = label;
		}

		# create node and link
		label2 = "totalpinjoint_" $2 "_point2";
		node_absolute_add(label2, j_totalpinjoint[$2, 2, 1], j_totalpinjoint[$2, 2, 2], j_totalpinjoint[$2, 2, 3], "totaljoint_node");

		edge_add("totalpinjoint_" $2, label1, label2, "totaljoint_edge");
	}
}

isvan == 0 && /^beam2:/ && show["beam"] {
	if (!exclude["beam", $2]) {
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

		label1 = $3;
		label2 = $7;

		if (beam2[$2, 1, 1] != 0. || beam2[$2, 1, 2] != 0. || beam2[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "beam_" $2 "_point1";
			node_add(label, beam2[$2, 1],
				beam2[$2, 1, 1],
				beam2[$2, 1, 2],
				beam2[$2, 1, 3],
				"beam_node");

			edge_add("beam_" $2 "_offset1", label1, label, "beam_offset");

			label1 = label;
		}

		if (beam2[$2, 2, 1] != 0. || beam2[$2, 2, 2] != 0. || beam2[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "beam_" $2 "_point2";
			node_add(label, beam2[$2, 2],
				beam2[$2, 2, 1],
				beam2[$2, 2, 2],
				beam2[$2, 2, 3],
				"beam_node");

			edge_add("beam_" $2 "_offset2", label2, label, "beam_offset");

			label2 = label;
		}

		edge_add("beam_" $2, label1, label2, "beam_edge");
	}
}

isvan == 0 && /^beam3:/ && show["beam"] {
	if (!exclude["beam", $2]) {
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

		label1 = $3;
		label2 = $7;
		label3 = $11;

		if (beam3[$2, 1, 1] != 0. || beam3[$2, 1, 2] != 0. || beam3[$2, 1, 3] != 0.) {
			# create offset node and link
			label = "beam_" $2 "_point1";
			node_add(label, beam3[$2, 1],
				beam3[$2, 1, 1],
				beam3[$2, 1, 2],
				beam3[$2, 1, 3],
				"beam_node");

			edge_add("beam_" $2 "_offset1", label1, label, "beam_offset");

			label1 = label;
		}

		if (beam3[$2, 2, 1] != 0. || beam3[$2, 2, 2] != 0. || beam3[$2, 2, 3] != 0.) {
			# create offset node and link
			label = "beam_" $2 "_point2";
			node_add(label, beam3[$2, 2],
				beam3[$2, 2, 1],
				beam3[$2, 2, 2],
				beam3[$2, 2, 3],
				"beam_node");

			edge_add("beam_" $2 "_offset2", label2, label, "beam_offset");

			label2 = label;
		}

		if (beam3[$2, 3, 1] != 0. || beam3[$2, 3, 2] != 0. || beam3[$2, 3, 3] != 0.) {
			# create offset node and link
			label = "beam_" $2 "_point3";
			node_add(label, beam3[$2, 3],
				beam3[$2, 3, 1],
				beam3[$2, 3, 2],
				beam3[$2, 3, 3],
				"beam_node");

			edge_add("beam_" $2 "_offset3", label3, label, "beam_offset");

			label3 = label;
		}

		edge_add("beam_" $2 "_1", label1, label2, "beam_edge");
		edge_add("beam_" $2 "_2", label2, label3, "beam_edge");
	}
}

isvan == 0 && (/^shell4:/ || /^membrane4:/) && show["shell"] {
	if (!exclude["shell", $2]) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by shell4("$2") as node 1 not found" > "/dev/stderr";
			exit;
		}
		if (!($4 in strnode)) {
			print "structural node("$4") requested by shell4("$2") as node 2 not found" > "/dev/stderr";
			exit;
		}
		if (!($5 in strnode)) {
			print "structural node("$5") requested by shell4("$2") as node 3 not found" > "/dev/stderr";
			exit;
		}
		if (!($6 in strnode)) {
			print "structural node("$6") requested by shell4("$2") as node 4 not found" > "/dev/stderr";
			exit;
		}
		shell4_label[shell4_num] = $2;
		shell4[$2, 1] = $3;
		shell4[$2, 2] = $4;
		shell4[$2, 3] = $5;
		shell4[$2, 4] = $6;
		shell4_num++;

		label1 = $3;
		label2 = $4;
		label3 = $5;
		label4 = $6;

		label = "shell_" $2;
		side[side_num] = label;
		side[side_num, "N"] = 4;
		side[side_num, 1] = label1;
		side[side_num, 2] = label2;
		side[side_num, 3] = label3;
		side[side_num, 4] = label4;
		side[side_num, "prop"] = "shell";
		side_num++;
	}
}

isvan == 0 && /^aero0:/ && show["aero"] {
	if (!exclude["aerodynamic element", $2]) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero0("$2") as node not found" > "/dev/stderr";
			exit;
		}

		#   2-----4
		#   |  1  |
		#   1-----3

		aero0_label[aero0_num] = $2;
		aero0[$2] = $3;
		aero0[$2, 1, 1] = $4;
		aero0[$2, 1, 2] = $5;
		aero0[$2, 1, 3] = $6;
		aero0[$2, 2, 1] = $7;
		aero0[$2, 2, 2] = $8;
		aero0[$2, 2, 3] = $9;
		aero0[$2, 3, 1] = $10;
		aero0[$2, 3, 2] = $11;
		aero0[$2, 3, 3] = $12;
		aero0[$2, 4, 1] = $13;
		aero0[$2, 4, 2] = $14;
		aero0[$2, 4, 3] = $15;
		aero0_num++;

		# create offset node and side
		label1 = "aero_" $2 "_point1";
		node_add(label1, aero0[$2],
			aero0[$2, 1, 1],
			aero0[$2, 1, 2],
			aero0[$2, 1, 3],
			"aero_node");

		label2 = "aero_" $2 "_point2";
		node_add(label2, aero0[$2],
			aero0[$2, 2, 1],
			aero0[$2, 2, 2],
			aero0[$2, 2, 3],
			"aero_node");

		label3 = "aero_" $2 "_point3";
		node_add(label3, aero0[$2],
			aero0[$2, 3, 1],
			aero0[$2, 3, 2],
			aero0[$2, 3, 3],
			"aero_node");

		label4 = "aero_" $2 "_point4";
		node_add(label4, aero0[$2],
			aero0[$2, 4, 1],
			aero0[$2, 4, 2],
			aero0[$2, 4, 3],
			"aero_node");

		label = "aero_" $2;
		side[side_num] = label;
		side[side_num, "N"] = 4;
		side[side_num, 1] = label2;
		side[side_num, 2] = label4;
		side[side_num, 3] = label3;
		side[side_num, 4] = label1;
		side[side_num, "prop"] = "aero";
		side_num++;
	}
}

isvan == 0 && /^aero2:/ && show["aero"] {
	if (!exclude["aerodynamic element", $2]) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero2("$2") as node 1 not found" > "/dev/stderr";
			exit;
		}
		if (!($10 in strnode)) {
			print "structural node("$10") requested by aero2("$2") as node 2 not found" > "/dev/stderr";
			exit;
		}

		#   2-----4
		# 1 |     | 2
		#   1-----3

		aero2_label[aero2_num] = $2;
		aero2[$2, 1] = $3;
		aero2[$2, 1, 1] = $4;
		aero2[$2, 1, 2] = $5;
		aero2[$2, 1, 3] = $6;
		aero2[$2, 2, 1] = $7;
		aero2[$2, 2, 2] = $8;
		aero2[$2, 2, 3] = $9;
		aero2[$2, 2] = $10;
		aero2[$2, 3, 1] = $11;
		aero2[$2, 3, 2] = $12;
		aero2[$2, 3, 3] = $13;
		aero2[$2, 4, 1] = $14;
		aero2[$2, 4, 2] = $15;
		aero2[$2, 4, 3] = $16;
		aero2_num++;

		# create offset node and side
		label1 = "aero_" $2 "_point1";
		node_add(label1, aero2[$2, 1],
			aero2[$2, 1, 1],
			aero2[$2, 1, 2],
			aero2[$2, 1, 3],
			"aero_node");

		label2 = "aero_" $2 "_point2";
		node_add(label2, aero2[$2, 1],
			aero2[$2, 2, 1],
			aero2[$2, 2, 2],
			aero2[$2, 2, 3],
			"aero_node");

		label3 = "aero_" $2 "_point3";
		node_add(label3, aero2[$2, 2],
			aero2[$2, 3, 1],
			aero2[$2, 3, 2],
			aero2[$2, 3, 3],
			"aero_node");

		label4 = "aero_" $2 "_point4";
		node_add(label4, aero2[$2, 2],
			aero2[$2, 4, 1],
			aero2[$2, 4, 2],
			aero2[$2, 4, 3],
			"aero_node");

		label = "aero_" $2;
		side[side_num] = label;
		side[side_num, "N"] = 4;
		side[side_num, 1] = label2;
		side[side_num, 2] = label4;
		side[side_num, 3] = label3;
		side[side_num, 4] = label1;
		side[side_num, "prop"] = "aero";
		side_num++;
	}
}

isvan == 0 && /^aero3:/ && show["aero"] {
	if (!exclude["aerodynamic element", $2]) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero3("$2") as node 1 not found" > "/dev/stderr";
			exit;
		}
		if (!($10 in strnode)) {
			print "structural node("$10") requested by aero3("$2") as node 2 not found" > "/dev/stderr";
			exit;
		}
		if (!($17 in strnode)) {
			print "structural node("$17") requested by aero3("$2") as node 3 not found" > "/dev/stderr";
			exit;
		}

		#   2-----4-----6
		# 1 |     | 2   | 3
		#   1-----3-----5

		aero3_label[aero3_num] = $2;
		aero3[$2, 1] = $3;
		aero3[$2, 1, 1] = $4;
		aero3[$2, 1, 2] = $5;
		aero3[$2, 1, 3] = $6;
		aero3[$2, 2, 1] = $7;
		aero3[$2, 2, 2] = $8;
		aero3[$2, 2, 3] = $9;
		aero3[$2, 2] = $10;
		aero3[$2, 3, 1] = $11;
		aero3[$2, 3, 2] = $12;
		aero3[$2, 3, 3] = $13;
		aero3[$2, 4, 1] = $14;
		aero3[$2, 4, 2] = $15;
		aero3[$2, 4, 3] = $16;
		aero3[$2, 3] = $17;
		aero3[$2, 5, 1] = $18;
		aero3[$2, 5, 2] = $19;
		aero3[$2, 5, 3] = $20;
		aero3[$2, 6, 1] = $21;
		aero3[$2, 6, 2] = $22;
		aero3[$2, 6, 3] = $23;
		aero3_num++;

		# create offset node and side
		label1 = "aero_" $2 "_point1";
		node_add(label1, aero3[$2, 1],
			aero3[$2, 1, 1],
			aero3[$2, 1, 2],
			aero3[$2, 1, 3], "aero_node");

		label2 = "aero_" $2 "_point2";
		node_add(label2, aero3[$2, 1],
			aero3[$2, 2, 1],
			aero3[$2, 2, 2],
			aero3[$2, 2, 3],
			"aero_node");

		label3 = "aero_" $2 "_point3";
		node_add(label3, aero3[$2, 2],
			aero3[$2, 3, 1],
			aero3[$2, 3, 2],
			aero3[$2, 3, 3],
			"aero_node");

		label4 = "aero_" $2 "_point4";
		node_add(label4, aero3[$2, 2],
			aero3[$2, 4, 1],
			aero3[$2, 4, 2],
			aero3[$2, 4, 3],
			"aero_node");

		label5 = "aero_" $2 "_point5";
		node_add(label5, aero3[$2, 3],
			aero3[$2, 5, 1],
			aero3[$2, 5, 2],
			aero3[$2, 5, 3],
			"aero_node");

		label6 = "aero_" $2 "_point6";
		node_add(label6, aero3[$2, 3],
			aero3[$2, 6, 1],
			aero3[$2, 6, 2],
			aero3[$2, 6, 3],
			"aero_node");

		label = "aero_" $2;
		side[side_num] = label;
		side[side_num, "N"] = 6;
		side[side_num, 1] = label2;
		side[side_num, 2] = label4;
		side[side_num, 3] = label6;
		side[side_num, 4] = label5;
		side[side_num, 5] = label3;
		side[side_num, 6] = label1;
		side[side_num, "prop"] = "aero";
		side_num++;
	}
}

function node_pos(i, X) {
	if (node[i, "relative"] >= 0) {
		label = node[i, "relative"];
		v2[1] = strnode[label, 1];
		v2[2] = strnode[label, 2];
		v2[3] = strnode[label, 3];

		if (strnode[label, "orientation"] == "phi") {
			phi2R(strnode[label, 4], strnode[label, 5], strnode[label, 6], R);

		} else if (strnode[label, "orientation"] == "euler313") {
			euler313_2R(strnode[label, 4], strnode[label, 5], strnode[label, 6], R);

		} else if (strnode[label, "orientation"] == "euler321") {
			euler321_2R(strnode[label, 4], strnode[label, 5], strnode[label, 6], R);

		} else if (strnode[label, "orientation"] == "euler123") {
			euler123_2R(strnode[label, 4], strnode[label, 5], strnode[label, 6], R);

		} else if (strnode[label, "orientation"] == "mat") {
			R[1, 1] = strnode[label, 4];
			R[1, 2] = strnode[label, 5];
			R[1, 3] = strnode[label, 6];
			R[2, 1] = strnode[label, 7];
			R[2, 2] = strnode[label, 8];
			R[2, 3] = strnode[label, 9];
			R[3, 1] = strnode[label, 10];
			R[3, 2] = strnode[label, 11];
			R[3, 3] = strnode[label, 12];
		}

		v1[1] = node[i, 1];
		v1[2] = node[i, 2];
		v1[3] = node[i, 3];

		mat3_mul_vec3_add_vec3(R, v1, v2, X);

	} else {
		label = node[i];
		X[1] = strnode[label, 1];
		X[2] = strnode[label, 2];
		X[3] = strnode[label, 3];
	}
}

isvan == 0 && /^###/ {
	printf("# this is a comment\n") >> volfile;

	printf("#  0 = Black\n") >> volfile;
	printf("#  1 = Blue\n") >> volfile;
	printf("#  2 = Green\n") >> volfile;
	printf("#  3 = Cyan\n") >> volfile;
	printf("#  4 = Red\n") >> volfile;
	printf("#  5 = Magenta\n") >> volfile;
	printf("#  6 = Brown\n") >> volfile;
	printf("#  7 = Light Gray\n") >> volfile;
	printf("#  8 = Dark Gray\n") >> volfile;
	printf("#  9 = Light Blue\n") >> volfile;
	printf("# 10 = Light Green\n") >> volfile;
	printf("# 11 = Light Cyan\n") >> volfile;
	printf("# 12 = Light Red\n") >> volfile;
	printf("# 13 = Light Magenta\n") >> volfile;
	printf("# 14 = Yellow\n") >> volfile;
	printf("# 15 = White\n") >> volfile;

	printf("# node properties\n") >> volfile;
	for (i = 0; i < nodeprop_num; i++) {
		printf("prop %s %s %f\n", nodeprop[i, "name"], nodeprop[i, "color"], nodeprop[i, "radius"]) >> volfile;
	}

	printf("# nodes\n") >> volfile;
	printf("%d\n", node_num) >> volfile;
	for (i = 0; i < node_num; i++) {
		node_pos(i, X);
		printf("%s %e %e %e %s\n", node[i], X[1], X[2], X[3], node[i, "prop"]) >> volfile;
	}

	printf("# edge properties\n") >> volfile;
	for (i = 0; i < edgeprop_num; i++) {
		printf("prop %s %s %f\n", edgeprop[i, "name"], edgeprop[i, "color"], edgeprop[i, "radius"]) >> volfile;
	}

	printf("# edges\n") >> volfile;
	printf("%d\n", edge_num) >> volfile;
	for (i = 0; i < edge_num; i++) {
		printf("%s %s %s %s\n", edge[i], edge[i, 1], edge[i, 2], edge[i, "prop"]) >> volfile;
	}

	printf("# side properties\n") >> volfile;
	for (i = 0; i < sideprop_num; i++) {
		printf("prop %s %s\n", sideprop[i, "name"], sideprop[i, "color"]) >> volfile;
	}

	printf("# sides\n") >> volfile;
	printf("%d\n", side_num) >> volfile;
	for (i = 0; i < side_num; i++) {
		s = side[i, 1];
		for (j = 2; j <= side[i, "N"]; j++) {
			s = s " " side[i, j];
		}
		printf("%s %s %s %s\n", side[i], side[i, "N"], s, side[i, "prop"]) >> volfile;
	}

	# flush .vol
	fflush(volfile);

	isvan = 1;
	i = 0;
	getline;
}

# all
isvan == 1 {
	strnode[$1, 1] = $2;
	strnode[$1, 2] = $3;
	strnode[$1, 3] = $4;
	strnode[$1, 4] = $5;
	strnode[$1, 5] = $6;
	strnode[$1, 6] = $7;
	if (strnode[$1, "orientation"] == "mat") {
		strnode[$1, 7] = $8;
		strnode[$1, 8] = $9;
		strnode[$1, 9] = $10;
		strnode[$1, 10] = $11;
		strnode[$1, 11] = $12;
		strnode[$1, 12] = $13;
	}

	if (++i == strnode_num) {
		# compute nodes
		for (i = 0; i < node_num; i++) {
			node_pos(i, X);

			printf("%e %e %e\n", X[1], X[2], X[3]) >> vanfile;
		}
		steps++;
		i = 0;

		if (every) {
			for (j = 0; j < every*strnode_num; j++) {
				getline;
			}
		}
	}
}

END {
	if (steps > 0) {
		printf "successfully converted %d steps\n", steps;
	}
}
