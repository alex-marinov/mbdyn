#!/bin/awk
# usage: awk -f GeoOut.awk <file>.mov
# -v INCR=<time step>		: time step (1.0)
# -v START=<initial time>	: initial time (0.0)
# -v SKIP=<skip>		: print every <skip> time steps
# -v STOP=<stop>		: stop after <stop> time steps

# Prints the labels
function print_preamble() {
	# triads
	printf("# Preamble\n");
	printf("object \"TriadLabels\" class array type int rank 1 shape 1 items %d data follows\n", strnode_num);
	for (i = 0; i < strnode_num; i++) {
		printf("\t%d\n", strnode_label[i]);
	}
	printf("attribute \"dep\" string \"positions\"\n");

	# beams
	beam_num = beam2_num + 2*beam3_num;
	if (beam_num > 0) {
		printf("object \"Beam2Connections\" class array type int shape 2 items %d data follows\n", beam_num);
		for (i = 0; i < beam2_num; i++) {
			l = beam2_label[i];
			printf("\t%d %d\n", beam2[l, 1], beam2[l, 2]);
		}
		for (i = 0; i < beam3_num; i++) {
			l = beam3_label[i];
			printf("\t%d %d\n", beam3[l, 1], beam3[l, 2]);
			printf("\t%d %d\n", beam3[l, 2], beam3[l, 3]);
		}
		printf("attribute \"element type\" string \"lines\"\n");
		printf("object \"Beam2Labels\" class array type int shape 1 items %d data follows\n", beam_num);
		for (i = 0; i < beam2_num; i++) {
			printf("\t%d\n", beam2_label[i]);
		}
		for (i = 0; i < beam3_num; i++) {
			printf("\t%d\n", beam3_label[i]);
			printf("\t%d\n", beam3_label[i]);
		}
		printf("attribute \"dep\" string \"connections\"\n");
		printf("object \"Beam2Offsets\" class array type float shape 6 items %d data follows\n", beam_num);
		for (i = 0; i < beam2_num; i++) {
			l = beam2_label[i];
			printf("\t%f %f %f %f %f %f\n", 
				beam2[l, 1, 1], beam2[l, 1, 2], beam2[l, 1, 3],
				beam2[l, 2, 1], beam2[l, 2, 2], beam2[l, 2, 3]);
		}
		for (i = 0; i < beam3_num; i++) {
			l = beam3_label[i];
			printf("\t%f %f %f %f %f %f\n", 
				beam3[l, 1, 1], beam3[l, 1, 2], beam3[l, 1, 3],
				beam3[l, 2, 1], beam3[l, 2, 2], beam3[l, 2, 3]);
			printf("\t%f %f %f %f %f %f\n", 
				beam3[l, 2, 1], beam3[l, 2, 2], beam3[l, 2, 3],
				beam3[l, 3, 1], beam3[l, 3, 2], beam3[l, 3, 3]);
		}
		printf("attribute \"dep\" string \"connections\"\n");
	}

	# aero
	aero_num = aero0_num + aero2_num + 2*aero3_num;
	if (aero_num > 0) {
		printf("object \"Aero4Connections\" class array type int shape 4 items %d data follows\n", aero_num);
		for (i = 0; i < aero0_num; i++) {
			l = aero0_label[i];
			printf("\t%d %d %d %d\n", aero0[l, 1], aero0[l, 1], aero0[l, 1], aero0[l, 1]);
		}
		for (i = 0; i < aero2_num; i++) {
			l = aero2_label[i];
			printf("\t%d %d %d %d\n", aero2[l, 1], aero2[l, 2], aero2[l, 1], aero2[l, 2]);
		}
		for (i = 0; i < aero3_num; i++) {
			l = aero3_label[i];
			printf("\t%d %d %d %d\n", aero3[l, 1], aero3[l, 2], aero3[l, 1], aero3[l, 2]);
			printf("\t%d %d %d %d\n", aero3[l, 2], aero3[l, 3], aero3[l, 2], aero3[l, 3]);
		}
		printf("attribute \"element type\" string \"quads\"\n");
		printf("object \"Aero4Labels\" class array type int shape 1 items %d data follows\n", aero_num);
		for (i = 0; i < aero0_num; i++) {
			printf("\t%d\n", aero0_label[i]);
		}
		for (i = 0; i < aero2_num; i++) {
			printf("\t%d\n", aero2_label[i]);
		}
		for (i = 0; i < aero3_num; i++) {
			printf("\t%d\n", aero3_label[i]);
			printf("\t%d\n", aero3_label[i]);
		}
		printf("attribute \"dep\" string \"connections\"\n");
		printf("object \"Aero4Offsets\" class array type float shape 12 items %d data follows\n", aero_num);
		for (i = 0; i < aero0_num; i++) {
			l = aero0_label[i];
			printf("\t%f %f %f %f %f %f %f %f %f %f %f %f\n", 
				aero0[l, 1, 1], aero0[l, 1, 2], aero0[l, 1, 3],
				aero0[l, 3, 1], aero0[l, 3, 2], aero0[l, 3, 3],
				aero0[l, 2, 1], aero0[l, 2, 2], aero0[l, 2, 3],
				aero0[l, 4, 1], aero0[l, 4, 2], aero0[l, 4, 3]);
		}
		for (i = 0; i < aero2_num; i++) {
			l = aero2_label[i];
			printf("\t%f %f %f %f %f %f %f %f %f %f %f %f\n", 
				aero2[l, 1, 1], aero2[l, 1, 2], aero2[l, 1, 3],
				aero2[l, 3, 1], aero2[l, 3, 2], aero2[l, 3, 3],
				aero2[l, 2, 1], aero2[l, 2, 2], aero2[l, 2, 3],
				aero2[l, 4, 1], aero2[l, 4, 2], aero2[l, 4, 3]);
		}
		for (i = 0; i < aero3_num; i++) {
			l = aero3_label[i];
			printf("\t%f %f %f %f %f %f %f %f %f %f %f %f\n", 
				aero3[l, 1, 1], aero3[l, 1, 2], aero3[l, 1, 3],
				aero3[l, 3, 1], aero3[l, 3, 2], aero3[l, 3, 3],
				aero3[l, 2, 1], aero3[l, 2, 2], aero3[l, 2, 3],
				aero3[l, 4, 1], aero3[l, 4, 2], aero3[l, 4, 3]);
			printf("\t%f %f %f %f %f %f %f %f %f %f %f %f\n", 
				aero3[l, 3, 1], aero3[l, 3, 2], aero3[l, 3, 3],
				aero3[l, 5, 1], aero3[l, 5, 2], aero3[l, 5, 3],
				aero3[l, 4, 1], aero3[l, 4, 2], aero3[l, 4, 3],
				aero3[l, 6, 1], aero3[l, 6, 2], aero3[l, 6, 3]);
		}
		printf("attribute \"dep\" string \"connections\"\n");
	}
}

# Prints the objects at one step
function print_step() {
	t = Start + Step*Skip*Incr;
	printf("# Step %f\n", t);
	printf("object \"TriadPositions%f\" class array type float rank 1 shape 3 items %d data follows\n", t, strnode_num);
	for (i = 0; i < strnode_num; i++) {
		l = strnode_label[i];
		printf("\t%f %f %f\n", 
			strnode_pos[l, 1], strnode_pos[l, 2], strnode_pos[l, 3]);
	}
	printf("object \"TriadCosines%f\" class array type float rank 1 shape 9 items %d data follows\n", t, strnode_num);
	for (i = 0; i < strnode_num; i++) {
		l = strnode_label[i];
		printf("\t%f %f %f %f %f %f %f %f %f\n", 
			strnode_cos[l, 1], strnode_cos[l, 2], strnode_cos[l, 3], 
			strnode_cos[l, 4], strnode_cos[l, 5], strnode_cos[l, 6], 
			strnode_cos[l, 7], strnode_cos[l, 8], strnode_cos[l, 9]);
	}
	printf("attribute \"dep\" string \"positions\"\n");
	printf("object \"OrientedTriads%f\" class field\n", t);
	printf("component \"positions\" \"TriadPositions%f\"\n", t);
	printf("component \"cosines\" \"TriadCosines%f\"\n", t);
	printf("component \"labels\" \"TriadLabels\"\n");
	if (beam_num > 0) {
		printf("object \"Beams2%f\" class field\n", t);
		printf("component \"positions\" \"TriadPositions%f\"\n", t);
		printf("component \"connections\" \"Beam2Connections\"\n");
		printf("component \"cosines\" \"TriadCosines%f\"\n", t);
		printf("component \"labels\" \"Beam2Labels\"\n");
		printf("component \"nodelabels\" \"TriadLabels\"\n");
		printf("component \"offsets\" \"Beam2Offsets\"\n");
	}
	if (aero_num > 0) {
		printf("object \"AeroQuads%f\" class field\n", t);
		printf("component \"positions\" \"TriadPositions%f\"\n", t);
		printf("component \"connections\" \"Aero4Connections\"\n");
		printf("component \"cosines\" \"TriadCosines%f\"\n", t);
		printf("component \"labels\" \"Aero4Labels\"\n");
		printf("component \"nodelabels\" \"TriadLabels\"\n");
		printf("component \"offsets\" \"Aero4Offsets\"\n");
	}
	printf("object \"MBDynSym%f\" class group\n", t);
	printf("member \"OrientedTriads\" value \"OrientedTriads%f\"\n", t);
	if (beam_num > 0) {
		printf("member \"Beams2\" value \"Beams2%f\"\n", t);
	}
	if (aero_num > 0) {
		printf("member \"AeroQuads\" value \"AeroQuads%f\"\n", t);
	}
}

# Prints the members of the group at the end
function print_members() {
	printf("# Series\n");
	printf("object \"series\" class series\n");
	for (i = 0; i <= Step; i++) {
		t = Start + i*Skip*Incr;
		printf("member %d position %f value \"MBDynSym%f\"\n", i, t, t);
	}
	printf("end\n");
}

# Customize rule to filter out undesired nodes/beams
function skip_node_label(l) {
	return 0;
}

function skip_beam_label(l) {
	return 0;
}

function skip_aero_label(l) {
	return 0;
}

# Initialize vars
BEGIN {
	FirstLabel = -1;
	FirstStep = 1;

	deg2rad = 1.74532925199433e-02;
	rad2deg = 5.72957795130823e+01;
	AngleScale = deg2rad;

	# time step
	Incr = 1.0;
	if (INCR != 0.0) {
		Incr = INCR;
	}

	# initial time
	Start = 0.0;
	if (START != 0.0) {
		Start = START;
	}

	# print every SKIP steps
	Skip = 1;
	if (SKIP != 0) {
		Skip = SKIP;
	}

        # stop after STOP steps 
        Stop = 10000000000000000000000000000000;
        if (STOP != 0) {
                Stop = STOP;
        }

	# Counters
	strnode_num = 0;
	beam2_num = 0;
	beam3_num = 0;
	aero0_num = 0;
	aero2_num = 0;
	aero3_num = 0;

	Preamble = 1;
	Node = 0;
}

# End of preamble 
/^======/ {
	print_preamble();
	Preamble = 0;
	getline;
}

/^structural node:/ {
	if (!Preamble) {
		print "got \"structural node\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_node_label($3)) {
		strnode_label[strnode_num] = $3;
		strnode[$3] = strnode_num;
		strnode_pos[$3, 1] = $4;
		strnode_pos[$3, 2] = $5;
		strnode_pos[$3, 3] = $6;
		strnode_num++;
	}
}

/^beam2:/ {
	if (!Preamble) {
		print "got \"beam2\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_beam_label($2)) {
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
}

/^beam3:/ {
	if (!Preamble) {
		print "got \"beam3\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_beam_label($2)) {
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
}

/^aero0:/ {
	if (!Preamble) {
		print "got \"aero0\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_aero_label($2)) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero0("$2") as node not found" > "/dev/stderr";
			exit;
		}

		aero0_label[aero0_num] = $2;
		aero0[$2, 1] = $3;
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
	}
}

/^aero2:/ {
	if (!Preamble) {
		print "got \"aero2\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_aero_label($2)) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero2("$2") as node not found" > "/dev/stderr";
			exit;
		}

		if (!($10 in strnode)) {
			print "structural node("$10") requested by aero2("$2") as node not found" > "/dev/stderr";
			exit;
		}

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
	}
}

/^aero3:/ {
	if (!Preamble) {
		print "got \"aero3\" out of preamble" > "/dev/stderr";
		exit;
	}

	if (!skip_aero_label($2)) {
		if (!($3 in strnode)) {
			print "structural node("$3") requested by aero3("$2") as node not found" > "/dev/stderr";
			exit;
		}

		if (!($10 in strnode)) {
			print "structural node("$10") requested by aero3("$2") as node not found" > "/dev/stderr";
			exit;
		}

		if (!($17 in strnode)) {
			print "structural node("$17") requested by aero3("$2") as node not found" > "/dev/stderr";
			exit;
		}

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
	}
}

# General rule
{
	if (!Preamble) {
		if ($1 == FirstLabel) {
			if (FirstStep) {
				FirstStep = 0;
			}
			if (!(blocks % Skip)) {
				print_step();
				Step++;
			}
			if (Step >= Stop) {
				print_members();
				exit;
			}
			blocks++;
		}
		if (blocks % Skip) {
			next;
		}

		if (FirstLabel == -1) {
			FirstLabel = $1;
		}

		if (!skip_node_label($1)) {
			strnode_pos[$1, 1] = $2;
			strnode_pos[$1, 2] = $3;
			strnode_pos[$1, 3] = $4;
	
			dCosAlpha = cos($5*AngleScale);
			dSinAlpha = sin($5*AngleScale);
			dCosBeta = cos($6*AngleScale);
			dSinBeta = sin($6*AngleScale);
			dCosGamma = cos($7*AngleScale);
			dSinGamma = sin($7*AngleScale);
			strnode_cos[$1, 1] = dCosBeta*dCosGamma;
			strnode_cos[$1, 2] = dCosAlpha*dSinGamma+dSinAlpha*dSinBeta*dCosGamma;
			strnode_cos[$1, 3] = dSinAlpha*dSinGamma-dCosAlpha*dSinBeta*dCosGamma;
			strnode_cos[$1, 4] = -dCosBeta*dSinGamma;
			strnode_cos[$1, 5] = dCosAlpha*dCosGamma-dSinAlpha*dSinBeta*dSinGamma;
			strnode_cos[$1, 6] = dSinAlpha*dCosGamma+dCosAlpha*dSinBeta*dSinGamma;
			strnode_cos[$1, 7] = dSinBeta;
			strnode_cos[$1, 8] = -dSinAlpha*dCosBeta;
			strnode_cos[$1, 9] = dCosAlpha*dCosBeta;
		}
	}
}

# Print last step and members
END {
	print_step();
	print_members();
}

