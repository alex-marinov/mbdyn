BEGIN {
	DY=0.01;
	DZ=0.01;

	edgeprop_add("edge_y", 4, 1.);
	edgeprop_add("edge_z", 2, 1.);

	for (i = 0; i <= 20; i++) {
		node_add(i "_y", i, 0., DY, 0., "hide");
		node_add(i "_z", i, 0., 0., DZ, "hide");

		edge_add(i "_y", i, i "_y", "edge_y");
		edge_add(i "_z", i, i "_z", "edge_z");
	}
}
